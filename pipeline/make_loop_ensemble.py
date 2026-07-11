#!/usr/bin/env python
"""
make_loop_ensemble.py — generate a loop-conformer ensemble for a membrane
protein target, holding the transmembrane + cytoplasmic core fixed and letting
the extracellular loops sample their accessible conformational space.

WHY (design decision, documented in DESIGN_DECISIONS.md):
  The extracellular loops of NKG7 carry no disulfides and are conformationally
  mobile, yet AlphaFold gives a single static snapshot. Designing binders
  against one frozen loop conformation risks shape-complementarity to a
  conformation the loop rarely adopts. We instead generate a small physically
  grounded ensemble (restrained MD: core frozen, loops free at 310 K) and
  distribute the design campaign across conformers, preferring binders robust
  across the ensemble.

METHOD: implicit-solvent (GBn2) restrained MD. Heavy atoms of the core
  (everything outside the loop windows) are position-restrained; loop atoms
  are free. Snapshots are harvested at fixed intervals after equilibration.

REUSABLE: driven entirely by CLI args — target PDB, loop residue windows,
  restraint strength, number of conformers. Not NKG7-specific.

Usage:
  python make_loop_ensemble.py \
      --pdb targets/NKG7_human_Q16617_AF.pdb \
      --loops 30-60 113-132 \
      --n-conformers 12 --out targets/ensembles/human
"""
import argparse, os, sys
import numpy as np


def parse_windows(specs):
    w = []
    for s in specs:
        a, b = s.split("-")
        w.append((int(a), int(b)))
    return w


def in_windows(resid, windows):
    return any(a <= resid <= b for a, b in windows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--loops", nargs="+", required=True,
                    help="loop residue windows, e.g. 30-60 113-132 (1-based, inclusive)")
    ap.add_argument("--n-conformers", type=int, default=12)
    ap.add_argument("--restraint-k", type=float, default=10.0,
                    help="core heavy-atom restraint force constant, kcal/mol/A^2")
    ap.add_argument("--temp", type=float, default=310.0, help="K")
    ap.add_argument("--equil-steps", type=int, default=25000)
    ap.add_argument("--sample-interval", type=int, default=10000,
                    help="MD steps between harvested conformers")
    ap.add_argument("--pad-core", type=int, default=1,
                    help="residues at each loop edge to also restrain (anchor the loop bases)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--platform", default="auto",
                    help="OpenMM platform: auto|CUDA|OpenCL|CPU|Reference (auto picks fastest available)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import openmm as mm
    import openmm.app as app
    import openmm.unit as u
    from pdbfixer import PDBFixer

    os.makedirs(args.out, exist_ok=True)
    windows = parse_windows(args.loops)
    # shrink loop windows by pad-core so the loop *bases* stay restrained (anchored)
    free_windows = [(a + args.pad_core, b - args.pad_core) for a, b in windows]
    print(f"[ensemble] target={args.pdb}  free loop windows={free_windows}", flush=True)

    # --- prepare structure: add missing atoms + hydrogens ---
    fixer = PDBFixer(filename=args.pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    prepared = os.path.join(args.out, "prepared.pdb")
    with open(prepared, "w") as fh:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, fh)

    # --- system: Amber14 + GBn2 implicit solvent ---
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    modeller = app.Modeller(fixer.topology, fixer.positions)
    system = ff.createSystem(modeller.topology,
                             nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.5 * u.nanometer,
                             constraints=app.HBonds)

    # --- position restraints on core heavy atoms ---
    restraint = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", args.restraint_k * u.kilocalories_per_mole / u.angstrom**2)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    positions = modeller.positions
    n_restrained = 0
    for atom in modeller.topology.atoms():
        resid = atom.residue.id
        try:
            resid = int(resid)
        except ValueError:
            continue
        if atom.element is not None and atom.element.symbol == "H":
            continue
        if not in_windows(resid, free_windows):  # core -> restrain
            x, y, z = positions[atom.index].value_in_unit(u.nanometer)
            restraint.addParticle(atom.index, [x, y, z])
            n_restrained += 1
    system.addForce(restraint)
    print(f"[ensemble] restrained {n_restrained} core heavy atoms; loops free", flush=True)

    integrator = mm.LangevinMiddleIntegrator(args.temp * u.kelvin,
                                             1.0 / u.picosecond,
                                             2.0 * u.femtosecond)
    integrator.setRandomNumberSeed(args.seed)

    def pick_platform(name):
        if name != "auto":
            return mm.Platform.getPlatformByName(name)
        avail = {mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())}
        for pref in ("CUDA", "OpenCL", "CPU", "Reference"):
            if pref in avail:
                return mm.Platform.getPlatformByName(pref)
        return mm.Platform.getPlatform(0)

    plat = pick_platform(args.platform)
    print(f"[ensemble] OpenMM platform: {plat.getName()}", flush=True)
    sim = app.Simulation(modeller.topology, system, integrator, plat)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(args.temp * u.kelvin, args.seed)
    sim.step(args.equil_steps)

    # --- harvest conformers ---
    conformers = []
    for i in range(args.n_conformers):
        sim.step(args.sample_interval)
        state = sim.context.getState(getPositions=True)
        pos = state.getPositions()
        out_pdb = os.path.join(args.out, f"conformer_{i:02d}.pdb")
        with open(out_pdb, "w") as fh:
            app.PDBFile.writeFile(sim.topology, pos, fh)
        conformers.append(out_pdb)
        print(f"[ensemble] conformer {i:02d} -> {out_pdb}", flush=True)

    print(f"[ensemble] DONE: {len(conformers)} conformers in {args.out}", flush=True)


if __name__ == "__main__":
    main()
