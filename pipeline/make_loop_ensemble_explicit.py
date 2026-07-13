#!/usr/bin/env python
"""
make_loop_ensemble_explicit.py — explicit-solvent restrained-MD loop ensemble.

ESCALATION of make_loop_ensemble.py (implicit GBn2) triggered when the cheap
multi-seed implicit-solvent basin-convergence check shows DIVERGENCE. This
variant swaps the physics model: TIP3P explicit water in a periodic box + PME
electrostatics + Monte-Carlo barostat (NPT), which is the physically faithful
test of whether the multi-basin loop behaviour is real or a GBn2 artifact.

Same restraint logic as the implicit script: core heavy atoms position-restrained,
loop atoms free; conformers harvested at fixed intervals after equilibration.

Key differences vs implicit:
  - ForceField: amber14-all.xml + amber14/tip3p.xml (explicit water)
  - modeller.addSolvent(padding, neutralising ions)
  - nonbondedMethod=PME, periodic, 1.0 nm cutoff
  - MonteCarloBarostat (1 bar) for NPT density equilibration
  - longer: 50 ns production default (25,000,000 steps @ 2 fs)
  - conformers stripped back to protein-only on write (water discarded)

Usage:
  python make_loop_ensemble_explicit.py \
      --pdb targets/NKG7_human_Q16617_AF.pdb --loops 30-60 113-132 \
      --n-conformers 20 --production-ns 50 --seed 1 --out runs/md_explicit/human
"""
import argparse, os
import numpy as np


def parse_windows(specs):
    return [(int(s.split("-")[0]), int(s.split("-")[1])) for s in specs]


def in_windows(resid, windows):
    return any(a <= resid <= b for a, b in windows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--loops", nargs="+", required=True)
    ap.add_argument("--n-conformers", type=int, default=20)
    ap.add_argument("--restraint-k", type=float, default=10.0)
    ap.add_argument("--temp", type=float, default=310.0)
    ap.add_argument("--production-ns", type=float, default=50.0,
                    help="production length in ns (conformers harvested evenly across it)")
    ap.add_argument("--equil-ns", type=float, default=0.5, help="NPT equilibration length ns")
    ap.add_argument("--pad-core", type=int, default=1)
    ap.add_argument("--padding", type=float, default=1.0, help="water box padding, nm")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--platform", default="auto")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import openmm as mm
    import openmm.app as app
    import openmm.unit as u
    from pdbfixer import PDBFixer

    os.makedirs(args.out, exist_ok=True)
    windows = parse_windows(args.loops)
    free_windows = [(a + args.pad_core, b - args.pad_core) for a, b in windows]
    dt = 2.0 * u.femtosecond
    steps_per_ns = int(round((1.0 * u.nanosecond) / dt))
    equil_steps = int(args.equil_ns * steps_per_ns)
    prod_steps = int(args.production_ns * steps_per_ns)
    interval = max(1, prod_steps // args.n_conformers)
    print(f"[explicit] target={args.pdb} free loops={free_windows}", flush=True)
    print(f"[explicit] {args.production_ns} ns production = {prod_steps} steps, "
          f"harvest every {interval} steps -> {args.n_conformers} conformers", flush=True)

    # --- prepare: add missing atoms + hydrogens ---
    fixer = PDBFixer(filename=args.pdb)
    fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    # --- explicit-solvent system: TIP3P water box + ions ---
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(ff, padding=args.padding * u.nanometer,
                        model="tip3p", ionicStrength=0.15 * u.molar,
                        neutralize=True)
    n_atoms = modeller.topology.getNumAtoms()
    print(f"[explicit] solvated system: {n_atoms} atoms", flush=True)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=1.0 * u.nanometer, constraints=app.HBonds)
    # NPT barostat
    system.addForce(mm.MonteCarloBarostat(1.0 * u.bar, args.temp * u.kelvin, 25))

    # --- position restraints on protein core heavy atoms (periodic form) ---
    restraint = mm.CustomExternalForce(
        "0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    restraint.addGlobalParameter("k", args.restraint_k * u.kilocalories_per_mole / u.angstrom**2)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    pos = modeller.positions
    n_restrained = 0
    protein_res = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
                   "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
                   "HID","HIE","HIP","CYX"}
    for atom in modeller.topology.atoms():
        if atom.residue.name not in protein_res:  # skip water/ions
            continue
        if atom.element is not None and atom.element.symbol == "H":
            continue
        try:
            resid = int(atom.residue.id)
        except ValueError:
            continue
        if not in_windows(resid, free_windows):  # core -> restrain
            x, y, z = pos[atom.index].value_in_unit(u.nanometer)
            restraint.addParticle(atom.index, [x, y, z])
            n_restrained += 1
    system.addForce(restraint)
    print(f"[explicit] restrained {n_restrained} protein-core heavy atoms", flush=True)

    integrator = mm.LangevinMiddleIntegrator(args.temp * u.kelvin, 1.0 / u.picosecond, dt)
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
    print(f"[explicit] OpenMM platform: {plat.getName()}", flush=True)
    sim = app.Simulation(modeller.topology, system, integrator, plat)
    sim.context.setPositions(modeller.positions)
    print("[explicit] minimising...", flush=True)
    sim.minimizeEnergy(maxIterations=5000)
    sim.context.setVelocitiesToTemperature(args.temp * u.kelvin, args.seed)
    print(f"[explicit] NPT equilibration {args.equil_ns} ns...", flush=True)
    sim.step(equil_steps)

    # protein atom indices for protein-only conformer output
    protein_idx = [a.index for a in modeller.topology.atoms() if a.residue.name in protein_res]
    prot_top = app.Modeller(modeller.topology, modeller.positions)
    to_delete = [a for a in prot_top.topology.atoms() if a.residue.name not in protein_res]
    prot_top.delete(to_delete)

    conformers = []
    for i in range(args.n_conformers):
        sim.step(interval)
        state = sim.context.getState(getPositions=True)
        allpos = state.getPositions(asNumpy=True)
        # fancy-index the nm Quantity array (keeps units so writeFile converts nm->Angstrom)
        prot_pos = allpos[protein_idx]
        out_pdb = os.path.join(args.out, f"conformer_{i:02d}.pdb")
        with open(out_pdb, "w") as fh:
            app.PDBFile.writeFile(prot_top.topology, prot_pos, fh)
        conformers.append(out_pdb)
        print(f"[explicit] conformer {i:02d} -> {out_pdb}", flush=True)

    print(f"[explicit] DONE: {len(conformers)} conformers in {args.out}", flush=True)


if __name__ == "__main__":
    main()
