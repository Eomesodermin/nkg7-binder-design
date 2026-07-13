# AlphaFlow predict.py patches (torch 1.12.1 compatibility)

AlphaFlow (github.com/bjing2016/alphaflow) pins torch 1.12.1+cu113 for the
OpenFold CUDA-kernel compile, but its `predict.py` assumes torch >= 2.x in two
places. Apply these after cloning, before inference:

1. **line ~76** — remove the `weights_only` kwarg (torch >= 1.13 only; 1.12.1's
   default already matches `weights_only=False`):
   ```
   -  ckpt = torch.load(args.weights, map_location='cpu', weights_only=False)
   +  ckpt = torch.load(args.weights, map_location='cpu')
   ```
2. **line ~77** — the checkpoint's `hyper_parameters` dict already contains a
   `training` key, which collides with the explicit `training=False` kwarg. Pop
   it first:
   ```
      ckpt['hyper_parameters'].pop('training', None)
      model = model_class(**ckpt['hyper_parameters'], training=False)
   ```

Two further non-code fixes (missing runtime deps the repo imports but does not
pin): `pip install requests tqdm wandb` into the alphaflow env.

Full build recipe + dependency pins: see the `protein-loop-flexibility` skill,
`references/dl_ensemble.md`.
