import ctypes.util
import pathlib as _pl
import os as _os
import sys as _sys

# Preload the active environment's libstdc++ (e.g. .pixi/envs/<env>/lib) to
# avoid ABI conflicts with the system library. sys.prefix is the running env,
# so this works for both the default (Stardist) and instanseg pixi envs.
_libstdcpp = _pl.Path(_sys.prefix) / "lib" / "libstdc++.so.6"
if _libstdcpp.exists():
    ctypes.CDLL(str(_libstdcpp))

# Force TF to use legacy Keras 2 via the standalone tf_keras package.
# TF >=2.16 defaults to Keras 3 which compiles Stardist's model with
# XLA JIT, triggering a cuDNN autotuner failure on 1x1 convolutions
# with CUDA 12/cuDNN 9 on HPC. Legacy Keras 2 uses the non-XLA cuDNN
# path and retains GPU acceleration. Requires tf_keras pip package.
# Harmless in the instanseg env, which has no TensorFlow.
_os.environ["TF_USE_LEGACY_KERAS"] = "1"

from cellsurvey.cli import main

if __name__ == '__main__':
    main()
