import ctypes.util
import pathlib as _pl
import os as _os

# Preload the pixi environment's libstdc++ to avoid ABI conflicts with the
# system library. Locate it relative to this script's pixi env directory.
_pixi_env_lib = _pl.Path(__file__).resolve().parent / ".pixi" / "envs" / "default" / "lib"
_libstdcpp = _pixi_env_lib / "libstdc++.so.6"
if _libstdcpp.exists():
    ctypes.CDLL(str(_libstdcpp))

# Force TF to use legacy Keras 2 via the standalone tf_keras package.
# TF >=2.16 defaults to Keras 3 which compiles Stardist's model with
# XLA JIT, triggering a cuDNN autotuner failure on 1x1 convolutions
# with CUDA 12/cuDNN 9 on HPC. Legacy Keras 2 uses the non-XLA cuDNN
# path and retains GPU acceleration. Requires tf_keras pip package.
_os.environ["TF_USE_LEGACY_KERAS"] = "1"

from cellsurvey.cli import main

if __name__ == '__main__':
    main()
