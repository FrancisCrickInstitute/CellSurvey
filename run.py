import ctypes.util
import pathlib as _pl
import os as _os

# Preload the pixi environment's libstdc++ to avoid ABI conflicts with the
# system library. Locate it relative to this script's pixi env directory.
_pixi_env_lib = _pl.Path(__file__).resolve().parent / ".pixi" / "envs" / "default" / "lib"
_libstdcpp = _pixi_env_lib / "libstdc++.so.6"
if _libstdcpp.exists():
    ctypes.CDLL(str(_libstdcpp))

# Force POT (Optimal Transport library used by MuSpAn) to skip TensorFlow
# backend import. TF 2.10 is compiled against numpy 1.x and crashes on
# numpy 2.x import, but we only need TF for Stardist segmentation later.
_os.environ["POT_BACKEND"] = "numpy"

from cellsurvey.cli import main

if __name__ == '__main__':
    main()
