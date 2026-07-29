import ctypes.util
import pathlib as _pl

# Preload the pixi environment's libstdc++ to avoid ABI conflicts with the
# system library. Locate it relative to this script's pixi env directory.
_pixi_env_lib = _pl.Path(__file__).resolve().parent / ".pixi" / "envs" / "default" / "lib"
_libstdcpp = _pixi_env_lib / "libstdc++.so.6"
if _libstdcpp.exists():
    ctypes.CDLL(str(_libstdcpp))

from cellsurvey.cli import main

if __name__ == '__main__':
    main()
