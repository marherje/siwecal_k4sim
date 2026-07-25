# 1. Source key4hep (release pinned in .key4hep-release; override with KEY4HEP_RELEASE).
#    Delegated to init_key4hep.sh so the pin is resolved in exactly one place;
#    "|| exit 1" because a sourced file can only 'return', which would otherwise
#    let the build carry on against whatever stack happens to be loaded.
_K4_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )"
source "${_K4_ROOT}/init_key4hep.sh" || exit 1

# 2. Build plugin
mkdir -p build && cd build
cmake -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH" \
      -DCMAKE_INSTALL_PREFIX=../install \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make install -j16
cd ..

# 3. Expose plugins to DD4hep & Gaudi algorithms to python
export LD_LIBRARY_PATH=$PWD/install/lib64:$PWD/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/install/lib64:$PWD/install/lib:$PWD/install/python:$PYTHONPATH

