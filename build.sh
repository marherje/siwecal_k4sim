# 1. Source key4hep (release pinned in .key4hep-release; override with KEY4HEP_RELEASE)
_K4_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )"
KEY4HEP_RELEASE="${KEY4HEP_RELEASE:-$(cat "${_K4_ROOT}/.key4hep-release" 2>/dev/null || echo 2026-04-08)}"
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "${KEY4HEP_RELEASE}"

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

