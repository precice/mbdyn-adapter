MBDyn-preCICE adapter
----------------------------------------------------

## Installing the package
After cloning this repository, execute the following command:
```
pip3 install --user .
```
Alternativly the package can also be used by including the source code folder in the python path environment:
```
export PYTHONPATH=$PYTHONPATH:<path-to-adapter-repository>/src/
```

### Dependencies:
*    **[MBDyn](https://public.gitlab.polimi.it/DAER/mbdyn)** with Python interface (*develop*-branch for Python 3 support)
*    **[preCICE](https://github.com/precice/precice)** and **[preCICE/python-bindings](https://github.com/precice/python-bindings)** (Version 2)

### Installing MBDyn
Clone the MBDyn repository and switch to *develop*-branch, configure and build with the Python interface enabled:
```
git clone https://public.gitlab.polimi.it/DAER/mbdyn.git
cd mbdyn && git checkout develop
./bootstrap.sh
CPPFLAGS=-I/usr/include/suitesparse PYTHON_VERSION=3 ./configure --enable-python=yes
make
sudo make install
```

Dependencies:
*    UMFPACK (Part of libsuitesparse-dev in Ubuntu repository.)

### Installing preCICE
Install preCICE via an appropriate debian package from the [release page](https://github.com/precice/precice/releases) and [python-bindings from PyPI](https://pypi.org/project/pyprecice/) using pip3:
```
sudo apt install ./libprecice2_<version>_<codename>.deb
pip3 install --user pyprecice
```
## Usage

Add mbdyn executable path to $PATH
```
export PATH=$PATH:/usr/local/mbdyn/bin
```

Add the Python interface path to $PYTHONPATH
```
export PYTHONPATH=$PYTHONPATH:/usr/local/mbdyn/libexec/mbpy
```

Add adapter to $PYTHONPATH
```
export PYTHONPATH=$PYTHONPATH:<path-to-adapter-root>/mbdynAdapter
```

Requires, see example folder:
*   pyhton script for calling adapter
*   config file and mesh file currently only support gmsh-files version 2 with ascii format

###Warning:
Meshes with large number of nodes or  can lead to desync of socket stream and mbdyn api failing. Workaround implemented by setting sleep timer before reading in mbdyn script, set in [input.py  `__force_coupling_str`](https://github.com/Hag3nL/mbdyn-adapter/blob/master/mbdynAdapter/input.py#L290).

## Acknowledgements
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 642682 (AWESCO).
