# MBDyn-preCICE adapter

## Installing the package

### Using pip3

After cloning this repository, execute the following command:

```
pip3 install --user .
```

### Not installing

Alternativly the package can also be used by including the source code folder in the python path environment variable, the advantage here is that any changes to the source code can immediatly be utilized:

```
export PYTHONPATH=$PYTHONPATH:<path-to-adapter-repository>/src/
```

### Dependencies

* **[MBDyn](https://public.gitlab.polimi.it/DAER/mbdyn)** with Python interface (use *develop*-branch)
* **[preCICE](https://github.com/precice/precice)** and **[preCICE/python-bindings](https://github.com/precice/python-bindings)** (Version 2)

Note: Since MBDyn releases became unfrequent use repository instead, as of June 2021 the only branch the garuantees support for recent versions of Numpy and Python 3 is the *develop*-branch.

## Installing MBDyn

Clone the MBDyn repository and switch to *develop*-branch, configure and build with the Python interface enabled:

```
git clone https://public.gitlab.polimi.it/DAER/mbdyn.git
cd mbdyn && git checkout develop
./bootstrap.sh
CPPFLAGS=-I/usr/include/suitesparse PYTHON_VERSION=3 ./configure --enable-python=yes
make
```

Dependencies:
* UMFPACK (Part of libsuitesparse-dev in Ubuntu repository.)

Optionally install MBDyn with `sudo make install`. By default it will be installed in `/usr/local/mbdyn/`, the installation path can be set by configureing with the option `--prefix=/path`. See MBDyn [Software Installation Page](https://www.mbdyn.org/?Software_Installation) for more details.

## Installing preCICE

See [preCICE](https://precice.org/installation-overview.html) and [Python bindings](https://precice.org/installation-bindings-python.html) installation instructions.

## Usage

Before the adapter can be used the mbdyn binary needs to be added to `PATH` and the Python API to `PYTHONPATH`, e.g.:

```
export PATH=$PATH:<path-to-mbdyn>/bin
export PYTHONPATH=$PYTHONPATH:<path-to-mbdyn>/libexec/mbpy
```

The MBDyn input script is created automatically by defining the structure in a gmsh file and setting parameters in a configuration file, see example folder. A Python script is then used to start the adapter.

### Warning

Meshes with large number of nodes can lead to desync of the socket stream and mbdyn API failing. Workaround implemented by setting sleep timer before recieving in the mbdyn input script, set in input.py `__force_coupling_str`.

## Acknowledgements

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 642682 (AWESCO).
