# MBDyn-preCICE adapter

<a style="text-decoration: none" href="https://github.com/precice/fenics-adapter/blob/master/LICENSE" target="_blank">
    <img src="https://img.shields.io/github/license/precice/mbdyn-adapter.svg" alt="GNU LGPL license">
</a>

:heart: **This adapter needs a maintainer. [Read more here](https://github.com/precice/mbdyn-adapter/issues/3).** :heart:

## Installation ##
```
sudo python setup.py install
```

Dependencies:
*    MBDyn with Python interface
*    preCICE with Python interface
*    OpenFOAM-adapter for preCICE

## MBDyn installation ##
```
wget https://www.mbdyn.org/userfiles/downloads/mbdyn-1.7.3.tar.gz
tar -xf mbdyn-1.7.3.tar.gz
cd mbdyn-1.7.3
CPPFLAGS=-I/usr/include/suitesparse ./configure --enable-python=yes
make 
make install
```

Add mbdyn executable path to $PATH
```
export PATH=$PATH:/usr/local/mbdyn/bin
```

Add the Python interface path to $PYTHONPATH
```
export PYTHONPATH=$PYTHONPATH:/usr/local/mbdyn/libexec/mbpy
```

Dependencies:
*    UMFPACK (Part of libsuitesparse-dev in Ubuntu repository.)

## Acknowledgements
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 642682 (AWESCO).
