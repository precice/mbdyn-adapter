# preCICE adapter for MBDyn #

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
