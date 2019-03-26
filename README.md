# preCICE adapter for MBDyn #

## Installation ##
```
sudo python setup.py install
```

Requirements:
*    MBDyn with Python interface
*    preCICE with Python interface

Additionally, the example tutorial requires [FOAM-FSI](https://github.com/davidsblom/FOAM-FSI) 
## MBDyn installation ##
```
wget https://www.mbdyn.org/userfiles/downloads/mbdyn-1.7.3.tar.gz
tar -xf mbdyn-1.7.3.tar.gz
cd mbdyn-1.7.3
CPPFLAGS=-I/usr/include/suitesparse ./configure
make 
make install
```

Add the Python interface path to $PYTHONPATH
```
export PYTHONPATH=$PYTHONPATH:/usr/local/mbdyn/libexec/mbpy
```

Requirements:
*    UMFPACK (Part of libsuitesparse-dev in Ubuntu repository.)
