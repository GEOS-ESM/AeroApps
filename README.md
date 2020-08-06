# AeroApps
Aerosol applications (needs improved description)


## To clone

1.  Load GMAO tool environment
```
$ module load GEOSenv
```

2. Clone AeroApps and its dependencies  (note side branch for now)
```
$ mepo clone git@github.com:GEOS-ESM/AeroApps.git -b feature/use_cmake
```

3. Create build directory
```
$ cd AeroApps
$ mkdir build
$ cd build
```

4. Load compiler environment
```
$ ../@env/g5_modules
```

5. Build and install.
Note that the default install path in is in the ```./build/install```
```
$ cmake .. -DBASEDIR=$BASEDIR/Linux  [-DCMAKE_INSTALL_PREFIX=<install-path>]
$ make -j 4 install
```
