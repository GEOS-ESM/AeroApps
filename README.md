# Aeroapps

## How to build AeroApps

### Preliminary Steps

#### Load Build Modules

In your `.bashrc` or `.tcshrc` or other rc file add a line:

##### NCCS

```
module use -a /discover/swdev/gmao_SIteam/modulefiles-SLES15
```

##### NAS
```
module use -a /nobackup/gmao_SIteam/modulefiles
```

##### GMAO Desktops

On the GMAO desktops, the SI Team modulefiles should automatically be
part of running `module avail` but if not, they are in:

```
module use -a /ford1/share/gmao_SIteam/modulefiles
```

Also do this in any interactive window you have. This allows you to get module files needed to correctly checkout and build the model.

Now load the `GEOSenv` module:
```
module load GEOSenv
```
which obtains the latest `git`, `CMake`, etc. modules needed to build.

#### Obtain the Model

```
git clone git@github.com:GEOS-ESM/AeroApps.git
```

---
#### Use mepo to checkout the model
```
cd AeroApps
mepo clone
```

##### Slow clones

If you notice your clone is taking a while, we recommend running:

```
mepo config set clone.partial blobless
```

This is a one-time command that tells mepo to use blobless clones for all future clones. Blobless clones are much faster than the default clone method, especially for repositories with a large history like MAPL.

#### Build the Model

##### Load Compiler, MPI Stack, and Baselibs

On tcsh:
```
source env@/g5_modules
```
or on bash:
```
source env@/g5_modules.sh
```

##### Create Build Directory

##### Run CMake

CMake generates the Makefiles needed to build the model.

On bash you can run:

```
cmake -B build -S . --install-prefix=$(pwd)/install
```

and on tcsh:

```
cmake -B build -S . --install-prefix=`pwd`/install
```

This will build the code in `build/` and will install in `install`.

NOTE: You can choose any directory for your build and install. That
said, what you pass to `--install-prefix=` *MUST* be a full path.

###### Building with Debugging Flags

To build with debugging flags add:

```
-DCMAKE_BUILD_TYPE=Debug
```
to the cmake line.

##### Build and Install with Make

```
cmake --build build --target install -j 6
```

If you put your build in a directory other than `build` use that in the
above command.

## Contributing

Please check out our [contributing guidelines](CONTRIBUTING.md).

## License

All files are currently licensed under the Apache-2.0 license, see [`LICENSE`](LICENSE).

Previously, the code was licensed under the [NASA Open Source Agreement, Version 1.3](LICENSE-NOSA).
