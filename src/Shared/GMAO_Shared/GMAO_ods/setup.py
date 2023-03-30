#!/usr/bin/env python3
#
# This can be made more robust...
#
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    import os

    baselib = os.environ.get("BASELIB")
    esmalib = os.environ.get("ESMALIB")
    odslib  = os.environ.get("ODSLIB")

    if (baselib is None) or (esmalib is None) or (odslib is None):
        raise ValueError("BASELIB, ESMALIB or ODSLIB environment variables not set; run this from GNUmakefile")

    print("Using BASELIB = ", baselib)
    print("Using ESMALIB = ", esmalib)
    print("Using ODSLIB = ",  odslib)

    config = Configuration('pyods',package_path='pyods')
    config.add_extension('pyods_',['pyods_.F90',], \
                          library_dirs=[odslib,esmalib,baselib],
                          libraries=['GMAO_ods','GMAO_eu',
                                     'netcdf','hdf5_hl','hdf5','mfhdf','df','jpeg','curl','z','sz',
                                     'imf','m','dl','guide','stdc++','gcc_s.1','mpi',])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
