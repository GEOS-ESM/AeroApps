#!/usr/bin/env python
"""
add oci_wav to rsr files
"""

from pyhdf.SD import SD, SDC

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    inFile = 'OCI_RSR_2.5nm_v7.hdf'
    outFile = 'OCI_RSR_2.5nm_v7_amended.hdf'

    hdf = SD(inFile,SDC.READ)
    wav_rsr = hdf.select('rsrwave')[:]
    wav_cen = hdf.select('centerwavelength')[:]
    band = hdf.select('bandwidth')[:]
    rsr = hdf.select('RSR')[:]
    nrsrwave,nwave = rsr.shape
    hdf.end()

    wave = np.arange(310,889.5,2.5)
    swir = np.array([940.615, 1038.97, 1250.15, 1378.04, 1615.98, 2130.4, 2260.64])
    wave = np.append(wave,swir)


    hdf = SD(outFile,SDC.WRITE|SDC.CREATE)
    hdf.Title = "OCI relative spectral response (RSR)"

    var = hdf.create('wave',SDC.FLOAT32,(nwave))
    var.long_name = 'wavelengths'
    var.units = 'nm'
    dim1 = var.dim(0)
    dim1.setname('nwave')
    var[:] = wave.astype(np.float32)
    var.endaccess()

    var = hdf.create('rsrwave',SDC.FLOAT32,(nrsrwave))
    var.long_name = 'rsrwavelengths'
    var.units = 'nm'
    dim1 = var.dim(0)
    dim1.setname('nrsrwave')
    var[:] = wav_rsr
    var.endaccess()
   
    var = hdf.create('centerwavelength',SDC.FLOAT32,(nwave))
    var.long_name = 'centerwavelengths'
    var.units = 'nm'
    dim1 = var.dim(0)
    dim1.setname('nwave')
    var[:] = wav_cen
    var.endaccess()     

    var = hdf.create('bandwidth',SDC.FLOAT32,(nwave))
    var.long_name = 'bandwidths'
    var.units = 'nm'
    dim1 = var.dim(0)
    dim1.setname('nwave')
    var[:] = band
    var.endaccess()

    var = hdf.create('RSR',SDC.FLOAT64,(nrsrwave,nwave))
    var.long_name = 'RSR'
    var.units = 'dimensionless'
    dim1 = var.dim(0)
    dim1.setname('nrsrwave')
    dim2 = var.dim(1)
    dim2.setname('nwave')
    var[:] = rsr
    var.endaccess()

    hdf.end()
