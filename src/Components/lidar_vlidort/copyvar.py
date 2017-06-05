#----
def _copyVar(ncIn,ncOut,name,dtype='f4',zlib=False,verbose=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.variables[name]
    if verbose:
        print 'copy variable ',name,x.dimensions
    y = ncOut.createVariable(name,dtype,x.dimensions,zlib=zlib)
    if hasattr(x,'long_name'): y.long_name = x.long_name
    if hasattr(x,'units'): y.units = x.units 
    try:
        y.missing_value = x.missing_value
    except:
        pass
    rank = len(x.shape)

    if rank == 1:
        y[:] = x[:]
    elif rank == 2:
        y[:,:] = x[:,:]
    elif rank == 3:
        y[:,:,:] = x[:,:,:]
    else:
        raise ValueError, "invalid rank of <%s>: %d"%(name,rank)