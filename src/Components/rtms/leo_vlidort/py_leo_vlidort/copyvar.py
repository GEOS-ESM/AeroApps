#----
def _copyVar(ncIn,ncOut,name,dtype='f4',zlib=False,verbose=False,rename=None):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.variables[name]
    if verbose:
        print('copy variable ',name,x.dimensions)
    if rename is None:
        outname = name
    else:
        outname = rename
    y = ncOut.createVariable(outname,dtype,x.dimensions,zlib=zlib)
    if hasattr(x,'long_name'): y.long_name = x.long_name
    if hasattr(x,'units'): y.units = x.units 
    try:
        y.missing_value = x.missing_value
    except:
        pass
    data = x[:]
    y[:] = data
