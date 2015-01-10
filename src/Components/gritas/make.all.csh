#!/bin/csh -f
#
#  Builds gritas/grisas at 1x1, 2x2.5 and 4x5 resolution
#

make distclean
./configure

foreach res ( 1x1 2x2.5 4x5 )

    rm -f grids.h
    ln -s grids_$res.h grids.h

    make clean
    make all
    if ( $status ) exit 1

    mv gritas.x gritas_$res.exe
    mv grisas.x grisas_$res.exe

end

foreach exe ( *.exe )
    set bin = {$exe:r}.x
    mv $exe $bin
end

    echo " "
    echo "Created: " `ls -1 *.x`
    echo " "

    exit 0
