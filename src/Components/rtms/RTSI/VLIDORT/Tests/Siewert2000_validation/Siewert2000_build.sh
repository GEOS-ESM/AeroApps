#  store original PARS file, replace with specific file

mv ../../includes/VLIDORT.PARS ../../includes/VLIDORT.PARS_hold
cp VLIDORT.PARS_Siewert2000 ../../includes/VLIDORT.PARS

#  make the executable

make clean
make

#  restore original PARS file

rm ../../includes/VLIDORT.PARS
mv ../../includes/VLIDORT.PARS_hold ../../includes/VLIDORT.PARS
