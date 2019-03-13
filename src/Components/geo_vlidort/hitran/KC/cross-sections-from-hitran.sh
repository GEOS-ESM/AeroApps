pgf90 -O3 -o cross-sections-from-hitran.x cross-sections-from-hitran.f90
#
# test for IR calculations of solar absorption
#
cross-sections-from-hitran.x <<EOF
hitran08.part
cross-section-test.out
1, 4000., 0.02, 50001., 0.500, 300., 500, 0.50, 1
EOF
