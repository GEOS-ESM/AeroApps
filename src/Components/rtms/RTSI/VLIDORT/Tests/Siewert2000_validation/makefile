
# general notes
###############

# Name of the compiler and the LINK.f statement should be supplied.

# Ordering of modules follows that required for GNU - the linker
# needs the modules in a certain order (a called module must be
# already compiled and linked before).

#  Compilers
#  =========

#  --------- g77

#VLIDORT_COMPILE = g77 -c -Wimplicit -Wall -ffortran-bounds-check -O2
#LAPACK_COMPILE  = g77 -c -ffortran-bounds-check -O2
#LAPACK_NOPT_COMPILE = g77 -c -ffortran-bounds-check
#F77 = g77

#  --------- Ifort

#VLIDORT_COMPILE = ifort -c     -warn all -zero   -CB -traceback  -check all
#LAPACK_COMPILE  = ifort -c     -warn all -zero     -CB -traceback  -check all
#LAPACK_NOPT_COMPILE = ifort -c -warn all -zero -CB -traceback  -check all
#VLIDORT_COMPILE = ifort -c -O3 -warn all -check bounds -zero  
#LAPACK_COMPILE  = ifort -c -O3  -zero 
#LAPACK_NOPT_COMPILE = ifort -c -O3 -zero
#F77 = ifort

#  ---------- gfortran
#   WARNING - must use some level of optimization (-O1, -O2 or -O3)

#VLIDORT_COMPILE = gfortran -c -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
#LAPACK_COMPILE  = gfortran -c -g -C -Wall  -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
#LAPACK_NOPT_COMPILE = gfortran -c -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
VLIDORT_COMPILE = gfortran -c -O3 -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
LAPACK_COMPILE  = gfortran -c -O3 -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
LAPACK_NOPT_COMPILE = gfortran -c -O3 -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
F77 = gfortran 

#  ---------- pgf90

#VLIDORT_COMPILE = pgf90 -c #-g -C -Minform,inform
#LAPACK_COMPILE  = pgf90 -c
#LAPACK_NOPT_COMPILE = pgf90 -c 
#F77 = pgf90


# Link definition
#################

LINK.f = $(F77)

# Path variables (User must change all of these)
################

VLIDORT_PATH = ../..
SPATH_S = $(VLIDORT_PATH)/sourcecode_str/
IPATH_S = $(VLIDORT_PATH)/includes/
LPATH_S = $(VLIDORT_PATH)/LAPACK_directory/
OBJ     = $(VLIDORT_PATH)/OBJECTS

#  OBJECT MODULES
#################

#  Masters sets

OBJECTS_VLIDORT_MASTERS1   = $(OBJ)/vlidort_masters_lambertian.o
OBJECTS_VLIDORT_MASTERS2   = $(OBJ)/vlidort_masters_brdf.o

OBJECTS_VLIDORT_L_MASTERS1 = $(OBJ)/vlidort_l_masters_lambertian.o
OBJECTS_VLIDORT_L_MASTERS2 = $(OBJ)/vlidort_l_masters_brdf.o

#  Basic set for Stokes vectors

OBJECTS_VLIDORT_SOLUTIONS = $(OBJ)/vlidort_solutions_str.o   \
                            $(OBJ)/vlidort_multipliers_str.o \
                            $(OBJ)/vlidort_bvproblem_str.o

OBJECTS_VLIDORT_INTENSITY = $(OBJ)/vlidort_intensity_str.o   \
                            $(OBJ)/vlidort_qm_intensity.o    \
                            $(OBJ)/vlidort_corrections_str.o \
                            $(OBJ)/vlidort_converge_str.o

OBJECTS_VLIDORT_THERMAL   = $(OBJ)/vlidort_thermalsup.o

OBJECTS_VLIDORT_SETUPS    = $(OBJ)/vlidort_miscsetups_str.o \
                            $(OBJ)/vlidort_inputs.o 

OBJECTS_VLIDORT_BRDF      = $(OBJ)/vlidort_brdf_supplement.o \
                            $(OBJ)/vlidort_brdf_kernels.o

OBJECTS_VLIDORT_AUX       = $(OBJ)/vlidort_aux.o          \
                            $(OBJ)/vlidort_writemodules_str.o

# linearized set

OBJECTS_VLIDORT_L_SOLUTIONS = $(OBJ)/vlidort_l_solutions_str.o   \
                              $(OBJ)/vlidort_l_multipliers_str.o \
                              $(OBJ)/vlidort_l_bvproblem_str.o

OBJECTS_VLIDORT_L_JACOBIANS = $(OBJ)/vlidort_l_wfatmos_str.o     \
                              $(OBJ)/vlidort_l_wfsurface_str.o   \
                              $(OBJ)/vlidort_qm_jacobians.o      \
                              $(OBJ)/vlidort_l_corrections_str.o \
                              $(OBJ)/vlidort_l_converge_str.o

OBJECTS_VLIDORT_L_THERMAL   = $(OBJ)/vlidort_l_thermalsup.o

OBJECTS_VLIDORT_L_SETUPS    = $(OBJ)/vlidort_l_miscsetups_str.o \
                              $(OBJ)/vlidort_l_inputs.o 

OBJECTS_VLIDORT_L_BRDF      = $(OBJ)/vlidort_l_brdf_supplement.o \
                              $(OBJ)/vlidort_l_brdf_kernels.o

OBJECTS_VLIDORT_L_AUX       = $(OBJ)/vlidort_l_writemodules_str.o

# environment & interface modules external to the main package

OBJECTS_VLIDORT_V2P4RT1  = $(OBJ)/Siewert2000_validation.o

#  Lapack source code

OBJECTS_LAPACK            = $(OBJ)/LAPACK-BLAS-tools.o  \
                            $(OBJ)/LAPACK-drivers.o     \
                            $(OBJ)/LAPACK-ev.o          \
                            $(OBJ)/LAPACK-dlamch.o      \
                            $(OBJ)/LAPACK-tools.o       \
                            $(OBJ)/LAPACK_Z16-tools.o

# EXECUTABLES
#############

APPLICATION1      = Siewert2000_validation.exe 

all: $(APPLICATION1)

$(APPLICATION1):   $(OBJECTS_VLIDORT_V2P4RT1)    \
                  $(OBJECTS_VLIDORT_MASTERS1)    \
                  $(OBJECTS_VLIDORT_SOLUTIONS)   \
                  $(OBJECTS_VLIDORT_INTENSITY)   \
                  $(OBJECTS_VLIDORT_THERMAL)     \
                  $(OBJECTS_VLIDORT_SETUPS)      \
                  $(OBJECTS_VLIDORT_AUX)         \
                  $(OBJECTS_LAPACK)
	$(LINK.f) -o $(APPLICATION1)             \
                  $(OBJECTS_VLIDORT_V2P4RT1)     \
                  $(OBJECTS_VLIDORT_MASTERS1)    \
                  $(OBJECTS_VLIDORT_SOLUTIONS)   \
                  $(OBJECTS_VLIDORT_INTENSITY)   \
                  $(OBJECTS_VLIDORT_THERMAL)     \
                  $(OBJECTS_VLIDORT_SETUPS)      \
                  $(OBJECTS_VLIDORT_AUX)         \
                  $(OBJECTS_LAPACK)              \
                  $(FLIBRARIES)

# COMPILATIONS
##############

clean:
	rm -f *.exe ../../OBJECTS/*.o

# -------------------------------------------------
#--------------------------Environment modules-----
# -------------------------------------------------

$(OBJ)/Siewert2000_validation.o: Siewert2000_validation.f \
	 $(IPATH_S)VLIDORT.PARS            \
	 $(IPATH_S)VLIDORT_INPUTS.VARS     \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS   \
	 $(IPATH_S)VLIDORT_RESULTS.VARS
	$(VLIDORT_COMPILE) Siewert2000_validation.f -o $(OBJ)/Siewert2000_validation.o

# ----------------------------------------------------
#----------------------VLIDORT master modules --------
# ----------------------------------------------------

$(OBJ)/vlidort_masters_lambertian.o: $(SPATH_S)vlidort_masters_lambertian.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_THERMALSUP.VARS     \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_masters_lambertian.f -o $(OBJ)/vlidort_masters_lambertian.o

$(OBJ)/vlidort_masters_brdf.o: $(SPATH_S)vlidort_masters_brdf.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_THERMALSUP.VARS     \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_masters_brdf.f -o $(OBJ)/vlidort_masters_brdf.o

# ----------------------------------------------------
#----------------------VLIDORT L_master modules ------
# ----------------------------------------------------

$(OBJ)/vlidort_l_masters_lambertian.o: $(SPATH_S)vlidort_l_masters_lambertian.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_THERMALSUP.VARS     \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS 
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_masters_lambertian.f -o $(OBJ)/vlidort_l_masters_lambertian.o

$(OBJ)/vlidort_l_masters_brdf.o: $(SPATH_S)vlidort_l_masters_brdf.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_THERMALSUP.VARS     \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_masters_brdf.f -o $(OBJ)/vlidort_l_masters_brdf.o

# ----------------------------------------------------
#----------------------VLIDORT Radiance modules ------
# ----------------------------------------------------

#--------------------------- VLIDORT solutions modules

$(OBJ)/vlidort_solutions_str.o: $(SPATH_S)vlidort_solutions_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_solutions_str.f -o $(OBJ)/vlidort_solutions_str.o

$(OBJ)/vlidort_multipliers_str.o: $(SPATH_S)vlidort_multipliers_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_multipliers_str.f -o $(OBJ)/vlidort_multipliers_str.o

$(OBJ)/vlidort_bvproblem_str.o: $(SPATH_S)vlidort_bvproblem_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_REFLECTANCE.VARS    \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_bvproblem_str.f -o $(OBJ)/vlidort_bvproblem_str.o

#------------------------------ VLIDORT Intensity modules

$(OBJ)/vlidort_intensity_str.o: $(SPATH_S)vlidort_intensity_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_THERMALSUP.VARS     \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_REFLECTANCE.VARS    \
	 $(IPATH_S)VLIDORT_RESULTS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_intensity_str.f -o $(OBJ)/vlidort_intensity_str.o

$(OBJ)/vlidort_qm_intensity.o: $(SPATH_S)vlidort_qm_intensity.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_RESULTS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_qm_intensity.f -o $(OBJ)/vlidort_qm_intensity.o

$(OBJ)/vlidort_converge_str.o: $(SPATH_S)vlidort_converge_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SINGSCAT.VARS       \
	 $(IPATH_S)VLIDORT_RESULTS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_converge_str.f -o $(OBJ)/vlidort_converge_str.o

$(OBJ)/vlidort_corrections_str.o: $(SPATH_S)vlidort_corrections_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_SINGSCAT.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_corrections_str.f -o $(OBJ)/vlidort_corrections_str.o

#-------------------------- VLIDORT Setups modules

$(OBJ)/vlidort_miscsetups_str.o: $(SPATH_S)vlidort_miscsetups_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_REFLECTANCE.VARS    \
	 $(IPATH_S)VLIDORT_SETUPS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_miscsetups_str.f -o $(OBJ)/vlidort_miscsetups_str.o

$(OBJ)/vlidort_inputs.o: $(SPATH_S)vlidort_inputs.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_inputs.f -o $(OBJ)/vlidort_inputs.o
#---------------------------- VLIDORT Thermal supplement

$(OBJ)/vlidort_thermalsup.o: $(SPATH_S)vlidort_thermalsup.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_THERMALSUP.VARS     \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_thermalsup.f -o $(OBJ)/vlidort_thermalsup.o


#---------------------------- VLIDORT BRDF modules

$(OBJ)/vlidort_brdf_supplement.o: $(SPATH_S)vlidort_brdf_supplement.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_SINGSCAT.VARS       \
	 $(IPATH_S)VLIDORT_REFLECTANCE.VARS    \
	 $(IPATH_S)VLIDORT_BRDF.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_brdf_supplement.f -o $(OBJ)/vlidort_brdf_supplement.o

$(OBJ)/vlidort_brdf_kernels.o: $(SPATH_S)vlidort_brdf_kernels.f    \
	 $(IPATH_S)VLIDORT.PARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_brdf_kernels.f -o $(OBJ)/vlidort_brdf_kernels.o

#------------------------------- VLIDORT Auxiliary modules

$(OBJ)/vlidort_aux.o: $(SPATH_S)vlidort_aux.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS
	$(LAPACK_COMPILE) $(SPATH_S)vlidort_aux.f -o $(OBJ)/vlidort_aux.o

$(OBJ)/vlidort_writemodules_str.o: $(SPATH_S)vlidort_writemodules_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_RESULTS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_writemodules_str.f -o $(OBJ)/vlidort_writemodules_str.o


#------------------------------------------- LAPACK modules

$(OBJ)/LAPACK-BLAS-tools.o: $(LPATH_S)LAPACK-BLAS-tools.f
	$(LAPACK_COMPILE) $(LPATH_S)LAPACK-BLAS-tools.f -o $(OBJ)/LAPACK-BLAS-tools.o

$(OBJ)/LAPACK-drivers.o: $(LPATH_S)LAPACK-drivers.f
	$(LAPACK_COMPILE) $(LPATH_S)LAPACK-drivers.f -o $(OBJ)/LAPACK-drivers.o

$(OBJ)/LAPACK-dlamch.o: $(LPATH_S)LAPACK-dlamch.f
	$(LAPACK_NOPT_COMPILE) $(LPATH_S)LAPACK-dlamch.f -o $(OBJ)/LAPACK-dlamch.o

$(OBJ)/LAPACK-ev.o: $(LPATH_S)LAPACK-ev.f
	$(LAPACK_COMPILE) $(LPATH_S)LAPACK-ev.f -o $(OBJ)/LAPACK-ev.o

$(OBJ)/LAPACK-tools.o: $(LPATH_S)LAPACK-tools.f
	$(LAPACK_COMPILE) $(LPATH_S)LAPACK-tools.f -o $(OBJ)/LAPACK-tools.o

$(OBJ)/LAPACK_Z16-tools.o: $(LPATH_S)LAPACK_Z16-tools.f
	$(LAPACK_COMPILE) $(LPATH_S)LAPACK_Z16-tools.f -o $(OBJ)/LAPACK_Z16-tools.o

# ---------------------------------------------------
#----------------------VLIDORT Jacobian modules ------
# ---------------------------------------------------

#---------------VLIDORT linearized solutions modules

$(OBJ)/vlidort_l_solutions_str.o: $(SPATH_S)vlidort_l_solutions_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_L_SETUPS.VARS       \
	 $(IPATH_S)VLIDORT_L_MULTIPLIERS.VARS  \
	 $(IPATH_S)VLIDORT_L_SOLUTION.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_solutions_str.f -o $(OBJ)/vlidort_l_solutions_str.o

$(OBJ)/vlidort_l_multipliers_str.o: $(SPATH_S)vlidort_l_multipliers_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_L_SETUPS.VARS       \
	 $(IPATH_S)VLIDORT_L_SOLUTION.VARS     \
	 $(IPATH_S)VLIDORT_L_MULTIPLIERS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_multipliers_str.f -o $(OBJ)/vlidort_l_multipliers_str.o

$(OBJ)/vlidort_l_bvproblem_str.o: $(SPATH_S)vlidort_l_bvproblem_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_REFLECTANCE.VARS    \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_L_SETUPS.VARS       \
	 $(IPATH_S)VLIDORT_L_THERMALSUP.VARS   \
	 $(IPATH_S)VLIDORT_L_SOLUTION.VARS     \
	 $(IPATH_S)VLIDORT_L_MULTIPLIERS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_bvproblem_str.f -o $(OBJ)/vlidort_l_bvproblem_str.o

#---------------------------- VLIDORT Linearized Thermal supplement

$(OBJ)/vlidort_l_thermalsup.o: $(SPATH_S)vlidort_l_thermalsup.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_THERMALSUP.VARS     \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_L_SETUPS.VARS       \
	 $(IPATH_S)VLIDORT_L_SOLUTION.VARS     \
	 $(IPATH_S)VLIDORT_L_THERMALSUP.VARS   \
	 $(IPATH_S)VLIDORT_L_MULTIPLIERS.VARS 
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_thermalsup.f -o $(OBJ)/vlidort_l_thermalsup.o

#------------------------- VLIDORT weighting function modules

$(OBJ)/vlidort_l_wfatmos_str.o: $(SPATH_S)vlidort_l_wfatmos_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_THERMALSUP.VARS     \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_REFLECTANCE.VARS    \
	 $(IPATH_S)VLIDORT_RESULTS.VARS        \
	 $(IPATH_S)VLIDORT_L_SETUPS.VARS       \
	 $(IPATH_S)VLIDORT_L_SOLUTION.VARS     \
	 $(IPATH_S)VLIDORT_L_THERMALSUP.VARS   \
	 $(IPATH_S)VLIDORT_L_MULTIPLIERS.VARS  \
	 $(IPATH_S)VLIDORT_L_RESULTS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_wfatmos_str.f -o $(OBJ)/vlidort_l_wfatmos_str.o

$(OBJ)/vlidort_qm_jacobians.o: $(SPATH_S)vlidort_qm_jacobians.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_L_SETUPS.VARS       \
	 $(IPATH_S)VLIDORT_L_SOLUTION.VARS     \
	 $(IPATH_S)VLIDORT_L_MULTIPLIERS.VARS  \
	 $(IPATH_S)VLIDORT_L_RESULTS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_qm_jacobians.f -o $(OBJ)/vlidort_qm_jacobians.o

$(OBJ)/vlidort_l_wfsurface_str.o: $(SPATH_S)vlidort_l_wfsurface_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_REFLECTANCE.VARS    \
	 $(IPATH_S)VLIDORT_L_SETUPS.VARS       \
	 $(IPATH_S)VLIDORT_L_SOLUTION.VARS     \
	 $(IPATH_S)VLIDORT_L_MULTIPLIERS.VARS  \
	 $(IPATH_S)VLIDORT_L_RESULTS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_wfsurface_str.f -o $(OBJ)/vlidort_l_wfsurface_str.o

$(OBJ)/vlidort_l_converge_str.o: $(SPATH_S)vlidort_l_converge_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_L_SINGSCAT.VARS     \
	 $(IPATH_S)VLIDORT_L_RESULTS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_converge_str.f -o $(OBJ)/vlidort_l_converge_str.o

$(OBJ)/vlidort_l_corrections_str.o: $(SPATH_S)vlidort_l_corrections_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_MULTIPLIERS.VARS    \
	 $(IPATH_S)VLIDORT_SINGSCAT.VARS       \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_L_SINGSCAT.VARS     \
	 $(IPATH_S)VLIDORT_L_SETUPS.VARS       \
	 $(IPATH_S)VLIDORT_L_SOLUTION.VARS     \
	 $(IPATH_S)VLIDORT_L_MULTIPLIERS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_corrections_str.f -o $(OBJ)/vlidort_l_corrections_str.o

#--------------------- VLIDORT Linearized Setups modules

$(OBJ)/vlidort_l_miscsetups_str.o: $(SPATH_S)vlidort_l_miscsetups_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_L_SETUPS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_miscsetups_str.f -o $(OBJ)/vlidort_l_miscsetups_str.o

$(OBJ)/vlidort_l_inputs.o: $(SPATH_S)vlidort_l_inputs.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_inputs.f -o $(OBJ)/vlidort_l_inputs.o

#---------------------- VLIDORT Linearized BRDF modules

$(OBJ)/vlidort_l_brdf_supplement.o: $(SPATH_S)vlidort_l_brdf_supplement.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_SETUPS.VARS         \
	 $(IPATH_S)VLIDORT_REFLECTANCE.VARS    \
	 $(IPATH_S)VLIDORT_BRDF.VARS           \
	 $(IPATH_S)VLIDORT_SOLUTION.VARS       \
	 $(IPATH_S)VLIDORT_SINGSCAT.VARS       \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_L_SINGSCAT.VARS     \
	 $(IPATH_S)VLIDORT_L_SETUPS.VARS       \
	 $(IPATH_S)VLIDORT_L_SOLUTION.VARS     \
	 $(IPATH_S)VLIDORT_L_BRDF.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_brdf_supplement.f -o $(OBJ)/vlidort_l_brdf_supplement.o

$(OBJ)/vlidort_l_brdf_kernels.o: $(SPATH_S)vlidort_l_brdf_kernels.f    \
	 $(IPATH_S)VLIDORT.PARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_brdf_kernels.f -o $(OBJ)/vlidort_l_brdf_kernels.o

#----------------------- VLIDORT Linearized Auxiliary modules

$(OBJ)/vlidort_l_writemodules_str.o: $(SPATH_S)vlidort_l_writemodules_str.f    \
	 $(IPATH_S)VLIDORT.PARS                \
	 $(IPATH_S)VLIDORT_INPUTS.VARS         \
	 $(IPATH_S)VLIDORT_BOOKKEEP.VARS       \
	 $(IPATH_S)VLIDORT_RESULTS.VARS        \
	 $(IPATH_S)VLIDORT_L_INPUTS.VARS       \
	 $(IPATH_S)VLIDORT_L_RESULTS.VARS
	$(VLIDORT_COMPILE) $(SPATH_S)vlidort_l_writemodules_str.f -o $(OBJ)/vlidort_l_writemodules_str.o

