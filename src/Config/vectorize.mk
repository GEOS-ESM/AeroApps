#
# System dependent vectorization flags
#


  ifeq ($(ESMA_FC), ifort)
  ifneq ("$(BOPT)", "g")

     # ----
     # FOPT
     # ----

     # For vectorizing on Sandybridge and higher processors with AVX and AVX2 instructions
     # NOTE: No guarantee of zero-diff between *bridge and *well/*lake processors
     #FOPT := $(FOPT3) -axAVX,CORE-AVX2 -fma -qopt-report0 -ftz -align all -fno-alias -align array32byte

     # For vectorizing on Haswell and higher processors with AVX2 instructions
     FOPT := $(FOPT3) -xCORE-AVX2 -fma -qopt-report0 -ftz -align all -fno-alias -align array32byte

     # For vectorizing on Skylake and higher processors with CORE-AVX512 instructions
     #FOPT := $(FOPT3) -xCORE-AVX512 -qopt-zmm-usage=high -fma -qopt-report0 -ftz -align all -fno-alias -align array64byte

     # Add common FOPT flags
     FOPT += -traceback -assume realloc_lhs

     # ---
     # FPE
     # ---

     # For lower precision, but (usually) better performance from AVX instructions, enable this
     # Allows for MPI layout regression in testing
     FPE := -fpe3 -fp-model consistent -g -assume noold_maxminloc

     # For lower precision, but better performance from AVX instructions, enable this
     #FPE := -fpe3 -fp-model fast=2 -no-prec-div -g -assume noold_maxminloc

  endif
  endif
