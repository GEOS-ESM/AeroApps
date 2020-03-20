NOTES on linearized 2-stream model
==================================

	27 March 2009.

	   R. Spurr

	   RT Solutions Inc.
           9 Channing Street, Cambridge, MA 02138

General Notes
-------------

2stream code works equally with either stsream value (0.5, 0.57335).

2stream source code
-------------------

There are no include files and no dependencies. 

Note that this version fo the code has no "exact single scatter" 
corrections. This SS code would be the same as that in LIDORT anyway.

All routines are contained in the directory "2stream_sourcecode". 
The complete list of Fortran modules is as follows

    2stream_aux_lapack.f
    2stream_bvproblem.f
    2stream_inputs.f
    2stream_intensity.f
    2stream_l_bvproblem.f
    2stream_l_jacobians.f
    2stream_l_masters_basic.f
    2stream_l_miscsetups.f
    2stream_l_solutions.f
    2stream_masters_basic.f
    2stream_miscsetups.f
    2stream_solutions.f

The subdirectory "SS_declarations" contains some older routines which
still have the single scatter options (DO_SSFULL, DO_SSCORR_OUTGOING etc.)
in the argument lists. No code was ever written to use these options.

The two masters are

    2stream_masters_basic.f		intensity only
    2stream_l_masters_basic.f	intensity + Jacobians

The main routines for intensity calculations are

    2stream_bvproblem.f
    2stream_intensity.f
    2stream_miscsetups.f
    2stream_solutions.f

The main routines for Jacobian calculations are

    2stream_l_bvproblem.f
    2stream_l_jacobians.f
    2stream_l_miscsetups.f
    2stream_l_solutions.f

The following routine contains some input checks and an error routine

    2stream_inputs.f

The following routine contains auxiliary software from LAPACK, and should
be used when you are running 2stream alone (without LIDORT). If you are
using LIDORT as well, then you can leave out this routine from the
compilation, since LIDORT has its own LAPACK auxiliary code
(which is of course the same).

    2stream_aux_lapack.f

Testing directory
-----------------

   fort.87				  prepared 60-level atmosphere data

Makefile with ifort, generates 3 exectubles (no path change needed)

   makefile				  

These files are the 2stream Jacobian tests

   column_surface_wfs_SV1.out_save	  saved output from "test_2sj"
   column_surface_wfs_SV2.out_save        
   profile_surface_wfs_SV1.out_save
   profile_surface_wfs_SV2.out_save

SV1 = stream value 1 = 0.5
SV2 = stream value 2 = 0.57735 = SQRT(1/3)

The 3 calling programs are

   test_2sj.f			3-layer   2stream Jacobian testing
   test_2sjl.f			3-layer	  LIDORT calculation
   test_lidort_2stream.f	60-layer  LIDORT + 2stream comparison

The LIDORT input files are

   test_2sjl.inp			  Input file for test_2sjl.f
   test_lidort_2stream.inp                Input file for test_lidort_2stream.f

These routines are not used

   read_prop.f
   test_lidort_2stream.f_slow

The testing using "test_2sj.exe" will generate one of four files with
suffix ".out" - compare this with the same file with suffix ".out_save".
You can change from stream 1 to stream 2 by going inside the calling
program and doing it by hand. Simililary you can change from profile
to column weighting functions by hand - look in the calling program
for guidance (fairly easy). The linearized 2stream code will run with
PROFILE or COLUMN but not both Jacobians. The surface Jacobian is
separate from the atmospheric Jacobians. [LIDORT is not being run here].
This usage is similar to LIDORT and VLIDORT.

The testing using "test_2sjl.exe" is really just a validation using LIDORT.
This uses the same artifical 3-layer set up that is in "test_2sj". This is
a single call to LIDORT [2stream is not being run here].

The 3rd exectuable is "test_lidort_2stream.exe", which is a joint validation
of LIDORT and 2stream without the LInearization. This uses a pre-prepared
60-level atmosphere, and just compares Intensity at TOA and BOA. This is not
changed from the earlier package I sent in early February.

Optical Property setups
-----------------------

The optical property setups for the linearized 2-stream model are similar
to those for LIDORT, the main difference being that it is not necessary to
use normalized definitions of linearized IOPs. In the 3-layer case, we have--

TAU   = COL * ABS + SCA
OMEGA = SCA / TAU

(Here, COL is a scaling that is set to 1.0)

For the profile Jacobians, I have considered two weighting functions, one
with respect to ABS, the other w.r.t SCA. Thus---

    dTAU             dOMG      OMG
    ---- = COL       ---- = -  --- . COL
    dABS             dABS      TAU

    dTAU             dOMG      (1.0-SSA)
    ---- = 1.0       ---- = +  ---------
    dSCA             dSCA         TAU

If you use these inputs, you get the derivatives immediately.

I have used the following forms in the 2stream calculation, but this is
only for convenience of finite difference calculations - it is not
necessary to use them. If you do use these, then the final Jacobian
results must be divided by ABS and SCA to get a pure derivative.

       dTAU            dOMG
   ABS.----        ABS.--- 
       dABS            dABS

       dTAU            dOMG
   SCA.----        SCA.----
       dSCA            dSCA

The equivalent LIDORT definitions are the fully-normalized forms:

   ABS dTAU        ABS dOMG
   ---.----        ---.--- 
   TAU dABS        OMG dABS

   SCA dTAU        SCA dOMG
   ---.----        ---.--- 
   TAU dSCA        OMG dSCA

Similar considerations apply to the total column weighting functions.
This time, we differentiate w.r.t the variable COL to get

    dTAU             dOMG      OMG
    ---- = ABS       ---- = -  --- . ABS
    dCOL             dCOL      TAU

