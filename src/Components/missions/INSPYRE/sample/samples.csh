#!/bin/csh
#SBATCH --time=1:00:00
#SBATCH --nodes=1 --ntasks-per-node=125
#SBATCH --job-name=inspyre_sample_traj
#SBATCH --account=s3339

  cd /home/pcolarco/INSPYRE_25/pcolarco/sample
  setenv AEROAPPS $NOBACKUP/AeroApps
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH

# run sampling scripts
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-CORE_GV_20260727_RA.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-CORE_GV_20260729_RA.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-CORE_GV_20260803_RA.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-CORE_GV_20260806_RC_L1.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-CORE_GV_20260806_RC_L2.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-CORE_GV_20260807_RC.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-CORE_GV_20260808_RC.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-CORE_GV_20260812_RC.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-MetNav_ER2_20260727_RA.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-MetNav_ER2_20260729_RA.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-MetNav_ER2_20260801_RA.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-MetNav_ER2_20260807_RA.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-MetNav_ER2_20260808_RA.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-MetNav_ER2_20260811_RA.ict
# ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-MetNav_ER2_20260812_RA.ict
 ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-CORE_GV_20260817_RC.ict
 ./sample_trajectory.py fp inst3_3d_aer_Nv ../../data/INSPYRE-MetNav_ER2_20260817_RA.ict

# plot curtain extinction
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-CORE_GV_20260727_RA.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-CORE_GV_20260729_RA.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-CORE_GV_20260803_RA.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-CORE_GV_20260806_RC_L1.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-CORE_GV_20260806_RC_L2.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-CORE_GV_20260807_RC.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-CORE_GV_20260808_RC.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-CORE_GV_20260812_RC.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-MetNav_ER2_20260727_RA.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-MetNav_ER2_20260729_RA.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-MetNav_ER2_20260801_RA.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-MetNav_ER2_20260807_RA.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-MetNav_ER2_20260808_RA.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-MetNav_ER2_20260811_RA.ict fp
# ./plot_trj_ext_curtain.py ../../data/INSPYRE-MetNav_ER2_20260812_RA.ict fp
 ./plot_trj_ext_curtain.py ../../data/INSPYRE-CORE_GV_20260817_RC.ict fp
 ./plot_trj_ext_curtain.py ../../data/INSPYRE-MetNav_ER2_20260817_RA.ict fp

# plot curtain mass
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-CORE_GV_20260727_RA.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-CORE_GV_20260729_RA.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-CORE_GV_20260803_RA.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-CORE_GV_20260806_RC_L1.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-CORE_GV_20260806_RC_L2.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-CORE_GV_20260807_RC.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-CORE_GV_20260808_RC.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-CORE_GV_20260812_RC.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-MetNav_ER2_20260727_RA.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-MetNav_ER2_20260729_RA.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-MetNav_ER2_20260801_RA.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-MetNav_ER2_20260807_RA.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-MetNav_ER2_20260808_RA.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-MetNav_ER2_20260811_RA.ict fp
# ./plot_trj_mass_curtain.py ../../data/INSPYRE-MetNav_ER2_20260812_RA.ict fp
 ./plot_trj_mass_curtain.py ../../data/INSPYRE-CORE_GV_20260817_RC.ict fp
 ./plot_trj_mass_curtain.py ../../data/INSPYRE-MetNav_ER2_20260817_RA.ict fp

# plot extinction along compared to CAPS
# ./plot_trj_ext_along.py ../../data/INSPYRE-CORE_GV_20260727_RA.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-CORE_GV_20260729_RA.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-CORE_GV_20260803_RA.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-CORE_GV_20260806_RC_L1.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-CORE_GV_20260806_RC_L2.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-CORE_GV_20260807_RC.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-CORE_GV_20260808_RC.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-CORE_GV_20260812_RC.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-MetNav_ER2_20260727_RA.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-MetNav_ER2_20260729_RA.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-MetNav_ER2_20260801_RA.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-MetNav_ER2_20260807_RA.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-MetNav_ER2_20260808_RA.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-MetNav_ER2_20260811_RA.ict fp
# ./plot_trj_ext_along.py ../../data/INSPYRE-MetNav_ER2_20260812_RA.ict fp
 ./plot_trj_ext_along.py ../../data/INSPYRE-CORE_GV_20260817_RC.ict fp
 ./plot_trj_ext_along.py ../../data/INSPYRE-MetNav_ER2_20260817_RA.ict fp

# plot mass along compared to SP2
# ./plot_trj_mass_along.py ../../data/INSPYRE-CORE_GV_20260727_RA.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-CORE_GV_20260729_RA.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-CORE_GV_20260803_RA.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-CORE_GV_20260806_RC_L1.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-CORE_GV_20260806_RC_L2.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-CORE_GV_20260807_RC.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-CORE_GV_20260808_RC.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-CORE_GV_20260812_RC.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-MetNav_ER2_20260727_RA.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-MetNav_ER2_20260729_RA.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-MetNav_ER2_20260801_RA.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-MetNav_ER2_20260807_RA.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-MetNav_ER2_20260808_RA.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-MetNav_ER2_20260811_RA.ict fp
# ./plot_trj_mass_along.py ../../data/INSPYRE-MetNav_ER2_20260812_RA.ict fp
 ./plot_trj_mass_along.py ../../data/INSPYRE-CORE_GV_20260817_RC.ict fp
 ./plot_trj_mass_along.py ../../data/INSPYRE-MetNav_ER2_20260817_RA.ict fp
