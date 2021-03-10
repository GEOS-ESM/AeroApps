!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_odsmeta --- Defines the ODS metadata conventions
!
! !INTERFACE:

      MODULE m_odsmeta

! !USES:

      Implicit none

! !DESCRIPTION:
! \label{MODS:Meta}
!  Defines ODS metadata conventions.
!
! !TO DO:
!
!  1. Define history and exclusion marks for superobbing
!  2. Implement bitwise encoding of history mark
!  3. Add kx/kt parameters
!
! !REVISION HISTORY:
!
!  01Feb99 (Dee)      See v1.1 revision history for previous revisions
!  13Dec99 (RT)       Added preprocessing flags defined by Meta and Yelena
!  20Dec99 (Dee/RT)   Removed restriction on preprocessing and online
!                      flag values; made it into a module.
!  23dec99  Todling/
!           da Silva  Created from qc_flag.h, added GEOS-4 specific flags.
!  28Dec99  Todling   Adding kt info (use to live in kt_max.h).
!  09Feb00  Dee       Added X_NOSTATS, X_PSAS
!  09Feb00  Todling   Defined kxmax here.
!  06Jun00  Todling   Updated kt list.
!  04Oct00  Todling   Increased kxmax to 512.
!  13May02  Dee       Modified some SQC flags and removed obsolete ones
!  03Jun02  Dee       Updated kt list; added X_BADLAYER
!  23Apr03  Sienkiewicz   Add history & exclusion marks for preprocessing
!  04Jun04  Todling   Added kxmod (per Jing Guo's m_Sndx)
!
!EOP
!-------------------------------------------------------------------------
    
      include 'odsmeta.h'

      end MODULE m_odsmeta


