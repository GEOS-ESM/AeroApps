* file: ods_worksp.h -
*
*     Include file for ODS 
*
*................................................................

* REVISION HISTORY
*     28Sep1998  C. Redder   Redefined I_Min, I_Max, I4_Min
*                            and I4_Max to prevent overflows
*                            on some machines during compilation
*     29Sep1998  C. Redder   Correction for defining I4_Min
*     13Sep2002  C. Redder   Changed declaration of I1_Val from an array
*                            of one character strings to a long character
*                            string.
*     27Feb2009  D. Nadeau   Fixed variable type mismatch (char/int1)
*
*     notes: A hyperslab is a multidimensional block (or subset)
*            of data within a NetCDF file that are to be written
*            or read.
*
*            The variable, start, is an array of NetCDF file
*            indicies specifying the location of the corner of
*            the hyperslab where the first of the data values
*            are be written.
*
*            The variable, count, is an array of edge lengths
*            for the NetCDF file hyperslab.


*     Declaration for variables ( including work space ) accessed
*     by the modules calling the routines, ODS_DefWSp, ODS_ValWSp
*     and ODS_NCVWSp
*     -----------------------------------------------------------

*     Maximum amount of available work space ...
*     ------------------------------------------
      integer     Max_WorkSp              ! ... for two and four byte
      parameter ( Max_WorkSp  = 100000 )  !     integers and real numbers
*      parameter ( Max_WorkSp = 20 )
      integer     Max_CWorkSp             ! ... for one byte integers and
      parameter ( Max_CWorkSp = 30000  )  !     character strings (size may
                                          !     be machine dependent)

*     work (or scratch) space in native format for ...
*     ------------------------------------------------
      integer            I_Val  ( Max_WorkSp ) ! integers
      real               R_Val  ( Max_WorkSp ) ! real numbers
      character *               ( Max_CWorkSp )
     .                   C_Val                 ! character data

*     work (or scratch) space for ( machine dependent ) ...
*     -----------------------------------------------------
      character *               ( Max_CWorkSp )
     .                   IC_Val                ! one  byte character
      integer   *   1    I1_Val ( Max_WorkSp ) ! one  byte integers
      integer   *   2    I2_Val ( Max_WorkSp ) ! two  byte integers
      integer   *   4    I4_Val ( Max_WorkSp ) ! four byte integers
      real      *   4    R4_Val ( Max_WorkSp ) ! four byte real numbers

*     Index values used to define the block
*     of values or characters from 1-D array
*     --------------------------------------
      integer     ival_beg, ichar_beg
      integer     ival_end, ichar_end

*     Start and count values used to read/write the block of
*     values in the work space using the NetCDF routines
*     ------------------------------------------------------
      integer     start_block ( MAXNCDIM )
      integer     count_block ( MAXNCDIM )

*     common storage
*     --------------
      common / ods_worksp / start_block,
     .                      count_block,
     .                      R_Val,
     .                      C_Val,
     .                      I_Val,
     .                      R4_Val,
     .                      I4_Val,
     .                      I2_Val,
     .                      I1_Val


*     Declaration for variables accessed within the routines,
*     ODS_DefWSp, ODS_ValWSp and ODS_NCVWSp and used to
*     partitioned the data and allocate the work space

*     note: All blocks (or in other words subslabs) are
*           grouped and processed in cycles according to
*           the mathematical expression:
*
*               iCycle = ( iBlock - 1 ) / NBlocks_Cycle + 1
*
*           where iCycle and iBlock are the cycle and block
*           index numbers and NBlocks_Cycle is defined below.
*           All blocks in every cycle will have a size of
*           Max_BlockSz (defined below), except possibly the
*           last block to be processed in the cycle.
*           ------------------------------------------------
      integer     NDim_SubSlab  ! The number of NetCDF dimensions
                                !   for a sub-hyperslab which
                                !   contains an entire block of
                                !   partitioned data.
      integer     Max_SlabLen   ! maximum edge length of a subslab
                                !   for the dimension, NDim_SubSlab
      integer     Max_BlockSz   ! The maximum number of values
                                !   in a one-dimension block of 
                                !   data.
      integer     NBlocks_Cycle ! The number of data blocks per
                                !   cycle.  The first block in
                                !   a cycle corresponds to a
                                !   NetCDF dimension size of
                                !   1 for the dimension,
                                !   NDim_SubSlab.  The last block
                                !   in the cycle may have a size
                                !   of less than Max_BlockSz.
      integer     NVal_Cycle    ! The number of values corre-
                                !   sponding to NBlocks_Cycle

*     Initial values for start and count
*     ----------------------------------
      integer     start_HSlab  ( MAXNCDIM )
      integer     count_HSlab  ( MAXNCDIM )

*     common storage
*     --------------
      common / ods_wspdef / NDim_SubSlab,
     .                      Max_SlabLen,
     .                      Max_BlockSz,
     .                      NBlocks_Cycle,
     .                      NVal_Cycle,
     .                      start_HSlab,
     .                      count_HSlab

*     Parameter describing machine error ( machine dependent ) for ...
*     ----------------------------------------------------------------
      real             R_Error          ! single precision numbers
      parameter      ( R_Error  = 1.0e-6  )
      double precision D_Error          ! double precision numbers
      parameter      ( D_Error  = 1.0d-12 )

*     Maximum and minimum values for some
*     variable types in native format ( machine dependent )
*     -----------------------------------------------------
      integer     I_Min,   I_Max     ! integers
      double precision 
     .            I_MinD,  I_MaxD    ! in double precision format
      parameter ( I_Min    = -( 2 ** 30 - 1 + 2 ** 30 ) )
      parameter ( I_Max    =    2 ** 30 - 1 + 2 ** 30 )  ! Numbers
                                     !   split to prevent overflows
      parameter ( I_MinD   =  I_Min * ( 1.0d0 + D_Error ) )
      parameter ( I_MaxD   =  I_Max * ( 1.0d0 + D_Error ) )

      real        R_Min,   R_Max     ! integers
      parameter ( R_Min    =  -3.4e38 * ( 1.0  + R_Error ) )
      parameter ( R_Max    =   3.4e38 * ( 1.0  + R_Error ) )

*     Maximum and minimum values for each
*     recognized NetCDF variable type
*     -----------------------------------
      integer     I1_Min,  I1_Max    ! one-byte  integers
      real        I1_MinR, I1_MaxR   ! in floating point format 
      parameter ( I1_Min   =    0 )
      parameter ( I1_Max   =  255 )
      parameter ( I1_MinR  = - R_Error )
      parameter ( I1_MaxR  =  I1_Max * ( 1.0   + R_Error ) )

      integer     I2_Min,  I2_Max    ! two-byte  integers
      real        I2_MinR, I2_MaxR   ! in floating point format
      parameter ( I2_Min   = -( 2 ** 15 - 1 ) )
      parameter ( I2_Max   =    2 ** 15 - 1   )
      parameter ( I2_MinR  =  I2_Min * ( 1.0   + R_Error ) )
      parameter ( I2_MaxR  =  I2_Max * ( 1.0   + R_Error ) )

      integer     I4_Min,  I4_Max    ! four-byte integers
      double precision 
     .            I4_MinD, I4_MaxD   ! in double precision format
      parameter ( I4_Min   = -( 2 ** 30 - 1 + 2 ** 30 ) )
      parameter ( I4_Max   =    2 ** 30 - 1 + 2 ** 30 ) ! Numbers
                                     !   split to prevent overflows
      parameter ( I4_MinD  =  I4_Min * ( 1.0d0 + D_Error ) )
      parameter ( I4_MaxD  =  I4_Max * ( 1.0d0 + D_Error ) )

      real        R4_Min,  R4_Max    ! four-byte real numbers
      parameter ( R4_Min   =  -3.4e38 * ( 1.0  + R_Error ) )
      parameter ( R4_Max   =   3.4e38 * ( 1.0  + R_Error ) )


