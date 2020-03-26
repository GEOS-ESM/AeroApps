!	(qea.h)

!-----------------------------------------------------------------------
!///////////////////////////////////////////////////////////////////////
!-----------------------------------------------------------------------
! 18Jul95 - Jing G.	- Separated from gridxx.h, for Q.E.A. grid
!			  settings only.
!.......................................................................

!-----------------------------------------------------------------------
! Data block associated with the EA grid of "Ainc"
!-----------------------------------------------------------------------
!  Q.E.A parameters

      integer npole
      parameter ( npole = 4 )                 ! No. grid points at poles
                                              ! NOTE: In reality, at the
                                              !       poles all possible
                                              !       longitudes are 
                                              !       possible. Leave
                                              !       npole alone.
                                                
      real eaytresh                           ! treshold latitude for
                                              ! quasi-equal area grid



!     Indices of quasi-equal area  (q.e.a.) and regular lat/lon grid
!     --------------------------------------------------------------
      integer  lea_beg       ! pointer to begining q.e.a. gridpoint
                             !  for a given latitude    
      integer  lea_len       ! how many zonal gridpoints at that 
                             !  longitude

      integer  j_north   
      integer  j_south   ! from j_south to j_north lat/lon gridpoints
                         !  coincide with q.e.a. grid points


!     Longitude of q.e.q. gridpoints (lat not needed right now)
!     ---------------------------------------------------------
      real     ea_lon     ! longitude


      common / gridprm / eaytresh, j_south, j_north

      common / gridarr / lea_beg(jdimx),lea_len(jdimx),
     &	ea_lon(idimx*jdimx)
     
