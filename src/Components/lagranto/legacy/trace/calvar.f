c     -------------------------------------------------------------------------
c     Potential temperature (TH)
c     -------------------------------------------------------------------------

      subroutine calc_TH (pt,t,p)
 
      implicit none
      
c     Argument declaration
      real  pt     ! Potential temperature [K]
      real  t      ! Temperature [either in C or in K]
      real  p      ! Pressure [hPa]

c     Physical parameters
      real      rdcp,tzero
      data      rdcp,tzero /0.286,273.15/
 
c     Calculation - distinction between temperature in C or in K
      if (t.lt.100.) then
         pt = (t+tzero) * ( (1000./p)**rdcp ) 
      else
         pt = t * ( (1000./p)**rdcp )
      endif

      end

c     -------------------------------------------------------------------------
c     Density (RHO)
c     -------------------------------------------------------------------------

      subroutine calc_RHO (rho,t,p)
 
      implicit none

c     Argument declaration
      real  rho    ! Density [kg/m^3]
      real  t      ! Temperature [either in C or in K]
      real  p      ! Pressure [hPa]

c     Physical parameters
      real      rd,tzero
      data      rd,tzero /287.05,273.15/

c     Auxiliary variables
      real      tk
 
c     Calculation - distinction between temperature in C or in K
      if (t.lt.100.) then
         tk = t + tzero
      else
         tk = t
      endif

      rho = 100.*p/( tk * rd ) 

      end

c     -------------------------------------------------------------------------
c     Relative humidity (RH)
c     -------------------------------------------------------------------------

      subroutine calc_RH (rh,t,p,q)
      
      implicit none

c     Argument declaration
      real  rh     ! Relative humidity [%]
      real  t      ! Temperature [either in C or in K]
      real  p      ! Pressure [hPa]
      real  q      ! Specific humidity [kg/kg]

c     Physical parameters
      real      rdcp,tzero
      data      rdcp,tzero /0.286,273.15/
      real      b1,b2w,b3,b4w,r,rd
      data      b1,b2w,b3,b4w,r,rd /6.1078, 17.2693882, 273.16, 35.86,
     &                              287.05, 461.51/

c     Auxiliary variables
      real      ge
      real      gqd
      real      tc
      real      pp,qk

c     Calculation - distinction between temperature in C or in K
      if (t.gt.100.) then
         tc = t - tzero
      else
         tc = t
      endif
      qk = q

      ge  = b1*exp(b2w*tc/(tc+b3-b4w))
      gqd = r/rd*ge/(p-(1.-r/rd)*ge)
      rh  = 100.*qk/gqd

      end

c     -------------------------------------------------------------------------
c     Equivalent potential temperature (THE)
c     -------------------------------------------------------------------------

      subroutine calc_THE (the,t,p,q)

      implicit none

c     Argument declaration
      real  the    ! Equivalent potential temperature [K]
      real  t      ! Temperature [either in C or in K]
      real  p      ! Pressure [hPa]
      real  q      ! Specific humidity [kg/kg]

c     Physical parameters 
      real     rdcp,tzero
      data     rdcp,tzero /0.286,273.15/

c     Auxiliary variables 
      real     tk,qk
 
c     Calculation - distinction between temperature in C or in K
      if (t.lt.100.) then
         tk = t + tzero
      else
         tk = t
      endif
      qk = q
      
      the = tk*(1000./p)
     +      **(0.2854*(1.0-0.28*qk))*exp(
     +      (3.376/(2840.0/(3.5*alog(tk)-alog(
     +      100.*p*max(1.0E-10,qk)/(0.622+0.378*
     +      q))-0.1998)+55.0)-0.00254)*1.0E3*
     +      max(1.0E-10,qk)*(1.0+0.81*qk))

      end

c     -------------------------------------------------------------------------
c     Latent heating rate (LHR)
c     -------------------------------------------------------------------------

      subroutine calc_LHR (lhr,t,p,q,omega,rh)

      implicit none 

c     Argument declaration
      real  lhr    ! Latent heating rate [K/6h]
      real  t      ! Temperature [either in C or in K]
      real  p      ! Pressure [hPa]
      real  q      ! Specific humidity [kg/kg]
      real  omega  ! Vertical velocity [Pa/s]
      real  rh     ! Relative humidity [%] 

c     Physical parameters 
      real      p0,kappa,tzero
      data      p0,kappa,tzero /1000.,0.286,273.15/
      real      blog10,cp,r,lw,eps
      data      blog10,cp,r,lw,eps /.08006,1004.,287.,2.5e+6,0.622/

c     Auxiliary variables
      real      tk
      real      qk
      real      tt
      real      esat,c

c     Calculation - distinction between temperature in C or in K
      if (t.lt.100.) then
         tk = t + tzero
      else
         tk = t
      endif
      qk = q

      if (rh.lt.80.) then 
         lhr = 0.    
      else if (omega.gt.0.) then 
         lhr = 0.
      else
         c   = lw/cp*eps*blog10*esat(tk)/p
         tt  = (tk*(p0/p)**kappa) 
         lhr = 21600.*     
     >        (1.-exp(.2*(80.-rh))) 
     >        *(-c*kappa*tt*omega/(100.*p))/(1.+c)   
      
      endif
      
      end

c     -------------------------------------------------------------------------
c     Wind speed (VEL)
c     -------------------------------------------------------------------------

      subroutine calc_VEL (vel,u,v)

      implicit none
 
c     Argument declaration
      real  vel    ! Wind speed [m/s]
      real  u      ! Zonal wind component [m/s]
      real  v      ! Meridional wind component [m/s]

      vel = sqrt ( u*u + v*v )

      end

c     -------------------------------------------------------------------------
c     Wind direction (DIR)
c     -------------------------------------------------------------------------

      subroutine calc_DIR (dir,u,v)

      implicit none
 
c     Argument declaration
      real  dir    ! Wind direction [deg]
      real  u      ! Zonal wind component [m/s]
      real  v      ! Meridional wind component [m/s]

      call getangle(1.,0.,u,v,dir)

      end

c     -------------------------------------------------------------------------
c     Zonal derivative of U (DUDX)
c     -------------------------------------------------------------------------

      subroutine calc_DUDX (dudx,u1,u0,lat)

      implicit none
  
c     Argument declaration
      real  dudx   ! Derivative of U in zonal direction [s^-1]
      real  u1     ! U @ LON + 1 DLON [m/s]
      real  u0     ! U @ LON - 1 DLON [m/s]
      real  lat    ! Latitude [deg]

c     Physical parameters
      real      pi180  
      parameter (pi180=31.14159/180.)
      real      deltay
      parameter (deltay =1.11e5)
      
      dudx = (u1-u0) / ( 2. * deltay * cos(pi180 * lat) )

      end

c     -------------------------------------------------------------------------
c     Zonal derivative of V (DVDX)
c     -------------------------------------------------------------------------

      subroutine calc_DVDX (dvdx,v1,v0,lat)
 
c     Argument declaration
      real  dvdx   ! Derivative of V in zonal direction [s^-1]
      real  v1     ! V @ LON + 1 DLON [m/s]
      real  v0     ! V @ LON - 1 DLON [m/s]
      real  lat    ! Latitude [deg]

c     Physical parameters
      real      pi180  
      parameter (pi180=3.14159/180.)
      real      deltay
      parameter (deltay =1.11e5)

      dvdx = (v1-v0) / ( 2. * deltay * cos(pi180 * lat) )

      end

c     -------------------------------------------------------------------------
c     Zonal derivative of T (DTDX)
c     -------------------------------------------------------------------------

      subroutine calc_DTDX (dtdx,t1,t0,lat)

      implicit none
  
c     Argument declaration
      real  dtdx   ! Derivative of T in zonal direction [K/m]
      real  t1     ! T @ LON + 1 DLON [K]
      real  t0     ! T @ LON - 1 DLON [K]
      real  lat    ! Latitude [deg]

c     Physical parameters
      real      pi180  
      parameter (pi180=3.14159/180.)
      real      deltay
      parameter (deltay =1.11e5)
      
      dtdx = (t1-t0) / ( 2. * deltay * cos(pi180 * lat) )

      end

c     -------------------------------------------------------------------------
c     Zonal derivative of TH (DTHDX)
c     -------------------------------------------------------------------------

      subroutine calc_DTHDX (dthdx,t1,t0,p,lat)

      implicit none
  
c     Argument declaration
      real  dthdx  ! Derivative of TH in zonal direction [K/m]
      real  t1     ! T @ LON + 1 DLON [K]
      real  t0     ! T @ LON - 1 DLON [K]
      real  p      ! P [hPa]
      real  lat    ! Latitude [deg]

c     Physical parameters
      real      pi180  
      parameter (pi180=3.14159/180.)
      real      deltay
      parameter (deltay =1.11e5)
      real      rdcp,pref
      data      rdcp,pref /0.286,1000./

      dthdx = (pref/p)**rdcp * 
     >        (t1-t0) / ( 2. * deltay * cos(pi180 * lat) )

      end

c     -------------------------------------------------------------------------
c     Meridional derivative of U (DUDY)
c     -------------------------------------------------------------------------

      subroutine calc_DUDY (dudy,u1,u0)

      implicit none
  
c     Argument declaration
      real  dudy   ! Derivative of U in meridional direction [s^-1]
      real  u1     ! U @ LAT + 1 DLAT [m/s]
      real  u0     ! U @ LAT - 1 DLAT [m/s]

c     Physical parameters
      real      deltay
      parameter (deltay =1.11e5)
      
      dudy = (u1-u0) / ( 2. * deltay )

      end

c     -------------------------------------------------------------------------
c     Meridional derivative of V (DVDY)
c     -------------------------------------------------------------------------

      subroutine calc_DVDY (dvdy,v1,v0)

      implicit none
  
c     Argument declaration
      real  dvdy   ! Derivative of V in meridional direction [s^-1]
      real  v1     ! V @ LAT + 1 DLAT [m/s]
      real  v0     ! V @ LAT - 1 DLAT [m/s]

c     Physical parameters
      real      deltay
      parameter (deltay =1.11e5)
      
      dvdy = (v1-v0) / ( 2. * deltay )

      end

c     -------------------------------------------------------------------------
c     Meridional derivative of T (DTDY)
c     -------------------------------------------------------------------------

      subroutine calc_DTDY (dtdy,t1,t0)

      implicit none
  
c     Argument declaration
      real  dtdy   ! Derivative of T in meridional direction [K/m]
      real  t1     ! T @ LAT + 1 DLAT [K]
      real  t0     ! T @ LAT - 1 DLAT [K]

c     Physical parameters
      real      deltay
      parameter (deltay =1.11e5)
      
      dtdy = (t1-t0) / ( 2. * deltay )

      end

c     -------------------------------------------------------------------------
c     Meridional derivative of TH (DTHDY)
c     -------------------------------------------------------------------------

      subroutine calc_DTHDY (dthdy,t1,t0,p)

      implicit none
  
c     Argument declaration
      real  dthdy  ! Derivative of TH in meridional direction [K/m]
      real  t1     ! TH @ LAT + 1 DLAT [K]
      real  t0     ! TH @ LAT - 1 DLAT [K]
      real  p      ! P [hPa]

c     Physical parameters
      real      deltay
      parameter (deltay =1.11e5)
      real      rdcp,pref
      data      rdcp,pref /0.286,1000./

      dthdy = (pref/p)**rdcp * (t1-t0) / ( 2. * deltay )

      end

c     -------------------------------------------------------------------------
c     Wind shear of U (DUDP)
c     -------------------------------------------------------------------------

      subroutine calc_DUDP (dudp,u1,u0,p1,p0)

      implicit none
 
c     Argument declaration
      real  dudp   ! Wind shear [m/s per Pa]
      real  u1     ! U @ P + 1 DP [m/s]
      real  u0     ! U @ P - 1 DP [m/s]
      real  p1     ! P + 1 DP [hPa]
      real  p0     ! P - 1 DP [hPa]

      dudp = 0.01 * (u1-u0) / (p1-p0)

      end

c     -------------------------------------------------------------------------
c     Wind shear of V (DVDP)
c     -------------------------------------------------------------------------

      subroutine calc_DVDP (dvdp,v1,v0,p1,p0)

      implicit none

c     Argument declaration
      real  dvdp   ! Wind shear [m/s per Pa]
      real  v1     ! V @ P + 1 DP [m/s]
      real  v0     ! V @ P - 1 DP [m/s]
      real  p1     ! P + 1 DP [hPa]
      real  p0     ! P - 1 DP [hPa]

      dvdp = 0.01 * (v1-v0) / (p1-p0)

      end

c     -------------------------------------------------------------------------
c     Vertical derivative of T (DTDP)
c     -------------------------------------------------------------------------

      subroutine calc_DTDP (dtdp,t1,t0,p1,p0)
 
      implicit none

c     Argument declaration
      real  dtdp   ! Vertical derivative of T [K/Pa]
      real  t1     ! T @ P + 1 DP [K]
      real  t0     ! T @ P - 1 DP [K]
      real  p1     ! P + 1 DP [hPa]
      real  p0     ! P - 1 DP [hPa]

      dtdp = 0.01 * (t1-t0) / (p1-p0)

      end

c     -------------------------------------------------------------------------
c     Vertical derivative of TH (DTHDP)
c     -------------------------------------------------------------------------

      subroutine calc_DTHDP (dthdp,t1,t0,p1,p0,p,t)

      implicit none
 
c     Argument declaration
      real  dthdp  ! Vertical derivative of TH [K/Pa]
      real  t1     ! T @ P + 1 DP [K]
      real  t0     ! T @ P - 1 DP [K]
      real  t      ! T [K]
      real  p1     ! P + 1 DP [hPa]
      real  p0     ! P - 1 DP [hPa]
      real  p      ! P [hPa]

c     Physical parameters
      real      rdcp,tzero,pref
      data      rdcp,tzero,pref /0.286,273.15,1000./
 
c     Auxiliary variables 
      real      tk1,tk0,tk

      if (t0.lt.100.) then
         tk0 = t0 + tzero
      endif
      if (t1.lt.100.) then
         tk1 = t1 + tzero
      endif      
      if (t.lt.100.) then
         tk  = t  + tzero
      endif      

      dthdp = 0.01*(pref/p)**rdcp * 
     >      ( (tk1-tk0)/(p1-p0) - rdcp * tk/p ) 

      end

c     -------------------------------------------------------------------------
c     Squared Brunt-Vaisäla frequency (NSQ)
c     -------------------------------------------------------------------------

      subroutine calc_NSQ (nsq,dthdp,th,rho)

      implicit none
 
c     Argument declaration
      real  nsq    ! Squared Brunt-Vaisäla frequency [s^-1]
      real  dthdp  ! D(TH)/DP [K/Pa]
      real  th     ! K
      real  rho    ! Density [kg m^-3] 

c     Physical parameters
      real      g 
      parameter (g=9.81)

      nsq = -g**2/th * rho * dthdp

      end

c     -------------------------------------------------------------------------
c     Relative vorticity (RELVORT)
c     -------------------------------------------------------------------------

      subroutine calc_RELVORT (relvort,dudy,dvdx,u,lat)

      implicit none
 
c     Argument declaration
      real  relvort    ! Relative vorticity [s^-1]
      real  u          ! Zonal wind component [m/s]
      real  dudy       ! du/dy [s^-1]
      real  dvdx       ! dv/dx [s^-1]
      real  lat        ! Latitude [deg]

c     Physical parameters
      real      pi180
      parameter (pi180=3.14159/180.)
      real      deltay
      parameter (deltay =1.11e5)

      relvort = dvdx - dudy + u * pi180/deltay * tan(pi180 * lat) 
      
      end

c     -------------------------------------------------------------------------
c     Absolute vorticity (ABSVORT)
c     -------------------------------------------------------------------------

      subroutine calc_ABSVORT (absvort,dudy,dvdx,u,lat)

      implicit none

c     Argument declaration
      real  absvort    ! Absolute vorticity [s^-1]
      real  u          ! Zonal wind component [m/s]
      real  dudy       ! du/dy [s^-1]
      real  dvdx       ! dv/dx [s^-1]
      real  lat        ! Latitude [deg]

c     Physical parameters
      real      pi180
      parameter (pi180=3.14159/180.)
      real      deltay
      parameter (deltay =1.11e5)
      real      omega
      parameter (omega=7.292e-5)

      absvort = dvdx - dudy + u * pi180/deltay * tan(pi180 * lat) +
     >          2. * omega * sin(pi180 * lat) 
      
      end

c     -------------------------------------------------------------------------
c     Divergence (DIV)
c     -------------------------------------------------------------------------

      subroutine calc_DIV (div,dudx,dvdy,v,lat)

      implicit none
 
c     Argument declaration
      real  div        ! Divergence [s^-1]
      real  v          ! Meridional wind component [m/s]
      real  dudx       ! du/dx [s^-1]
      real  dvdy       ! dv/dy [s^-1]
      real  lat        ! Latitude [deg]

c     Physical parameters
      real      pi180
      parameter (pi180=3.14159/180.)
      real      deltay
      parameter (deltay =1.11e5)

      div = dudx + dvdy - v * pi180/deltay * tan(pi180 * lat) 
      
      end

c     -------------------------------------------------------------------------
c     Deformation (DEF)
c     -------------------------------------------------------------------------

      subroutine calc_DEF (def,dudx,dvdx,dudy,dvdy)

      implicit none

c     Argument declaration
      real  def        ! Deformation [s^-1]
      real  dudx       ! du/dx [s^-1]
      real  dvdx       ! dv/dy [s^-1]
      real  dudy       ! du/dx [s^-1]
      real  dvdy       ! dv/dy [s^-1]

c     Physical parameters
      real      pi180
      parameter (pi180=3.14159/180.)

      def = sqrt( (dvdx+dudy)**2 + (dudx-dvdy)**2 )
      
      end

c     -------------------------------------------------------------------------
c     Potential Vorticity (PV)
c     -------------------------------------------------------------------------

      subroutine calc_PV (pv,absvort,dthdp,dudp,dvdp,dthdx,dthdy)

      implicit none

c     Argument declaration
      real  pv         ! Ertel-PV [PVU]
      real  absvort    ! Absolute vorticity [s^-1]
      real  dthdp      ! dth/dp [K/Pa]
      real  dudp       ! du/dp [m/s per Pa]
      real  dvdp       ! dv/dp [m/s per Pa]
      real  dthdx      ! dth/dx [K/m]
      real  dthdy      ! dth/dy [K/m]

c     Physical and numerical parameters
      real      scale
      parameter (scale=1.E6)
      real      g
      parameter (g=9.80665)

      pv = -scale * g * ( absvort * dthdp + dudp * dthdy - dvdp * dthdx)

      end

c     -------------------------------------------------------------------------
c     Richardson number (RI)
c     -------------------------------------------------------------------------

      subroutine calc_RI (ri,dudp,dvdp,nsq,rho)

      implicit none

c     Argument declaration
      real  ri         ! Richardson number
      real  dudp       ! Du/Dp [m/s per Pa]
      real  dvdp       ! Dv/Dp [m/s per Pa]
      real  nsq        ! Squared Brunt-Vailälä frequency [s^-1]
      real  rho        ! Density [kg/m^3]

c     Physical and numerical parameters
      real      g
      parameter (g=9.80665)

      ri = nsq / ( dudp**2 + dvdp**2 ) / ( rho * g )**2

      end

c     -------------------------------------------------------------------------
c     Ellrod and Knapp's turbulence indicator (TI)
c     -------------------------------------------------------------------------

      subroutine calc_TI (ti,def,dudp,dvdp,rho)
      
      implicit none

c     Argument declaration
      real  ti         ! Turbulence idicator
      real  def        ! Deformation [s^-1]
      real  dudp       ! Du/Dp [m/s per Pa]
      real  dvdp       ! Dv/Dp [m/s per Pa]
      real  rho        ! Density [kg/m^3]

c     Physical and numerical parameters
      real      g
      parameter (g=9.80665)

      ti = def * sqrt ( dudp**2 + dvdp**2 ) * ( rho * g )

      end

c     -------------------------------------------------------------------------
c     Distance from starting position
c     -------------------------------------------------------------------------

      subroutine calc_DIST0 (dist0,lon0,lat0,lon1,lat1)

      implicit none

c     Argument declaration
      real  dist0      ! Distance from starting position [km]
      real  lon0,lat0  ! Starting position
      real  lon1,lat1  ! New position 

c     Externals 
      real     sdis
      external sdis

      dist0 = sdis(lon0,lat0,lon1,lat1)
      
      end

c     -------------------------------------------------------------------------
c     Heading of the trajectory (HEAD)
c     -------------------------------------------------------------------------

      subroutine calc_HEAD (head,lon0,lat0,lon1,lat1)
      
      implicit none

c     Argument declaration
      real  head       ! Heading angle (in deg) relativ to zonal direction
      real  lon0,lat0  ! Starting position
      real  lon1,lat1  ! New position 

c     Physical parameters
      real      pi180
      parameter (pi180=3.14159/180.)

c     Auixiliary variables
      real     dx,dy

      dx = (lon1-lon0) * cos(pi180*0.5*(lat0+lat1))
      dy = lat1-lat0

      call getangle(1.,0.,dx,dy,head)

      end

c 
c     *************************************************************************
c     Auxiliary subroutines and functions
c     *************************************************************************
     
c     -------------------------------------------------------------------------
c     Saturation vapor pressure over water
c     -------------------------------------------------------------------------

      real function esat(t)
 
C     This function returns the saturation vapor pressure over water
c     (mb) given the temperature (Kelvin).
C     The algorithm is due to Nordquist, W. S. ,1973: "Numerical
C     Approximations of Selected Meteorological Parameters for Cloud
C     Physics Problems" ECOM-5475, Atmospheric Sciences Laboratory,
c     U. S. Army Electronics Command, White Sands Missile Range,
c     New Mexico 88002.
 
      real p1,p2,c1,t
 
      p1=11.344-0.0303998*t
      p2=3.49149-1302.8844/t
      c1=23.832241-5.02808*log10(t)
      esat=10.**(c1-1.3816e-7*10.**p1+8.1328e-3*10.**p2-2949.076/t)

      end

c     --------------------------------------------------------------------------
c     Angle between two vectors
c     --------------------------------------------------------------------------

      SUBROUTINE getangle (ux1,uy1,ux2,uy2,angle)

c     Given two vectors <ux1,uy1> and <ux2,uy2>, determine the angle (in deg)
c     between the two vectors.

      implicit none

c     Declaration of subroutine parameters
      real ux1,uy1
      real ux2,uy2
      real angle

c     Auxiliary variables and parameters
      real len1,len2,len3
      real val1,val2,val3
      real vx1,vy1
      real vx2,vy2
      real pi
      parameter (pi=3.14159265359)

      vx1 = ux1
      vx2 = ux2
      vy1 = uy1
      vy2 = uy2

      len1=sqrt(vx1*vx1+vy1*vy1)
      len2=sqrt(vx2*vx2+vy2*vy2)

      if ((len1.gt.0.).and.(len2.gt.0.)) then
         vx1=vx1/len1
         vy1=vy1/len1
         vx2=vx2/len2
         vy2=vy2/len2

         val1=vx1*vx2+vy1*vy2
         val2=-vy1*vx2+vx1*vy2

         len3=sqrt(val1*val1+val2*val2)

         if ( (val1.ge.0.).and.(val2.ge.0.) ) then
            val3=acos(val1/len3)
         else if ( (val1.lt.0.).and.(val2.ge.0.) ) then
            val3=pi-acos(abs(val1)/len3)
         else if ( (val1.ge.0.).and.(val2.le.0.) ) then
            val3=-acos(val1/len3)
         else if ( (val1.lt.0.).and.(val2.le.0.) ) then
            val3=-pi+acos(abs(val1)/len3)
         endif
      else
         val3=0.
      endif

      angle=180./pi*val3

      END


c     --------------------------------------------------------------------------
c     Spherical distance between lat/lon points                                                          
c     --------------------------------------------------------------------------

      real function sdis(xp,yp,xq,yq)
c
c     calculates spherical distance (in km) between two points given
c     by their spherical coordinates (xp,yp) and (xq,yq), respectively.
c
      real      re
      parameter (re=6370.)
      real      pi180
      parameter (pi180=3.14159/180.)
      real      xp,yp,xq,yq,arg

      arg=sin(pi180*yp)*sin(pi180*yq)+
     >    cos(pi180*yp)*cos(pi180*yq)*cos(pi180*(xp-xq))
      if (arg.lt.-1.) arg=-1.
      if (arg.gt.1.) arg=1.

      sdis=re*acos(arg)

      end

