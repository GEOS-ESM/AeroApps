!
! This include file is meant to associate the main PlumeRise object with flat
! arrays/scalars in the original thread unsafe implementation.
!

! The PlumeRise object
! --------------------
  Type(PlumeRise), target    ::  this

!  The flattened variables
!  -----------------------
   real,pointer,dimension(:) ::  w,t,qv,qc,qh,qi,sc,  &  ! blob
                                 vth,vti,rho,txs,  &
                                 est,qsat,qpas,qtotal,vel_p,rad_p
   
   real,pointer,dimension(:) ::  wc,wt,tt,qvt,qct,qht,qit,sct,vel_t,rad_t
   real,pointer,dimension(:) ::  dzm,dzt,zm,zt,vctr1,vctr2, &
                                 vt3dc,vt3df,vt3dk,vt3dg,scr1

!
                                                    ! environment at plume grid
   real,pointer,dimension(:) ::  pke,the,thve,thee,pe,te,qvenv,rhe,dne,sce,upe,vpe,vel_e
                                                    ! environment at RAMS  grid
   real,pointer,dimension(:) ::  ucon,vcon,wcon,thtcon ,rvcon,picon,tmpcon,&
                                 dncon,prcon,zcon,zzcon,scon 

   real,pointer     :: DZ,DQSDZ,VISC(:),VISCOSITY,TSTPF   
   integer, pointer :: N,NM1,L
   
   real,pointer     :: ADVW,ADVT,ADVV,ADVC,ADVH,ADVI,CVH(:),CVI(:),ADIABAT,&
                       WBAR,ALAST(:),VHREL,VIREL  ! advection

   real,pointer     :: ZSURF,ZBASE,ZTOP
   integer, pointer :: LBASE

   real,pointer     :: AREA,RSURF,ALPHA,RADIUS(:)  ! entrain

   real,pointer     :: HEATING(:),FMOIST,BLOAD   ! heating

   real,pointer     :: DT,TIME,TDUR
   integer,pointer  :: MINTIME,MDUR,MAXTIME
   
   real,pointer     :: ztop_(:)

!  Associate flat arrays/scalar pointers with object
!  -------------------------------------------------

     w => this%w
     t => this%t
     qv => this%qv
     qc => this%qc
     qh => this%qh
     qi => this%qi
     sc => this%sc
     vth => this%vth
     vti => this%vti
     rho => this%rho
     txs => this%txs
     est => this%est
     qsat => this%qsat
     qpas => this%qpas
     qtotal => this%qtotal
     wc => this%wc
     wt => this%wt
     tt => this%tt
     qvt => this%qvt
     qct => this%qct
     qht => this%qht
     qit => this%qit
     sct => this%sct
     dzm => this%dzm
     dzt => this%dzt
     zm => this%zm
     zt => this%zt
     vctr1 => this%vctr1
     vctr2 => this%vctr2
     vt3dc => this%vt3dc
     vt3df => this%vt3df
     vt3dk => this%vt3dk
     vt3dg => this%vt3dg
     scr1 => this%scr1
     pke => this%pke
     the => this%the
     thve => this%thve
     thee => this%thee
     pe => this%pe
     te => this%te
     qvenv => this%qvenv
     rhe => this%rhe
     dne => this%dne
     sce => this%sce
     ucon => this%ucon
     vcon => this%vcon
     wcon => this%wcon
     thtcon => this%thtcon         !1-km
     rvcon => this%rvcon           !1-km
     picon => this%picon           !1-km
     tmpcon => this%tmpcon
     dncon => this%dncon           !1-km
     prcon => this%prcon           !1-km
     zcon => this%zcon             !1-km
     zzcon => this%zzcon           !1-m1
     scon => this%scon
     dz => this%dz
     dqsdz => this%dqsdz
     visc => this%visc
     viscosity => this%viscosity
     tstpf => this%tstpf
     n => this%n
     nm1 => this%nm1
     l => this%l
     advw => this%advw
     advt => this%advt
     advv => this%advv
     advc => this%advc
     advh => this%advh
     advi => this%advi
     cvh => this%cvh
     cvi => this%cvi
     adiabat => this%adiabat
     wbar => this%wbar
     alast => this%alast
     vhrel => this%vhrel
     virel => this%virel
     zsurf => this%zsurf
     zbase => this%zbase
     ztop => this%ztop
     lbase => this%lbase
     area => this%area
     rsurf => this%rsurf
     alpha => this%alpha
     radius => this%radius
     heating => this%heating
     fmoist => this%fmoist
     bload => this%bload
     dt => this%dt
     time => this%time
     tdur => this%tdur
     mintime => this%mintime
     mdur => this%mdur
     maxtime => this%maxtime
     ztop_ => this%ztop_
     
     
     !srf-04may2009-AWE
     upe => this%upe
     vpe => this%vpe
     vel_e=>this%vel_e
     vel_p=>this%vel_p
     rad_p=>this%rad_p
     vel_t=>this%vel_t
     rad_t=>this%rad_t
     !-srf-end
