32c32
< ! #  Release Date :   August 2014    (2.7)                      #
---
> ! #  Release Date :   June 2014      (2.7)                      #
41,42c41,42
< ! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
< ! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
---
> ! #       SURFACE-LEAVING Mark 2, LBBF Jacobians       (2.7)    #
> ! #       VBRDF upgrade / TAYLOR Series / ThreadSafe   (2.7)    #
128a129,135
> ! A Sayer, 04 Nov 2014
> ! Apply a normalisation factor of 1/(pi/mu0) to output water-leaving reflectance, to
> ! bring things in line with results from e.g. 6S simulations and expected behaviour.
> ! Think this is a subtlety related to reflectance vs. normalised radiance treatment,
> ! although it is very obvious if you don't do it. Correction applied at end of the
> ! subroutine.
> 
453a461,482
>       ! Correction of normalisation factor: divide by (pi/mu0). A Sayer 04 Nov 2014.
>       ! Have kept this outside the above loops so it is more obvious.
>       WLeaving_Iso(J)=WLeaving_Iso(J)/(pi/cos(szas(J)*dtr))
>       do i = 1, nstreams
>          WLeaving_SD(J,I) = WLeaving_SD(J,I)/(pi/cos(szas(J)*dtr))
>       enddo
>       do i = 1, nvzas
>          WLeaving_SV(J,I) = WLeaving_SV(J,I)/(pi/cos(szas(J)*dtr))
>       enddo
>       ! And I think we need to do the same for the gradients, too.
>       dC_WLeaving_Iso(J)=dC_WLeaving_Iso(J)/(pi/cos(szas(J)*dtr))
>       dW_WLeaving_Iso(J)=dW_WLeaving_Iso(J)/(pi/cos(szas(J)*dtr))
>       do i = 1, nstreams
>          dC_WLeaving_SD(J,I) = dC_WLeaving_SD(J,I)/(pi/cos(szas(J)*dtr))
>          dW_WLeaving_SD(J,I) = dW_WLeaving_SD(J,I)/(pi/cos(szas(J)*dtr))
>       enddo
>       do i = 1, nvzas
>          dC_WLeaving_SV(J,I) = dC_WLeaving_SV(J,I)/(pi/cos(szas(J)*dtr))
>          dW_WLeaving_SV(J,I) = dW_WLeaving_SV(J,I)/(pi/cos(szas(J)*dtr))
>       enddo
> 
> 
576a606,619
> ! Updated A Sayer November 03 2014:
> ! Extended functionality down to 200 nm. Achieved by:
> ! - Extended data arrays down to 200 nm (another 40 elements).
> ! - Changed logic check for contribution to 0.2-0.9 microns from 0.4-0.9 microns, and started table lookup calculation from 0.2 microns instead of 0.4 microns.
> ! Note, this is based on a simple extension of the published optical model for vis wavelengths. Possible that other scatterers/absorbers
> ! which are neglected in this model may be important at UV wavelengths.
> ! Do linear interpolation of optical property LUTs, rather than nearest neighbour, to remove discontinuities. Achieved by:
> ! - Replicated final element of LUTs to avoid potential for extrapolation errors.
> ! - Replace nint() call with floor() call to correctly get lower bound
> ! - Define variable dwl, fractional distance along the 5 nm LUT grid
> ! - Implement the interpolation using dwl rather than direct lookup of nearest value.
> ! Also:
> ! - Corrected Prieur and Sathyendranath, Limnol. Oceanogr. reference year to 1981 instead of 1983.
> ! - Corrected typo in water scattering coefficient at 410 nm: 0.0068 was written instead of 0.0061. Removes artificial spike at 410 nm in calculated reflectance.
615c658
<       real water_abs_coeff(101),water_scat_coeff(101),abs_coeff_chl(101)
---
>       real water_abs_coeff(142),water_scat_coeff(142),abs_coeff_chl(142)
619c662
<       real a_wat,b_wat,b_wat_all,a_chl,f,a_tot,b_tot,a_ph,a_cdom,v,bp,bbp,eta,x,z
---
>       real a_wat,b_wat,b_wat_all,a_chl,f,a_tot,b_tot,a_ph,a_cdom,v,bp,bbp,eta,x,z,dwl
631a675,682
>        3.0700,2.5300,1.9900,1.6500,1.3100,&
>        1.1185,0.9270,0.8235,0.7200,0.6395,&
>        0.5590,0.5080,0.4570,0.4150,0.3730,&
>        0.3305,0.2880,0.2515,0.2150,0.1780,&
>        0.1410,0.1230,0.1050,0.0907,0.0844,&
>        0.0761,0.0678,0.0620,0.0561,0.0512,&
>        0.0463,0.0421,0.0379,0.0340,0.0300,&
>        0.0260,0.0220,0.0206,0.0191,0.0181,&
652c703
<        6.7858/
---
>        6.7858,6.7858/
661c712,720
<        0.0076,0.0072,0.0086,0.0064,0.0061,&
---
>        0.1510,0.1350,0.1190,0.1093,0.0995,&
>        0.0908,0.0820,0.0753,0.0685,0.0630,&
>        0.0575,0.0530,0.0485,0.0450,0.0415,&
>        0.0384,0.0353,0.0329,0.0305,0.0284,&
>        0.0262,0.0246,0.0229,0.0215,0.0200,&
>        0.0188,0.0175,0.0164,0.0153,0.0144,&
>        0.0134,0.0127,0.0120,0.0113,0.0106,&
>        0.0100,0.0094,0.0089,0.0084,0.0080,&
>        0.0076,0.0072,0.0068,0.0064,0.0061,&
681c740
<        0.0002/
---
>        0.0002,0.0002/
685a745,753
>       data abs_coeff_chl/&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.053,0.123,0.195,0.264,&
>        0.335,0.405,0.476,0.546,0.617,&
698,706c766,774
<        0.215,0.0,0.0,0.0,0.0,&
<        0.0,0.0,0.0,0.0,0.0,&
<        0.0,0.0,0.0,0.0,0.0,&
<        0.0,0.0,0.0,0.0,0.0,&
<        0.0,0.0,0.0,0.0,0.0,&
<        0.0,0.0,0.0,0.0,0.0,&
<        0.0,0.0,0.0,0.0,0.0,&
<        0.0,0.0,0.0,0.0,0.0,&
<        0.0/
---
>        0.215,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000,0.000,0.000,0.000,&
>        0.000,0.000/
716c784
<       if (wl.lt.0.400.or.wl.gt.0.900)then
---
>       if (wl.lt.0.200.or.wl.gt.0.900)then
721,724c789,797
<       iwl=1+nint((wl-0.400)/0.005)
<       a_wat=water_abs_coeff(iwl)
<       b_wat_all=water_scat_coeff(iwl)
<       a_chl=abs_coeff_chl(iwl)
---
>       iwl=1+floor((wl-0.200)/0.005)
> ! 03 Nov 2014 A Sayer now linear interpolation rather than nearest neighbour
> !      a_wat=water_abs_coeff(iwl)
> !      b_wat_all=water_scat_coeff(iwl)
> !      a_chl=abs_coeff_chl(iwl)
>       dwl=(wl-0.200)/0.005-floor((wl-0.200)/0.005)
>       a_wat=water_abs_coeff(iwl)+dwl*(water_abs_coeff(iwl+1)-water_abs_coeff(iwl))
>       b_wat_all=water_scat_coeff(iwl)+dwl*(water_scat_coeff(iwl+1)-water_scat_coeff(iwl))
>       a_chl=abs_coeff_chl(iwl)+dwl*(abs_coeff_chl(iwl+1)-abs_coeff_chl(iwl))
