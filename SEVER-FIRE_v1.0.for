

c     This code is the description of the global fire model SEVER-FIRE v1.0
c     Reference: Venevsky et al., 2018, Geoscientific Model Development
c     https://doi.org/10.5194/gmd-2018-178

c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE SEVER_FIRE
c     Describes separately lightning and human-induced fires
c..............................................................................      

      subroutine sever_fire(pftpar,dtmax,dtmin,dtemp,litter_ag,
     * acflux_fire,afire_frac,dcflux_fire,lm_ind,rm_ind,sm_ind,hm_ind,
     * nind,dw1,present,tree,fdi,pconv,elev,fpc_grid,humfires_old,
     * lightfires_old,dist_city,annum_fire_hum,annum_fire_nat,
     * annum_fire,an_areafires_hum,an_areafires_nat,an_areafires,
     * popdens,rurratio,wealth,lat,timing_rural,timing_urban,lon,year,
     * nmfl,nmfh,area_hum,area_nat,wind_d)
        
      implicit none
c     PARAMETERS
      integer npft,npftpar,nsoilpar
        parameter (npft=10,npftpar=35,nsoilpar=7)
      real pi
        parameter (pi=3.14159265)
      real minfuel
        parameter (minfuel=100.0)  !fuel threshold to carry a fire (gC/m2)
      real pixelarea

c     ARGUMENTS
      real pftpar(1:npft,1:npftpar)
      real dtemp(1:365)
	  real dtmax(1:365),dtmin(1:365)
      real litter_ag(1:npft)
	  real lat
      real acflux_fire
      real lm_ind(1:npft),rm_ind(1:npft)
      real sm_ind(1:npft),hm_ind(1:npft)
      real nind(1:npft)
      real dw1(1:365)
	  real wind_d(1:365)
      real afire_frac
      logical present(1:npft),tree(1:npft)
     
c     LOCAL VARIABLES
      integer pft,d
      real fire_length,fuel,fire_prob(1:365)
	  real fdi(1:365)                         !daily Nesterov Index (Venevsky et.al., 2002)
	  real light_ign (1:365)                  !potential daily lightning ignitions per grid cell
	  real human_ign (1:365)                  !potential daily human ignitions per grid cell
	  real num_fire_nat (1:365)               !natural number of fires per grid cell
	  real num_fire_hum (1:365)               !human number of fires per grid cell
	  real num_fire (1:365)                   !number of fires per mln. ha
	  real area_burnt_hum(1:365)
	  real area_burnt_nat(1:365)
	  real area_burnt (1:365)                 !potential daily lightning ignitions
	  real dfire_frac (1:365)                 !potential daily lightning ignitions
      real pconv (1:365)					  !convective precipitation mm/day
	  real elev								  !elevation in m
      real fire_index,disturb
      real resist(1:npft)
      real moistfactor,flam,litter_ag_total
      real fire_term
      real dcflux_fire(1:365)
	  real dens_fuel,depth_fuel
      real fpc_grid(1:npft)
      real annum_fire
	  real an_areafires
      real U_front
	  real humfires_old		                  !number of ongoing fires from the previous day human
	  real lightfires_old                     !number of ongoing fires from the previous day lightning
      real annum_fire_hum,annum_fire_nat,
     * an_areafires_hum,an_areafires_nat
	  real dist_city
	  real timing_rural(1:365),timing_urban(1:365)
	  integer year
	  real lon
	  real nmfl
	  real nmfh
	  real dcflux_fire_nat(1:365),dcflux_fire_hum(1:365)
      real afire_frac_hum,afire_frac_nat,
     *     acflux_fire_nat,acflux_fire_hum
	  real dfire_frac_nat(1:365),dfire_frac_hum(1:365)
      real emiss
         parameter(emiss=0.55)               !emissitivity of flames unitless Venevsky, et.al. 2002

      real fire_durat,wind_speed
         parameter (wind_speed=3.7)
      integer j
	  real popdens							 !population density person/km2
      real rurratio						     !percenatge of rural population in a baseline year
      real wealth                            !ignition potential for one person, depending on wealth
	  integer vsp
	  real fire_durat_old
	  real area_hum,area_nat

c    calculate bulk density of fuel
      dens_fuel=0.01
      do j=1,npft
         dens_fuel=dens_fuel+fpc_grid(j)*pftpar(j,35)
      enddo


c     Calculate total above-ground litter
      litter_ag_total=0.0
      do pft=1,npft
        litter_ag_total=litter_ag_total+litter_ag(pft)
      enddo
      if (dens_fuel.le.0.0) then 
	    depth_fuel=0.0
	  else
        depth_fuel=(litter_ag_total/dens_fuel*0.1)
c     depth of fuel is in cantimeters
      endif


c     Calculate moisture of extinction by weightning of PFTs
      moistfactor=0.0
      do pft=1,npft
        flam=pftpar(pft,6)
        if (litter_ag_total.gt.0.0) then
          moistfactor=moistfactor+(litter_ag(pft)/litter_ag_total)*flam
        else
          moistfactor=0.0
        endif   
      enddo


c     First calculate the annual fraction of the grid cell affected by fire
c     Initialise
      fire_length=0.0
      fuel=0.0
      acflux_fire=0.0
      annum_fire=0.0
      an_areafires=0.0
      U_front=0.0
      annum_fire_hum=0.0
      annum_fire_nat=0.0
      an_areafires_hum=0.0
      an_areafires_nat=0.0

      do d=1,365
        dcflux_fire(d)=0.0
        num_fire(d)=0.0
        area_burnt(d)=0.0
        dfire_frac(d)=0.0
        human_ign(d)=0.0
	    light_ign(d)=0.0
      enddo


c     Assign a minimum fire fraction (for presentational purposes)
      afire_frac=0.001
      afire_frac_nat=0.0005
	  afire_frac_hum=0.0005


c     Assign PFT resistance to fire
      do pft=1,npft
        resist(pft)=pftpar(pft,8)
      enddo


c     Calculate the length of the fire season (units=days)
      do d=1,365

c     Calculate today's fire probability, fire_prob
c     Assume fire is only possible when temperature is above zero
        if (dtemp(d).gt.0.0.and.moistfactor.gt.0.0) then
          fire_prob(d)=EXP((-pi)*(dw1(d)/moistfactor)**2)
        else 
          fire_prob(d)=0.0
        endif
        fire_length=fire_length+fire_prob(d)
      enddo

      call firedanger(dtmax,dtmin,fire_prob,fpc_grid,fdi,wind_d)
	  call natignition(pconv,elev,light_ign,dtemp,dw1,
     *     dens_fuel,depth_fuel,lat,nmfl,lon,year)
      call humignition(human_ign,popdens, rurratio,
     * 	dist_city,wealth,timing_rural,timing_urban,lat,dens_fuel)


c     Make ignitions inpossible when there is no fuel
      do pft=1,npft
        fuel=fuel+litter_ag(pft)
      enddo
      if (fuel.lt.minfuel) then
         do d=1,365
	     light_ign(d)=0.0
	     human_ign(d)=0.0
	     enddo
      endif
      fire_durat=2*(1-exp(-0.001*dist_city))
	if (fire_durat.lt.0.20) fire_durat=0.20


c     allow at least six hours of burn for lighning fires
      if (dist_city.lt.0.0) fire_durat=0.0
c     city - no fire possible
	  if (fire_durat.gt.2.0) fire_durat=2.0


c     fire duration in days
c     salmon-chavill forest fire data	 
c     areas burnt and number of fires are calculated  
      do d=1,365

       if (fire_durat.gt.1.0) then
	       num_fire_hum(d)=fdi(d)*human_ign(d)
	       num_fire_nat(d)=fdi(d)*light_ign(d)


           U_front=(1+wind_d(d))*(3*emiss)/(dens_fuel*
     *     (16.0+100.0*dw1(d)))

c     frontal speed is in m per second, density of fuel in kg per m3
c     area burnt is in m2
		   area_burnt_hum(d)=pi/8.0*(60.0*60.0*24.0)**2.0*
     *     U_front**2.0*num_fire_hum(d)
     *     + pi/8.0*(60.0*60.0*24.0)**2.0*fire_durat_old**2.0*
     *     U_front**2.0*humfires_old

c 6 hours of active state and other time of a day is smoldering state - why? where did I read this??
		   area_burnt_nat(d)=pi/8.0*(60.0*60.0*24.0)**2.0*
     *     U_front**2.0*num_fire_nat(d)
     *     + pi/8.0*(60.0*60.0*24.0)**2.0*fire_durat_old**2.0*
     *     U_front**2.0*lightfires_old

	       humfires_old=num_fire_hum(d)
           lightfires_old=num_fire_nat(d)
	       fire_durat_old=fire_durat-1.0


c      It is assumed that fires do not continue longer than two days
        else
	     num_fire_hum(d)=fdi(d)*human_ign(d)
	     num_fire_nat(d)=fdi(d)*light_ign(d)

           U_front=(1+wind_d(d))*(3*emiss)/(dens_fuel*
     *     (16.0+100.0*dw1(d)))

           area_burnt_hum(d)=pi/8.0*(60.0*60.0*24.0)**2.0*fire_durat**2.0*
     *     U_front**2.0*num_fire_hum(d)
     *     + pi/8.0*(60.0*60.0*24.0)**2.0*fire_durat_old**2.0*
     *     U_front**2.0*humfires_old

		   area_burnt_nat(d)=pi/8.0*(60.0*60.0*24.0)**2.0*fire_durat**2.0*
     *     U_front**2.0*num_fire_nat(d)
     *     + pi/8.0*(60.0*60.0*24.0)**2.0*fire_durat_old**2.0*
     *     U_front**2.0*lightfires_old

	       humfires_old=0.0
           lightfires_old=0.0
	       fire_durat_old=0.0
       endif

c      Obtain annual number of fires,area burnt, and burnt fraction
c      Number of fires
		  num_fire(d)=num_fire_hum(d)+num_fire_nat(d)
          annum_fire_hum=annum_fire_hum+num_fire_hum(d)
          annum_fire_nat=annum_fire_nat+num_fire_nat(d)
          annum_fire=annum_fire+num_fire(d)
c      Area burnt
          area_burnt(d)=area_burnt_hum(d)+area_burnt_nat(d)
          an_areafires_hum=an_areafires_hum+area_burnt_hum(d)
          an_areafires_nat=an_areafires_nat+area_burnt_nat(d)
          an_areafires=an_areafires_hum+an_areafires_nat
c      burnt fraction
          dfire_frac(d)=area_burnt(d)/(pixelarea(lat))
          dfire_frac_nat(d)=area_burnt_nat(d)/(pixelarea(lat))
          dfire_frac_hum(d)=area_burnt_hum(d)/(pixelarea(lat))
          afire_frac=afire_frac+dfire_frac(d)
		  afire_frac_nat=afire_frac_nat+dfire_frac_nat(d)
		  afire_frac_hum=afire_frac_hum+dfire_frac_hum(d)
      enddo

      
c     Calculate the available fuel (above-ground litter) to carry the fire
      do pft=1,npft
        fuel=fuel+litter_ag(pft)
      enddo


c     Reduce fraction of grid cell affected by fire when fuel
c     becomes limiting (reduced carrying capacity)
      if (fuel.lt.minfuel) then
        afire_frac=0.0
        afire_frac_nat=0.0
	    afire_frac_hum=0.0
      endif


      if (afire_frac.lt.0.001) then
	   afire_frac=0.001
	   afire_frac_nat=0.0005
	   afire_frac_hum=0.0005
	  else
	  endif

      if (afire_frac.ge.1.0) afire_frac=1.0
      if (afire_frac_hum.ge.1.0) then
	    afire_frac_hum=0.95555
	    afire_frac_nat=0.04445
      else
	  endif

	  if (afire_frac_nat.ge.1.0) then
        afire_frac_nat=0.95555
	    afire_frac_hum=0.04445
      else
	  endif

      if ((afire_frac_nat+afire_frac_hum).gt.1.0) then      
      afire_frac_nat=afire_frac_nat/(afire_frac_nat+afire_frac_hum)
      afire_frac_hum=afire_frac_hum/(afire_frac_nat+afire_frac_hum)
      else
	  endif


c     Implement the effect of the fire on vegetation structure and litter
c     in the disturbed fraction.  

c     Each PFT is assigned a resistance to fire, representing the fraction of
c     the PFT which survives a fire. Grasses assumed already to have completed
c     their life cycle and thus are not affected by fire, giving them
c     a competitive advantage against woody PFTs.

      do pft=1,npft
        if (present(pft).and.tree(pft)) then

c         Calculate the fraction of individuals in grid cell which die
          disturb=(1.0-resist(pft))*afire_frac

c         Calculate carbon flux to atmosphere (gC/m2) due to burnt biomass
          acflux_fire=acflux_fire+disturb*(nind(pft)*
     *      (lm_ind(pft)+sm_ind(pft)+hm_ind(pft)+rm_ind(pft)))

c         Update the individual density
          nind(pft)=nind(pft)*(1.0-disturb)
        endif


c       Add combusted litter to carbon flux to atmosphere term
        acflux_fire=acflux_fire+(afire_frac*litter_ag(pft))
        acflux_fire_nat=acflux_fire*
     *	  afire_frac_nat/(afire_frac_nat+afire_frac_hum)
        acflux_fire_hum=acflux_fire*
     *	  afire_frac_hum/(afire_frac_nat+afire_frac_hum)


c       Update the above ground litter term
        litter_ag(pft)=(1.0-afire_frac)*litter_ag(pft)
      enddo


      do d=1,365
	   if (fire_length.eq.0) then
        dcflux_fire(d)=0.0
        dcflux_fire_nat(d)=0.0
        dcflux_fire_hum(d)=0.0
       else
        dcflux_fire(d)=acflux_fire*fire_prob(d)/fire_length
        dcflux_fire_nat(d)=acflux_fire_nat*fire_prob(d)/fire_length
        dcflux_fire_hum(d)=acflux_fire_hum*fire_prob(d)/fire_length
       endif
      enddo

      nmfh=annum_fire_hum
      nmfl=annum_fire_nat
c     area in thousand hectares
      area_hum=an_areafires_hum*0.001*0.0001
	  area_nat=an_areafires_nat*0.001*0.0001

      return
      end



c----------------------------------------------------------------------
c     FIRE DANGER INDEX
c     Calculation of the daily fire danger index based on the Nesterov Index
c     (equations (4 and 7) Venevsky et.al., 2002)
     
     
      subroutine firedanger(dtmax,dtmin,fire_prob,fpc_grid,fdi,wind_d)
      implicit none
      
      real dtmax(1:365),dtmin(1:365),fire_prob(1:365)
      real fdi(1:365)
	  real wind_d(1:365)
c     daily fire danger index (range 0-1))     
      integer d
      integer j
	  integer npft
	    parameter (npft=10)
      real alpha
        parameter (alpha=0.000337) !coefficient of fire risk function

c     LOCAL VARIABLES
      real fpc_grid(1:npft)
      real d_NI
      real fpc_grass
      real alpha_new
c     daily Nesterov Index (range from 0 to 4000) 
c     it is assumed that flammability of grass higher, i.e. alpha two times higher -Alexander
      
c     calculate relevant monthly minimum and maximum temperatures     
      alpha_new=alpha
	  fpc_grass=fpc_grid(npft-1)+fpc_grid(npft)

	  if (fpc_grass. gt.0.5) alpha_new=2.0*alpha

      d_NI=0.0

      do d=1,365
       
       if (fire_prob(d).le.0.0001.or.dtmin(d).le.0.0) then
	    d_NI=0.0
	   else
c     wind speed dependence is taken from Canadian ISI and scaled to average global wind)
	    d_NI=d_NI+fire_prob(d)*exp(wind_d(d)/3.2*3.6*0.05039)*
     * ((dtmax(d)**2.0 - dtmin(d)**2.0)/4 + 2.0 * (dtmax(d) + dtmin(d)))      
       endif

      fdi(d)=1.0 - exp(-1.0*alpha_new*d_NI)

      enddo
      return
      end



c----------------------------------------------------------------------
c     POTENTIAL NUMBER OF LIGHTNING IGNITIONS

      subroutine natignition(pconv,elev,light_ign,dtemp,dw1,
     *         dens_fuel,depth_fuel,lat,nmfl,lon,year)

	  real light_ign (1:365)            !potential daily lightning ignitions per mln. ha
      real pconv (1:365)				!convective precipitation mm/day
      real dtemp(1:365)
	  real elev						    !elevation in m
      real a,b,c,dc,e                   !coefficients of polinomial approximating number of lightning strokes
	   parameter (a=0.0375,b=-0.0476,c=0.00541,dc=0.000321,e=-0.00000293)
      real ar30
        parameter (ar30=66689.4)
c     area 2.5x2 at 30 degrees Pickering
	  real timethund						  !average time of active lightning in minutes in a cell
	    parameter (timethund=80.)
c     number of flashes is scaled from area 2,5x 2.0 to area 0.5x0.5 and multiplied to 24*60
c     reference Blyth Half a day is taken durnal activity Uman and Krider average time of thunderstorm is 40 minutes, are 100 km2

      real numgcfl(1:365)				!daily number of ground-cloud flashes per grid cell
	  real positfl(1:365)				!number of positive large peak flashes
	  real negatfl(1:365)				!number of neagtive large peak flashes
      real positign(1:365)				!number of ignitions cased by positive large peak flashes
	  real negatign(1:365)				!number of ignitions cased by neagtive large peak flashes
      real dw1(1:365)
	  real lat
	  real lon
	  integer d
	  real nflash_ann
	  real nmfl
	  real numgcfl_month(1:12)
	  integer year
      real nflash_ann_new

      nflash_ann=0.0
   
      do d=1,365

c     convective precipitation is 6 hours triangle distibution is assumed in our case
c     and rectangle in Alan and Pickering
	   pconv(d)=1.65*pconv(d)
       numgcfl(d)=a+pconv(d)*(b+pconv(d)*(c+pconv(d)*(dc+e*pconv(d))))
c
	   if (pconv(d).lt.5.0) numgcfl(d)=0.0
c
       pie = 4. * ATAN(1.)
       dip = pie/180.
	   numgcfl(d)=cos(lat*dip)/0.866*numgcfl(d)
       nflash_ann=nflash_ann+numgcfl(d)
	  					
      enddo

       call monthly(numgcfl_month,numgcfl,.true.)
       call daily(numgcfl_month,numgcfl,.false.)

       do d=1,365

        positfl(d)=numgcfl(d)*0.1*0.075*dens_fuel/4.
        negatfl(d)=numgcfl(d)*0.9*0.025*dens_fuel/4.

c      scaling of fuel efficiency is approximated from Latham
c      Percent of occurence negative 90% and positive 10%

        if (dtemp(d).gt.0.) then
        positign(d)=(1/(1+exp(5.5*(1./1.5)**(((16.-dens_fuel)/16.)*5)*
     *	     1.25-1.2*0.5**((16.-dens_fuel)/16.*5.+0.1)*depth_fuel)))*
     *          positfl(d) 
        negatign(d)=(1/(1+exp(5.5*(1./1.5)**(((16.-dens_fuel)/16.)*5)
     *	             -1.2*0.5**((16.-dens_fuel)/16.*5.)*depth_fuel)))*
     *                 negatfl(d)
        light_ign(d)= (positign(d)+negatign(d))*(1-dw1(d)*0.1)*0.15
     *	            *(1+0.0001*(elev-1000.))**4

        else
        positign(d)=0.
	    negatign(d)=0.
        light_ign(d)=0.
	    endif
      enddo
      nflash_ann=nflash_ann/pixelarea(lat)*1000000.0
      nmfl=nflash_ann

      do d=1,365
	   nflash_d(d)=numgcfl(d)
	  enddo
	            
      return
      end



c----------------------------------------------------------------------
c     POTENTIAL NUMBER OF HUMAN IGNITIONS

      subroutine humignition(human_ign,popdens, rurratio,
     * 	dist_city,wealth,timing_rural,timing_urban,lat,dens_fuel) 
      real human_ign(1:365)
	  real human_ign_rur(1:365)
	  real human_ign_urb(1:365)
      integer d
	  real popdens
      real wealth
	  real rurratio
	  real lat
	  real dens_fuel
	  real timing_rural(1:365),timing_urban(1:365)
c     timing_rural and timing_urban distribution of exposure of population to natural landscape 
c     in a year in realtive units (six typical curves), depending on the region
      
	  if (dist_city.lt.0) dist_city=0.0

      do d=1,365


       human_ign_rur(d)= (exp(wealth*(-7.65*10**(-2.))*(0.0001)*dens_fuel/8.*
     *       rurratio*timing_rural(d)*6.8*popdens**(0.43)*
     *       (pixelarea(lat)*0.000001)
c    wealth is per mln ha, but almost constant 

       human_ign_urb(d)= (exp(wealth*(-7.65*10**(-2.))*(0.0001)*dens_fuel/8.*
     *	     (1-rurratio)*timing_urban(d)*6.8*popdens**(0.43)*
     *       (pixelarea(lat)*0.000001)

       human_ign(d)=human_ign_rur(d)+human_ign_urb(d)
	  enddo

      return
	  end


c----------------------------------------------------------------------
c     PREPARE HUMAN VARIABLES IN HUMAN IGNITIONS

      subroutine humpop_status(year,lat,popdensm,rurratiom,dist_citym,
    *           popdens,rurratio,dist_city,timing_rural,timing_urban)
      real popdensm(1:95)
      real rurratiom(1:95)
      real dist_citym(1:95)
      real popdens            !population density
      real rurratio           !ratio of rural population to urban
      real dist_city          !distance from the nearest city(km)
      real timing_rural(1:365),timing_urban(1:365)
      real lat
      integer d
      integer year,actualyear


      if (lat.ge.0) then
c     Northern Hemisphere: rural population goes in fields at spring and autumn,
c     urban population has contact with natural ecosystems at Easter, Christmas and summer vacation

c     initialisation
       do d=1,365
        timing_rural(d)=0.5
        timing_urban(d)=0.0
       enddo
c     rural population spring
       do d=60,151
        timing_rural(d)=1.0
       enddo
c     rural population autaumn
       do d=244,334
        timing_rural(d)=1.0
       enddo
c     urban population Easter
       do d=92,151
        timing_urban(d)=1.0
       enddo
c    urban population Christmas
       do d=335,365
        timing_urban(d)=1.0
       enddo
c    urban population summer vacations
       do d=182,243
        timing_urban(d)=1.0
       enddo

      else
c     Southern Hemisphere: rural population goes in fields at spring and autumn,
c     urban population has contact with natural ecosystems at Winter vacation time and at Easter
       do d=1,365
        timing_rural(d)=0.75
c     rural population activity higher then in Northern Hemisphere
        timing_urban(d)=0.0
       enddo

c     rural population spring
       do d=60,151
        timing_rural(d)=1.0
       enddo
c     rural population autaumn
       do d=244,334
        timing_rural(d)=1.0
       enddo

c     urban population Easter
       do d=92,151
        timing_urban(d)=1.0
       enddo

c     urban population summer vacations
       do d=1,59
        timing_urban(d)=1.0
       enddo

       do d=304,365
        timing_urban(d)=1.0
       enddo
      endif

c     spinup case the status of population is taken for the year 1957 (spinup 1012 years)
      if (year.le.1012) then

       popdens=0.0
       rurratio=0.0
       dist_city=5000.0
c   we consider that population is absent in spinup years (here is an example from the period 1957-2002)
      else
       actualyear=year-1012+1+2
       popdens=popdensm(actualyear)
       rurratio=rurratiom(actualyear)
       dist_city=dist_citym(actualyear)
      endif

      return
      end


----------------------------------------------------------------------
c     GRID AREA
      function pixelarea(lat)
      implicit none

      real pixelarea,lat
      real pie,dip

      pie = 4. * ATAN(1.)
      dip = pie/180.
c     0.5 degrees at the equator is 111km
      pixelarea=((111.0*10**3*0.5)**2)*cos(lat*dip)

      return
      end






