* switch_planck: version with planck-sized disk masses, switch on
* planck: version with planck-sized disk masses
* 20241107: smearing of all geometry variables
* big_G: more geometry inputs, start name simplification
* G_2_ctr_cyl_geom: newtonian, cylinder and cylidrical counterweight on detector; centered on source; loop over geomteries
* G_2_ctr_cyl: newtonian, cylinder and cylidrical counterweight on detector; centered on source
* G_2_cylinder: newtonian, cylinder and cylidrical counterweight on detector
* G_1_cylinder: newtonian, single cylinder on detector
* newtonian_with_cylinder_2: attempt add cylindrical extrusion on source
* newtonian_with_cylinder: attempt add cylindrical extrusion
* newtonian_sampling_20230720: insert newtonian potential; set importance sampling flag off (0) 
* vpp_point_point_20230214: jcl version of vpp potential
* vpp_point_point_20230101: vpp point-point approx code from A Desai.
* vps_point-point_20221111: vps interaction, point-point geomtery, importance sampling flag,
* saves geometry information
* yukawa_simple_20220523
* based on Yukawa_simple_20200423: added variables for spin vectors in detector and source; gap in
* terms of closest approach; no "Nflag"
* based on yukawa_simple_20200422: essential elements only 
*	  
      program planck
*
      implicit none
*
* variables
*
      real*4 Results(2,100,100) ! results stored in Results(i,j,k)
				! k runs over games, j runs over phase points
				! Results(1,j,k)=W: computed virtual work
				! Results(2,j,k)=points: number of inside points
*
      real*4 Stats(3,100) ! statistics stored in Stats(i,j)
	  		  ! j runs over phase points
	  		  ! Stats(1,j)=omega: source mass phase
			  ! Stats(2,j)=Wmean: mean work at this phase
			  ! Stats(3,j)=Wdev: stan dev of mean
*	  						 
      real*4 range ! range of force
*
      integer*4 throws ! number of trial points per game
      integer*4 games  ! number of games per phase point
      integer*4 phpnts ! number of phase points
*
      integer*4 points ! number of points inside in a game
*
      real*4 W     ! virtual work
      real*4 Wampl ! Fourier amplitude of virtual work
*
      real*4 bounds1(6) ! limits for point in mass 1,
			! cartesian coordinates, order is
			! x1 upper, x1 lower, y1 upper, y1 lower
			! z1 upper, z1 lower
*
      real*4 bounds2(6) ! limits for point in mass 2
			! spherical coords about z axis, order is
			! R upper, R lower, theta upper, theta lower,
			! phi upper, phi lower
*
      real*4 shape1(36)	! shape of mass 1 (detector) in its cartesian
	  	        ! includes detector mode shape
      real*4 shape2(36)	! shape of mass 2 (source) in its cartesians
	  		! includes source mode shape
*
      real*4 gap            ! resting gap
      real*4 xmax,ymax,zmax ! maximum separations of points in x, y, z for bounding box 
      real*4 amp	    ! source mass amplitude
      real*4 omega	    ! source mass phase
      real*4 asinom         ! amp*sin(omega)
      real*4 D(3)	    ! source position relative to detector
	  	            ! in detector cartesians
      real*4 dmin           ! closest approach (fixed) between test masses
*
      integer*4 i,j,k ! index
*
      real*4 volume    ! bounds box volume
      real*4 point1(3) ! point in box 1; x, y, z
      real*4 point2(3) ! point in box 2; scaled R, theta, phi
      real*4 modez     ! z mode shape function
*
      real*4 snth, csth   ! sin(theta), cos(theta)
      real*4 snphi, csphi ! sin(phi), cos(phi)
      real*4 R		  ! radial coord, mass 2
*
      logical*4 inside ! .true. if both points are inside masses
*	  
      real*4 scaler ! scaler factor in integrand
      real*4 dot    ! dot product in integral
*
* spin-dependent variables
*  
      real*4 theta1, phi1   ! polar, azimuthal angle of detector spin in source corrdinates
      real*4 snth1, csth1   ! sin(theta1), cos(theta1)
      real*4 snphi1, csphi1 ! sin(phi1), cos(phi1)
*
      real*4 theta2, phi2   ! polar, azimuthal angle of source spin
      real*4 snth2, csth2   ! sin(theta2), cos(theta2)
      real*4 snphi2, csphi2 ! sin(phi2), cos(phi2)
*
* factors for Vpp
*
      real*4  sdots            ! spin dot product
      real*4  s1dotr,s2dotr    ! spin-displacement dot products
      real*4  dsdotr1,dsdotr2  ! spin-displacement derivatives
      real*4  r1,r2            ! range factors
*
      integer*4 sflag ! importance sampling flag
*
      real*4 Pi ! c/d of a circle
*
      integer*4 u              ! config counter
*
* test mass relative densities
*
      real*4 d1, d2
*
* uncertainties for pinned variables
*
      real*4 ushape120,ushape121,ushape124,ushape125,ushape218,ushape219         !detector G mass x-position uncertainty
      real*4 ud1,ud2
*
*     flag for switching on disk masses
*
      integer*4 pflag
*
* initialize
*
      Pi=3.1415926535 ! c/d of a circle
*
      sflag=0                   ! set to 1 for importance sampling
      pflag = 0                 ! set to 1 for disk masses
*	
      throws=10000000
      games=80
      phpnts=9
*
* input file
* 
      open(20,file='big_G_input.txt',status='unknown')
*
      do 999 u = 1,1
         read(20,*) amp,gap,shape1(1),shape1(2),shape1(3),shape2(1)
     >        ,shape2(2),shape2(3),shape1(18),shape1(19),ushape120,
     >	      ushape121,shape1(22),shape1(23),ushape124,ushape125,
     >        shape2(16),shape2(17),ushape218,ushape219,ud1,ud2
         dmin=gap-amp
         print *,'gap = ',gap
*		
         bounds2(3)=Pi/2.0      ! theta upper
         bounds2(4)=0.0         ! theta lower
         bounds2(5)=2.0*Pi      ! phi upper
         bounds2(6)=0.0         ! phi lower
*
*
         bounds1(1)= 1.0        ! x upper	
         bounds1(2)=-(shape1(1)+1) ! x lower
         bounds1(3)= (shape1(2)/2+1) ! y upper
         bounds1(4)=-(shape1(2)/2+1) ! y lower
*
         shape1(4)=0            ! detector surface curvature polynomial:
         shape1(5)=0            ! s1(4) + s1(5)x + s1(6)y
         shape1(6)=0            ! + s1(7)xy + s1(8)x**2 + s1(9)y**2
         shape1(7)=0            ! all curvature and mode polynomials
         shape1(8)=0            ! have this form and all are
         shape1(9)=0            ! in microns of displacement
*	  
         shape1(10)=0           ! detector z mode shape polynomial
         shape1(11)=0 
         shape1(12)=2/shape1(2)  
         shape1(13)=0 
         shape1(14)=0 
         shape1(15)=0 
*	  
         shape1(16)=0           ! detector x virtual displacement
         shape1(17)=0           ! detector y virtual displacement
*
         shape1(20)=-shape1(1)/2.0 + ushape120 ! cylinder x-position
         shape1(21)=shape1(2)/2.0-shape1(18)+ ushape121 ! cylinder y-position
*
         shape1(24)=-shape1(1)/2.0+ ushape124 ! counterweight x-position
         shape1(25)=-shape1(2)/2.0+shape1(22)+ ushape125 ! counterweight y-position
*	  
*	  
         shape2(4)=0            ! source surface curvature polynomial
         shape2(5)=0 
         shape2(6)=0 
         shape2(7)=0 
         shape2(8)=0 
         shape2(9)=0 
*	  
         shape2(10)=1.0         ! source z mode shape polynomial
         shape2(11)=-1.0/shape2(1) ! 
         shape2(12)=0           ! conventional normalization is 
         shape2(13)=0           ! one micron at cartesian origin
         shape2(14)=0           ! 
         shape2(15)=0           ! 
*
         shape2(18)=5080.0/2.0 + ushape218 ! source cylinder x-position (rel. to source zero)
         shape2(19)=0.0 + ushape219 ! source cylinder y-position (rel. to source zero)
*	  
         D(1)=-shape1(1) - 1000.0  ! source displacement x
         D(2)=0.0               ! y
         D(3)=gap+(shape1(3)/2.0)+(shape2(3)/2.0) ! z
         if(pflag .gt. 0)then
            D(3)=gap+(shape1(3)/2.0)+(shape2(3)/2.0)+shape1(19)
     >       +shape2(17)         ! z
         endif
*
         xmax=shape1(1)+shape2(1)+abs(D(1)) ! max distance between any 2 points in x
         ymax=shape1(2)+shape2(2)+abs(D(2)) ! max distance between any 2 points in y
         zmax=shape1(3)+shape2(3)+abs(D(3))+amp ! max distance between any 2 points in z
*
         bounds1(5)= (shape1(3)/2.0+shape1(19)) ! z upper
         bounds1(6)=-(shape1(3)/2.0+shape1(23)) ! z lower
         bounds2(1)=sqrt(xmax**2+ymax**2+zmax**2) ! R upper, based on max between any points
         bounds2(2)=dmin-1      ! R lower, closest approach between source & detector
*
* detector spins: along z
*
         theta1=0.0             ! detector spin polar angle
         phi1  =0.0             ! detector spin azimuth angle
         csth1=cos(theta1)
         snth1=sin(theta1)
         csphi1=cos(phi1)
         snphi1=sin(phi1)
*
* source spins : along z
*
         theta2 =0.0            ! source spin polar angle
         phi2 =0.0              ! source spin azimuth angle
         csth2=cos(theta2)
         snth2=sin(theta2)
         csphi2=cos(phi2)
         snphi2=sin(phi2)
*
         sdots = snth1*csphi1*snth2*csphi2+snth1*snphi1*snth2*snphi2
     >        +csth1*csth2      ! spin dot product    
*
* densities
*
         d1 = 1.0
         d2 = 1.0
*
* replace R limits with q (scaled R) limits (importance sampling)
*
         if (sflag .gt. 0) then 
            bounds2(1)=exp(-bounds2(1)/range)
            bounds2(2)=exp(-bounds2(2)/range)
         end if
* 
* output files 
*
         open(9,file='big_G_W_out.txt',status='unknown')
         open(10,file='big_G_W_phase.txt',status='unknown')
         open(11,file='big_G_inside_points.txt',status='unknown')
         open(12,file='big_G_geometry.txt',status='unknown')
*
* output geometry
*
         write(12,35) shape1(1),shape1(2),shape1(3),shape2(1),shape2(2),
     >        shape2(3),theta1,phi1,theta2,phi2,range,gap,amp

 35      format(e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,
     >        e13.5,e13.5,e13.5,e13.5)
*                
* compute volume
*
         volume=(bounds1(1)-bounds1(2))*(bounds1(3)-bounds1(4))*
     >        (bounds1(5)-bounds1(6))*(bounds2(1)-bounds2(2))*
     >        (bounds2(3)-bounds2(4))*(bounds2(5)-bounds2(6))
*
* loop over source mass phase values
*
         do 90 j=1,phpnts
            print *,'start phase loop',j
            omega=6.2832*(j-1)/phpnts
            asinom=amp*sin(omega)
            Stats(1,j)=omega
*
* loop over games
*
            do 80 k=1,games
               points=0
               W=0.0
*
* loop over throws
*
               do 50 i=1,throws	
                  call dice(point1,point2,bounds1,bounds2) ! get a point
*     
                  snth=sin(point2(2))
                  csth=cos(point2(2))
                  snphi=sin(point2(3))
                  csphi=cos(point2(3))
                  R=point2(1)
*
                  call Plates(shape1,shape2,D,point1,point2,asinom,
     >                 R,snth,csth,snphi,csphi,inside,d1,d2,pflag) ! is the point inside?
*
                  if (inside) then
                     points=points+1 ! bump points count
*
                     scaler = 1
                     dot = -modez(point1,shape1)*snth*csth
*
                     if (u .eq. 1) then
                        write(11,45) point1(1),
     >                       point1(2),point1(3),point2(1),point2(2)
     >                       ,point2(3),j,d1,d2
 45                     format(e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,i4
     >                       ,e13.5,e13.5)
                     endif
                     W=W+d1*d2*scaler*dot
                  end if	  
 50            continue
*
* normalize and record results
*
               W=W*volume/throws
               Results(1,j,k)=W
               Results(2,j,k)=points
*
* report progress
*
*     		write(9,65) j,k
*65    		format (i4,i4)
*	  
 80         continue            ! close loop over games		
*
 90      continue               ! close loop over phase points
*
* compute means
*
         do 85 j=1,phpnts
            Stats(2,j)=0
            do 82 k=1,games
               Stats(2,j)=Stats(2,j)+Results(1,j,k)	
 82         continue
            Stats(2,j)=Stats(2,j)/games	
 85      continue
*
* standard deviation of mean
*
         do 95 j=1,phpnts
            Stats(3,j)=0
            do 92 k=1,games
               Stats(3,j)=Stats(3,j)+(Results(1,j,k)-Stats(2,j))**2  
 92         continue
            Stats(3,j)=(Stats(3,j)**0.5)/games	
 95      continue
*
* write Results
*
*      write(9,125)
*125   format(/,'phase  game   W   points',/)
*
*      do 135 j=1,phpnts
*      	do 136 k=1,games
*      	write(9,100) j,k,Results(1,j,k),Results(2,j,k)
*100    format(i5,i5,e12.5,e12.5)	  
*136 	continue		
*135   continue
*
* write Stats
*
*      write(9,225)
         write(10,225)
 225     format('  omega        Wmean        Wdev')
         do 235 j=1,phpnts
*     	write(9,200) Stats(1,j),Stats(2,j),Stats(3,j)
            write(10,200) Stats(1,j),Stats(2,j),Stats(3,j)
 200        format(e13.5,e13.5,e13.5)	  		
 235     continue
*
* compute amplitude of virtual work and write
*
         Wampl=0.0
         do 245 j=1,phpnts
            omega=6.2832*(j-1)/phpnts
            Wampl=Wampl+sin(omega)*Stats(2,j) ! Fourier amplitude		
 245     continue
         Wampl=Wampl*2.0/phpnts ! scale Fourier amplitude
*      write(9,255)
*255   format(/,'Wampl=',/)
         write(9,260) Wampl
 260     format(e12.5)
*
         print *, "done"
*
 999  continue                  !close loop over configs
      end
*
* functions and subroutines
*
*** modez ***
*
*  mode shape function
*
      function modez(point,shape)
      real*4 modez,point(3),shape(36)
      modez = shape(10) + shape(11)*point(1) + shape(12)*point(2)
     >     + shape(13)*point(1)*point(2)
     >     + shape(14)*point(1)**2 + shape(15)*point(2)**2
      end
*	  
*** curve ***
*
*  curvature function
*
      function curve(point,shape)
      real*4 curve,point(3),shape(36)
      curve = shape(4) + shape(5)*point(1) + shape(6)*point(2)
     >     + shape(7)*point(1)*point(2)
     >     + shape(8)*point(1)**2 + shape(9)*point(2)**2
      end
*
*** Plates.f ***
*
*  returns inside=.true. iff point1 is inside mass 1 according to shape1
*  and point2 is inside mass 2 according to shape2.  This version is for
*  our planar torsional oscillator.
*
      subroutine Plates(shape1,shape2,D,point1,point2,asinom,
     >     R,snth,csth,snphi,csphi,inside,d1,d2,pflag)
      implicit none
*
* variables passed
*
      real*4 shape1(36)         ! shape of mass 1 (detector) in its cartesian
                                ! includes detector mode shape
      real*4 shape2(36)         ! shape of mass 2 (source) in its cartesians
                                ! includes source mode shape
      real*4 D(3)               ! source position relative to detector
	  		        ! in detector cartesians
      real*4 point1(3)          ! point in box 1; x, y, z
      real*4 point2(3)          ! point in box 2; scaled R, theta, phi
      real*4 asinom             ! amp*sin(omega)
      real*4 R                  ! radial coord, mass 2
      real*4 snth, csth         ! sin(theta), cos(theta)
      real*4 snphi, csphi       ! sin(phi), cos(phi)

      logical*4 inside          ! .true. if both points are inside masses
*
      real*4 d1, d2             ! test mass relative density
      integer*4 pflag
*
* variables
*	  
      real*4  modez             ! mode shape function
      real*4  curve             ! curvature function
      real*4  sdelta            ! source deflection
      real*4  ddelta            ! detector deflection
      real*4  pn2(3)            ! cartesian coords of point 2 in mass 2 cartesians
*
* compute point2 in source mass cartesians
*
      pn2(1)=point1(1)+R*snth*csphi-D(1)
      pn2(2)=point1(2)+R*snth*snphi-D(2)
      pn2(3)=point1(3)+R*csth-D(3)
*
*  deflections
*
      sdelta=asinom*modez(pn2,shape2)+curve(pn2,shape2)
      ddelta=curve(point1,shape1)
*
* return if outside
*
      inside=.false.
      if(pflag .gt. 0)then
*     print *,'plates called; start plates',sdelta,ddelta
*     
*     source cylinder
*     
         if (pn2(3) .ge. -shape2(3)/2.0-shape2(17)+sdelta .and.
     >        pn2(3) .le. -shape2(3)/2.0+sdelta .and.
     >        sqrt((pn2(1)-shape2(18))**2+(pn2(2)-shape2(19))**2) .le.
     >        shape2(16)) then
            d2 = 8.3913
*     d2 = 1.0
            
         else if (pn2(3) .le. shape2(3)/2.0+sdelta .and. ! source z max
     >           pn2(3) .ge. -(shape2(3)/2.0)+sdelta .and. ! source z min
     >           pn2(1) .le. shape2(1) .and. ! source x max
     >           pn2(1) .ge. 0.0 .and. ! source x min
     >           pn2(2) .le. shape2(2)/2.0 .and. ! source y max
     >           pn2(2) .ge. -shape2(2)/2.0) then ! source y min
            d2 = 1.0
            
         else
            return
         end if
*     
*     detector cylinder
*     
         if (point1(3) .le. (shape1(3)/2.0)+shape1(19)+ddelta .and.
     >        point1(3) .ge. shape1(3)/2.0+ddelta .and.
     >        sqrt((point1(1)-shape1(20))**2+(point1(2)-shape1(21))**2)
     >        .le. shape1(18)) then
            d1 = 8.3913
*     d1 = 1.0
*     
*     detector counterweight
*     
         else if (point1(3) .ge. (-shape1(3)/2.0)-shape1(23)+ddelta
     >           .and. point1(3) .le. -shape1(3)/2.0+ddelta .and.
     >           sqrt((point1(1)-shape1(24))**2+(point1(2)
     >           -shape1(25))**2) .le. shape1(22)) then
            d1 = 8.3913
*     d1 = 1.0
*     
         else if (point1(3) .le. (shape1(3)/2.0)+ddelta .and. ! detector z max
     >           point1(3) .ge. -(shape1(3)/2.0)+ddelta .and. ! detector z min
     >           point1(1) .le. 0.0 .and. ! detector x max
     >           point1(1) .ge. -shape1(1) .and. ! detector x min
     >           point1(2) .le. shape1(2)/2.0 .and. ! detector y max
     >           point1(2) .ge. -shape1(2)/2.0) then ! detector y min
            d1 = 1.0
            
         else
            return
         end if
*
      else
         if (pn2(3) .le. shape2(3)/2.0+sdelta .and. ! source z max
     >        pn2(3) .ge. -(shape2(3)/2.0)+sdelta .and. ! source z min
     >        pn2(1) .le. shape2(1) .and. ! source x max
     >        pn2(1) .ge. 0.0 .and. ! source x min
     >        pn2(2) .le. shape2(2)/2.0 .and. ! source y max
     >        pn2(2) .ge. -shape2(2)/2.0) then ! source y min
            d2 = 1.0
            
         else
            return
         end if
         if (point1(3) .le. (shape1(3)/2.0)+ddelta .and. ! detector z max
     >        point1(3) .ge. -(shape1(3)/2.0)+ddelta .and. ! detector z min
     >        point1(1) .le. 0.0 .and. ! detector x max
     >        point1(1) .ge. -shape1(1) .and. ! detector x min
     >        point1(2) .le. shape1(2)/2.0 .and. ! detector y max
     >        point1(2) .ge. -shape1(2)/2.0) then ! detector y min
            d1 = 1.0
            
         else
            return
         end if
      end if
*         
      inside=.true.             ! both points inside
*     
      return
*     
      end
*	
*** dice.f ***
*
* returns point1 inside the box defined by bounds1 and point2 inside
* the box defined by bounds2.
*
      subroutine dice(point1,point2,bounds1,bounds2)
      implicit none
*
* variables passed
*
      real*4 point1(3)          ! point in box 1: x, y, z
      real*4 point2(3)          ! point in box 2: R, theta, phi
      real*4 bounds1(6)         ! box around mass 1, order is:
	         		!    x1 upper, x1 lower, y1 upper, y1 lower
		        	!    z1 upper, z1 lower
      real*4 bounds2(6)         ! box around polar coordinates, order is:
			        !    R upper, R lower, theta upper, theta lower
			        !    phi upper, phi lower
*
* variables
*
      integer*4 randi           ! function for getting random integers
      integer*4 rand1,rand2     ! scratch random numbers	
      integer*4 i               ! index
*
* throw the dice
*
      do 20 i=1,3
*      print *,'dice called; start dice loop' 
         rand1=randi()
         rand2=randi()
*     	print *, rand1,rand2
         point1(i)=((bounds1(2*i-1)+bounds1(2*i))+
     >        (bounds1(2*i-1)-bounds1(2*i))*(rand1/32767.0))/2.0
         point2(i)=((bounds2(2*i-1)+bounds2(2*i))+
     >        (bounds2(2*i-1)-bounds2(2*i))*(rand2/32767.0))/2.0
 20   continue
      return
      end
*	  
***  rand.f  ***
*
*  function returns a random integer [-32768,+32767]
*  derived from Absoft example random.f
*
      FUNCTION randi()  
      implicit none
      REAL*4 rando, rand
      integer*4 randi
*
      rando=rand(0)
      randi = -32768 + INT(65536 * rando)
*
      END
*
