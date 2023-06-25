c the large-scale field configuration of ADAF

 
	 parameter(n=200,np=200,n_b=100,n_z=150,n_zh=20)
	 implicit double precision(a-g,o-z)

	 dimension p_ijd(np+1,np+1),ad(np+1,np+1)
	 dimension pd_ijd(np+1,np+1)
	 dimension pminus_ijd(np+1,np+1)
	 dimension r_id(n+1),r_ibd(n_b+1),bd(n+1),rb_infd(n+1),ajsd(n+1)


	 dimension r_xd(n+1),theta_d(n+1)

	 dimension q_ijd(np+1,np+1),qminus_ijd(np+1,np+1)
	 dimension rvarphi_d(n+1),rvarphi_dd(n+1),rvarphi_infd(n+1)
	 dimension rvarphif_d(n+1),rvarphif_dd(n+1)
	 dimension rvarphix_d(n+1)

	 dimension b_s_d(np+1,n_z+1),z_d(np+1,n_z+1)

	 dimension sigma_d(n+1),vr_d(n+1),ah_d(n+1),akappa_d(n+1)
	 dimension bz_d(n+1),omega_d(n+1),zs_d(n+1)
	 dimension r_ad(n+1),t_m_d(n+1)
	 dimension vr_fd(n+1)
	
       dimension rvar_arr_d(n+1,n+1)


c magnetic field strengths at the disk surface 

	 dimension bz_sd(n_b+1),br_sd(n_b+1),bz_0d(n_b+1)



c vr_d: radial velocity
c sigma_d: surface density (=1, at r_out)
c ah_d: H/R
c akappa_d: B_z/B_r at z=H
c bz_d: z-component magnetic field strength at z=H
c omega_d: \Omega
c zs_d: z_s/H, z_s: sonic point
c r_ad: Alfven radius along the field line threading the disk at r_i
c t_m_d: magnetic torque: T_m/\Sigma R\Omega_K^2



c b_s_d(i,j): magnetic field shape
c b_s_d(i,j)=r_i at z=z_j
c i: r_id(i), j: z_j, j=1,z=0.0


	 integer indx(n+1)

	 pi=2.0*asin(1.0)
	 eta_i=tanh(1.0)

	 phi=0.5*pi

c err_b: error control for B-field shape

	 err_b=1.0d-6

c err_alf: err control for the Alfven point

	 err_alf=5.0d-3

c  q_d(n+1,n+1): q_ij
c  qminus_d(n+1,n+1): q_ij^-1
c  rvarphi_dd(n+1): r\varphi_d(potential of the disk)
c  rvarphi_infd(n+1): r\varphi_inf(potential of the current at infinity)

c disk parameters:
c theta: c_s/R\Omega



c	 theta=0.4
c       theta=0.09


c viscosity: alpha
c alpha=0.1 case, stronger fields obtained


       alpha=0.1
c       alpha=0.2

c eps_mad: the radial velocity of a MAD is: eps_mad*v_k

c	 eps_mad=1.0e-3
	 eps_mad=1.0e-2
	 

c ah: the disk half-thickness	H/R

c	    ah=theta
c	    c_1=0.5
c          c_1=0.16

c	 ah=1.0*theta

c p_m: eta/(0.1*c_s H)
c Prandtle number

c        p_m=1.0
       p_m=0.5
c       p_m=0.2

c         p_m=1.5

c contribution of the outflows on the angular momentum transfer in the disk
c f_m=0.0: no magnetic torque

c       f_m=2.0
c       f_m=1.0
c       f_m=0.0
        f_m=9.232002E-01


c al_in: the angular momentum at r_in, in units of Keplerian value

c	  al_in=0.5

c	  delta_q=1.0
	  delta_q=0.5
c	  delta_q=0.3
c	  delta_q=0.2
c	  delta_q=0.1
c	  delta_q=0.05

	    
	   alambda=0.5
c	   alambda=0.1
c	   alambda=0.25


c bz0: dimensionless magnetic field strength at t=0,
c bz0=1.0: equipartition with gas pressure
c intial field lines are assumed to vertically thread the disk 



c *x.dat
  
       bz0=0.5


c 19 \le i_try \le 188


c disk extended from r_in to r_out
c the z-range of the magnetic field, z=0 to z_max         

	 r_in=4.0
c	 r_out=1.0d5
c	 r_out=1000.0
c	 r_out=180.0
c t=t_1
c	 r_out=r_in*102.2
c t=t_2
	 r_out=r_in*32.831490
c	 r_out=r_in*6.5


c r_mad: size of MAD

c        r_mad=r_in*1.18	
c         r_mad=r_in*1.16	

        r_mad=r_in*1.1701
c         r_mad=0.0

c only consider the winds from r_d_in to r_d_out

c	 r_d_in=1.5
c	 r_d_out=75.0

c	 r_d_in=1.0
c	 r_d_out=101.0

c	 r_min=1.0

c	 r_min=1.35

	 r_max=r_out


	 z_min=0.0
	 z_max=200.0

c the magetic flux in the region of r<r_mad, or r<r_isco if no mad state

       phi_b=0.0 


c *****w1p0.dat    total torque/viscous torque=1.0  i.e., no  magnetic winds
c ********b.dat    p_m=0.5   	 eps_mad=1.0e-3
c ********c.dat    p_m=0.5     eps_mad=1.0e-2
c ********a.dat    p_m=1.0     eps_mad=1.0e-2



 	 open(1,file='mjs_pm_br32p8rm1p1701w0p92c.dat',status='unknown')
c 	 open(2,file='rvarphi_r_pm1x.dat',status='unknown')
c 	 open(3,file='rvarphi_z_pm1x.dat',status='unknown')
c 	 open(4,file='rvarphi_var_pm1x.dat',status='unknown')
c 	 open(5,file='disk_h.dat',status='unknown')

c 	 open(10,file='disk_h_pm1x.dat',status='unknown')


	  open(11,file='ma_brz_br32p8rm1p1701w0p92c.dat',status='unknown')
c	  open(12,file='adaf_ah2_aa.dat',status='unknown')

	  open(13,file='mad_bz_br32p8rm1p1701w0p92c.dat',status='unknown')
	  open(16,file='ma_phi_br32p8rm1p1701w0p92c.dat',status='unknown')
	  open(15,file='mi_brz_br32p8rm1p1701w0p92c.dat',status='unknown')



c 	 open(21,file='adaf_vr_b.dat',status='old')
c 	 open(22,file='adaf_theta_b.dat',status='old')
c 	 open(23,file='adaf_ah_b.dat',status='old')
	 


  

c n_int: steps in the integral over phi

c	 n_int=1000
	 n_int=500

	 dlnr=(log(r_out)-log(r_in))/n

	   do i=1,n+1

c	    t_m_d(i)=0.0

	    r_id(i)=exp(log(r_in)+(i-1.0)*dlnr)
c	    read(21,*)r_id(i),vr_d(i)

c test self-similar model
c	    vr_d(i)=-0.105
		
	    rb_infd(i)=-r_id(i)*bz0

	    rvarphi_infd(i)=0.5*r_id(i)**2*bz0
	    rvarphi_dd(i)=0.0

	    rvarphi_d(i)=rvarphi_dd(i)+rvarphi_infd(i)

c	     read(22,*)ry,theta_d(i)

	 

	
c	    read(23,*)rx,ah_d(i)

c	    write(5,*)rx,ah_d(i)*rx


c 		ah_d(i)=ah

c v_r: radial velocity (/R\Omega)

c	    vr_d(i)=-c_1*alpha

c	    vr_d(i)=-alpha*theta**2
c	1	/(sqrt(1.0-theta**2)-al_in/sqrt(r_id(i)))

c	    vr_d(i)=-alpha*theta**2
c	1	/(sqrt(1.0-2.5*theta**2)-al_in/sqrt(r_id(i)))

c test self-similar model
c	    vr_d(i)=-0.10465
c	    ah_d(i)=0.5906
c          theta_d(i)=ah_d(i)

c for thin accretion disk

c	    ah_d(i)=0.1

c for an ADAF

c	    ah_d(i)=0.7
c	    ah_d(i)=0.8

	    ah_d(i)=0.4
          theta_d(i)=ah_d(i)
		 
c	    vr_d(i)=-1.0*(3.0/2.0)*alpha*ah_d(i)**2
c	    vr_d(i)=-2.0*(3.0/2.0)*alpha*ah_d(i)**2
c	    vr_d(i)=-5.0*(3.0/2.0)*alpha*ah_d(i)**2

	    if(r_id(i).gt.r_mad)then
	    vr_d(i)=-(1.0+f_m)*(3.0/2.0)*alpha*theta_d(i)*2.0*ah_d(i)
	    else
	    vr_d(i)=-eps_mad
	    endif

	  


	write(*,*)r_id(i),vr_d(i)
	
	   enddo 

	   r_out=r_id(n+1)
	   r_in=r_id(1)

	   dlnr=(log(r_out)-log(r_in))/n

	   r_max=r_out
	   r_min=r_in

c  changed on 2022-04-18
c         r_min=r_in*0.5     


	   vr_fd=vr_d

2	  continue

	

c+++++++++++++++++++++++++++++++++++++++++++++++++

	  do i=1,np+1
	     do j=1,np+1  

	     r_i=r_id(i)
	     ah=ah_d(i)
	     r_j=r_id(j)
	     z_i=0.0

	       
	      if(i.ne.j)then

c
		  call subp_ij(n_int,n_zh,r_i,z_i,ah,r_j,p_ij)

	      p_ijd(i,j)=p_ij
	      pd_ijd(i,j)=p_ij

	      elseif(i.eq.j)then

	      r_j=exp(log(r_id(j))+alambda*dlnr)

	      call subp_ij(n_int,n_zh,r_i,z_i,ah,r_j,p_ij1)

	      r_j=exp(log(r_id(j))-alambda*dlnr)

	      call subp_ij(n_int,n_zh,r_i,z_i,ah,r_j,p_ij2)

	      p_ijd(i,j)=0.5*(p_ij1+p_ij2)
	1	  +2.0*pi*(p_m*alpha)*theta_d(i)*r_i/(vr_d(i)*dlnr)
	      pd_ijd(i,j)=0.5*(p_ij1+p_ij2)
	      endif

c	   write(*,*)'i   j   p_ij',i,j,p_ijd(i,j)


	   enddo
	  enddo
c+++++++++++++++++++++++++++++++++++++++++++++++++
	
c         pause

	   ad=p_ijd
	   bd=rb_infd


	   call ludcmp_dp(ad,n+1,np+1,indx,d)

         pminus_ijd=ad

c	   pause

	   ad=pminus_ijd

	   call lubksb_dp(ad,n+1,np+1,indx,bd)
	
         ajsd=bd/dlnr

	do i=1,n+1
	write(*,*)i,ajsd(i)
	write(1,*)r_id(i),ajsd(i)
c	write(10,*)r_id(i),ah*r_id(i)

	enddo


c	   pause

c++++++++++++++++++++++++++++++++++++++
c check the solution!

	 do i=1,n+1

	 x=0.0 
	 do j=1,n+1
	 	x=x+p_ijd(i,j)*ajsd(j)*dlnr


	 enddo

	 write(*,*)ajsd(i),abs(x-rb_infd(i))

	enddo
	  
	pause
c++++++++++++++++++++++++++++++++++++++++++

c to calculate the magnetic field strengths at the disk surface
c z=H

c	  delta_r=1.0d-2
c	  delta_z=1.0d-2


	 dlnrb=(log(r_max)-log(r_min))/n_b
	  do i=1,n_b+1
	   r_ibd(i)=exp(log(r_min)+(i-1.0)*dlnrb)
	  enddo


	  do i=1,n_b+1

	    ah=1.0*ah_d(2*(i-1)+1)
c	    ah=1.0*ah_d(4*(i-1)+1)

c	    ah=ah_d(i)

	    r_i=r_ibd(i)


c	   do j=1,n
	    
c	    if(i.eq.1.or.i.eq.(n_b+1))then
           
c c      
c		elseif(r_i.gt.r_id(j).and.r_i.le.r_id(j+1))then
c	     sl_ah=(ah_d(j+1)-ah_d(j))/(log(r_id(j+1))-log (r_id(j)))
c	      ahx=ah_d(j)+sl_ah*(log(r_i)-log(r_id(j)))
c	    elseif(r_i.lt.r_id(1))then
c	      ahx=ah_d(1)
c	    endif
 
c	   enddo

	    z_i=ah*r_ibd(i)


	 call subbrz(n,n_zh,r_i,z_i,ah_d,r_id,ajsd,
	1 bz0,b_r,b_z)

		br_sd(i)=b_r
		bz_sd(i)=b_z   


c+++++++++++++++++++++++++++++++++++++
c b_z at z=0


	    z_i=0.0

	 call subbrz(n,n_zh,r_i,z_i,ah_d,r_id,ajsd,
	1 bz0,b_r,b_z)

		bz_0d(i)=b_z



c++++++++++++++++++++++++++++++++++++++++




	  write(*,*)'br_sd      bz_sd',r_id(i),br_sd(i),bz_sd(i)
	  write(*,*)r_ibd(i),br_sd(i),bz_sd(i)
c	write(12,*)r_i,ah

c	pause

	  write(11,*)r_ibd(i),br_sd(i),bz_sd(i)	  
	  write(13,*)r_ibd(i),bz_0d(i),b_r

c	  pause
	  enddo



c+++++++++++++++++++++++++++++++++++++


	 do i=1,n_b+1

	 if(i.ge.2.and.i.le.n_b)then
	
	 xi_bz=-(log(bz_sd(i+1))-log(bz_sd(i-1)))
	1 /(log(r_ibd(i+1))-log(r_ibd(i-1)))

	 xi_br=-(log(br_sd(i+1))-log(br_sd(i-1)))
	1 /(log(r_ibd(i+1))-log(r_ibd(i-1)))
	 
	 elseif(i.eq.1)then
	 
	 
	 xi_bz=-(log(bz_sd(i+1))-log(bz_sd(i)))
	1 /(log(r_ibd(i+1))-log(r_ibd(i)))


	 xi_br=-(log(br_sd(i+1))-log(br_sd(i)))
	1 /(log(r_ibd(i+1))-log(r_ibd(i)))

	 elseif(i.eq.(n_b+1))then

	 xi_bz=-(log(bz_sd(i))-log(bz_sd(i-1)))
	1 /(log(r_ibd(i))-log(r_ibd(i-1)))

	 xi_br=-(log(br_sd(i))-log(br_sd(i-1)))
	1 /(log(r_ibd(i))-log(r_ibd(i-1)))

	 endif

	 write(*,*)r_ibd(i),xi_br,xi_bz
	 write(15,*)r_ibd(i),xi_br,xi_bz

	 enddo


c++++++++++++++++++++++++++++++++++++++++


c to calculate the array of r\varphi(r_i,z_i)


c	 dr_i=(r_max-r_min)/n
c	 dz_i=(z_max-z_min)/n



c	  do i_r=1,n
c	   r_i=r_min+(i_r-0.5)*dr_i

c	    do i_z=1,n
c	    z_i=z_min+(i_z-0.5)*dz_i


c	 call subrvarphi_drz(n,n_zh,r_i,z_i,ah,dlnr,r_id,ajsd,
c	1 rvarphi_drz)

c	    write(2,*)r_i
c	    write(3,*)z_i


c 	    write(*,*)'r  z   rvarphi',r_i,z_i,rvarphi_drz+0.5*r_i**2*bz0


c	     write(4,*)rvarphi_drz+0.5*r_i**2*bz0

c	    if(i_z.eq.1)then
c	    rvarphix_d(i_r)=rvarphi_drz+0.5*r_i**2*bz0
c	    r_xd(i_r)=r_i
c	    endif


c	    enddo
	   
c	   write(*,*)'r_i',r_i


c	   enddo




	phi_b=pi*r_in**2*bz_0d(1)

	  write(*,*)'phi_b=',phi_b

	 

	do i=1,n_b+1
	  if(r_mad.gt.r_in.and.r_ibd(i).le.r_mad)then
	  phi_b=phi_b+2.0*pi*r_ibd(i)*bz_0d(i)*r_ibd(i)*dlnr
	  write(*,*)r_mad,r_in,r_ibd(i)
	  write(*,*)i,r_ibd(i),phi_b
c	  pause
	  endif
	enddo



c      pause

c	  do i_r=1,n
c	   write(*,*)r_xd(i_r),rvarphix_d(i_r)
c	  enddo

	  write(16,*)'alpha=',alpha
	  write(16,*)'eps_mad=',eps_mad
	  write(16,*)'p_m=',p_m
	  write(16,*)'f_m=',f_m
	  write(16,*)'r_in=',r_in
	  write(16,*)'r_tr/r_in=',r_out/r_in
	  write(16,*)'r_mad/r_in',r_mad/r_in
	  write(16,*)'phi_b=',phi_b

	  write(*,*)'phi_b=',phi_b


	  pause





	  close(1)
	  close(2)
	  close(3)
	  close(4)




	      stop
	      end



c calculating the sonic point z_s





	 subroutine subz_s(akappa,r_i,ah,omega,theta,z_s,dphi_eff)
	 implicit double precision(a-g,o-z)

  	 eta_i=tanh(1.0)

       z_min=0.0

	  z=z_min

	  i=0

10	  continue

	  i=i+1

	  r=r_i*ah*sqrt(1.0-eta_i**2+eta_i**2*z**2)
	1  /(akappa*eta_i**2)
	1  -r_i*ah*sqrt(1.0-eta_i**2)/(akappa*eta_i**2)
	1  +r_i

	  drdz=r_i*ah*z/(akappa*sqrt(1.0-eta_i**2+eta_i**2*z**2))

	  f_slope=drdz/(r**2*sqrt(1.0+(z*ah)**2))
	1  +z*ah**2/(r*(1.0+(z*ah)**2)**1.5)
	1  -omega**2*r*drdz/r_i**3

c	  write(*,*)'z   f_slope',z,f_slope
c	  pause


	  
	   if(i.gt.2.and.f_slope0*f_slope.le.0.0)then

	    z_max=z
	    z_min=z0
	    f_slope_ma=f_slope
	    f_slope_mi=f_slope0

	    goto 50

	   endif
	    

	  
	  if(i.eq.1)then
	  z0=z
	  z=z+0.01
	  else
	  f_slope0=f_slope
	  z0=z
	  z=z*1.5
	  endif


	  goto 10

50	  continue

	  z=0.5*(z_min+z_max)

	  r=r_i*ah*sqrt(1.0-eta_i**2+eta_i**2*z**2)
	1  /(akappa*eta_i**2)
	1  -r_i*ah*sqrt(1.0-eta_i**2)/(akappa*eta_i**2)
	1  +r_i

	  drdz=r_i*ah*z/(akappa*sqrt(1.0-eta_i**2+eta_i**2*z**2))

	  f_slope=drdz/(r**2*sqrt(1.0+(z*ah)**2))
	1  +z*ah**2/(r*(1.0+(z*ah)**2)**1.5)
	1  -omega**2*r*drdz/r_i**3


	  if(abs(z_max-z_min).lt.1.0d-3.or.abs(f_slope).lt.1.0d-8)then
	  z_s=z
	  goto 100
	  endif


c	  write(*,*)'z_min   z_max',z_min,z_max
c	  write(*,*)'slope',f_slope_mi,f_slope_ma

c	  pause

	  if(f_slope*f_slope_mi.gt.0.0)then
	   z_min=z
	   f_slope_mi=f_slope
	   elseif(f_slope*f_slope_ma.gt.0.0)then
	   z_max=z
	   f_slope_ma=f_slope
	  endif

	  goto 50

100     continue

	  z=0

	  r=r_i*ah*sqrt(1.0-eta_i**2+eta_i**2*z**2)
	1  /(akappa*eta_i**2)
	1  -r_i*ah*sqrt(1.0-eta_i**2)/(akappa*eta_i**2)
	1  +r_i
	  pot_0=-1.0/(r*sqrt(1.0+(z*ah)**2))
	1  -0.5*omega**2*r**2/r_i**3

	  z=z_s

	  r=r_i*ah*sqrt(1.0-eta_i**2+eta_i**2*z**2)
	1  /(akappa*eta_i**2)
	1  -r_i*ah*sqrt(1.0-eta_i**2)/(akappa*eta_i**2)
	1  +r_i
	  pot_s=-1.0/(r*sqrt(1.0+(z*ah)**2))
	1  -0.5*omega**2*r**2/r_i**3


	  dphi_eff=r_i*(pot_s-pot_0)/theta**2

	  return
	  end



c++++++++++++++++++++++++++++
c calculate B_r(r,z) and B_z(r,z) of the disk for given ajsd (currents)


	 
	 subroutine subbrz(n,n_zh,r_i,z_i,ah_d,r_id,ajsd,
	1 bz0,b_r,b_z)
	 implicit double precision(a-g,o-z)
	 
	 dimension r_id(n+1),ajsd(n+1),ah_d(n+1)

	 pi=2.0*asin(1.0)

c	 phi=0.5*pi

	  n_phi=180

	  dphi=pi/n_phi
	 	
c	 alambda=0.5
	 alambda=0.0

	 b_r=0.0
	 b_z=0.0

	 dlnr=log(r_id(2))-log(r_id(1))
	  
	 do j=1,n+1


	     r_j=r_id(j)
	     ah=ah_d(j)

	   dz=2.0*ah*r_j/n_zh

	   do i_z=1,n_zh

	   z_h=-ah*r_j+(i_z-0.5)*dz


		   do j_phi=1,n_phi

               phi=(j_phi-0.5)*dphi


	     if(abs(r_i/r_j-1.0).ge.1.0d-8)then


	     f_b=r_i**2+r_j**2+(z_i-z_h)**2
	1    -2.0*r_i*r_j*cos(phi)
	    
	      b_z=b_z
	1    +r_j**2*dlnr*cos(phi)*dphi*ajsd(j)*dz
	1    /(f_b**0.5*ah*r_j*r_i)
	1    +r_j**2*dlnr*cos(phi)*dphi*ajsd(j)*(r_j*cos(phi)-r_i)*dz
	1    /(f_b**1.5*ah*r_j)
	
	 
	      b_r=b_r
	1    +r_j**2*dlnr*cos(phi)*dphi*ajsd(j)*(z_i-z_h)*dz
	1    /(f_b**1.5*ah*r_j)



           elseif(abs(r_i/r_j-1.0).lt.1.0d-8)then
		
c		elseif(sqrt((r_i-r_j)**2+(z_i-z_h)**2).lt.0.0)then
	        
c		elseif(sqrt((r_i-r_j)**2+(z_i-z_h)**2).lt.10.0*dz)then
	 	
	        r_j=exp(log(r_id(j))+alambda*dlnr)

	     f_b=r_i**2+r_j**2+(z_i-z_h)**2
	1    -2.0*r_i*r_j*cos(phi)
	    
	      db_z1=r_j**2*dlnr*cos(phi)*dphi*ajsd(j)*dz
	1    /(f_b**0.5*ah*r_j*r_i)
	1    +r_j**2*dlnr*cos(phi)*dphi*ajsd(j)*(r_j*cos(phi)-r_i)*dz
	1    /(f_b**1.5*ah*r_j)
	
	 
	      db_r1=r_j**2*dlnr*cos(phi)*dphi*ajsd(j)*(z_i-z_h)*dz
	1    /(f_b**1.5*ah*r_j)

	        r_j=exp(log(r_id(j))-alambda*dlnr)


	     f_b=r_i**2+r_j**2+(z_i-z_h)**2
	1    -2.0*r_i*r_j*cos(phi)
	    
	      db_z2=r_j**2*dlnr*cos(phi)*dphi*ajsd(j)*dz
	1    /(f_b**0.5*ah*r_j*r_i)
	1    +r_j**2*dlnr*cos(phi)*dphi*ajsd(j)*(r_j*cos(phi)-r_i)*dz
	1    /(f_b**1.5*ah*r_j)
	
	 
	      db_r2=r_j**2*dlnr*cos(phi)*dphi*ajsd(j)*(z_i-z_h)*dz
	1    /(f_b**1.5*ah*r_j)

	        b_z=b_z+0.5*(db_z1+db_z2)
	        b_r=b_r+0.5*(db_r1+db_r2)

	     endif

	 enddo


	 enddo


	 enddo


	 b_z=b_z+bz0

	 return
	 end







c++++++++++++++++++++++++++++++++++++++++++++++
      



c++++++++++++++++++++++++++++
c calculate rvarphi(r,z) of the disk for given ajsd (currents)


	 
	 subroutine subrvarphi_drz(n,n_zh,r_i,z_i,ah_d,dlnr,r_id,ajsd,
	1 rvarphi_drz)
	 implicit double precision(a-g,o-z)
	 
	 dimension r_id(n+1),ajsd(n+1),ah_d(n+1)

	 pi=2.0*asin(1.0)

	 phi=0.5*pi
	 	
	 alambda=0.5

	 rvarphi_drz=0.0
	  
	   do j=1,n+1


	     r_j=r_id(j)
	     ah=ah_d(j)


	   q_ij=0.0
	   dz=2.0*ah*r_j/n_zh

	   do i_z=1,n_zh

	   z_h=-ah*r_j+(i_z-0.5)*dz


	   ak=sqrt((4.0*r_i*r_j)/(r_i**2+r_j**2+(z_i-z_h)**2+2.0*r_i*r_j))


	     if(abs(r_i/r_j-1.0).ge.1.0d-8)then

c		if(sqrt((r_i-r_j)**2+(z_i-z_h)**2).ge.1.0*dz)then

	     fun_K=ellf_dp(phi,ak)
	     fun_E=elle_dp(phi,ak)

	     
		q_ij=q_ij+
	1	r_j*4.0*(r_i*r_j/sqrt(r_i**2+r_j**2+(z_i-z_h)**2+2.0*r_i*r_j))
	1	*((2.0-ak**2)*fun_K-2.0*fun_E)/ak**2

           elseif(abs(r_i/r_j-1.0).lt.1.0d-8)then
		
c		elseif(sqrt((r_i-r_j)**2+(z_i-z_h)**2).lt.0.0)then
	        
c		elseif(sqrt((r_i-r_j)**2+(z_i-z_h)**2).lt.10.0*dz)then
	 	
	        r_j=exp(log(r_id(j))+alambda*dlnr)
 	        ak=sqrt((4.0*r_i*r_j)
	1	    /(r_i**2+r_j**2+(z_i-z_h)**2+2.0*r_i*r_j))

	        fun_K=ellf_dp(phi,ak)
	        fun_E=elle_dp(phi,ak)

		    q_ij1=r_j*4.0*(r_i*r_j/
	1		sqrt(r_i**2+r_j**2+(z_i-z_h)**2+2.0*r_i*r_j))
	1	    *((2.0-ak**2)*fun_K-2.0*fun_E)/ak**2

	 
	        r_j=exp(log(r_id(j))-alambda*dlnr)
 	        ak=sqrt((4.0*r_i*r_j)
	1	    /(r_i**2+r_j**2+(z_i-z_h)**2+2.0*r_i*r_j))

	        fun_K=ellf_dp(phi,ak)
	        fun_E=elle_dp(phi,ak)

		    q_ij2=r_j*4.0*(r_i*r_j/
	1		sqrt(r_i**2+r_j**2+(z_i-z_h)**2+2.0*r_i*r_j))
	1	    *((2.0-ak**2)*fun_K-2.0*fun_E)/ak**2

	        q_ij=q_ij+0.5*(q_ij1+q_ij2)

	     endif

	 enddo

	  q_ij=q_ij*dz/(2.0*ah*r_j)

	  rvarphi_drz=rvarphi_drz+q_ij*ajsd(j)*dlnr


	 enddo

	 return
	 end







c++++++++++++++++++++++++++++++++++++++++++++++
c calculate the P_ij
	 


	  subroutine subp_ij(n_int,n_z,r_i,z_i,ah,r_j,p_ij)
	  implicit double precision(a-g,o-z)

	  n=n_int
	  

	  pi=2.0*asin(1.0)

	  dphi=2.0*pi/n

	  p_ij=0.0

	   dz=2.0*ah*r_j/n_z

	   do k=1,n

		 do i=1,n_z

		 z_h=-ah*r_j+(i-0.5)*dz

	    phi=(k-0.5)*dphi
	    p_ij=p_ij+(r_j**2+(z_i-z_h)**2-r_i*r_j*cos(phi))*cos(phi)
	1	/sqrt((r_j**2+r_i**2+(z_i-z_h)**2-2.0*r_i*r_j*cos(phi))**3)

		  enddo

	   enddo

	  p_ij=p_ij*r_j**2*dphi*dz/(2.0*ah*r_j)

	  return
	  end

c++++++++++++++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++++++++++++++++++
c calculate the P_ij(z)
	 


	  subroutine subpz_ij(n_int,r_i,r_j,z,p_ij)
	  implicit double precision(a-g,o-z)

	  n=n_int
	  pi=2.0*asin(1.0)

	  dphi=2.0*pi/n

	  p_ij=0.0

	   do k=1,n
	    phi=(k-0.5)*dphi
	    p_ij=p_ij+(r_j**2+z**2-r_i*r_j*cos(phi))*cos(phi)
	1	/sqrt((r_j**2+r_i**2+z**2-2.0*r_i*r_j*cos(phi))**3)
	   enddo

	  p_ij=p_ij*r_j**2*dphi

	  return
	  end

c++++++++++++++++++++++++++++++++++++++++++++++





c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      SUBROUTINE lubksb_dp(a,n,np,indx,b)
	implicit double precision(a-g,o-z)
      INTEGER n,np,indx(n)
      dimension a(np,np),b(n)
      INTEGER i,ii,j,ll
c      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END



      SUBROUTINE ludcmp_dp(a,n,np,indx,d)
	implicit double precision(a-g,o-z)
      INTEGER n,np,indx(n),NMAX
      dimension a(np,np)
      PARAMETER (NMAX=2000,TINY=1.0d-20)
      INTEGER i,imax,j,k
      dimension vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END


 
 
 
 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c elliptal functions (double precision)

c K(ak), if phi=pi/2.0




      FUNCTION ellf_dp(phi,ak)
	   implicit double precision(a-g,o-z)
c      REAL ellf,ak,phi
CU    USES rf
c      REAL s,rf
      s=sin(phi)
      ellf_dp=s*rf_dp(cos(phi)**2,(1.-s*ak)*(1.+s*ak),1.d0)
      return
      END



c E(ak), if phi=pi/2.0



      FUNCTION elle_dp(phi,ak)
       implicit double precision(a-g,o-z)
c      REAL elle,ak,phi
CU    USES rd,rf
c      REAL cc,q,s,rd,rf
      s=sin(phi)
      cc=cos(phi)**2
      q=(1.-s*ak)*(1.+s*ak)
      elle_dp=s*(rf_dp(cc,q,1.d0)-((s*ak)**2)*rd_dp(cc,q,1.d0)/3.)
      return
      END


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++






      FUNCTION rf_dp(x,y,z)
	   implicit double precision(a-g,o-z)
c      REAL rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08,TINY=1.5d-38,BIG=3.E37,THIRD=1./3.,
     *C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
c      REAL alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf_dp=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)

      return
      END




      FUNCTION rd_dp(x,y,z)
	   implicit double precision(a-g,o-z)
c      REAL rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
      PARAMETER (ERRTOL=.05,TINY=1.e-25,BIG=4.5E21,C1=3./14.,C2=1./6.,
     *C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4)
c      REAL alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
c     *sqrtz,sum,xt,yt,zt


      if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rd'
      xt=x
      yt=y
      zt=z
      sum=0.
      fac=1.
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=.25*fac
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=.2*(xt+yt+3.*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.*eb
      ee=ed+ec+ec
      rd_dp=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)+delz
	1*(C2*ee+delz*(-C3*
     *ec+delz*C4*ea)))/(ave*sqrt(ave))
      return
      END



c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


             



