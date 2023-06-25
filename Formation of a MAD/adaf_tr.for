c using X-ray data to calculate the transition radius r_tr
c the mdot(t)\propto exp(-t/tau)

 
	 parameter(n=42)
	 implicit double precision(a-g,o-z)
	 dimension date_d(n),alum_c_d(n),alum_d_d(n),alum_tot_d(n),
	1 amdot_d(n),r_tr_d(n),alum_adaf_d(n),p_w_d(n)
	  real t,arate_d,r_tr,p_w,f_m

	 bhmass=8.5
	 alum_edd=1.2501e+38*bhmass
	 tau=57.61				 
c	 p_w=0.563
	 r_isco=4.0

 	 open(1,file='lum_x.dat',status='unknown')
 	 open(10,file='r_tr.dat',status='unknown')
 	 open(11,file='r_tr2.dat',status='unknown')

  	 open(15,file='lum_adaf.dat',status='unknown')


	 do i=1,42
	  read(1,*)date_d(i),alum_c_d(i),alum_d_d(i),alum_tot_d(i)
	  amdot_d(1)=alum_tot_d(1)/alum_edd
	  amdot_d(i)=amdot_d(1)*exp(-(date_d(i)-date_d(1))/tau)
c r_tr_d: r_tr/r_isco
	  r_tr_d(i)=amdot_d(i)*alum_tot_d(1)/amdot_d(1)/alum_d_d(i)
	  p_w_d(i)=log(alum_c_d(i)/alum_edd/(amdot_d(i)
	1  -alum_d_d(i)/alum_edd))/log(2.25/r_tr_d(i))
	  t=date_d(i)
	  arate_d=amdot_d(i)
	  r_tr=r_tr_d(i)
	  p_w=p_w_d(i)

	  c_f_m=(alum_c_d(i)/alum_edd)
	1  /(amdot_d(i)-alum_d_d(i)/alum_edd)

	  alpha=0.1

	  

	  f_m=1.0/c_f_m-1.0
c	  write(10,*)t,arate_d,r_tr,p_w,f_m
	  write(10,*)t,arate_d,r_tr,p_w,f_m,
	1  arate_d**(6.0/5.0)*(r_tr*4.0)**(1.0/20.0)*alpha**(-21.0/10.0)
	1  *5.61e-2*8.5**(-1.0/10.0),
	1  59.97*alpha**(-2)*arate_d**2*r_tr**(-1),
	1  1.64e-4*alpha**(-1.0/8.0)*8.5**(-1.0/8.0)
     1 *(r_tr_d(i)*4)**(21.0/16.0)

 	  write(11,*)t,arate_d,r_tr
         write(10,*)'f_m     p_w'


	  f_m=1.0
        p_w_2=log((alum_c_d(i)/alum_edd)*(1.0+f_m)
	1  /(amdot_d(i)-(alum_d_d(i)/alum_edd)))
	1  /log(2.25/r_tr_d(i))
         write(10,*)f_m,p_w_2

	  f_m=2.0
        p_w_2=log((alum_c_d(i)/alum_edd)*(1.0+f_m)
	1  /(amdot_d(i)-(alum_d_d(i)/alum_edd)))
	1  /log(2.25/r_tr_d(i))
         write(10,*)f_m,p_w_2

	   p_w_3=0.3
	   f_m=(amdot_d(i)-alum_d_d(i)/alum_edd)
	1  *(2.25/r_tr_d(i))**p_w_3
     1  /(alum_c_d(i)/alum_edd)-1.0

	   write(10,*)f_m,p_w_3

         write(10,*)'     '




c	  write(10,*)date_d(i),amdot_d(i),r_tr_d(i),p_w_d(i)
	  alum_adaf_d(i)=(amdot_d(i)-alum_d_d(i)/alum_edd)
	1  *(2.25/r_tr_d(i))**p_w_d(i)
	  write(15,*)date_d(i),alum_adaf_d(i)
	 enddo

	 pause

	 stop
	 end



