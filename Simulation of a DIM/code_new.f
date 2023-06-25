        program AccretionDisk        
        parameter (n1=4,n2=4,n=5000)
        parameter (ni=5000)
        implicit real*8 (a-h,o-z)
        real*8 m1,m2,mv,mb,mi,mu30
        real*8 M_sun,m_p,L_Edd,MdotEdd
        character*80 filename
        logical init/.true./,rate/.false./,noprint,nograph,coldstart
     1                                                  ,debut_alpha
        dimension xii(ni,8),xv(ni),yv(ni)
        dimension xf(n,n1+n2),xi(n,n1+n2),beta0(n1+n2)
        dimension deltax(n,n1+n2),dmax(n1+n2)
        dimension x2_p(n),x3_p(n),r_p(n)
        dimension dx2(n),dx3(n)
        dimension Reva(100),Evamdot(100)
        dimension DataX(701),DataY(701)
        common/sortie/noprint,nograph
        common/solution_prec/r_p,x2_p,x3_p,dx2,dx3,deltat
        common/parametres/m1,xmp16,alpha_c,alpha_h,rs,rstar,c_tid,p_orb
     1             ,eta_impact, size_impact,r_tid,mu30,fill,width
        common/alf/debut_alpha
        common/trecord/trecord,f_over,g_over
        common/Revap/Reva,Evamdot
        common/DataXray/DataX,DataY,fill_2
        data beta0/8*0.5/
        debut_alpha=.true.

        m1=8.5
        m2=0.6  
        p_orb=16.4518
        c_tid=1e-5
        c_tid=0
        width=10
        if(c_tid.ge.1d-20) width=0
        eta_impact=.0
        size_impact=10.
        mu30=2
        fill=1e-3
        alpha_c=0.02
        alpha_h=0.3
        xmp16_f=5.5
        f_over=0
        g_over=0    
 
        sep=3.53*(m1+m2)**0.3333333*p_orb**0.6666666
        rs=rcirc(m1/m2)*sep
        rstar=0.05
        r_tid=0.351*sep
        tol_1=0.04
        tol_2=0.04
        tol_3=0.04
        tolerance=1e-4
        noprint=.false.
        noprint=.true.
        nograph=.false.
C        nograph=.true.
        coldstart=.false.
        coldstart=.true.
        dt_write=1d44
        deltat_min=1
        deltat_max=1e6
        deltat=1e4
        trecord=0

        Gravatiy=6.67259*1e-11
        M_sun=1.9891*1e30
        c=299792458.0
        m_p=1.6726231*1e-27
        sigma_T=6.65245*1e-29 
        p_=0.1
        pi=3.141592653589793
        Rg_ = Gravatiy*m1*M_sun/c**2*100.0/1e10 
        Rs_=2*Rg_*1e10
        L_Edd = 4*pi*Gravatiy*m1*M_sun*m_p*c/sigma_T
        MdotEdd=L_Edd/p_/c**2*1000      
        
        epsilon_=40
        Reva_down=1e8
        Reva_up=1e11
        Reva(1)=Reva_down
        Evamdot(1)=0.008*MdotEdd*( (Reva(1)/Rs_)**(1.0/4.0) + 
     1                       epsilon_* (Reva(1)/(800*Rs_))**2)**(-1)
        do k=2,100
            Reva(k)=Reva_down*(10.0**(3.0/99.0))**(k-1)
            Evamdot(k)=0.008*MdotEdd*( (Reva(k)/Rs_)**(1.0/4.0) + 
     1                       epsilon_* (Reva(k)/(800*Rs_))**2)**(-1)
        end do

       a_ = 1.
       b_ = 10.63
       c_ = 8.914
       do i=1,701
          DataX(i) = 0+0.1*(i-1)
          DataY(i)= a_*exp(-((DataX(i)-b_)/c_)**2)
       end do

        if(coldstart) then
           temps=0
           deltat0=deltat
           deltat=1d10
           xmp16=100
           xmp160=xmp16
           alpha_c0=alpha_c
           alpha_c=alpha_h
           debut_alpha=.true.
           do k=1,n
              xi(k,1)=log(1e10*rstar)+(k-1.)**2/(n-1.)**2*
     1                  log(0.90*r_tid/rstar)
              r10=1d-10*exp(xi(k,1))
              f=(1.0000001d0-sqrt(rstar/r10))**0.25
        
              sig=7*alpha_h**(-0.8)*xmp16**0.7*r10**(-0.75)*f**2.8
              tc=3.1*alpha_h**(-0.2)*xmp16**0.3*r10**(-0.75)*f**1.2
              xi(k,2)=max(log(sig),-3d0)
              xi(k,3)=max(log(tc),1.5d0)
              xi(k,4)=5.31*xmp16
              xi(k,5)=1
              xi(k,6)=0.90*r_tid
              xi(k,7)=xmp16
              xi(k,8)=log(r_tid/rstar)
           end do
           do k=1,n
              r_p(k)=xi(k,1)
              x2_p(k)=xi(k,2)
              x3_p(k)=xi(k,3)
           end do
           call spline(r_p,x2_p,n,1d31,1d31,dx2)
           call spline(r_p,x3_p,n,1d31,1d31,dx3)
           pas=1d0/(n-1)
           imax0=8000
           imax=imax0
           call solution(xi,beta0,pas,tolerance,imax)
           if(imax.ge.imax0) then
              print *,' No steady solution found, imax = ', imax
              stop
           endif
           deltat=deltat0
           alpha_c=alpha_c0
           debut_alpha=.true.
           niter=imax
        else
           xmp160=xmp16_f
           print '(''initial time step : '',$)'
           deltat=1e4
           print '(''input file : '',$)'
           filename='infile3'
           open(12,file=filename,status='old')
           read(12,*) xii
           read(12,*) temps
           close(12)
           do k=1,ni
              xv(k)=(k-1d0)/(ni-1d0)
           end do
        
           do l=1,8
              do k=1,ni
                 yv(k)=xii(k,l)
              end do
              do k=1,n
                 xx=(k-1d0)/(n-1d0)
                 call inter(xv,yv,ni,xx,xi(k,l))
              end do
           end do
           niter=0
        endif

        open(10,file='outfile',status='unknown')
        iwrt15=0
        temps_prec=temps

        write(*,'('' binary parameters: m1='',f7.3,'' m2='',f8.3,
     1          ''Porb='',f8.3)') m1,m2,p_orb
        if(c_tid.le.1e20) then
           write(*,'('' exponential truncation at tidal radius r_tid='',
     1             f8.3,'' with width'',e10.3)'),r_tid,width
        else
           write(*,'('' r**5 tidal torque with coefficient'',e10.3)')
     1          ,c_tid
        endif
        write(*,'('' eta_impact='',f8.3,'' size_impact='',f8.3)') 
     1    eta_impact,size_impact
        write(*,'('' irradiation factor: fill='',e10.3)') fill
        write(*,'('' mu_30='',f8.3)') mu30
        write(*,'('' alpha_c='',f8.3,'' alpha_h='',f8.3)')
     1          alpha_c,alpha_h
        if(coldstart) then
          write(*,'('' initial condition: steady state with
     1          high Mdot'')')
        else
           write(*,'('' input file: '',a)') filename
        endif
        write(*,'('' deltat='',e10.3,'' min: '', e10.3,'' max :''
     1     ,e10.3)')  deltat,deltat_min,deltat_max
        write(*,800)
800     format(//,3x,'deltat',9x,'t',12x,'err',8x,'mdisk',9x,'rin',9x,
     1          'rout',8x,'Mdot,tr',6x,'Mdot,in',
     1  6x,'mb',6x,'mv',6x,'mi',4x,'rtrans1',4x,'rtrans2',2x,'niter')
        do k=1,n
           r_p(k)=xi(k,1)
           x2_p(k)=xi(k,2)
           x3_p(k)=xi(k,3)
        end do

        call graphique(xi,xmass,xmpin,rin,rout,mb,mv,mi,rtrans,rtrans1)

        call spline(r_p,x2_p,n,1d31,1d31,dx2)
        call spline(r_p,x3_p,n,1d31,1d31,dx3)

        deltat_prec=1d30
        imax0=4000 
        iwrite=0
        errm=0
        call dzero(n*(n1+n2),deltax)

        write(*,'(e10.3,e15.7,e12.3,5e13.5,3f8.3,2e11.4,i5,e10.4)') 
     1  deltat,
     1  temps/86400,errm-1,xmass,
     1  1e10*rin,1e10*rout,1e16*xmp16,xmpin,mb,mv,mi,1e10*rtrans1,
     1  1e10*rtrans,niter,fill

        do i_temps=1,200000
        if(temps.ge.200*86400.*1e4) stop
        if(rtrans.ge..3) then
           deltat_min=100
           deltat_min=10
        else
           deltat_min=1
        endif
        niter=0
        xmprec=xmass
        xmpinprec=xmpin
        xmp16=xmp16_f+(xmp160-xmp16_f)*max(1-1e-7*temps,0d0)

1       do m=1,n1+n2
           do k=1,n
              xf(k,m)=xi(k,m)+deltax(k,m)*deltat/deltat_prec
           end do
        end do

        pas=1d0/(n-1)
        imax=imax0
        call solution(xf,beta0,pas,tolerance,imax)
        niter=niter+imax
        if(imax.gt.imax0) then
           if(rate) deltat=deltat/2 
           if(deltat.lt.1e-5) then
            write(*,'("deltat taixiaole")')              
            stop 
           endif
           rate=.true.
           deltat_prec=1e30 
           goto 1
        endif
        
        rate=.false.
        
        call graphique(xf,xmass,xmpin,rin,rout,mb,mv,mi,rtrans,rtrans1)
        dxm=1e16*xmp16-0.5*(xmpin+xmpinprec)
        if(dxm.eq.0) dxm=1e-20
        errm=(xmass-xmprec)/dxm/deltat
        write(*,'(e10.3,e15.7,e12.3,5e13.5,3f8.3,2e11.4,i5,e10.3)')
     1   deltat,
     1  temps/86400,errm-1,xmass,
     1  1e10*rin,1e10*rout,1e16*xmp16,xmpin,mb,mv,mi,
     1  rtrans1,rtrans,niter,fill

        iwrite=iwrite+1
        if(iwrite.ge.10) then
           rewind(10)
           write(10,*) xf
           write(10,*) temps+deltat
           call flush(10)
           iwrite=0
        end if

        if(temps.ge.temps_prec+dt_write) then
           temps_prec=temps
           iwrt15=iwrt15+1
           write(filename,'(''outfile_'',i4.4)') iwrt15
           open(15,file=filename,status='unknown')
           write(15,*) xf
           write(15,*) temps+deltat
           close(15)
        endif

        do m=1,n1+n2
           dmax(m)=0
           do k=1,n
              deltax(k,m)=xf(k,m)-xi(k,m)
              dmax(m)=max(dmax(m),abs(deltax(k,m)))
              xi(k,m)=xf(k,m)
           end do
        end do

        do k=1,n
           r_p(k)=xi(k,1)
           x2_p(k)=xi(k,2)
           x3_p(k)=xi(k,3)
        end do

        if(deltat.ge.0.9*deltat_min) then
           call spline(r_p,x2_p,n,1d31,1d31,dx2)
           call spline(r_p,x3_p,n,1d31,1d31,dx3)
        else 
           call dzero(n,dx2)
           call dzero(n,dx3)
        endif
        imax0=100

        temps=temps+deltat
        trecord=temps         

        if(i_temps.ge.3) deltat_prec=deltat
        if(niter.ge.imax0/2) deltat_prec=1d60
        if((dmax(2).le.tol_2/2).and.(dmax(1).le.tol_1/2).and.
     1              (dmax(3).le.tol_3/2)) deltat=deltat*1.1
        if((dmax(2).gt.tol_2).or.(dmax(1).gt.tol_1).or.
     1                                       (dmax(3).gt.tol_3)) 
     1  deltat=deltat/1.2
        deltat=min(max(deltat,deltat_min),deltat_max)

        end do
        end program
        
        
        subroutine graphique(xf,xmass,xmpin,rin,rout,mb,mv,mi,rtrans,
     1             rtrans1)
        parameter (n1=4,n2=4,n=5000)
        implicit real*8 (a-h,o-z)
        logical debut/.true./,noprint,nograph
        integer pgopen
        real*8 m1,mv,mb,mi,mu30
        real*4 xg(n),xg1(n),yg1(n),yg2(n),yg3(n),yg4(n),yg5(n)
     1                              ,cons1,cons2,cons3,cons4
        character*80 filename2,ttext_     
        dimension t0(100),fb(100),fv(100),fi(100),alpha(n)
        dimension xf(n,n1+n2),fluxb(n),fluxv(n),fluxi(n)
        dimension sigmad(n),rd(n),sigma(n),t4(n),r10(n)
        dimension r_p(n),x2_p(n),x3_p(n)
        dimension dx2(n),dx3(n)
        dimension DataX(701),DataY(701)
        common/sortie/noprint,nograph
        common/parametres/m1,xmp16,alpha_c,alpha_h,rs,rstar,c_tid,p_orb,
     1        eta_impact,       size_impact,r_tid,mu30,fill,width
        common/solution_prec/r_p,x2_p,x3_p,dx2,dx3,deltat
        common/trecord/trecord,f_over,g_over
        common/DataXray/DataX,DataY,fill_2       
 
        data icnt/999/
        save
        if(debut) then
           debut=.false.
           open(25,file='bol_cor_t',status='old')
           do k=1,100
              read(25,*) t0(k),fb(k),fv(k),fi(k)
           end do
           close(25)
           alpha_m=(alpha_c+alpha_h)/2
           cons1=log10(13.4*m1**(-0.38)*alpha_c**(-0.83))

           if(.not.nograph) then
              write(filename2,'(''outputplot'',''.ps/CPS'')')
              id2=pgopen(filename2)
              if(id2.lt.1)stop
              call pgask(.false.)
C              id1=pgopen('/xw')
C              call pgask(.false.)
C              id2=pgopen('/xw')
C              call pgask(.false.)
            endif
        endif
            
        cons2=log10(8.3*m1**(-0.38)*alpha_h**(-0.77))
        i_r=0
        do k=1,n
           sigma(k)=exp(xf(k,2))
           r10(k)=exp(xf(k,1))/1e10
           r103=r10(k)**3
           t4(k)=exp(xf(k,3))
           rd(k)=1e10*r10(k)
           sigmad(k)=sigma(k)*2*3.1415926*rd(k)
           call setfill(fill,r10(1),xf(k,7),xf(k,6))

          if(trecord.lt.86400*15222)then
           till4=1.26d19*fill*min(xf(k,7),1.525359165863364d3)/r10(k)**2
          else
           till4_1=1.26d19*fill*min(xf(k,7),1.525359165863364d3)
     1                                                  /r10(k)**2
           till4_2=1.26d19*fill_2*17.9/r10(k)**2
           till4=max(till4_1,till4_2)
          endif
           ti=sqrt(sqrt(till4))

           alpha(k)=alfa(t4(k),ti,r10(k))
           yg5(k)=log10(ti)-4
           call fluxes(xf(k,1)-23.02585,xf(k,2),1d4*t4(k),xf(k,3)+9.210,
     1                  ti,te44,qplus)
           te44=(1d16*t4(k)**4-till4)*te44+1d-16*till4
           te44=max(te44,1d-16*till4)
           xlte4=0.25*log10(te44)
           te4=sqrt(sqrt(te44))

           call inter(t0,fb,100,xlte4,factb)
           call inter(t0,fv,100,xlte4,factv)
           call inter(t0,fi,100,xlte4,facti)
           fluxb(k)=10**factb*te44*r10(k)
           fluxv(k)=10**factv*te44*r10(k)
           fluxi(k)=10**facti*te44*r10(k)
        end do

        rout=r10(n)
        rin=r10(1)
        xmpin=0.1885*xf(1,4)*1e16
        k=n-1
        do while((alpha(k)-alpha_m)*(alpha(k+1)-alpha_m).ge.0 .and.
     1                                                         k.gt.1)
           k=k-1
        end do
        if(k.eq.1) then
           rtrans=0
        else
           rtrans=r10(k)+(r10(k+1)-r10(k))*(alpha(k)-alpha_m)/
     1          (alpha(k)-alpha(k+1))
        endif
        if(alpha(n).ge.alpha_m) rtrans=rout
        k=1
        do while((alpha(k)-alpha_m)*(alpha(k+1)-alpha_m).ge.0 .and.
     1                                                           k.lt.n)
           k=k+1
        end do
        if(k.eq.n) then
           rtrans1=0
        else
           rtrans1=r10(k)+(r10(k+1)-r10(k))*(alpha(k)-alpha_m)/
     1          (alpha(k)-alpha(k+1))
        endif
        if(alpha(1).ge.alpha_m) rtrans1=rin
        xmass=0
        do k=1,n-1
           xmass=xmass+(rd(k+1)-rd(k))*(sigmad(k+1)+sigmad(k))/2
        end do
        ftotb=0
        ftotv=0
        ftoti=0
        do k=1,n-1
           ftotb=ftotb+(r10(k+1)-r10(k))*(fluxb(k+1)+fluxb(k))/2
           ftotv=ftotv+(r10(k+1)-r10(k))*(fluxv(k+1)+fluxv(k))/2
           ftoti=ftoti+(r10(k+1)-r10(k))*(fluxi(k+1)+fluxi(k))/2
        end do
        mb=-2.5*log10(1.18e-3*ftotb)+12
        mv=-2.5*log10(1.18e-3*ftotv)+12
        mi=-2.5*log10(1.18e-3*ftoti)+12
        mb=max(-99d0,mb)
        mv=max(-99d0,mv)
        mi=max(-99d0,mi)

        write(ttext_,'(e15.7)') trecord/86400
        if(nograph) return
        if(trecord.lt.15200.*86400) return
        icnt=icnt+1        
        if(icnt.ge.10) then
           icnt=0
C           call pgslct(id1)
C           call pgenv(0.,n+1.,-4.,5.,0,0)
           call pgslct(id2)
           call pgsls(1)
           call pgslw(6)
           call pgenv(7.9,11.5,-4.,5.,0,0)
           call pgsci(3)
           call pgsls(2)
           call pgmove(8.,cons1-2.20)
           call pgdraw(12.,cons1+2.20)
           call pgmove(8.,cons2-2.22)
           call pgdraw(12.,cons2+2.22)
           call pgsci(1)
           call pgsls(1)
C        endif

          do k=1,n
            xg(k)=xf(k,1)*0.434294482
            xg1(k)=k
            yg1(k)=xf(k,3)*0.434294482
            yg2(k)=xf(k,2)*0.434294482
            yg3(k)=log10(abs(xf(k,4))+1e-30)
            yg4(k)=log10(alpha(k)+1e-30)
          end do
C        call pgslct(id1)
C        call pgline(n,xg1,yg1)
C        call pgsls(2)
C        call pgline(n,xg1,yg5)
C        call pgsls(1)
C        call pgsci(3)
C        call pgline(n,xg1,yg2)
Cc        call pgsci(4)
Cc        call pgline(n,xg1,yg3)
C        call pgsci(5)
C        call pgline(n,xg1,yg4)
C        call pgsci(15)
          call pgslct(id2)
          call pgslw(5)
          call pgline(n,xg,yg1)
          call pgsls(2)
          call pgline(n,xg,yg5)
          call pgsls(1)
          call pgsci(3)
          call pgline(n,xg,yg2)
C        call pgsci(4)
C        call pgline(n,xg,yg3)
          call pgsci(5)
          call pgline(n,xg,yg4)
          call pgsci(15)
          call pgsls(2)
          call pgsci(15)
          call pgtext(8.5,-2.5,ttext_)
          call pgtext(8.0,0,'.       ')
        endif
       
        return
        end        


        subroutine fg(xf,f,g,fp,c,debut)
        implicit real*8 (a-h,o-z)
        real*8 nu14,m1,kappa,nusig,mu30
        parameter (n1=4,n2=4,n=5000)
        logical debut
        dimension f(n,n1+n2),g(n,n1+n2),xf(n,n1+n2),fp(n),c(n)
        dimension r_p(n),x2_p(n),x3_p(n),x2_prec(n),x3_prec(n),
     1          alpha(n),nusig(n),coeff(n),rexp(n),ti(n),till4(n),
     1          cp(n),dsigmadt_1(n),bilan0(n),bilan(n),qimpact(n),cc(n)
        dimension sigma(n),t4(n),r10(n),sigma_prec(n),t4_prec(n),sr10(n)
        dimension dx2(n),dx3(n)
        dimension DataX(701),DataY(701)

        common/solution_prec/r_p,x2_p,x3_p,dx2,dx3,deltat
        common/parametres/m1,xmp16,alpha_c,alpha_h,rs,rstar,c_tid,p_orb,
     1          eta_impact,size_impact,r_tid,mu30,fill,width
        common/trecord/trecord,f_over,g_over
        common/DataXray/DataX,DataY,fill_2

        logical flag1,flag2,flag3,flag6,flag7
        data r10_av/1e10/,t4_av/1e10/,sigma_av/1e10/,x7_av/1e10/
        save
        
        flag1=.false.
        flag2=.false.
        flag3=.false.
        flag6=.false.
        flag7=.false.
        if(debut) sm1=sqrt(m1)

        if((xf(2,2).ne.sigma_av).or.debut) then
           do k=1,n
              sigma(k)=exp(xf(k,2))
           end do
           sigma_av=xf(2,2)
           flag2=.true.
        endif
        
        if((xf(2,3).ne.t4_av).or.debut) then
           do k=1,n
              t4(k)=exp(xf(k,3))
           end do
           t4_av=xf(2,3)
           flag3=.true.
        endif

        if((xf(2,1).ne.r10_av).or.debut) then
           i_r=0
           do k=1,n
              r10(k)=exp(xf(k,1))/1e10
              sr10(k)=sqrt(r10(k))
              if(c_tid.le.1d-20) then
                 rexp(k)=exp(width*(r10(k)-r_tid))
              else
                 rexp(k)=c_tid*r10(k)**5
              endif
              r=xf(k,1)
              i_r=i_r+1
              do while((r.lt.r_p(i_r)).and.(i_r.gt.1))
                 i_r=i_r-1
              end do
              if(i_r.ge.n) i_r=i_r-1
              do while((r.ge.r_p(i_r+1).and.(i_r).lt.n-1))
                 i_r=i_r+1
              end do
              h=r_p(i_r+1)-r_p(i_r)
              h=max(h,1d-30)
              a=max((r_p(i_r+1)-r)/h,0d0)
              a=min(a,1d0)
              b=min(1d0-a,1d0)
              b=1-a
              x2_prec(k)=a*x2_p(i_r)+b*x2_p(i_r+1)+((a**3-a)*dx2(i_r)
     1                         +(b**3-b)* dx2(i_r+1))*h**2/6
              x3_prec(k)=a*x3_p(i_r)+b*x3_p(i_r+1)+((a**3-a)*dx3(i_r)
     1                         +(b**3-b)* dx3(i_r+1))*h**2/6
              sigma_prec(k)=exp(x2_prec(k))
              t4_prec(k)=exp(x3_prec(k))
           end do
           r10_av=xf(2,1)
           flag1=.true.
        endif

        if(flag1.or.xf(6,1).ne.par6_av) then
           par6_av= xf(6,1)
           flag6=.true.
           do k=1,n
              rout=xf(k,6)
              if(c_tid.le.1d-20) then
                 rexpout=exp(width*(rout-r_tid))
              else
                 rexpout=c_tid*rout**5
              endif
              qimpact(k)=1.06e5*m1*xmp16*eta_impact*size_impact
     1              *r10(k)**2*exp(size_impact*(r10(k)/rout-1))/rout**3
           end do
        end if

        if(flag1.or.xf(2,7).ne.x7_av) then
           do k=1,n
          if(trecord.lt.86400*15222)then
             till4(k)=1.26d19*fill*min(xf(k,7),1.525359165863364d3)
     1                                                     /r10(k)**2
          else
           till4_1=1.26d19*fill*min(xf(k,7),1.525359165863364d3)
     1                                                  /r10(k)**2
           till4_2=1.26d19*fill_2*17.9/r10(k)**2
           till4(k)=max(till4_1,till4_2)
          endif
              ti(k)=sqrt(sqrt(till4(k)))
           end do
           x7_av=xf(2,7)
           flag7=.true.
        end if 

        if(flag1.or.flag3) then
           do k=1,n

          if(trecord.lt.86400*15222)then
           till4(k)=1.26d19*fill*min(xf(k,7),1.525359165863364d3)
     1                                                     /r10(k)**2
          else
           till4_1=1.26d19*fill*min(xf(k,7),1.525359165863364d3)
     1                                                  /r10(k)**2
           till4_2=1.26d19*fill_2*17.9/r10(k)**2
           till4(k)=max(till4_1,till4_2)
          endif
          
              ti(k)=sqrt(sqrt(till4(k)))
              alpha(k)=alfa(t4(k),ti(k),r10(k))
           end do
        endif
        
        if(flag1.or.flag2.or.flag3.or.flag7) then
           do k=1,n
              call fluxes(xf(k,1)-23.02585,xf(k,2),1d4*t4(k),
     1          xf(k,3)+9.210,ti(k),te44,qplus)
              t42=t4(k)**2
              t44=t42*t42
              te44=(1d16*t44-till4(k))*te44
              qplus=qplus*alpha(k)
              r102=r10(k)**2
              r103=r102*r10(k)
              om2=1.334e-4*m1/r103
              nu14=0.8888888889e-14*qplus/om2/sigma(k)
              nusig(k)=nu14*sigma(k)
              a=1.89e-3*t44*r103/sigma(k)/m1
              b=1.013*t4(k)*r103/m1
              h8=a+sqrt(a*a+b)
              rho=5e-9*sigma(k)/h8
              call param(rho,1e4*t4(k),log(rho),xf(k,3)+9.21034,p,kappa,
     1          gradad,q,1)
              cp(k)=1e-4*p*q/(rho*t4(k)*gradad)
              dt4dt=0.5d0*(t4(k)+t4_prec(k))*(xf(k,3)-x3_prec(k))/deltat
              dsigmadt=0.5d0*(sigma(k)+sigma_prec(k))*(xf(k,2)
     1                 -x2_prec(k))/deltat
              dsigmadt=(sigma(k)-sigma_prec(k))/deltat
              dsigmadt_1(k)=3.3333333e5*r102*dsigmadt
              cc(k)=3e-12*p/rho*sigma(k)
              dtidal=2.015e4*m1*nusig(k)*rexp(k)*
     1                  sr10(k)*(1-0.151*r10(k)*sr10(k)/sm1/p_orb)
              bilan0(k)=.01*r102*cp(k)*sigma(k)*dt4dt+1.134e6*r102*te44-
     1                  3e4*nusig(k)*m1/r10(k)-dtidal
              coeff(k)=1e8*r10(k)/(3d-4*nusig(k)*cp(k)*t4(k))
           end do
        endif

        if(flag1.or.flag2.or.flag3.or.flag6.or.flag7) then
           do k=1,n
              bilan(k)=bilan0(k)-qimpact(k)
           end do
        endif
        do k=1,n
           c(k)=cc(k)
           dlogtdlogr=xf(k,5)*coeff(k)
           f12=0.05*(xf(k,4)/nusig(k))**2+1
     1             +10*(width*rexp(k)/rexp(n))**2+0.1*dlogtdlogr**2
           dlogrdq=xf(k,8)/sqrt(f12)
           fp(k)=xf(k,4)/sigma(k)
           f(k,1)=1e4*r10(k)*xf(k,5)
           f(k,2)=xf(k,4)
           f(k,3)=xf(k,1)
           f(k,4)=xf(k,3)
           f(k,5)=nusig(k)*sr10(k)
           f(k,6)=xf(k,6)
           f(k,7)=xf(k,7)
           f(k,8)=xf(k,8)
           g(k,1)=bilan(k)-3e-8*t4(k)*xf(k,4)*cp(k)*dlogtdlogr
           g(k,1)=g(k,1)*dlogrdq
C          g(k,2)=dsigmadt_1(k)*dlogrdq
           call setf_over
           f_=f_over*max(0.,(alpha(k)-alpha_c)/(alpha_h-alpha_c))
     1                                       /3.3333333e5/r10(k) 
           g(k,2)=(dsigmadt_1(k)+
     1             f_*sigma(k)*(r10(k)**2)*3.3333333e5)*dlogrdq
           g(k,3)=dlogrdq
           g(k,4)=dlogtdlogr*dlogrdq
C          g(k,5)=(xf(k,4)*sr10(k)-r10(k)**2*rexp(k)*nusig(k))*dlogrdq
           g(k,5)=(xf(k,4)*sr10(k)-2./3.*f_*sigma(k)*(r10(k)**2) 
     1         *sr10(k)*1e6*g_over*3.3333333e5*r10(k)
     1         -r10(k)**2*rexp(k)*nusig(k) )*dlogrdq
           g(k,6)=0
           g(k,7)=0
           g(k,8)=0
        end do
        end


        subroutine h1(h,par)
        implicit real*8 (a-z)
        integer n1,n2,n
        parameter (n1=4,n2=4,n=5000)
        dimension h(4),par(8)
        dimension x2_p(n),x3_p(n),r_p(n),dx2(n),dx4(n),dx3(n)
        common/parametres/m1,xmp16,alpha_c,alpha_h,rs,rstar,c_tid,p_orb,
     1          eta_impact,size_impact,r_tid,mu30,fill,width
        common/solution_prec/r_p,x2_p,x3_p,dx2,dx3,deltat

        xmpin=0.1885*par(4)
        xmp18=xmpin/100
        rin=rstar*1d10
        if(mu30.gt.0) then
           rin=setrmin(xmp18*1e18)
        endif
        h(1)=par(2)+3
        h(2)=par(1)-log(rin)
        h(3)=par(5)
        h(4)=par(7)-0.1885*par(4)
        return
        end


        subroutine h2(h,par)
        implicit real*8 (a-z)
        integer n1,n2,n
        parameter (n1=4,n2=4,n=5000)
        integer i_r
        dimension h(4),par(8)
        dimension x2_p(n),x3_p(n),r_p(n)
        dimension dx2(n),dx3(n)
        dimension DataX(701),DataY(701)
        common/solution_prec/r_p,x2_p,x3_p,dx2,dx3,deltat
        common/parametres/m1,xmp16,alpha_c,alpha_h,rs,rstar,c_tid,p_orb,
     1          eta_impact,size_impact,r_tid,mu30,fill,width
        common/DataXray/DataX,DataY,fill_2

        r=par(1)
        r10=exp(r)/1e10
        r10p=exp(r_p(n))/1e10
        sigma=exp(par(2))
        t4=exp(par(3))
        sigmap=exp(x2_p(n))
        if(r10.ge.r10p) then
           sigmapp=sigmap
        else
           i_r=n-1
           do while((r.lt.r_p(i_r)).and.(i_r.gt.1))
              i_r=i_r-1
           end do
           hh=r_p(i_r+1)-r_p(i_r)
           a=max((r_p(i_r+1)-r)/hh,0d0)
           a=min(a,1d0)
           b=min(1d0-a,1d0)
           x2_prec=a*x2_p(i_r)+b*x2_p(i_r+1)
           sigmapp=exp(x2_prec)
        endif
       
        if(trecord.lt.86400*15222)then
          till4=1.26d19*fill*min(par(7),1.525359165863364d3)/r10**2
        else
           till4_1=1.26d19*fill*min(par(7),1.525359165863364d3)
     1                                                  /r10**2
           till4_2=1.26d19*fill_2*17.9/r10**2
           till4=max(till4_1,till4_2)
        endif

        ti=sqrt(sqrt(till4))
        call fluxes(par(1)-23.02585,par(2),1d4*t4,par(3)+9.210,ti,te44,
     1               qplus)
        alpha=alfa(t4,ti,r10)
        qplus=qplus*alpha
        om2=1.334e-4*m1/r10**3
        nu14=0.8888888889e-14*qplus/om2/sigma
        rpoin2=(r10**2-r10p**2)/deltat

        h(1)=0.1885*par(4)/xmp16-1+3.1415926e4*rpoin2*(0.500*sigmap+
     1          0.500*sigmapp)/xmp16
        h(2)=xmp16*(1-sqrt(rs/r10))-0.0942*nu14*sigma
        h(3)=par(5)
        h(4)=par(6)-r10
        return
        end


        function alfa(t,tt,r)
        implicit real*8 (a-z)
        logical debut
        common/parametres/m1,xmp16,alpha_c,alpha_h,rs,rstar,c_tid,p_orb,
     1          eta_impact,size_impact,r_tid,mu30,fill,width
        common/alf/debut
        save

        lc=log10(alpha_c)
        lh=log10(alpha_h)
        epillw=(tt/1e4)**2.0
        logTcmaxw=-0.1*lc
     1              -0.05*epillw*log10(r)
        Tcmaxw4=10700*exp(2.302585093*logTcmaxw)

        Tcminepill=20900.0-11300.0*epillw
        logTcminw=-0.22*lh
     1       -0.01*log10(m1)+(0.05-0.12*epillw)*log10(r)
        Tcminw4=Tcminepill*exp(2.302585093*logTcminw)

        T4_crit=(Tcmaxw4+Tcminw4)/2.0
        T4_crit=max(T4_crit,0.0)
        T4_crit8=T4_crit**8./1e32

        t2=t*t
        t4=t2*t2
        t8=t4*t4
        f=1/(1+T4_crit8/t8)
        alfa=(lh-lc)*f+lc
        alfa=exp(2.302585093*alfa)
        return
        end


        subroutine param(rho,t,lrho1,lt1,p,kappa,gradad,q,n)
        implicit real*8 (a-h,k-z)
        integer nn1,nn2,n,k,n1,n2
        parameter (n1=217,n2=108)
        real*4 logrho(n1),logt(n2),tp(n1,n2),tkappa(n1,n2),
     1          tgradad(n1,n2),tq(n1,n2)
        dimension rho(n),t(n),lrho1(n),lt1(n),p(n),kappa(n),gradad(n),
     1          q(n)
        logical first/.true./
        save

        if(first) then
           open(11,file='p30.dat',status='old',form='unformatted')
           read(11) nn1,nn2
           read(11) logrho,logt
           read(11) tp,tkappa,tgradad,tq
           close(11) 
           first=.false.
           deltalogrho=(logrho(n1)-logrho(1))/(n1-1d0)
           deltalogt=(logt(n2)-logt(1))/(n2-1d0)
           deltalogrho1=1/deltalogrho
           deltalogt1=1/deltalogt
           lr0=logrho(1)
           lt0=logt(1)
        end if

        do k=1,n
        lt=lt1(k)*.434294482d0
        lrho=lrho1(k)*.434294482d0
        i1=max(min(int((lrho-lr0)*deltalogrho1+1),n1-1),1)
        j1=max(min(int((lt-lt0)*deltalogt1+1),n2-1),1)
        i2=i1+1
        j2=j1+1
        dr=(lrho-logrho(i1))*deltalogrho1
        dt=(lt-logt(j1))*deltalogt1

        p1=(1-dr)*tp(i1,j1)+dr*tp(i2,j1)
        p2=(1-dr)*tp(i1,j2)+dr*tp(i2,j2)
        p(k)=(1-dt)*p1+dt*p2

        gradad1=(1-dr)*tgradad(i1,j1)+dr*tgradad(i2,j1)
        gradad2=(1-dr)*tgradad(i1,j2)+dr*tgradad(i2,j2)
        gradad(k)=(1-dt)*gradad1+dt*gradad2

        q1=(1-dr)*tq(i1,j1)+dr*tq(i2,j1)
        q2=(1-dr)*tq(i1,j2)+dr*tq(i2,j2)
        q(k)=(1-dt)*q1+dt*q2

        p(k)=p(k)*rho(k)*t(k)*1e8+2.523e-15*t(k)**4
        end do
        return
        end


        subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)
        implicit real*8 (a-h,o-z)
        dimension abd(lda,n),ipvt(n)
        m=ml+mu+1
        ju=0
        do k=1,n-1
           lm=min0(ml,n-k)
           l=0
           smax=abs(abd(m,k))
           do kk=1,lm
              stest=abs(abd(m+kk,k))
              if(stest.ge.smax) then
                 l=kk
                 smax=stest
              endif
           end do
           ipvt(k)=l+k
           l=l+m
           t=abd(l,k)
           abd(l,k)=abd(m,k)
           abd(m,k)=t
           t=-1d0/abd(m,k)
           do kk=1,lm
              abd(m+kk,k)=t*abd(m+kk,k)
           end do
           ju=min0(max0(ju,mu+ipvt(k)),n) 
           mm=m
           do j=k+1,ju
              l=l-1
              mm=mm-1
              t=abd(l,j)
              abd(l,j)=abd(mm,j)
              abd(mm,j)=t
              do kk=1,lm
                abd(mm+kk,j)=abd(mm+kk,j)+t*abd(m+kk,k)
              end do
           end do
        end do
        ipvt(n)=n
        return
        end


        subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
        integer lda,n,ml,mu,ipvt(n),job
        integer k,kb,l,lb,lm,m
        real*8 abd(lda,n),b(n)
        real*8 t
        m=mu+ml+1
        do k=1,n-1
           lm=min0(ml,n-k)
           l=ipvt(k)
           t=b(l)
           b(l)=b(k)
           b(k)=t
           do kk=1,lm
              b(k+kk)=b(k+kk)+t*abd(m+kk,k)
           end do
        end do
        do kb=1,n
           k=n+1-kb
           b(k)=b(k)/abd(m,k)
           lm=min0(k,m)-1
           lb=k-lm
           t=-b(k)
           do kk=lb,lb+lm-1
              b(kk)=b(kk)+t*abd(m-k+kk,k)
           end do
        end do
        return
        end


        subroutine dzero(n,x)
        real*8 x(*)
        do 10 i=1,n
           x(i)=0d0
10      continue
        return
        end


        integer function idamax(n,sx,incx) 
        real*8 sx(*),smax,stest
        integer i,incx,n
        idamax=0
        if(n.lt.1) return
        idamax=1
        if(n.eq.1) return
        smax=abs(sx(1))
        do i=2,n
           stest=abs(sx(i))
           if(stest.gt.smax) then
              idamax=i
              smax=stest
           endif
        end do
        return
        end


        subroutine solution(xf,beta,pas,tolerance,nmax)
        implicit real*8 (a-h,o-z)
        parameter (n1=4,n2=4,n=5000)
        dimension xf(n,n1+n2),beta(n1+n2)
        dimension dte(n*(n1+n2)),xl(n1+n2)
        logical flag,noprint,nograph
        real*8 lambda,lambda_max,lambda_prec,lambda_test
        common/sortie/noprint,nograph
        save
        alpha1=0.15
        alpha2=0.15
        lambda_max=1
        n_total=0
        errmax=1e37
        err_prec=1e37
        lambda_prec=1
        unsurn=1d0/n
        fdp=1e-8
        flag=.false.
        do while((abs(errmax).gt.tolerance).and.(n_total.lt.nmax))
                n_total=n_total+1
           call increment(xf,beta,pas,dte,xl,dtemax,fdp)
           imax=idamax(n*(n1+n2),dte,1)
           kmax=(imax-1)/(n1+n2)+1
           lmax=imax-(kmax-1)*(n1+n2)
           errmax=dte(imax)

           if(n_total.ge.(nmax/1.5/lambda_max)) lambda_max=
     1                                                0.75d0*lambda_max
           lambda_test=1
           if(err_prec.ne.0) lambda_test=(1+lambda+errmax/err_prec)/2
           if(lambda_test.le.0.2 .or. lambda_test.ge.1.5 
     1          .or. imax.ne.imax_prec) flag=.false.
           err_prec=errmax
           lambda_prec=lambda
           imax_prec=imax
           lambda=lambda_max
           do m=1,n1+n2
              do k=1,n
                 i=(n1+n2)*(k-1)+m
                 lambda=min(lambda,(alpha1*abs(xf(k,m))+alpha2*xl(m))/
     1                  (dabs(dte(i))+1d-30))
              end do
           end do
           
           if(lambda.le.1e-6) n_total=nmax+1
           if(flag.and.lambda.ge.0.5) lambda=lambda_test
           flag=.true.
           if(.not.noprint) print '(''errmax = '',e10.4,'' lambda = '',
     1                  f10.8,2x,i4,1x,i3,
     1          3(1x,e13.7))',errmax,lambda,kmax,lmax,xf(kmax,lmax)
           do m=1,n1+n2
              do k=1,n
                 i=m+(k-1)*(n1+n2)
                 xf(k,m)=xf(k,m)+lambda*dte(i)
              end do
           end do
        end do
        nmax=n_total
        if(abs(errmax).ge.tolerance) nmax=nmax+1
        return
        end


        subroutine inter(x,y,n,t,val)
        implicit real*8 (a-h,o-z)
        dimension x(n),y(n)
        data i/1/
        save
        i=min(n-1,i)
        if(x(1).le.x(2)) then
           do while((t.lt.x(i)).and.(i.gt.1))
              i=i-1
           end do
           do while((t.ge.x(i+1)).and.(i.lt.n-1))
              i=i+1
           end do
        else
           do while((t.ge.x(i)).and.(i.gt.1))
              i=i-1
           end do
           do while((t.lt.x(i+1)).and.(i.lt.n-1))
              i=i+1
           end do
        endif
        if(abs(x(i+1)-x(i)).le.1e-15) then
           val=y(i)
        else
           val=y(i)+(y(i+1)-y(i))/(x(i+1)-x(i))*(t-x(i))
        endif
        return
        end


        subroutine spline(x,y,n,yp1,ypn,y2)
        implicit real*8 (a-h,o-z)
        parameter (nmax=20000)
        dimension x(n),y(n),y2(n),u(nmax)
        if (yp1.gt..99e30) then
           y2(1)=0
           u(1)=0
        else
           y2(1)=-0.5d0
           u(1)=(3/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif
        do i=2,n-1
           sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
           p=sig*y2(i-1)+2
           y2(i)=(sig-1)/p
           u(i)=(6*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1          /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        end do
        if (ypn.gt..99e30) then
           qn=0
           un=0
        else
           qn=0.5d0
           un=(3/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif
        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1)
        do k=n-1,1,-1
           y2(k)=y2(k)*y2(k+1)+u(k)
        end do
        return
        end


        subroutine increment(xf,beta,pas,dte,xl,dtemax,fdp,debut1)
        implicit real*8 (a-h,o-z)
        parameter (n1=4,n2=4,n=5000)
        dimension xf(n,n1+n2),beta(n1+n2),b1(n1+n2),b2(n1+n2)
        dimension f1(n1+n2),f2(n1+n2)
        dimension xf1(n,n1+n2,n1+n2),xf2(n,n1+n2,n1+n2),dpv(n,n1+n2)
        dimension ff1(n,n1+n2,n1+n2),gg1(n,n1+n2,n1+n2),
     1  ff2(n,n1+n2,n1+n2),
     1  gg2(n,n1+n2,n1+n2),dpv1(n,n1+n2),fp1(n,n1+n2),fp2(n,n1+n2),
     1  c1(n,n1+n2),c2(n,n1+n2),derfp(n,n1+n2),derc(n,n1+n2),
     1  c0(n),fp0(n)
        dimension derf(n,n1+n2,n1+n2),derg(n,n1+n2,n1+n2)
        dimension tab(5*n1+4*n2-2,n*(n1+n2)),par1(n1+n2),par2(n1+n2),
     1  dte(n*(n1+n2)),xl(n1+n2),ipvt((n1+n2)*n),xxl(n,n1+n2)
        logical debut,debut1,init/.true./,noprint,nograph
        common/sortie/noprint,nograph
        data nav/10/
        save
        if(init) then
           init=.false.
           do m=1,n2+n1
              b1(m)=pas*beta(m)
              b2(m)=pas*(1-beta(m))
           end do
           unsurn=1d0/n
           do k=1,n
              derf(k,n1+n2,n1+n2)=1
           end do
        endif
        debut=.true.
        call dzero(n*(n1+n2)*(5*n1+4*n2-2),tab)

        do m=1,n1+n2
           xl(m)=0
           do k=1,n
              xl(m)=xl(m)+abs(xf(k,m))*unsurn*0.2d0
           end do
           if(xl(m).le.1d-6) xl(m)=1
c          if(xl(m).le.1d-9) xl(m)=1
        end do

        do m=1,n1+n2
           do k=1,n
              xxl(k,m)=0
              k1=max(k-nav/2,1)
              k2=min(k1+nav,n)
              if(k2.eq.n) k1=k2-nav
              do l=1,nav
                xxl(k,m)=xxl(k,m)+abs(xf(l,m))
              end do
              xxl(k,m)=xxl(k,m)/nav
           end do
        end do

        do i=1,n1+n2
           par1(i)=xf(1,i)
           par2(i)=xf(1,i)
        end do
        call h1(f1,par1)
        do l=1,n1
           dte(l)=-f1(l)
        end do
        do m=1,n1+n2
           dp=fdp*(abs(par1(m))+1d-0*xxl(1,m))
           par1(m)=par1(m)-dp
           par2(m)=par2(m)+dp
           call h1(f1,par1)
           call h1(f2,par2)
           par1(m)=par1(m)+dp
           par2(m)=par2(m)-dp
           do l=1,n1
              id=3*n1+3*n2+l-m-1
              der=0.5d0*(f2(l)-f1(l))/dp
              tab(id,m)=der
           end do
        end do

        call fg(xf,ff1,gg1,fp0,c0,debut)
        debut=.false.
        do l=1,n1+n2
           do k=2,n
              iligne=l+n1+(n1+n2)*(k-2)
              dte(iligne)=-ff1(k,l,1)+ff1(k-1,l,1)+b1(l)*gg1(k,l,1)+
     1                          b2(l)*gg1(k-1,l,1)
           end do
        end do
        do k=2,n
           iligne=1+n1+(n1+n2)*(k-2)
           dte(iligne)=dte(iligne)+0.5d0*(c0(k-1)+c0(k))*(fp0(k-1)-
     1                                                       fp0(k))
        end do

        do m=1,n1+n2-1
           do k=1,n
              dpv(k,m)=fdp*(abs(xf(k,m))+1d-0*xxl(k,m))
              dpv1(k,m)=1/dpv(k,m)
           end do
        end do
        do m=1,n1+n2-1
           do l=1,n1+n2
              do k=1,n
                 xf1(k,l,m)=xf(k,l)
                 xf2(k,l,m)=xf(k,l)
              end do
           end do
        end do
        do m=1,n1+n2-1
           do k=1,n
              xf1(k,m,m)=xf(k,m)-dpv(k,m)
              xf2(k,m,m)=xf(k,m)+dpv(k,m)
           end do
        end do
        tmp=1/xf(1,n1+n2)
        do l=1,n1+n2
           do k=1,n
              derg(k,l,n1+n2)=gg1(k,l,1)*tmp
           end do
        end do
        do m=n1+n2-1,1,-1
           call fg(xf1(1,1,m),ff1(1,1,m),gg1(1,1,m),fp1(1,m),c1(1,m),
     1                                                            debut)
           call fg(xf2(1,1,m),ff2(1,1,m),gg2(1,1,m),fp2(1,m),c2(1,m),
     1                                                            debut)
           do l=1,n1+n2
              do k=1,n
                 derf(k,l,m)=0.5d0*(ff2(k,l,m)-ff1(k,l,m))*dpv1(k,m)
                 derg(k,l,m)=0.5d0*(gg2(k,l,m)-gg1(k,l,m))*dpv1(k,m)
              end do
           end do
           do k=1,n
              derfp(k,m)=0.5d0*(fp2(k,m)-fp1(k,m))*dpv1(k,m)
              derc(k,m)=0.5d0*(c2(k,m)-c1(k,m))*dpv1(k,m)
          end do
        end do

        call fill(tab,derf,derg,derfp,derc,fp0,c0,b1,b2)

        do i=1,n1+n2
           par1(i)=xf(n,i)
           par2(i)=xf(n,i)
        end do
        call h2(f1,par1)
        do l=1,n2
           iligne=(n1+n2)*(n-1)+n1+l
           dte(iligne)=-f1(l)
        end do
        do m=1,n1+n2
           dp=fdp*(abs(par1(m))+1d-0*xxl(n,m))
           par1(m)=par1(m)-dp
           par2(m)=par2(m)+dp
           call h2(f1,par1)
           call h2(f2,par2)
           par1(m)=par1(m)+dp
           par2(m)=par2(m)-dp
           icol=(n1+n2)*(n-1)+m
           do l=1,n2
              der=(f2(l)-f1(l))/2/dp
              id=4*n1+3*n2-1-m+l
              tab(id,icol)=der
           end do
        end do

        if(.not.noprint) then
           i=idamax(n*(n1+n2),dte,1)
           dtemax=abs(dte(i))
           kmax=(i-n1-1)/(n1+n2)+1
           lmax=i-n1-(kmax-1)*(n1+n2)
           print '('' dte max:'',e10.3,i6,i3,4x,$)',dtemax,
     1                  kmax,lmax
        endif

        call dgbfa(tab,5*n1+4*n2-2,(n1+n2)*n,2*n1+n2-1,n1+2*n2-1,
     1                                                       ipvt,ier)
        call dgbsl(tab,5*n1+4*n2-2,(n1+n2)*n,2*n1+n2-1,n1+2*n2-1,ipvt,
     1                                                          dte,0)

        return
        end
        
        
        subroutine fill(tab,derf,derg,derfp,derc,fp0,c0,b1,b2)
        implicit real*8 (a-h,o-z)
        parameter (n1=4,n2=4,n=5000)
        dimension b1(n1+n2),b2(n1+n2)
        dimension derfp(n,n1+n2),derc(n,n1+n2),c0(n),fp0(n)
        dimension derf(n,n1+n2,n1+n2),derg(n,n1+n2,n1+n2)
        dimension tab(5*n1+4*n2-2,n*(n1+n2))

        do m=1,n1+n2
           do l=1,n1+n2
              id=4*n1+3*n2-1-m+l
              do k=1,n-1
                 icol=m+(k-1)*(n1+n2)
                 tab(id,icol)=-derf(k,l,m)-b2(l)*derg(k,l,m)
              end do
           end do
           do l=1,n1+n2
              id=4*n1+3*n2-1-m+l
              do k=2,n
                 icol=m+(k-2)*(n1+n2)
                 tab(id-n1-n2,icol+n1+n2)=derf(k,l,m)-b1(l)*derg(k,l,m)
              end do
           end do
           id=4*n1+3*n2-m
           do k=2,n
              icol=m+(k-2)*(n1+n2)
              tab(id,icol)=tab(id,icol)+0.5d0*derc(k-1,m)*(fp0(k)
     1                    -fp0(k-1))-0.5d0*(c0(k)+c0(k-1))*derfp(k-1,m)
              tab(id-n1-n2,icol+n1+n2)=tab(id-n1-n2,icol+n1+n2)+
     1               0.5d0*derc(k,m)*(fp0(k)-fp0(k-1))+0.5d0*(c0(k)+
     1               c0(k-1))*derfp(k,m)
           end do
        end do
        return
        end


        function rcirc(q)
        implicit real*8 (a-h,o-z)
        dimension x(14),y(14)
        data x/0.0667,0.15,0.3,0.5,0.75,1.,1.6667,3.,4.5,5.5,7.,8.5,
     1        10.,15./
        data y/0.0404,0.0492,0.0589,0.0683,0.0780,0.0865,0.1058,
     1        0.1363,0.1631,0.1781,0.1978,0.2157,0.2297,0.2698/
        call inter(x,y,14,q,rc)
        rcirc=rc
        return
        end


        subroutine fluxes(lrn,lsigman,tc,ltcn,ti,te44,qplus)
        implicit real*8 (a-h,l,o-z)
        parameter (n=9,nr=70,ns=120,nt=190)
        real*8 m1,mu30
        real*4 te(nt,ns,nr,n),fm(nt,ns,nr,n)
        real*8 xx(n),yy1(n),yy2(n),yp(n)
        character*80 filename
        common/parametres/m1,xmp16,alpha_c,alpha_h,rs,rstar,c_tid,p_orb
     1             ,eta_impact, size_impact,r_tid,mu30,fill,width
        logical debut2/.true./
        save
        if(debut2) then
           debut2=.false.
           do i=1,n
              write(filename,'(''teff_m=1.2msun.ill_'',i1,''8.bin'')')
     1                                                               i-1
              open(90,file=filename,status='old',
     1                                             form='unformatted')
              read(90) r0,s0,t0,dr0,ds0,dt0
              do ir=1,70
                 read(90) ((te(it,is,ir,i),it=1,nt),is=1,ns)
                 read(90) ((fm(it,is,ir,i),it=1,nt),is=1,ns)
              end do
              close(90)
           end do
           r0=r0+(log10(m1)-0.07918124605d0)/3d0
           dr1=1d0/dr0
           ds1=1d0/ds0
           dt1=1d0/dt0
        endif
        lr=lrn*0.434294482
        lsigma=lsigman*0.434294482
        ltc=ltcn*0.434294482

        ir1=max(min(int((lr-r0)*dr1+1),nr-1),1)
        is1=max(min(int((lsigma-s0)*ds1+1),ns-1),1)
        it1=max(min(int(1+(ltc-t0)*dt1),nt-1),1)
        
        ir2=ir1+1
        is2=is1+1
        it2=it1+1

        ds=(lsigma-s0-ds0*(is1-1))*ds1
        dt=(ltc-t0-dt0*(it1-1))*dt1
        dr=(lr-r0-dr0*(ir1-1))*dr1

        do ix=1,n
           xx(ix)=(ix-1d0)/(n-1d0)*tc
           f1=(1-dt)*fm(it1,is1,ir1,ix)+dt*fm(it2,is1,ir1,ix)
           f2=(1-dt)*fm(it1,is2,ir1,ix)+dt*fm(it2,is2,ir1,ix)
           fa=(1-ds)*f1+ds*f2

           f1=(1-dt)*fm(it1,is1,ir2,ix)+dt*fm(it2,is1,ir2,ix)
           f2=(1-dt)*fm(it1,is2,ir2,ix)+dt*fm(it2,is2,ir2,ix)
           fb=(1-ds)*f1+ds*f2

           yy2(ix)=(1-dr)*fa+dr*fb

           t1=(1-dt)*te(it1,is1,ir1,ix)+dt*te(it2,is1,ir1,ix)
           t2=(1-dt)*te(it1,is2,ir1,ix)+dt*te(it2,is2,ir1,ix)
           ta=(1-ds)*t1+ds*t2

           t1=(1-dt)*te(it1,is1,ir2,ix)+dt*te(it2,is1,ir2,ix)
           t2=(1-dt)*te(it1,is2,ir2,ix)+dt*te(it2,is2,ir2,ix)
           tb=(1-ds)*t1+ds*t2

           yy1(ix)=(1-dr)*ta+dr*tb
        end do
        a=1+(n-1)*min(ti/tc,1d0)
        ix=min(a,n-0.999d0)
        b=a-ix
        a=1-b
        a3a=a*(a*a-1)
        b3b=b*(b*b-1)
        h26=(tc/(n-1))**2/6
        call splines(xx,yy1,n,yp)
        te44=a*yy1(ix)+b*yy1(ix+1)+(a3a*yp(ix)+b3b*yp(ix+1))*h26
        call splines(xx,yy2,n,yp)
        qplus=a*yy2(ix)+b*yy2(ix+1)+(a3a*yp(ix)+b3b*yp(ix+1))*h26

        te44=1d-16*exp(min(te44,40d0))
        qplus=exp(min(qplus,90d0))
        return
        print '(''File '',A,'' does not exist'')',filename
        stop
        end


        subroutine splines(x,y,n,y2)
        implicit real*8 (a-h,o-z)
        parameter (nmax=2000)
        dimension x(n),y(n),y2(n),u(nmax)
        y2(1)=0
        u(1)=0
        h2f6=6/(x(2)-x(1))**2
        do i=2,n-1
           y2(i)=-1d0/(4d0+y2(i-1))
           u(i)=-(h2f6*(y(i+1)+y(i-1)-2d0*y(i))-u(i-1))*y2(i)
        end do
        y2(n)=0
        do k=n-1,1,-1
           y2(k)=y2(k)*y2(k+1)+u(k)
        end do
        return
        end


        subroutine setfill(fillq,rin_,xmp_16,rout_)
        implicit real*8 (a-h,o-z)
        real*8 fillX,m1,mu30,fill_2
        dimension DataX(701),DataY(701)
        common/parametres/m1,xmp16,alpha_c,alpha_h,rs,rstar,c_tid,p_orb,
     1          eta_impact,size_impact,r_tid,mu30,fill,width
        common/trecord/trecord,f_over,g_over
        common/DataXray/DataX,DataY,fill_2

       if(rin_.gt.rstar) then
            if(rin_.lt.5.0e9)then
                fill_w=1e-3*(rin_/rstar)**(-2)
            else
                fill_w=1e-3*(xmp_16)**1.32
            endif
       else
            fill_w=1e-3
       endif
        fill=fill_w   

        if(trecord.lt.86400*(15222+70)) then
            tX_=(trecord-15222*86400)/86400
            call inter(DataX,DataY,701,tX_,fillX)
            fill_2=fillX*1e-3
        else
            fill_2=0.  
        endif
        return
        end


        subroutine setf_over
        implicit real*8 (a-h,o-z)
        real*8 m1,mu30
        common/parametres/m1,xmp16,alpha_c,alpha_h,rs,rstar,c_tid,p_orb,
     1          eta_impact,size_impact,r_tid,mu30,fill,width
        common/trecord/trecord,f_over,g_over
        f_over=0.
        g_over=9.5e-5
        if((trecord.gt.(15222+11)*86400).and.
     1                   (trecord.lt.(15222+15)*86400))then
            f_over=.9e-2/4*(trecord-(15222+11)*86400)/86400
        elseif((trecord.gt.(15222+15)*86400).and.
     1                   (trecord.lt.(15222+73)*86400))then
            f_over=.9e-2
        endif
        end


        function setrmin(Mdotin)
        implicit real*8 (a-h,o-z)
        real*8 r_evap,Mdotin
        dimension Reva(100),Evamdot(100)
        common/Revap/Reva,Evamdot
        common/trecord/trecord,f_over,g_over
        if(Mdotin.gt.1.525359165863364e+16) then
            r_evap=5e8
        else            
            call inter(Evamdot,Reva,100,Mdotin,r_evap)
                 if(Mdotin.lt.3.736069672544239e+14)then
                     r_evap=5.0e9
                 endif        
        endif
        setrmin=r_evap
        return
        end
