subroutine set_geojac !{{{
use prmtr
use variable
implicit none
double precision dxdxi,drdxi,dxdeta,drdeta
double precision det
integer i,j

do j=1,nj-1
   do i=1,ni
      dxdxi=x(i,j)-x(i-1,j)
      drdxi=r(i,j)-r(i-1,j)
      
      dxdeta=(x(i,j+1)+x(i-1,j+1)&
             -x(i,j-1)-x(i-1,j-1))*0.25d0
      drdeta=(r(i,j+1)+r(i-1,j+1)&
             -r(i,j-1)-r(i-1,j-1))*0.25d0
      
      det=dxdxi*drdeta-drdxi*dxdeta
      
      !geojaci(1,1,i,j)=dxdxi
      !geojaci(1,2,i,j)=dxdeta
      !geojaci(2,1,i,j)=drdxi
      !geojaci(2,2,i,j)=drdeta

      geojaci(1,1,i,j)= drdeta/det
      geojaci(1,2,i,j)=-dxdeta/det
      geojaci(2,1,i,j)=-drdxi /det
      geojaci(2,2,i,j)= dxdxi /det
   end do
end do

do j=1,nj
   do i=1,ni-1
      dxdeta=x(i,j)-x(i,j-1)
      drdeta=r(i,j)-r(i,j-1)

      dxdxi=(x(i+1,j)+x(i+1,j-1)&
            -x(i-1,j)-x(i-1,j-1))*0.25d0
      drdxi=(r(i+1,j)+r(i+1,j-1)&
            -r(i-1,j)-r(i-1,j-1))*0.25d0
      
      det=dxdxi*drdeta-drdxi*dxdeta
      
      !geojacj(1,1,i,j)=dxdxi
      !geojacj(1,2,i,j)=dxdeta
      !geojacj(2,1,i,j)=drdxi
      !geojacj(2,2,i,j)=drdeta

      geojacj(1,1,i,j)= drdeta/det
      geojacj(1,2,i,j)=-dxdeta/det
      geojacj(2,1,i,j)=-drdxi /det
      geojacj(2,2,i,j)= dxdxi /det
   end do
end do
end subroutine set_geojac  !}}}
subroutine set_viscous
use prmtr
use variable
implicit none
double precision, dimension(6) :: Ev, Fv
double precision, dimension(6) :: dwdxi,dwdata,dwdx,dwdr
double precision, dimension(0:ni+1,0:nj+1) :: Tdeg 
double precision temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
double precision k_heat
double precision dudxi,dudeta
double precision dudx,dudr
double precision dvdxi,dvdeta
double precision dvdx,dvdr
double precision dTdxi,dTdeta
double precision dTdx,dTdr
double precision dTurenedxi,dTurenedeta
double precision dTurenedx ,dTurenedr
double precision ddissipdxi,ddissipdeta
double precision ddissipdx ,ddissipdr
double precision rho,u,v,Temper,Turene,dissip
double precision Rgas
double precision tau_xx,tau_xr,tau_rr
double precision qx,qr
double precision divu

double precision D,nu,Rt,mu_t,nu_t,u_eta,non_y,f_mu,f_eta,R_stress1,R_stress2,R_stress3
double precision, parameter::Sc=1d0
double precision, parameter::Prentl=7d-1
double precision, parameter::C_mu      = 0.09d0
double precision, parameter::C_eta1    = 1.5d0
double precision, parameter::C_eta2    = 1.9d0
double precision, parameter::sigma_k   = 1.4d0
double precision, parameter::sigma_eta = 1.4d0
double precision temp10, temp11, temp12, temp13

integer i,j,k

!$omp parallel do shared(Tdeg,j),&
!$omp             private(i)
do j=0,nj
   do i=0,ni
      Tdeg(i,j)=w(4,i,j)/w(1,i,j)/gas_specific
   end do
end do
!$omp end parallel do

!set TGvi!{{{
!$omp parallel do shared(vis_i,j),&
!$omp             private(rho,u,v,Temper,&
!$omp             k_heat,dudxi,dvdxi,dTdxi,&
!$omp             dudeta,dvdeta,dTdeta,&
!$omp             dudx,dvdx,dTdx,dudr,dvdr,dTdr,&
!$omp             divu,tau_xx,tau_xr,tau_rr,qx,qr,&
!$omp             Ev,Fv,i)
do j=1,nj-1
   do i=1,ni
      rho    = (w(1 ,i,j)+w(1 ,i-1,j))*0.5d0
      u      = (w(2 ,i,j)+w(2 ,i-1,j))*0.5d0
      v      = (w(3 ,i,j)+w(3 ,i-1,j))*0.5d0
      Turene = (w(5 ,i,j)+w(5 ,i-1,j))*0.5d0
      dissip = (w(6 ,i,j)+w(6 ,i-1,j))*0.5d0
      Temper = (Tdeg(i,j)+Tdeg(i-1,j))*0.5d0 !total energy

      nu    = mu/rho
      k_heat= mu*gamma/(gamma-1d0)*gas_specific/Prentl               !k from mu by Prandtl number
      u_eta = (nu*dissip)**(1d0/4d0)
      Rt    = Turene**2/(nu*dissip)
      non_y = u_eta*length/nu

      f_mu  = (1d0-exp(-non_y/14d0 ))**2*(1d0+(5d0/(Rt**(3d0/4d0)))*exp(-(Rt/200d0)**2))
      f_eta = (1d0-exp(-non_y/3.1d0))**2*(1d0-                0.3d0*exp(-(Rt/6.5d0)**2))

      nu_t     = C_mu*f_mu*Turene**2/dissip
      mu_t     = nu_t*rho

      dudxi       =  w(2,i,j) -w(2,i-1,j)
      dvdxi       =  w(3,i,j) -w(3,i-1,j)
      dTdxi       = Tdeg(i,j)-Tdeg(i-1,j)
      dTurenedxi  =  w(5,i,j) -w(5,i-1,j)
      ddissipdxi  =  w(6,i,j) -w(6,i-1,j)

      dudeta      =  (w(2,i-1,j+1) +w(2,i,j+1)&
                     -w(2,i-1,j-1) -w(2,i,j-1))*0.25d0
      dvdeta      =  (w(3,i-1,j+1) +w(3,i,j+1)&
                     -w(3,i-1,j-1) -w(3,i,j-1))*0.25d0
      dTdeta      = (Tdeg(i-1,j+1)+Tdeg(i,j+1)&
                    -Tdeg(i-1,j-1)-Tdeg(i,j-1))*0.25d0
      dTurenedeta =  (w(5,i-1,j+1) +w(5,i,j+1)&
                     -w(5,i-1,j-1) -w(5,i,j-1))*0.25d0
      ddissipdeta =  (w(6,i-1,j+1) +w(6,i,j+1)&
                     -w(6,i-1,j-1) -w(6,i,j-1))*0.25d0

      dudx     =     dudxi*geojaci(1,1,i,j)+     dudeta*geojaci(2,1,i,j)
      dvdx     =     dvdxi*geojaci(1,1,i,j)+     dvdeta*geojaci(2,1,i,j)
      dTdx     =     dTdxi*geojaci(1,1,i,j)+     dTdeta*geojaci(2,1,i,j)
      dTurenedx=dTurenedxi*geojaci(1,1,i,j)+dTurenedeta*geojaci(2,1,i,j)
      ddissipdx=ddissipdxi*geojaci(1,1,i,j)+ddissipdeta*geojaci(2,1,i,j)

      dudr     =     dudxi*geojaci(1,2,i,j)+     dudeta*geojaci(2,2,i,j)
      dvdr     =     dvdxi*geojaci(1,2,i,j)+     dvdeta*geojaci(2,2,i,j)
      dTdr     =     dTdxi*geojaci(1,2,i,j)+     dTdeta*geojaci(2,2,i,j)
      dTurenedr=dTurenedxi*geojaci(1,2,i,j)+dTurenedeta*geojaci(2,2,i,j)
      ddissipdr=ddissipdxi*geojaci(1,2,i,j)+ddissipdeta*geojaci(2,2,i,j)

      divu=dudx+dvdr

      tau_xx=-2d0/3d0*mu*divu+2d0*mu*dudx
      tau_rr=-2d0/3d0*mu*divu+2d0*mu*dvdr
      tau_xr= mu*(dvdx+dudr)
      qx  = -k_heat*dTdx
      qr  = -k_heat*dTdr
   
      temp3 = (mu+mu_t/sigma_k)*dTurenedx
      temp4 = (mu+mu_t/sigma_k)*dTurenedr
      
      R_stress1 = (2d0/3d0*mu_t*divu+2d0*mu_t*dudx-2d0/3d0*rho*Turene)
      R_stress2 = (2d0/3d0*mu_t*divu+2d0*mu_t*dvdr-2d0/3d0*rho*Turene)
      R_stress3 = mu_t*(dudr+dvdx)

      !if(i==1.and.j==1)then
      !   print *,nu,k_heat,u_eta
      !   print *,Rt,non_y,f_mu
      !   print *,f_eta,nu_t
      !   print *,dissip,Turene
      !end if


      temp0 = mu+mu_t/sigma_k
      temp1 = mu+mu_t/sigma_eta
      temp2 = C_eta1*dissip/Turene

      Ev(1)=0d0
      Ev(2)=tau_xx+R_stress1
      Ev(3)=tau_xr+R_stress3
      Ev(4)=u*tau_xx+v*tau_xr-qx+temp3
      Ev(5)=temp0*dTurenedx
      Ev(6)=temp1*ddissipdx

      Fv(1)=0d0
      Fv(2)=tau_xr+R_stress3
      Fv(3)=tau_rr+R_stress2
      Fv(4)=u*tau_xr+v*tau_rr-qr+temp4
      Fv(5)=temp0*dTurenedr
      Fv(6)=temp1*ddissipdr

      vis_i(:,i,j)=Ev*nvi(1,i,j)+Fv*nvi(2,i,j)

      temp5 = R_stress1*dudx+R_stress3*dvdx
      temp6 = R_stress3*dudr+R_stress2*dvdr
      !temp7 = R_stress1*dudx+R_stress3*dvdx
      !temp8 = R_stress3*dudr+R_stress2*dvdr
      temp7 = temp5*nvi(1,i,j)+temp6*nvi(2,i,j)
      Sterm(5,i,j)=temp7-0.5d0*rho*dissip
                    
      Sterm(6,i,j)=temp7*rho*temp2-0.5d0*rho*C_eta2*f_eta*dissip**2/Turene
              
              
   end do
end do
!$omp end parallel do
!}}}

!set TGvj!{{{
!$omp parallel do shared(vis_j,j),&
!$omp             private(rho,u,v,Temper,&
!$omp             k_heat,dudxi,dvdxi,dTdxi,&
!$omp             dudeta,dvdeta,dTdeta,&
!$omp             dudx,dvdx,dTdx,dudr,dvdr,dTdr,&
!$omp             divu,tau_xx,tau_xr,tau_rr,qx,qr,&
!$omp             Ev,Fv,i)
do j=1,nj
   do i=1,ni-1

      rho    = (w(1    ,i,j-1)+w(1    ,i,j))*0.5d0
      u      = (w(2    ,i,j-1)+w(2    ,i,j))*0.5d0
      v      = (w(3    ,i,j-1)+w(3    ,i,j))*0.5d0
      Turene = (w(5    ,i,j-1)+w(5    ,i,j))*0.5d0
      dissip = (w(6    ,i,j-1)+w(6    ,i,j))*0.5d0
      Temper =    (Tdeg(i,j-1)+   Tdeg(i,j))*0.5d0

      nu    = mu/rho
      k_heat= mu*gamma/(gamma-1d0)*gas_specific/Prentl               !k from mu by Prandtl number
      u_eta = (nu*dissip)**(1d0/4d0)
      Rt    = Turene**2/(nu*dissip)
      non_y = u_eta*length/nu

      f_mu  = (1d0-exp(-non_y/14d0 ))**2*(1d0+(5d0/(Rt**(3d0/4d0)))*exp(-(Rt/200d0)**2))
      f_eta = (1d0-exp(-non_y/3.1d0))**2*(1d0-                0.3d0*exp(-(Rt/6.5d0)**2))

      nu_t     = C_mu*f_mu*Turene**2/dissip
      mu_t     = nu_t*rho


      dudxi      =  (w(2,i+1,j) +w(2,i+1,j-1)&
                    -w(2,i-1,j) -w(2,i-1,j-1))*0.25d0
      dvdxi      =  (w(3,i+1,j) +w(3,i+1,j-1)&
                    -w(3,i-1,j) -w(3,i-1,j-1))*0.25d0
      dTdxi      = (Tdeg(i+1,j)+Tdeg(i+1,j-1)&
                   -Tdeg(i-1,j)-Tdeg(i-1,j-1))*0.25d0
      dTurenedxi =  (w(5,i+1,j) +w(5,i+1,j-1)&
                    -w(5,i-1,j) -w(5,i-1,j-1))*0.25d0
      ddissipdxi =  (w(6,i+1,j) +w(6,i+1,j-1)&
                    -w(6,i-1,j) -w(6,i-1,j-1))*0.25d0


      dudeta       =  w(2,i,j) -w(2,i,j-1)
      dvdeta       =  w(3,i,j) -w(3,i,j-1)
      dTdeta       = Tdeg(i,j)-Tdeg(i,j-1)
      dTurenedeta  =  w(5,i,j) -w(5,i,j-1)
      ddissipdeta  =  w(6,i,j) -w(6,i,j-1)


      dudx     =     dudxi*geojacj(1,1,i,j)+     dudeta*geojacj(2,1,i,j)
      dvdx     =     dvdxi*geojacj(1,1,i,j)+     dvdeta*geojacj(2,1,i,j)
      dTdx     =     dTdxi*geojacj(1,1,i,j)+     dTdeta*geojacj(2,1,i,j)
      dTurenedx=dTurenedxi*geojacj(1,1,i,j)+dTurenedeta*geojacj(2,1,i,j)
      ddissipdx=ddissipdxi*geojacj(1,1,i,j)+ddissipdeta*geojacj(2,1,i,j)

      dudr=          dudxi*geojacj(1,2,i,j)+     dudeta*geojacj(2,2,i,j)
      dvdr=          dvdxi*geojacj(1,2,i,j)+     dvdeta*geojacj(2,2,i,j)
      dTdr=          dTdxi*geojacj(1,2,i,j)+     dTdeta*geojacj(2,2,i,j)
      dTurenedr=dTurenedxi*geojacj(1,2,i,j)+dTurenedeta*geojacj(2,2,i,j)
      ddissipdr=ddissipdxi*geojacj(1,2,i,j)+ddissipdeta*geojacj(2,2,i,j)

      !!2 dimentional-plane
      divu=dudx+dvdr

      tau_xx=-2d0/3d0*mu*divu+2d0*mu*dudx
      tau_rr=-2d0/3d0*mu*divu+2d0*mu*dvdr
      tau_xr= mu*(dvdx+dudr)
      qx  = -k_heat*dTdx
      qr  = -k_heat*dTdr

      temp3 = (mu+mu_t/sigma_k)*dTurenedx
      temp4 = (mu+mu_t/sigma_k)*dTurenedr

      R_stress1 = (2d0/3d0*mu_t*divu+2d0*mu_t*dudx-2d0/3d0*rho*Turene)
      R_stress2 = (2d0/3d0*mu_t*divu+2d0*mu_t*dvdr-2d0/3d0*rho*Turene)
      R_stress3 = mu_t*(dudr+dvdx)

      temp0 = mu+mu_t/sigma_k
      temp1 = mu+mu_t/sigma_eta
      temp2 = C_eta1*dissip/Turene

      Ev(1)=0d0
      Ev(2)=tau_xx+R_stress1
      Ev(3)=tau_xr+R_stress3
      Ev(4)=u*tau_xx+v*tau_xr-qx+temp3
      Ev(5)=temp0*dTurenedx
      Ev(6)=temp1*ddissipdx

      Fv(1)=0d0
      Fv(2)=tau_xr+R_stress3
      Fv(3)=tau_rr+R_stress2
      Fv(4)=u*tau_xr+v*tau_rr-qr+temp4
      Fv(5)=temp0*dTurenedr
      Fv(6)=temp1*ddissipdr

      vis_j(:,i,j)=Ev*nvj(1,i,j)+Fv*nvj(2,i,j)

      temp5 = R_stress1*dudx+R_stress3*dvdx
      temp6 = R_stress3*dudr+R_stress2*dvdr
      temp7 = temp5*nvj(1,i,j)+temp6*nvj(2,i,j)

      Sterm(5,i,j)= Sterm(5,i,j)+temp7-0.5d0*rho*dissip
                                 
      Sterm(6,i,j)= Sterm(6,i,j)+temp7*rho*temp2-0.5d0*rho*C_eta2*f_eta*dissip**2/Turene
   end do
end do
!$omp end parallel do
!}}}
!Sterm=Sterm*0.5d0

end subroutine set_viscous
