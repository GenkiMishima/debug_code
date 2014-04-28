subroutine set_jacobian
   use prmtr
   use variable
   implicit none
   integer i,j,m,n
   double precision rho, u, v, p, E, H, c
   double precision nuA, nuB, unA, unB
   double precision, parameter :: gamma_tilda=gamma-1.d0
   double precision, dimension(4,4) ::  A, B, Atilde, Btilde
   !$omp parallel do private(i,m,rho,u,v,p,E,H,c,&
   !$omp                     A,B)
   do j=0,nj
      do i=0,ni
         rho= w(1,i,j)
         u  = w(2,i,j)
         v  = w(3,i,j)
         p  = w(4,i,j)
         E  = p/gamma_tilda+0.5d0*rho*(u**2+v**2)
         H  = (E+p)/rho
         c  = sqrt(gamma*p/rho)

         unA=vnci(1,i,j)*u+vnci(2,i,j)*v
         unB=vncj(1,i,j)*u+vncj(2,i,j)*v
         nuA=abs(u)+c
         nuB=abs(v)+c
                            

         A(1,1)= 0.d0 
         A(2,1)= 0.5d0*(u**2+v**2)*gamma_tilda-u**2 
         A(3,1)= -u*v
         A(4,1)= 0.5d0*u*(-2.d0*H+(u**2+v**2)*gamma_tilda)

         A(1,2)= 1.d0
         A(2,2)= u*(3.d0-gamma)
         A(3,2)= v
         A(4,2)= H-u**2*gamma_tilda

         A(1,3)= 0.d0
         A(2,3)= -v*gamma_tilda
         A(3,3)= u
         A(4,3)= -u*v*gamma_tilda

         A(1,4)= 0.d0
         A(2,4)= gamma_tilda
         A(3,4)= 0.d0
         A(4,4)= u*gamma

         B(1,1)= 0.d0 
         B(2,1)= -u*v
         B(3,1)= 0.5d0*(u**2+v**2)*gamma_tilda-v**2 
         B(4,1)= 0.5d0*v*(-2.d0*H+(u**2+v**2)*gamma_tilda)

         B(1,2)= 0.d0
         B(2,2)= v
         B(3,2)= -u*gamma_tilda
         B(4,2)= -u*v*gamma_tilda

         B(1,3)= 1.d0
         B(2,3)= u
         B(3,3)= v*(3.d0-gamma)
         B(4,3)= H-v**2*gamma_tilda

         B(1,4)= 0.d0
         B(2,4)= 0.d0
         B(3,4)= gamma_tilda
         B(4,4)= v*gamma

         !set Atilde, Apm
         Atilde = vnci(1,i,j)*A+vnci(2,i,j)*B

         Ap(:,:,i,j)= Atilde*0.5d0
         Am(:,:,i,j)= Atilde*0.5d0

         do m=1,4
            Ap(m,m,i,j)=Ap(m,m,i,j)+nuA*0.5d0
            Am(m,m,i,j)=Am(m,m,i,j)-nuA*0.5d0
         end do

         !set Btilde, Bpm
         Btilde= vncj(1,i,j)*A+vncj(2,i,j)*B

         Bp(:,:,i,j)= Btilde*0.5d0
         Bm(:,:,i,j)= Btilde*0.5d0

         do m=1,4
            Bp(m,m,i,j)=Bp(m,m,i,j)+nuB*0.5d0
            Bm(m,m,i,j)=Bm(m,m,i,j)-nuB*0.5d0
         end do
         !set alpha
         temp0= Vol(i,j)/dt(i,j) + dsci(i,j)*nuA+dscj(i,j)*nuB
         alpha(i,j)=1d0/temp0
 

         !alpha(i,j)=1.d0+dt(i,j)/dsi(i,j)*nuA+dt(i,j)/dsj(i,j)*nuB
         !print *,alpha(i,j)
         !if(i==50.and.j==30)write(*,'(4es15.7)') nuA,nuB
         !if(i==50.and.j==30)write(*,'(4es15.7)')A(:,:)
         !if(i==50.and.j==30)write(*,'(4es15.7)')
         !if(i==50.and.j==30)write(*,'(4es15.7)')B(:,:)
      end do
   end do
   !$omp end parallel do
   !write(*,'(4es15.7)')
   !write(*,'(4es15.7)')Ap(:,:,50,30)
   !write(*,'(4es15.7)')
   !write(*,'(4es15.7)')Bp(:,:,50,30)
!   write(*,'(4es15.7)') Ap(:,:,50,41)
end subroutine set_jacobian
subroutine calc_next_step_imp
   use prmtr
   use variable
   implicit none
   integer i,j,m,n
   double precision, dimension(4) :: temp_lu
   double precision, dimension(4,0:ni  ,0:nj  )::q_plime,q_dplime
   double precision, dimension(4,4,0:ni  ,0:nj  )::Z,inv_Z

   do j=0,nj
      do i=0,ni
         do m=1,4
            Z(m,m,i,j)=1d0+dt(i,j)/dsi(i,j)*(Ap(m,m,i,j)-Am(m,m,i,j))+dt(i,j)/dsj(i,j)*(Bp(m,m,i,j)-Bm(m,m,i,j))
         end do
      end do
   end do
   inv_z=z**-1
   print *,z(:,:,10,20),inv_z(:,:,10,20)

   !$omp parallel do private(i,temp0)
   do j=1,nj-1
      do i=1,ni-1
         RHS(:,i,j)=dt(i,j)/area(i,j)*( dsj(i,j)*(X_Numerical(:,i  ,j  )-vis_i(:,i  ,j  ))&
                                       -dsj(i,j)*(X_Numerical(:,i+1,j  )-vis_i(:,i+1,j  ))&
                                       +dsi(i,j)*(Y_Numerical(:,i  ,j  )-vis_j(:,i  ,j  ))&
                                       -dsi(i,j)*(Y_Numerical(:,i  ,j+1)-vis_j(:,i  ,j+1)))
         !print *,RHS(:,i,j)
      end do
   end do
   !$omp end parallel do

   q_plime(:,0,:)=0.d0
   q_plime(:,:,0)=0.d0
   !$omp parallel do private(i,temp0)
   do j=1,nj-1
      do i=1,ni-1
         temp_lu=RHS(:,i,j)&
                +alpha(i,j)*dt(i,j)/dsi(i,j)*matmul(Ap(:,:,i-1,j),q_plime(:,i-1,j))&
                +alpha(i,j)*dt(i,j)/dsj(i,j)*matmul(Bp(:,:,i,j-1),q_plime(:,i,j-1))

         q_plime(:,i,j)=matmul(inv_Z(:,:,i,j),temp_lu(:))
                        
                        
      end do
   end do
   !$omp end parallel do

   do j=1,nj-1
      do i=1,ni-1
         q_dplime(:,i,j)=matmul(Z(:,:,i,j),q_plime(:,i,j))
      end do
   end do
                        
   q_imp(:,ni,:)=0.d0
   q_imp(:,:,nj)=0.d0
   !$omp parallel do private(i,m,temp0)
   do j=1,nj-1
      n=nj-j
      do i=1,ni-1
         m=ni-i
         temp_lu=q_dplime(:,m,n)&
              -alpha(m,n)*dt(i,j)/dsi(m,n)*matmul(Ap(:,:,m+1,n),q_imp(:,m+1,n))&
              -alpha(m,n)*dt(i,j)/dsj(m,n)*matmul(Bp(:,:,m,n+1),q_imp(:,m,n+1))
         q_imp(:,m,n)=matmul(inv_Z(:,:,i,j),temp_lu(:))
                      
                      
      end do
   end do
   !$omp end parallel do

   !$omp parallel do private(i,temp0)
   do j=1,nj-1
      do i=1,ni-1
               q(:,i,j)=q(:,i,j)+q_imp(:,i,j)
               !print *,q_imp
      end do
   end do
   !$omp end parallel do
end subroutine calc_next_step_imp

