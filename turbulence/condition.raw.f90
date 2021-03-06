subroutine set_IC
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
   !$omp parallel do shared(w,j), private(i)
   do j = 0, nj
      do i = 0, ni
            w(1,i,j)=1.d0
            w(2,i,j)=100.d0
            w(3,i,j)=0.d0
            w(4,i,j)=1.d5
            w(5,i,j)=1.d0
            w(6,i,j)=1.d0
      enddo
   enddo
   !$omp end parallel do
end subroutine set_IC
subroutine set_BC
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
   double precision sonic
   !Bottom Ceiling Boundary
   !$omp parallel do shared(w), private(i)
   do i = 1, ni-1
      temp0 = mu/w(1,i,1)
      !Bottom
      !if(i<100)then
      !   w(:,i,0   )= w(:,i,1   )
      !   w(2,i,0   )= w(2,i,1   )-2d0*(w(2,i,   1)*nvj(1,i,   1)+w(3,i,   1)*nvj(2,i,   1))*nvj(1,i,   1)
      !   w(3,i,0   )= w(3,i,1   )-2d0*(w(2,i,   1)*nvj(1,i,   1)+w(3,i,   1)*nvj(2,i,   1))*nvj(2,i,   1)
      !else
         w(:,i,0   )= w(:,i,1   )
         w(2,i,0   )=-w(2,i,1   )
         w(3,i,0   )=-w(3,i,1   )
         w(5,i,0   )= 0d0
         w(6,i,0   )= 2d0*temp0*w(5,i,0)/(1d0)**2
      !end if

      !Ceiling
      w(:,i,nj  )= w(:,i,nj-1)
      !w(2,i,nj  )= w(2,i,nj-1)-2d0*(w(2,i,nj-1)*nvj(1,i,nj-1)+w(3,i,nj-1)*nvj(2,i,nj-1))*nvj(1,i,nj-1)
      !w(3,i,nj  )= w(3,i,nj-1)-2d0*(w(2,i,nj-1)*nvj(1,i,nj-1)+w(3,i,nj-1)*nvj(2,i,nj-1))*nvj(2,i,nj-1)
      !w(1,i,nj)=1.d0
      !w(2,i,nj)=100.d0
      !w(3,i,nj)=0.d0
      !w(4,i,nj)=w(4,i,nj-1)
      !w(5,i,nj)=1.d-10
      !w(6,i,nj)=1.d-10
      if(i<89)then
         temp0 =  mu/w(1,i,200)
         w(:,i,199)= w(:,i,200)
         w(2,i,199)=-w(2,i,200)
         w(3,i,199)=-w(3,i,200)
         w(5,i,199)= 0d0
         w(6,i,199)= 2d0*temp0*w(5,i,200)/(1d0)**2
      end if
   enddo                  
   !$omp end parallel do
   !Left Right Boundary
   !Left
   !$omp parallel do shared(w), private(j)
   do j = 0, nj
      w(1,0  ,j)=1.d0
      w(2,0  ,j)=100.d0
      w(3,0  ,j)=0.d0
      w(4,0  ,j)=w(4,1,j)
      w(5,0  ,j)=1.d0
      w(6,0  ,j)=1.d0
      !w(:,0  ,j)=w(:,1,j)
   !$omp end parallel do
   !Right
   !$omp parallel do shared(w), private(j)
         w(:,ni  ,j  )= w(:,ni-1,j)
         w(4,ni  ,j  )= 1.d5
      if(j<200)then
         temp0 = mu/w(1,i,89)
         w(:,88,j)= w(:,89,j)
         w(2,88,j)=-w(2,89,j)
         w(3,88,j)=-w(3,89,j)
         w(5,88,j)= 0d0
         w(6,88,j)= 2d0*temp0*w(5,89,j)/(1d0)**2
      end if
   enddo
   !$omp end parallel do
end subroutine set_BC
subroutine set_w
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
   wa=w
   !$omp parallel do private(i,temp0)
   do j=1,nj-1
      do i=1,ni-1
         w(1,i,j)=q(1,i,j)
         temp0=1d0/w(1,i,j)
         w(2,i,j)=q(2,i,j)*temp0
         w(3,i,j)=q(3,i,j)*temp0
         w(4,i,j)=(gamma-1.d0)*(q(4,i,j)-0.5d0*w(1,i,j)*(w(2,i,j)**2+w(3,i,j)**2))
         w(5,i,j)=q(5,i,j)*temp0
         w(6,i,j)=q(6,i,j)*temp0
         if(w(4,i,j)<0.d0)then
            write(*,*) i, j, t
            write(*,'(4es15.7)') w(:,i,j)
            call exit(1)
         endif
      enddo
   enddo
   !$omp end parallel do
   wd = w-wa
end subroutine set_w
subroutine set_conservative_variable_vector
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
   !$omp parallel do private(i)
   do j=1,nj-1
      do i=1,ni-1
         q(1,i,j)=w(1,i,j)
         q(2,i,j)=w(2,i,j)*q(1,i,j)
         q(3,i,j)=w(3,i,j)*q(1,i,j)
         q(4,i,j)=w(4,i,j)/(gamma-1.d0)+0.5d0*w(1,i,j)*(w(2,i,j)**2+w(3,i,j)**2)
         q(5,i,j)=w(5,i,j)*q(1,i,j)
         q(6,i,j)=w(6,i,j)*q(1,i,j)
      enddo
   enddo
   !$omp end parallel do
end subroutine set_conservative_variable_vector
subroutine set_dt
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision rho,u,v,p,un,a,uds
   double precision temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
   temp0=1.d300
   !$omp parallel do reduction(min:temp0) private(i,temp1,temp2,temp3)
   do j=1,nj-1
      do i=1,ni-1

         temp3=sqrt(gamma*w(4,i,j)/w(1,i,j))
         temp1=abs(-nvj(2,i,j)*w(2,i,j)+nvj(1,i,j)*w(3,i,j))
         temp2=abs(-nvi(2,i,j)*w(2,i,j)+nvi(1,i,j)*w(3,i,j))
         temp1=dsi(i,j)/(temp1+temp3)
         temp2=dsj(i,j)/(temp2+temp3)
         temp0=min(temp0,temp1,temp2)

      enddo
   enddo
   !$omp end parallel do
   dt(:,:)=CFL*temp0
   !dt(:,:)=temp0
end subroutine set_dt
!subroutine set_separation
!   use prmtr
!   use variable
!end subroutine set_separation
