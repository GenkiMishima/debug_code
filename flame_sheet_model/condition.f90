subroutine set_IC
   use prmtr
   use variable
   implicit none
   integer i,j
   !$omp parallel do shared(w,j), private(i)
   do j = 0, nj
      do i = 0, ni
            w(1,i,j)=1.d0
            w(2,i,j)=0.d0
            w(3,i,j)=0.d0
            w(4,i,j)=1.d5
            w(5,i,j)=0.d0
      enddo
   enddo
   Chem_num(1,:,:)%moler_weight=4
   Chem_num(2,:,:)%moler_weight=16
   Chem_num(3,:,:)%moler_weight=11
   Chem_num(4,:,:)%moler_weight=9
   Chem_num(5,:,:)%moler_weight=15
   !$omp end parallel do
end subroutine set_IC
subroutine set_Chem
use prmtr
use variable
implicit none
double precision f_st
!print *,w(5,:,:)
f_st=

Chem_num(1,:,:)%mass_rate=
Chem_num(2,:,:)%mass_rate=
Chem_num(3,:,:)%mass_rate=
Chem_num(4,:,:)%mass_rate=
Chem_num(5,:,:)%mass_rate=

Chem_num(1,:,:)%rho=Chem_num(1,:,:)%mass_rate*w(1,:,:)
Chem_num(2,:,:)%rho=Chem_num(2,:,:)%mass_rate*w(1,:,:)
Chem_num(3,:,:)%rho=Chem_num(3,:,:)%mass_rate*w(1,:,:)
Chem_num(4,:,:)%rho=Chem_num(4,:,:)%mass_rate*w(1,:,:)
Chem_num(5,:,:)%rho=Chem_num(5,:,:)%mass_rate*w(1,:,:)
end subroutine set_Chem
subroutine set_BC
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision sonic
   !Bottom Ceiling Boundary
   !$omp parallel do shared(w), private(i)
   do i = 1, ni-1
      !Bottom
      if(i<55.and.i>105)then
         w(:,i,0   )= w(:,i,1   )
         w(4,i,0   )= 1d5
      else
         w(:,i,0   )= w(:,i,1   )
         w(2,i,0   )= w(2,i,1   )-2d0*(w(2,i,   1)*nvj(1,i,   1)+w(3,i,   1)*nvj(2,i,   1))*nvj(1,i,   1)
         w(3,i,0   )= w(3,i,1   )-2d0*(w(2,i,   1)*nvj(1,i,   1)+w(3,i,   1)*nvj(2,i,   1))*nvj(2,i,   1)

         !w(:,i,0   )= w(:,i,1   )
         !w(2,i,0   )=-w(2,i,1   )
         !w(3,i,0   )=-w(3,i,1   )
      end if

      !Ceiling
      if(i<55.and.i>105)then
         w(:,i,nj   )= w(:,i,nj-1)
         w(4,i,nj   )= 1d5
      else
         w(:,i,nj  )= w(:,i,nj-1)
         w(2,i,nj  )= w(2,i,nj-1)-2d0*(w(2,i,nj-1)*nvj(1,i,nj-1)+w(3,i,nj-1)*nvj(2,i,nj-1))*nvj(1,i,nj-1)
         w(3,i,nj  )= w(3,i,nj-1)-2d0*(w(2,i,nj-1)*nvj(1,i,nj-1)+w(3,i,nj-1)*nvj(2,i,nj-1))*nvj(2,i,nj-1)
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
      w(5,0  ,j)=0.5d0
   enddo
   !$omp end parallel do
   !Right
   !$omp parallel do shared(w), private(j)
   do j = 0, nj
      w(:,ni,j   )= w(:,ni-1,j)
      w(4,ni,j   )= 1d5
   enddo
   !$omp end parallel do
end subroutine set_BC
subroutine set_w
   use prmtr
   use variable
   implicit none
   integer i,j
   !$omp parallel do private(i,temp0)
   do j=1,nj-1
      do i=1,ni-1
         w(1,i,j)=q(1,i,j)
         temp0=1d0/w(1,i,j)
         w(2,i,j)=q(2,i,j)*temp0
         w(3,i,j)=q(3,i,j)*temp0
         w(4,i,j)=(gamma-1.d0)*(q(4,i,j)-0.5d0*w(1,i,j)*(w(2,i,j)**2+w(3,i,j)**2))
         w(5,i,j)=0.d0
         if(w(4,i,j)<0.d0)then
            write(*,*) i, j, t
            write(*,'(4es15.7)') w(:,i,j)
            call exit(1)
         endif
      enddo
   enddo
   !$omp end parallel do
end subroutine set_w
subroutine set_conservative_variable_vector
   use prmtr
   use variable
   implicit none
   integer i,j
   !$omp parallel do private(i)
   do j=1,nj-1
      do i=1,ni-1
         q(1,i,j)=w(1,i,j)
         q(2,i,j)=w(2,i,j)*q(1,i,j)
         q(3,i,j)=w(3,i,j)*q(1,i,j)
         q(4,i,j)=w(4,i,j)/(gamma-1.d0)+0.5d0*w(1,i,j)*(w(2,i,j)**2+w(3,i,j)**2)
         q(5,i,j)=0d0
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
