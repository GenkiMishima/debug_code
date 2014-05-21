subroutine set_IC
   use prmtr
   use material
   use variable
   implicit none
   integer i,j
   !$omp parallel do shared(w,j), private(i)
   do j = 0, nj
      do i = 0, ni
            w(1,i,j)=1.d0
            w(2,i,j)=100.d0
            w(3,i,j)=0.d0
            w(4,i,j)=1.d5
      enddo
   enddo
   !$omp end parallel do
end subroutine set_IC
subroutine set_material
   use prmtr
   use material
   use variable
   double precision temp1,temp2,temp_P
   double precision real_a,real_b,R_gas,T_wall,T_inf
   double precision P0
   R_gas = gas_specific
   
   real_a = 27.d0*(R_gas*T_cri)**2/(64.d0*P_cri)
   real_b = R_gas*T_cri/(8.d0*P_cri)

   do j = 0,nj
      do i = 0,ni
         temp_P = w(4,i,j)
         temp1 = 3.d0*real_a*real_b*(real_b*temp_P+R_gas*T_cri)-real_a**2
         temp2 = 2.d0*real_a**3+18.d0*real_a**2*real_b**2*temp_P-9.d0*real_a**2*real_b*R_gas*T_cri
         w(1,i,j) = -1.d0*2.d0**(1.d0/3.d0)*temp1/(3.d0*real_a*real_b*(temp2+dsqrt(4.d0*temp1**3&
                    +temp2**2))**(1.d0/3.d0))+(temp2+dsqrt(4.d0*temp1**3+temp2**2))**(1.d0/3.d0)&
                    /(3.d0*2.d0**(1.d0/3.d0)*real_a*real_b)+1.d0/(3.d0*real_b)
         P0 = 11.d6
         T_wall = 150.d0
         T_inf = 600.d0
         nonW(5,i,j) = dsqrt((R_gas*(temp_P+real_a*w(1,i,j)**2)+Cv_st*(temp_P+real_a*w(1,i,j)**2&
                       *(-1.d0+2.d0*real_b*w(1,i,j))))/(Cv_st*w(1,i,j)*(1.d0-real_b*w(1,i,j))))
         temp2=(R_gas**(temp_P+real_a*w(1,i,j)**2)+Cv_st*(temp_P+real_a*w(1,i,j)**2*(-1.d0+2.d0*real_b*w(1,i,j))))
      end do
   end do
   
   temp_out1=nonW(5,30,50)
   temp_out2=w(1,30,50)
end subroutine set_material
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
      if(i<55)then
         w(:,i,0   )= w(:,i,1   )
         w(2,i,0   )= w(2,i,1   )-2d0*(w(2,i,   1)*nvj(1,i,   1)+w(3,i,   1)*nvj(2,i,   1))*nvj(1,i,   1)
         w(3,i,0   )= w(3,i,1   )-2d0*(w(2,i,   1)*nvj(1,i,   1)+w(3,i,   1)*nvj(2,i,   1))*nvj(2,i,   1)
      else
         w(:,i,0   )= w(:,i,1   )
         w(2,i,0   )=-w(2,i,1   )
         w(3,i,0   )=-w(3,i,1   )
      end if

      !Ceiling
      w(:,i,nj  )= w(:,i,nj-1)
      w(2,i,nj  )= w(2,i,nj-1)-2d0*(w(2,i,nj-1)*nvj(1,i,nj-1)+w(3,i,nj-1)*nvj(2,i,nj-1))*nvj(1,i,nj-1)
      w(3,i,nj  )= w(3,i,nj-1)-2d0*(w(2,i,nj-1)*nvj(1,i,nj-1)+w(3,i,nj-1)*nvj(2,i,nj-1))*nvj(2,i,nj-1)
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
   enddo
   !$omp end parallel do
   !Right
   !$omp parallel do shared(w), private(j)
   do j = 0, nj
         w(:,ni  ,j  )= w(:,ni-1,j)
         w(4,ni  ,j  )= 1.d5
   enddo
   !$omp end parallel do
end subroutine set_BC
subroutine set_w
   use prmtr
   use variable
   implicit none
   integer i,j
   double precision temp0, temp1, temp2, temp3
   !$omp parallel do private(i,temp0)
   do j=1,nj-1
      do i=1,ni-1
         w(1,i,j)=q(1,i,j)
         temp0=1d0/w(1,i,j)
         w(2,i,j)=q(2,i,j)*temp0
         w(3,i,j)=q(3,i,j)*temp0
         w(4,i,j)=(gamma-1.d0)*(q(4,i,j)-0.5d0*w(1,i,j)*(w(2,i,j)**2+w(3,i,j)**2))
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
   double precision temp0, temp1, temp2, temp3
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
