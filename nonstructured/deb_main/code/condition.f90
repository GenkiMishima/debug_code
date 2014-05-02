subroutine set_IC
   use glbl_prmtr
   use form_format
   use grid_prmtr
   implicit none
   integer i,j
   !$omp parallel do shared(w,j), private(i)
   do i = 1,cell_counter-1
      if(i<10)then
         cells(i)%Val%w(1)=5.d0
         cells(i)%Val%w(2)=5.d0
         cells(i)%Val%w(3)=5.d0
         cells(i)%Val%w(4)=5.d0
      end if
   enddo
   !$omp end parallel do
end subroutine set_IC
subroutine set_BC_INPUT(boundary_condition,v)
   use glbl_prmtr
   use form_format
   use grid_prmtr
   implicit none
   integer i,j,v
   double precision,dimension(4)::boundary_condition

   boundary_condition(1)=1d2
   boundary_condition(2)=1d2
   boundary_condition(3)=1d2
   boundary_condition(4)=cells(v)%Val%w(4)

end subroutine set_BC_INPUT
subroutine set_BC_OUTPUT(boundary_condition,v)
   use glbl_prmtr
   use form_format
   use grid_prmtr
   implicit none
   integer i,j,v
   double precision,dimension(4)::boundary_condition

   boundary_condition(:)=cells(v)%Val%w(:)
   boundary_condition(4)=1d2

end subroutine set_BC_OUTPUT
subroutine set_BC_SLIP(boundary_condition,v,i)
   use glbl_prmtr
   use form_format
   use grid_prmtr
   implicit none
   integer i,j,v
   double precision,dimension(4)::boundary_condition

   boundary_condition(:)=cells(v)%Val%w(:)
   boundary_condition(2)=cells(v)%Val%w(2)-2d0*(cells(v)%Val%w(2)*lines(i)%Val%nv(1)&
                           +cells(v)%Val%w(3)*lines(i)%Val%nv(2))*lines(i)%Val%nv(1)
   boundary_condition(3)=cells(v)%Val%w(3)-2d0*(cells(v)%Val%w(2)*lines(i)%Val%nv(1)&
                           +cells(v)%Val%w(3)*lines(i)%Val%nv(2))*lines(i)%Val%nv(2)

end subroutine set_BC_SLIP
subroutine set_BC
   use glbl_prmtr
   use form_format
   use grid_prmtr
   implicit none
   integer i,j
   integer temp_int1,temp_int2,temp_int3
   double precision sonic

   do i=1,boundary_counter

   if(lines(i)%Val%bc.eq.8)then
      temp_int1=lines(i)%Val%c1
      temp_int2=lines(i)%Val%c2
      call set_BC_INPUT(cells(temp_int1)%Val%w,temp_int2)
   end if
   if(lines(i)%Val%bc.eq.9)then
      temp_int1=lines(i)%Val%c1
      temp_int2=lines(i)%Val%c2
      call set_BC_OUTPUT(cells(temp_int1)%Val%w,temp_int2)
   end if
   if(lines(i)%Val%bc.eq.10)then
      temp_int1=lines(i)%Val%c1
      temp_int2=lines(i)%Val%c2
      call set_BC_SLIP(cells(temp_int1)%Val%w,temp_int2,i)
   end if
   if(lines(i)%Val%bc.eq.11)then
      temp_int1=lines(i)%Val%c1
      temp_int2=lines(i)%Val%c2
      call set_BC_SLIP(cells(temp_int1)%Val%w,temp_int2,i)
   end if

   end do
      
end subroutine set_BC
subroutine set_w
   use glbl_prmtr
   use form_format
   use grid_prmtr
   implicit none
   integer i,j
   double precision temp0,temp1,temp2
   !$omp parallel do private(i,temp0)
   do i=1,cell_counter-1
      cells(i)%Val%w(1)=cells(i)%Val%q(1)
      temp0=1d0/cells(i)%Val%w(1)
      cells(i)%Val%w(2)=cells(i)%Val%q(2)*temp0
      cells(i)%Val%w(3)=cells(i)%Val%q(3)*temp0
      cells(i)%Val%w(4)=(gamma-1.d0)*(cells(i)%Val%q(4)-0.5d0*cells(i)%Val%w(1)*(cells(i)%Val%w(2)**2+cells(i)%Val%w(3)**2))
      if(cells(i)%Val%w(4)<0.d0)then
         write(*,*) i, t
         write(*,'(4es15.7)') cells(i)%Val%w(:)
         call exit(1)
      endif
   enddo
   !$omp end parallel do
end subroutine set_w
subroutine set_conservative_variable_vector
   use glbl_prmtr
   use form_format
   use grid_prmtr
   implicit none
   integer i,j
   !$omp parallel do private(i)
   do i=1,cell_counter-1
      cells(i)%Val%q(1)=cells(i)%Val%w(1)
      cells(i)%Val%q(2)=cells(i)%Val%w(2)*cells(i)%Val%q(1)
      cells(i)%Val%q(3)=cells(i)%Val%w(3)*cells(i)%Val%q(1)
      cells(i)%Val%q(4)=cells(i)%Val%w(4)/(gamma-1.d0)+0.5d0*cells(i)%Val%w(1)*(cells(i)%Val%w(2)**2+cells(i)%Val%w(3)**2)
   enddo
   !$omp end parallel do
end subroutine set_conservative_variable_vector
subroutine set_dt
   use glbl_prmtr
   use form_format
   use grid_prmtr
   implicit none
   integer i,j
   integer temp_int1,temp_int2,temp_int3
   double precision temp0,temp1,temp2,temp3
   double precision rho,u,v,p,un,a,uds
   temp0=1.d300
   !$omp parallel do reduction(min:temp0) private(i,temp1,temp2,temp3)
   do i=1,cell_counter-1
      temp_int1=cells(i)%Val%t_b(1)
      temp_int2=cells(i)%Val%t_b(2)
      !temp_int3=cells(i)%Val%t_b(3)
      print *,temp_int1,temp_int2

      temp3=sqrt(gamma*cells(i)%Val%w(4)/cells(i)%Val%w(1))

      !Boundary1
      temp1=abs(-lines(temp_int1)%Val%nv(2)*cells(i)%Val%w(2)+lines(temp_int1)%Val%nv(1)*cells(i)%Val%w(3))
      temp1=lines(temp_int1)%Val%ds/(temp1+temp3)
      !temp2=dsj(i,j)/(temp2+temp3)
      temp0=min(temp0,temp1)
      !Boundary2
      temp1=abs(-lines(temp_int2)%Val%nv(2)*cells(i)%Val%w(2)+lines(temp_int2)%Val%nv(1)*cells(i)%Val%w(3))
      temp1=lines(temp_int2)%Val%ds/(temp1+temp3)
      !temp2=dsj(i,j)/(temp2+temp3)
      temp0=min(temp0,temp1)

   enddo
   !$omp end parallel do
   dt=CFL*temp0
   !dt(:,:)=temp0
end subroutine set_dt
