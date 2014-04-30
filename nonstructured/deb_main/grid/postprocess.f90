program main
use form_format
use prmtr
implicit none
integer i,j,time
integer temp_int, temp_int1,temp_int2,temp_int3,temp_int4, bc
double precision temp, temp_x, temp_r,temp1, temp2,temp3
double precision,dimension(1:node_num)::rho_mat,u_mat,v_mat,p_mat,temp_mat
   !integer*4,external::access
   call grid
   call set_out_format(rho_mat,u_mat,v_mat,p_mat)
   !do time=1,500
   open(66,file='result.vtk')
   call vtk_grid
   call vtk_scalar(p_mat,'pressure1')
   close(66)
end program main

subroutine vtk_grid
use form_format
use prmtr
   implicit none
   integer i,j
   write(66,'(a26)') '# vtk DataFile Version 2.0'
   write(66,'(a40)') '2D Unstructured Grid of Linear Triangles'
   write(66,'(a5)') 'ASCII'
   write(66,*)
   write(66,'(a25)') 'DATASET UNSTRUCTURED_GRID'
   write(66,'(a7,i10,a6)') 'POINTS ', int(node_num),' float'
   do i=1,node_num
      !write(66,'(e14.6,e14.6,e14.6)') nodes(i)%x,nodes(i)%r,0
      write(66,*) nodes(i)%x,nodes(i)%r,0d0
   enddo
   write(66,*)
   write(66,'(a6,i10,i10)') 'CELLS ',cell_counter-1,4*(cell_counter-1)
   do i=1,cell_counter-1
      write(66,'(i10,i10,i10,i10)'),3,cells(i)%t_l1-1,cells(i)%t_l2-1,cells(i)%t_l3-1
   enddo
   write(66,*)
   write(66,'(a11,i10)') 'CELL_TYPES ', cell_counter-1
   do i=1,cell_counter-1
      write(66,'(i10)'),5
   enddo
   write(66,*)
   write(66,'(a10,i10)') 'POINT_DATA ',int(node_num)
end subroutine vtk_grid
subroutine vtk_scalar(arr,title)
use prmtr
   implicit none
   double precision, dimension(1:node_num), intent(in)::arr
   character(*),intent(in)::title

   integer i,j

   write(66,'(a)') 'SCALARS '//trim(title)//' float'
   write(66, '(a20)') 'LOOKUP_TABLE default'
   do i = 1, node_num
      write(66,'(e14.6e3)') arr(i)
   enddo
   write(66,*)
end subroutine vtk_scalar
!subroutine vtk_vector(arr1,arr2,title)
!use prmtr
!   implicit none
!   double precision,dimension(node_num),intent(in)::arr1,arr2
!   character(*),intent(in)::title
!
!   integer i,j
!
!   write(66,'(a22)') 'VECTORS '//trim(title)//' float'
!   do j = 1, node_num
!      write(66,'(e14.6e3,2(1x,e14.6e3))') arr1(j), arr2(j), 0.d0
!   end do
!   write(66,*)
!end subroutine vtk_vector
subroutine set_out_format(rho,u,v,p)
use prmtr
use form_format
implicit none
integer i,j
integer temp_int1
double precision temp1
double precision,dimension(1:node_num)::rho,u,v,p
!double precision,dimension(1:   node_num)::nodes
!double precision,dimension(1:element_num)::cells
do i= 1,node_num
   temp1=0d0
   do j=1,nodes(i)%node_num_count
      temp_int1=nodes(i)%cell_num(j)
      temp1=temp1+cells(temp_int1)%value1
   end do
   p(i)=temp1/dble(nodes(i)%node_num_count)
   !print *,nodes(i)%value1,temp1
end do


end subroutine set_out_format
