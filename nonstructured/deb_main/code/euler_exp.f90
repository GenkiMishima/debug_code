subroutine calc_next_step_exp
use glbl_prmtr
use form_format
use grid_prmtr
   implicit none
   integer i,j
   integer temp_int, temp_int1,temp_int2,temp_int3,temp_int4, bc
   double precision, dimension(4)::vis_term
   !temp0=10000.d0
   !temp_residual2=temp_residual1
   !temp1=0.d0
   !$omp parallel do private(i) shared(j,q)
   do i=1,cell_counter-1
      temp_int1=cells(i)%Val%t_b(1)
      temp_int2=cells(i)%Val%t_b(2)
      temp_int3=cells(i)%Val%t_b(3)
      cells(i)%Val%q(:)=cells(i)%Val%q(:)+dt/cells(i)%Val%area&
                       *(lines(temp_int1)%Val%ds*(lines(temp_int1)%Val%Flux(:))&
                        +lines(temp_int2)%Val%ds*(lines(temp_int2)%Val%Flux(:))&
                        +lines(temp_int3)%Val%ds*(lines(temp_int3)%Val%Flux(:)))
   enddo
   !$omp end parallel do

   !do j=1,nj-1
   !   do i=1,ni-1
   !      temp_residual1=temp_residual1+q(4,i,j) 
   !   enddo
   !enddo
   !!residual
   !residual=abs(temp_residual1-temp_residual2)

end subroutine calc_next_step_exp
