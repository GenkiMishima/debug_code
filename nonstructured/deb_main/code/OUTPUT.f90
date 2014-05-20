subroutine output
use glbl_prmtr
use form_format
use grid_prmtr
   implicit none
   integer i,j
   !if(time==1.or.(mod(time,out_time)== 0).or.temp_int==1) then
      !$omp parallel do private(i)
      !do j = 1, nj
      !   if(i==1.and.j==1)then
      !      wp(:,i,j)=w(:,i,j)
      !   elseif(i==ni.and.j==1)then
      !      wp(:,i,j)=w(:,i-1,j)
      !   elseif(i==1.and.j==nj)then
      !      wp(:,i,j)=w(:,i,j-1)
      !   elseif(i==ni.and.j==nj)then
      !      wp(:,i,j)=w(:,i-1,j-1)
      !   elseif(j==nj)then
      !      wp(:,i,j)=wp(:,i,j-1)
      !   elseif(j==0)then
      !      wp(:,i,j)=wp(:,i,j+1)
      !   !elseif((i>=step_forward.and.i<=step_backward).and.j<=step_height)then
      !   !   wp(:,i,j)=w(:,i,j)
      !   else
      !      wp(:,i,j)=0.25d0*(w(:,i-1,j-1)+w(:,i,j-1)+w(:,i-1,j)+w(:,i,j))
      !   endif
      !enddo
      !$omp end parallel do
      !$omp parallel do private(i)
      !do j=1,nj
      !   do i=1,ni
      !      Mach_number(i,j)=sqrt((wp(2,i,j)**2+wp(3,i,j)**2)*wp(1,i,j)/(gamma*wp(4,i,j)))
      !   enddo
      !enddo
      !$omp end parallel do

      write(tmpstring,'(i3.3)') int(time/out_time)
      open(10, file='data/density_'//trim(tmpstring)//'.d')
      open(11, file='data/U_Velocity_'//trim(tmpstring)//'.d')
      open(12, file='data/V_Velocity_'//trim(tmpstring)//'.d')
      open(13, file='data/pressure_'//trim(tmpstring)//'.d')
      !open(14, file='data/Mach_'//trim(tmpstring)//'.d')
      !open(15, file='data/Vectorx_'//trim(tmpstring)//'.d')
      !open(16, file='data/Vectorr_'//trim(tmpstring)//'.d')
!      !$omp parallel do private(i)
      do i=1,cell_counter-1
         write(10, *) cells(i)%Val%w(1)
         write(11, *) cells(i)%Val%w(2)
         write(12, *) cells(i)%Val%w(3)
         write(13, *) cells(i)%Val%w(4)
         !write(14, *) Mach_number(i,j)
         !write(15, *) nvj(1,i,j)
         !write(16, *) nvj(2,i,j)
         !y(1,j) = dy*0.5d0+dble(j-1)*dy
         !write(10, *) y(1,j), w(1,50,j)
         !write(11, *) y(1,j), w(2,50,j)
         !write(12, *) y(1,j), w(3,50,j)
         !write(13, *) y(1,j), w(4,50,j)
         !write(14, *) y(1,j), Mach_number(1,j)
      enddo
!      !$omp end parallel do
      close(10)
      close(11)
      close(12)
      close(13)
      !close(14)
      !close(15)
      !close(16)
      write(*, *) tmpString, t
      !open(15,file="restart.bin",form="unformatted")
      !   write(15) w
      !   write(15) time
      !   write(15) t
      !   write(15) dt
      !close(15)
      !counter1=0d0
      !do j=1,nj-1
      !   do i=1,ni-1
      !      temp0=q(4,i,j)-qq(4,i,j)
      !      temp1=temp0/q(4,i,j)
      !      if(abs(temp1)>0.01d0)then
      !      !print *,temp1
      !      counter1=counter1+1
      !      end if
      !   end do
      !end do
      !      print *,counter1
   !end if
end subroutine output
