!**************************************
!
!  READ THE .MSH FILES
!
!**************************************
subroutine grid
use glbl_prmtr
use form_format
use grid_prmtr
implicit none
integer i,j,k
integer temp_int, temp_int1,temp_int2,temp_int3,temp_int4, bc
double precision temp, temp_x, temp_r,temp1,temp2,temp3,temp4,delta_x,delta_r
integer node_temp1,node_temp2,node_temp3
integer, dimension(element_num) ::ele_type

   !Read Nodes{{{
   open(55, file='grid/test1.msh')
   !Header
   read(55,*) moji
   read(55,*) 
   read(55,*) moji
   !Nodes
   read(55,*) moji
   read(55,*) temp
   if(int(temp).ne.node_num)then
      print *,'node_num is ',temp
      call exit(1)
   end if
   do i = 1, node_num
      read(55,*) temp, nodes(i)%x,nodes(i)%r
   end do
   read(55,*) moji
   !Check Elements type
   read(55,*) moji
   read(55,*) temp
   if(temp.ne.element_num)then
      print *,'element_num is ',temp
      call exit(1)
   end if
   do i = 1, element_num
      read(55,*) temp, ele_type(i)
   end do
   close(55)
   !}}} 

   !Read elements{{{
   open(55, file='grid/test1.msh')
   !Header{{{
   read(55,*) 
   read(55,*) 
   read(55,*) 
   !Nodes
   read(55,*) 
   read(55,*) 
   do i = 1, node_num
      read(55,*) 
   end do
   read(55,*) moji
   !}}}
   !Elements
   read(55,*) moji
   read(55,*) temp
   if(temp.ne.element_num) call exit(2)
   !Clear_Counters
   line_counter=1
   cell_counter=1
   boundary_counter=0
   temp2=1d0/2d0
   temp3=1d0/3d0
   !boundary_counter=1
   !Read elements
   nodes%node_num_count=0
   do i = 1, element_num
      !Read points data{{{
      if(ele_type(i).eq.15)then
         read(55,*) 
      !}}}
      !Read pysical data{{{
      elseif(ele_type(i).eq.1)then
         boundary_counter=boundary_counter+1
         read(55,*) temp,temp,temp,bc,temp,node_temp1,node_temp2
         !print *, 'make',boundary_counter
         !lines(boundary_counter)%Val%line_number=boundary_counter
         lines(boundary_counter)%Val%bc=bc
         lines(boundary_counter)%Val%c1=1d6+bc
         lines(boundary_counter)%Val%centx=abs(nodes(node_temp1)%x+nodes(node_temp2)%x)*temp2
         lines(boundary_counter)%Val%centr=abs(nodes(node_temp1)%r+nodes(node_temp2)%r)*temp2
         delta_x=abs(nodes(node_temp1)%x-nodes(node_temp2)%x)
         delta_r=abs(nodes(node_temp1)%r-nodes(node_temp2)%r)
         temp1=sqrt(delta_x**2+delta_r**2)
         lines(boundary_counter)%Val%nv(1)= delta_r/temp1
         lines(boundary_counter)%Val%nv(2)=-delta_x/temp1
         lines(boundary_counter)%Val%n1=node_temp1
         lines(boundary_counter)%Val%n2=node_temp2
      !if(boundary_counter==9)then
      !do j = 1,boundary_counter
      !   write(*,'(i3,i3,i3,i3)') lines(j)%line_number,lines(j)%bc,lines(j)%c1,lines(j)%c2
      !end do
      !call exit(1)
      !end if
      !}}}
      !Read triangle data{{{
      elseif(ele_type(i).eq.2)then
         read(55,*) temp,temp,temp,temp,temp,node_temp1,node_temp2,node_temp3
         !cells(cell_counter)%Val%cell_number=cell_counter
         cells(cell_counter)%Val%t_l1=node_temp1
         cells(cell_counter)%Val%t_l2=node_temp2
         cells(cell_counter)%Val%t_l3=node_temp3
   
         temp_int1=nodes(node_temp1)%node_num_count
         temp_int2=nodes(node_temp2)%node_num_count
         temp_int3=nodes(node_temp3)%node_num_count
         
         nodes(node_temp1)%cell_num(temp_int1)=cell_counter
         nodes(node_temp2)%cell_num(temp_int2)=cell_counter
         nodes(node_temp3)%cell_num(temp_int3)=cell_counter
         nodes(node_temp1)%node_num_count=temp_int1+1
         nodes(node_temp2)%node_num_count=temp_int2+1
         nodes(node_temp3)%node_num_count=temp_int3+1

         cells(cell_counter)%Val%grax=(nodes(node_temp1)%x+nodes(node_temp2)%x&
                                 +nodes(node_temp3)%x)*temp3
         cells(cell_counter)%Val%grar=(nodes(node_temp1)%r+nodes(node_temp2)%r&
                                 +nodes(node_temp3)%r)*temp3

         cell_counter=cell_counter+1
      !}}}
      !Read else{{{
      else
         read(55,*)
      end if
      !}}}

   end do
   close(55)
   !}}}

   !!Print{{{
   !open(55,file="out.d")
   !!do i = 1,boundary_counter-1
   !!   write(55,'(6es15.7)') boundarys(i)
   !!end do
   !do i = 1,boundary_counter
   !   write(55,'(i3,i3,i3,i3,i3)') i, lines(i)%line_number,lines(i)%bc,lines(i)%c1,lines(i)%c2
   !end do
   !close(55)
   !}}}

!Cell
temp1=0d0
do i=1,cell_counter-1
   temp1=cells(i)%Val%grar
   do j = 1,i
      if(cells(j)%Val%grar>temp1)then
         temp_cell1=cells(j)%Val
         cells(j)%Val=cells(i)%Val
         !cells(i)%Val=temp_cell
         do k = j+1,i
            temp_cell2=cells(k)%Val
            cells(k)%Val=temp_cell1
            !cells(k-1)%Val=temp_cell
            temp_cell1=temp_cell2

            !cells(k)%Val%cell_number=cells(k)%Val%cell_number+1
         end do
         exit
      elseif(i==j) then
         exit
         !cells(i)%Val%cell_number=i
      end if
   end do
end do

temp2=1d0/2d0
temp3=1d0/3d0
cells(:)%Val%counter=0

do i=1,cell_counter-1

   node_temp1=cells(i)%Val%t_l1
   node_temp2=cells(i)%Val%t_l2
   node_temp3=cells(i)%Val%t_l3

   print *, 'node',node_temp1,node_temp2,node_temp3
   !{{{
   !!Middle point of Node1,2
   !delta_x=abs(nodes(node_temp1)%x-nodes(node_temp2)%x)
   !delta_r=abs(nodes(node_temp1)%r-nodes(node_temp2)%r)
   !temp1=sqrt(delta_x**2+delta_r**2)
   !temp_x=abs(nodes(node_temp1)%x+nodes(node_temp2)%x)*temp2
   !temp_r=abs(nodes(node_temp1)%r+nodes(node_temp2)%r)*temp2

   !!Ckeck that what elements face to this lines
   !do j = 1,boundary_counter
   !   temp_int3=lines(j)%Val%n1
   !   temp_int4=lines(j)%Val%n2
   !   temp_int=0
   !   !temp_int1=cells(i)%Val%counter
   !   if((node_temp1.eq.temp_int3.and.node_temp2.eq.temp_int4).or.(node_temp1.eq.temp_int4.and.node_temp2.eq.temp_int3)) then
   !      temp_int=1
   !      lines(j)%Val%c2=i
   !      cells(i)%Val%counter=cells(i)%Val%counter+1
   !      if(cells(i)%Val%counter.eq.3)print *,cells(i)%Val%t_b(:),j
   !      cells(i)%Val%t_b(cells(i)%Val%counter)=j
   !      !print *, 'exist', cell_counter, lines(j)%c1
   !      exit
   !   end if
   !end do
   !if(temp_int==0)then
   !   boundary_counter=boundary_counter+1
   !   !lines(j)%Val%line_number=boundary_counter
   !   lines(j)%Val%c1=i
   !   lines(j)%Val%centx=temp_x
   !   lines(j)%Val%centr=temp_r
   !   lines(j)%Val%nv(1)= delta_r/temp1
   !   lines(j)%Val%nv(2)=-delta_x/temp1
   !   lines(j)%Val%n1=node_temp1
   !   lines(j)%Val%n2=node_temp2
   !end if

   !!print *,'1-2',node_temp1,node_temp2,temp_int3,temp_int4


   !!Middle point of Node2,3
   !delta_x=abs(nodes(node_temp2)%x-nodes(node_temp3)%x)
   !delta_r=abs(nodes(node_temp2)%r-nodes(node_temp3)%r)
   !temp1=sqrt(delta_x**2+delta_r**2)
   !temp_x=abs(nodes(node_temp2)%x+nodes(node_temp3)%x)*temp2
   !temp_r=abs(nodes(node_temp2)%r+nodes(node_temp3)%r)*temp2

   !!Ckeck that what elements face to this lines
   !do j = 1,boundary_counter
   !   temp_int3=lines(j)%Val%n1
   !   temp_int4=lines(j)%Val%n2
   !   temp_int=0
   !   !temp_int1=cells(i)%Val%counter
   !   if((node_temp2.eq.temp_int3.and.node_temp3.eq.temp_int4).or.(node_temp2.eq.temp_int4.and.node_temp3.eq.temp_int3)) then
   !      temp_int=1
   !      lines(j)%Val%c2=i
   !      cells(i)%Val%counter=cells(i)%Val%counter+1
   !      if(cells(i)%Val%counter.eq.3)print *,cells(i)%Val%t_b(:),j
   !      cells(i)%Val%t_b(cells(i)%Val%counter)=j
   !      !print *, 'exist', cell_counter, lines(j)%c1
   !      exit
   !   end if
   !end do
   !if(temp_int==0)then
   !   boundary_counter=boundary_counter+1
   !   !lines(j)%Val%line_number=boundary_counter
   !   lines(j)%Val%c1=i
   !   lines(j)%Val%centx=temp_x
   !   lines(j)%Val%centr=temp_r
   !   lines(j)%Val%nv(1)= delta_r/temp1
   !   lines(j)%Val%nv(2)=-delta_x/temp1
   !   lines(j)%Val%n1=node_temp2
   !   lines(j)%Val%n2=node_temp3
   !end if

   !!print *,'2-3',node_temp2,node_temp3,temp_int3,temp_int4

   !!Middle point of Node3,1
   !delta_x=abs(nodes(node_temp3)%x-nodes(node_temp1)%x)
   !delta_r=abs(nodes(node_temp3)%r-nodes(node_temp1)%r)
   !temp1=sqrt(delta_x**2+delta_r**2)
   !temp_x=abs(nodes(node_temp3)%x+nodes(node_temp1)%x)*temp2
   !temp_r=abs(nodes(node_temp3)%r+nodes(node_temp1)%r)*temp2

   !!Ckeck that what elements face to this lines
   !do j = 1,boundary_counter
   !   temp_int3=lines(j)%Val%n1
   !   temp_int4=lines(j)%Val%n2
   !   temp_int=0
   !   !temp_int1=cells(i)%Val%counter
   !   if((node_temp3.eq.temp_int3.and.node_temp1.eq.temp_int4).or.(node_temp3.eq.temp_int4.and.node_temp1.eq.temp_int3)) then
   !      temp_int=1
   !      lines(j)%Val%c2=i
   !      cells(i)%Val%counter=cells(i)%Val%counter+1
   !      if(cells(i)%Val%counter.eq.3)print *,cells(i)%Val%t_b(:),j
   !      cells(i)%Val%t_b(cells(i)%Val%counter)=j
   !      !print *, 'exist', cell_counter, lines(j)%c1
   !      exit
   !   end if
   !end do
   !if(temp_int==0)then
   !   boundary_counter=boundary_counter+1
   !   !lines(j)%Val%line_number=boundary_counter
   !   lines(j)%Val%c1=i
   !   lines(j)%Val%centx=temp_x
   !   lines(j)%Val%centr=temp_r
   !   lines(j)%Val%nv(1)= delta_r/temp1
   !   lines(j)%Val%nv(2)=-delta_x/temp1
   !   lines(j)%Val%n1=node_temp3
   !   lines(j)%Val%n2=node_temp1
   !end if

   !!print *,'3-1',node_temp3,node_temp1,temp_int3,temp_int4
   !}}}
   
   !Middle point of Node1,2
   delta_x=abs(nodes(node_temp1)%x-nodes(node_temp2)%x)
   delta_r=abs(nodes(node_temp1)%r-nodes(node_temp2)%r)
   temp1=sqrt(delta_x**2+delta_r**2)
   temp_x=abs(nodes(node_temp1)%x+nodes(node_temp2)%x)*temp2
   temp_r=abs(nodes(node_temp1)%r+nodes(node_temp2)%r)*temp2

   !Ckeck that what elements face to this lines
   do j = 1,boundary_counter
      temp_int3=lines(j)%Val%n1
      temp_int4=lines(j)%Val%n2
      temp_int=0
      !temp_int1=cells(i)%Val%counter
      if((node_temp1.eq.temp_int3.and.node_temp2.eq.temp_int4).or.(node_temp1.eq.temp_int4.and.node_temp2.eq.temp_int3)) then
         temp_int=1
         lines(j)%Val%c2=i
         !cells(i)%Val%counter=cells(i)%Val%counter+1
         !if(cells(i)%Val%counter.eq.3)print *,cells(i)%Val%t_b(:),j
         !cells(i)%Val%t_b(cells(i)%Val%counter)=j
         !print *, 'exist', cell_counter, lines(j)%c1
         exit
      end if
   end do
   if(temp_int==0)then
      boundary_counter=boundary_counter+1
      !lines(j)%Val%line_number=boundary_counter
      lines(j)%Val%c1=i
      lines(j)%Val%centx=temp_x
      lines(j)%Val%centr=temp_r
      lines(j)%Val%nv(1)= delta_r/temp1
      lines(j)%Val%nv(2)=-delta_x/temp1
      lines(j)%Val%n1=node_temp1
      lines(j)%Val%n2=node_temp2
   end if

   !print *,'1-2',node_temp1,node_temp2,temp_int3,temp_int4


   !Middle point of Node2,3
   delta_x=abs(nodes(node_temp2)%x-nodes(node_temp3)%x)
   delta_r=abs(nodes(node_temp2)%r-nodes(node_temp3)%r)
   temp1=sqrt(delta_x**2+delta_r**2)
   temp_x=abs(nodes(node_temp2)%x+nodes(node_temp3)%x)*temp2
   temp_r=abs(nodes(node_temp2)%r+nodes(node_temp3)%r)*temp2

   !Ckeck that what elements face to this lines
   do j = 1,boundary_counter
      temp_int3=lines(j)%Val%n1
      temp_int4=lines(j)%Val%n2
      temp_int=0
      !temp_int1=cells(i)%Val%counter
      if((node_temp2.eq.temp_int3.and.node_temp3.eq.temp_int4).or.(node_temp2.eq.temp_int4.and.node_temp3.eq.temp_int3)) then
         temp_int=1
         lines(j)%Val%c2=i
         !cells(i)%Val%counter=cells(i)%Val%counter+1
         !if(cells(i)%Val%counter.eq.3)print *,cells(i)%Val%t_b(:),j
         !cells(i)%Val%t_b(cells(i)%Val%counter)=j
         !print *, 'exist', cell_counter, lines(j)%c1
         exit
      end if
   end do
   if(temp_int==0)then
      boundary_counter=boundary_counter+1
      !lines(j)%Val%line_number=boundary_counter
      lines(j)%Val%c1=i
      lines(j)%Val%centx=temp_x
      lines(j)%Val%centr=temp_r
      lines(j)%Val%nv(1)= delta_r/temp1
      lines(j)%Val%nv(2)=-delta_x/temp1
      lines(j)%Val%n1=node_temp2
      lines(j)%Val%n2=node_temp3
   end if

   !print *,'2-3',node_temp2,node_temp3,temp_int3,temp_int4

   !Middle point of Node3,1
   delta_x=abs(nodes(node_temp3)%x-nodes(node_temp1)%x)
   delta_r=abs(nodes(node_temp3)%r-nodes(node_temp1)%r)
   temp1=sqrt(delta_x**2+delta_r**2)
   temp_x=abs(nodes(node_temp3)%x+nodes(node_temp1)%x)*temp2
   temp_r=abs(nodes(node_temp3)%r+nodes(node_temp1)%r)*temp2

   !Ckeck that what elements face to this lines
   do j = 1,boundary_counter
      temp_int3=lines(j)%Val%n1
      temp_int4=lines(j)%Val%n2
      temp_int=0
      !temp_int1=cells(i)%Val%counter
      if((node_temp3.eq.temp_int3.and.node_temp1.eq.temp_int4).or.(node_temp3.eq.temp_int4.and.node_temp1.eq.temp_int3)) then
         temp_int=1
         lines(j)%Val%c2=i
         !cells(i)%Val%counter=cells(i)%Val%counter+1
         !if(cells(i)%Val%counter.eq.3)print *,cells(i)%Val%t_b(:),j
         !cells(i)%Val%t_b(cells(i)%Val%counter)=j
         !print *, 'exist', cell_counter, lines(j)%c1
         exit
      end if
   end do
   if(temp_int==0)then
      boundary_counter=boundary_counter+1
      !lines(j)%Val%line_number=boundary_counter
      lines(j)%Val%c1=i
      lines(j)%Val%centx=temp_x
      lines(j)%Val%centr=temp_r
      lines(j)%Val%nv(1)= delta_r/temp1
      lines(j)%Val%nv(2)=-delta_x/temp1
      lines(j)%Val%n1=node_temp3
      lines(j)%Val%n2=node_temp1
   end if

   !print *,'3-1',node_temp3,node_temp1,temp_int3,temp_int4

end do

do i=1,cell_counter-1

   node_temp1=cells(i)%Val%t_l1
   node_temp2=cells(i)%Val%t_l2
   node_temp3=cells(i)%Val%t_l3

   !middle point of node1,2
   do j = 1,boundary_counter
      temp_int3=lines(j)%val%n1
      temp_int4=lines(j)%val%n2
      !temp_int1=cells(i)%val%counter
      if((node_temp1.eq.temp_int3.and.node_temp2.eq.temp_int4).or.(node_temp1.eq.temp_int4.and.node_temp2.eq.temp_int3)) then
         cells(i)%val%counter=cells(i)%val%counter+1
         cells(i)%val%t_b(cells(i)%val%counter)=j
         if(cells(i)%val%counter.eq.3)print *,cells(i)%val%t_b(:),j
         !print *,'cells',i,cells(i)%val%t_b
         exit
      end if
   end do
   !middle point of node2,3
   do j = 1,boundary_counter
      temp_int3=lines(j)%val%n1
      temp_int4=lines(j)%val%n2
      !temp_int1=cells(i)%val%counter
      if((node_temp2.eq.temp_int3.and.node_temp3.eq.temp_int4).or.(node_temp2.eq.temp_int4.and.node_temp3.eq.temp_int3)) then
         cells(i)%val%counter=cells(i)%val%counter+1
         cells(i)%val%t_b(cells(i)%val%counter)=j
         if(cells(i)%val%counter.eq.3)print *,cells(i)%val%t_b(:),j
         !print *,'cells',i,cells(i)%val%t_b
         exit
      end if
   end do
   !middle point of node3,1
   do j = 1,boundary_counter
      temp_int3=lines(j)%val%n1
      temp_int4=lines(j)%val%n2
      !temp_int1=cells(i)%val%counter
      if((node_temp3.eq.temp_int3.and.node_temp1.eq.temp_int4).or.(node_temp3.eq.temp_int4.and.node_temp1.eq.temp_int3)) then
         cells(i)%val%counter=cells(i)%val%counter+1
         cells(i)%val%t_b(cells(i)%val%counter)=j
         if(cells(i)%val%counter.eq.3)print *,cells(i)%val%t_b(:),j
         !print *,'cells',i,cells(i)%val%t_b
         exit
      end if
   end do

end do
do i=1,cell_counter-1
   temp_int1=cells(i)%Val%t_l1
   temp_int2=cells(i)%Val%t_l2
   temp_int3=cells(i)%Val%t_l3

   temp1=nodes(temp_int1)%x*(nodes(temp_int2)%r-nodes(temp_int3)%r)
   temp2=nodes(temp_int2)%x*(nodes(temp_int3)%r-nodes(temp_int1)%r)
   temp3=nodes(temp_int3)%x*(nodes(temp_int1)%r-nodes(temp_int2)%r)

   cells(i)%Val%area=0.5d0*abs(temp1+temp2+temp3)

end do


!do i=1,boundary_counter
!   print *,i,lines(i)%Val%c1,lines(i)%Val%c2
!end do
!do i=1,cell_counter-1
!   print *,i,cells(i)%Val%grax,cells(i)%Val%grar
!end do
!do i=1,cell_counter-1
!   print *,i,cells(i)%Val%t_b
!end do
!do i=1,boundary_counter
!   print *,i,lines(i)%Val%n1,lines(i)%Val%n2
!end do
!do i=1,node_num
!   print *,i,nodes(i)%x,nodes(i)%r
!end do

end subroutine grid
