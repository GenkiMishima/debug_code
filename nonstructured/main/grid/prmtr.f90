module form_format
implicit none
type cell_comp
   !integer cell_number
   integer t_l1,t_l2,t_l3
   double precision grax,grar
   double precision value1
end type cell_comp
type Cell
   type (cell_comp) Val
end type Cell
type line_comp
   !integer line_number
   integer bc
   integer c1,c2
   double precision centx,centr
   double precision nv(2)
end type line_comp
type Line
   type (line_comp) Val
end type Line
type Node
   integer node_num_count
   integer,dimension(10)::cell_num
   double precision x,r
   double precision value1
end type Node
end module form_format
module prmtr
use form_format
implicit none
double precision, parameter ::node_num = 908
double precision, parameter ::line_num = 10
double precision, parameter ::triangle_num = 22
double precision, parameter ::element_num = 1814
integer line_counter,cell_counter,boundary_counter
!integer temp_int, temp_int1,temp_int2,temp_int3,temp_int4, bc
!double precision temp, temp_x, temp_r,temp1, temp2,temp3
!double precision node_temp1,node_temp2,node_temp3
!integer, dimension(element_num) ::ele_type
character moji*30
!type (triangles_list  )   t_lists( element_num+20)
!type (boundarys_list  )   b_lists( element_num+20)
!type (boundarys_format) boundarys( element_num+20)
!type (triangles_format) triangles(triangle_num+20)
type (Node) nodes (1:node_num   )
type (Line) lines (1:element_num+1000000)
type (Cell) cells (1:element_num)
type(cell_comp) temp_cell
type(line_comp) temp_line
end module prmtr
