module glbl_prmtr
implicit none
integer, parameter ::node_num = 18
integer, parameter ::line_num = 10
integer, parameter ::triangle_num = 22
integer, parameter ::element_num = 34
double precision, parameter ::gamma=1.4d0
double precision, parameter ::CFL=0.01d0
integer line_counter,cell_counter,boundary_counter
double precision t,dt
character moji*30

end module glbl_prmtr
module form_format
implicit none
type cell_comp
   integer counter
   integer t_l1,t_l2,t_l3
   integer t_b(2)
   double precision grax,grar
   double precision w(4), q(4)
end type cell_comp
type Cell
   type (cell_comp) Val
end type Cell
type line_comp
   !integer line_number
   integer bc
   integer(4) c1,c2
   double precision centx,centr
   double precision ds
   double precision nv(2)
end type line_comp
type Line
   type (line_comp) Val
end type Line
type Node
   integer node_num_count
   integer node_bound_count
   integer,dimension(0:9)::cell_num
   double precision x,r
   double precision value1
end type Node
!type boundary_condition
!   double precision w(4)
!end type boundary_condition
end module form_format
module grid_prmtr
use glbl_prmtr
use form_format
implicit none
type (Node) nodes (1:node_num   )
type (Line) lines (1:1d6+50)
type (Cell) cells (1:1d6+50)
type(cell_comp) temp_cell1,temp_cell2
type(line_comp) temp_line

end module grid_prmtr
