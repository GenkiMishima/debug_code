module variable
   use prmtr
   integer time
   double precision, dimension(5,-1:ni  ,-1:nj  ) :: w, w_left, w_right, w_down, w_up
   double precision, dimension(5, 1:ni  , 1:nj  ) :: wp,q,qq
   double precision, dimension(5, 1:ni  , 1:nj  ) :: X_Flux, X_Numerical,vis_i
   double precision, dimension(5, 1:ni  , 1:nj  ) :: Y_Flux, Y_Numerical,vis_j
   double precision, dimension(5, 0:ni  , 0:nj  ) :: RHS, q_imp
   double precision, dimension(5, 1:ni  , 1:nj  ) :: X_inst
   double precision, dimension(5, 1:ni  , 1:nj  ) :: Y_inst
   double precision, dimension(   0:ni  , 0:nj  ) :: Vol, Area
   double precision, dimension(   0:ni+1, 0:nj+1) :: x, r
   double precision, dimension(   0:ni  , 0:nj  ) :: grid_cen, dx, dr, dt
   double precision vnl, vnr, vtl, vtr
   double precision, dimension(   1:ni  , 1:nj  ) :: Mach_number
   double precision, dimension(   1:ni  , 1:nj  ) :: temp_matrix5,temp_matrix
   double precision, dimension(   1:ni  , 1:nj  ) :: xgrid, rgrid
   double precision, dimension(5)                 :: flux
   !double precision, dimension(4)                 :: corner_for, corner_back
   double precision, dimension(5)                 :: delta_x, delta_r
   double precision, dimension(   0:ni+1, 0:nj+1) :: dsi
   double precision, dimension(   0:ni+1, 0:nj+1) :: dsj
   double precision, parameter :: gamma_bar=gamma-1.d0
   double precision temp0, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
   double precision temp_residual1,temp_residual2
   double precision rhol, rhor, ul, ur, vl, vr, pl, pr, El, Er, Hl, Hr, t
   double precision uaver, vaver, Haver, aaver
   double precision AB, BC, CD, DA
   double precision residual
   integer          temp_int
   !double precision, dimension(     0:ni  , 0:nj  ) :: alpha
   !double precision, dimension(4,4, 0:ni  , 0:nj  ) :: Ap, Am, Bp, Bm

   double precision, dimension(2,   0:ni+1, 0:nj+1) :: nvi
   double precision, dimension(2,   0:ni+1, 0:nj+1) :: nvj
   double precision, dimension(2,2, 0:ni  , 0:nj  ) :: geojaci
   double precision, dimension(2,2, 0:ni  , 0:nj  ) :: geojacj

   type Chem_prmtr
      integer moler_weight
      double precision rho,mass_rate,Y!mass_fraction
   end type Chem_prmtr
   type (Chem_prmtr) Chem_num(nY,-1:ni,-1:nj)
   double precision temp_chem_mat (-1:ni,-1:nj)

   character*20 tmpstring
end module variable
