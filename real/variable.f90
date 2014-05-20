module variable
   use prmtr
   integer time
   !Fundamental_Variable_and_FLUX*********
   double precision, dimension(4,-1:ni  ,-1:nj  ) :: w,w_left,w_right,w_down,w_up !rho,u,v,P
   double precision, dimension(8,-1:ni  ,-1:nj  ) :: Nonw!rho,u,v,P,C
   double precision, dimension(4, 1:ni  , 1:nj  ) :: wp,q,qq
   double precision, dimension(   1:ni  , 1:nj  ) :: Mach_number
   double precision, dimension(4, 1:ni  , 1:nj  ) :: X_Flux, X_Numerical,vis_i
   double precision, dimension(4, 1:ni  , 1:nj  ) :: Y_Flux, Y_Numerical,vis_j

   !Geometry_Variable*******************
   double precision, dimension(   0:ni  , 0:nj  ) :: Vol, Area
   double precision, dimension(   0:ni+1, 0:nj+1) :: x, r
   double precision, dimension(   0:ni  , 0:nj  ) :: grid_cen, dx, dr, dt
   double precision, dimension(   1:ni  , 1:nj  ) :: xgrid, rgrid
   double precision, dimension(2, 0:ni+1, 0:nj+1) :: nvi
   double precision, dimension(2, 0:ni+1, 0:nj+1) :: nvj
   double precision vnl, vnr, vtl, vtr
   double precision, dimension(4)                 :: delta_x, delta_r
   double precision, dimension(   0:ni+1, 0:nj+1) :: dsi
   double precision, dimension(   0:ni+1, 0:nj+1) :: dsj
   double precision, dimension(2,2, 0:ni  , 0:nj  ) :: geojaci
   double precision, dimension(2,2, 0:ni  , 0:nj  ) :: geojacj

   !OTHER************************
   double precision, parameter :: gamma_bar=gamma-1.d0
   !double precision rhol, rhor, ul, ur, vl, vr, pl, pr, El, Er, Hl, Hr, t
   double precision AB, BC, CD, DA
   double precision residual,t
   integer          temp_int
   
   !LUSGS*************************
   !double precision, dimension(  4, 0:ni  , 0:nj  ) :: RHS, q_imp
   !double precision, dimension(     0:ni  , 0:nj  ) :: alpha
   !double precision, dimension(4,4, 0:ni  , 0:nj  ) :: Ap, Am, Bp, Bm

   character*20 tmpstring
end module variable
