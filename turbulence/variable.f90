module variable
use prmtr
integer time
double precision, dimension(6,-1:ni  ,-1:nj  ) :: w    !rho,u,v,P,Tenergy,dissipation factor
double precision, dimension(6, 1:ni  , 1:nj  ) :: wp,q,qq, Sterm
double precision, dimension(6, 1:ni  , 1:nj  ) :: X_Flux, X_Numerical,vis_i
double precision, dimension(6, 1:ni  , 1:nj  ) :: Y_Flux, Y_Numerical,vis_j
double precision, dimension(   0:ni  , 0:nj  ) :: Vol, Area
double precision, dimension(   0:ni  , 0:nj  ) :: grid_cen, dx, dr, dt
double precision, dimension(   0:ni+1, 0:nj+1) :: x, r
double precision vnl, vnr, vtl, vtr
double precision residual,t
double precision, dimension(   1:ni  , 1:nj  ) :: temp_matrix5,temp_matrix
double precision, dimension(4)                 :: corner_for, corner_back
double precision, dimension(   0:ni+1, 0:nj+1) :: dsi
double precision, dimension(   0:ni+1, 0:nj+1) :: dsj
double precision, parameter :: gamma_bar=gamma-1.d0
double precision temp_residual1,temp_residual2
double precision AB, BC, CD, DA
integer          temp_int
double precision, dimension(     0:ni  , 0:nj  ) :: alpha
double precision, dimension(4,4, 0:ni  , 0:nj  ) :: Ap, Am, Bp, Bm
double precision, dimension(2,   0:ni+1, 0:nj+1) :: nvi
double precision, dimension(2,   0:ni+1, 0:nj+1) :: nvj
double precision, dimension(2,2, 0:ni  , 0:nj  ) :: geojaci
double precision, dimension(2,2, 0:ni  , 0:nj  ) :: geojacj

character*20 tmpstring
end module variable
