module prmtr
   implicit none
   double precision, parameter :: pi=2.d0*acos(0.d0)
   double precision, parameter :: mu=1.8d-5 !viscous term
   double precision, parameter :: gamma=1.4d0
   double precision, parameter :: gas_constant=8.3144621d0
   double precision, parameter :: molar_mass=28.96d0
   double precision, parameter :: gas_specific=gas_constant/molar_mass*1.d3
   double precision, parameter::length=1d-6
   !Grid number
   integer         , parameter :: ni=224
   integer         , parameter :: nj=100
   !MUSCL INCREMENT
   double precision, parameter :: phi=1.d0
   double precision, parameter :: kappa=1.d0/3.d0
   !Division
   double precision, parameter :: CFL= 1d-1
   !STEP PARAMETER
   !integer, parameter :: step_forward =101
   integer, parameter :: step_backward=101
   integer, parameter :: step_height  =101
   !OUTPOT TIME
   integer, parameter :: out_time = 500
   integer, parameter :: time_max = 5000000
   !RESIDUAL
   double precision, parameter :: epsilon= 1d-2
end module prmtr
