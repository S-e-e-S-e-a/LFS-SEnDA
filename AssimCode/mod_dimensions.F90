module mod_dimensions
! LICOM model grid dimensions
#ifdef PACIFIC
   integer, parameter :: group=20
   integer, parameter :: nx=1440                ! i-dimension of model grid
   integer, parameter :: ny=720               ! j-dimension of model grid
   integer, parameter :: nz=55                 ! k-dimension of model grid
   integer, parameter :: startPx = 1
   integer, parameter :: startPy = 1
   integer, parameter :: subx1 = startPx, &
                         subx2 = startPx + nx - 1 , &
                         suby1 = startPy, &
                         suby2 = startPy + ny - 1
#endif
end module mod_dimensions
