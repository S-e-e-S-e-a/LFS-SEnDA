module m_get_mod_grid
contains 
subroutine get_mod_grid(modlon,modlat,depths,levels)
   use mod_dimensions
   implicit none
   real, dimension(nx,ny), intent(out) :: modlon,modlat,depths
   real*4, dimension(nz) :: levels
   character(len=8) tag8
   logical ex

!Add by fuww June 14, 2007
   integer, parameter :: licom_nx = 1440, licom_ny = 720
   real, dimension(licom_nx,licom_ny) :: licomlon, licomlat, licomdepths
   integer i,j,reclA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read position from model files
   inquire(file='newpos.uf',exist=ex)
   if (.not.ex) stop 'newpos.uf file does not exist'
   open(10,file='newpos.uf',form='formatted',status='old',action='read')
     read(10,*) licomlat, licomlon
   close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read depths from model files
   write (tag8,'(i4.4,a,i3.3)') licom_nx,'x', licom_ny
   inquire(file='depths'//tag8//'.uf',exist=ex)
   if (.not.ex) stop 'depths.uf file does not exist'
   open(10,file='depths'//tag8//'.uf',form='formatted',status='old',action='read')
      read(10,*)licomdepths,levels
   close(10)

   modlat(1:nx,1:ny) = licomlat(subx1:subx2,suby1:suby2)
   modlon(1:nx,1:ny) = licomlon(subx1:subx2,suby1:suby2)
!lzt
   depths(1:nx,1:ny) = licomdepths(subx1:subx2,suby1:suby2)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine  get_mod_grid
end module  m_get_mod_grid
