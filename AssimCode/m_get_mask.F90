module m_get_mask
contains 
subroutine get_mask(assim_mask)
   use mod_dimensions
   implicit none
   real, dimension(nx,ny), intent(out) :: assim_mask
   logical ex
   integer i,j,reclA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   inquire(file='mask.uf',exist=ex) 
   if (.not.ex) stop 'mask.uf file does not exist' !yuzp20200105
   open(10,file='mask.uf',form='formatted',status='old',action='read') !yuzp20200105
     read(10,*) assim_mask
   close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine  get_mask
end module  m_get_mask
