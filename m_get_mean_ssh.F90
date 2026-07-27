module m_get_mean_ssh
contains 
subroutine get_mean_ssh(mean_ssh)
   use mod_dimensions
   implicit none
   real, dimension(nx,ny), intent(out) :: mean_ssh
   logical ex

!Add by fuww June 14, 2007
!   integer, parameter :: licom_nx = 1440, licom_ny = 720
!   real, dimension(licom_nx,licom_ny) :: mean_ssh
   integer i,j,reclA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read position from model files
!   inquire(file='mean_ssh_obs.uf',exist=ex)
!   if (.not.ex) stop 'mean_ssh_obs.uf file does not exist'
!   open(10,file='mean_ssh_obs.uf',form='formatted',status='old',action='read')
   inquire(file='mdt_m0.uf',exist=ex) !zyz use mdt
   if (.not.ex) stop 'mdt.uf file does not exist' !yuzp20200105
   open(10,file='mdt_m0.uf',form='formatted',status='old',action='read') !yuzp20200105
     read(10,*) mean_ssh
   close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine  get_mean_ssh
end module  m_get_mean_ssh
