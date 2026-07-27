module m_get_obs
contains
integer function get_nrobs_d(filename,ex)
! Reads the observations to be used for assimilation from the file
! observation.uf. Each element is of type measurement and will be 
! stored in the vector d defined below.
   use mod_measurement
   implicit none

   integer reclO             ! Record length for an observation type
   type(measurement) Obs     ! measurements

   integer j
   logical ex
   character(len=256) :: filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find nb of  measurements stored
   inquire(file=trim(filename),exist=ex)
   if (.not.ex) stop 'File "observations.uf" does not exist'
   inquire(iolength=reclO)Obs
   open(10,file=trim(filename),form='unformatted',access='direct',recl=reclO)
   do j=1,4000000
      read(10,rec=j,err=200)Obs
   enddo
   200 get_nrobs_d=j-1
end function get_nrobs_d

subroutine get_obs_d(obs,nrobs,filename)
   use mod_measurement
   implicit none

   integer,           intent(in)  :: nrobs           ! Number of measurements
   type(measurement), intent(out) :: obs(nrobs)      ! measurements

   integer reclO,j
   character(len=256) :: filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read measurements and store in obs
   inquire(iolength=reclO)obs(1)
   open(10,file=trim(filename),form='unformatted',access='direct',recl=reclO)
   do j=1,nrobs
      read(10,rec=j)obs(j)
   enddo
   close(10)

end subroutine get_obs_d
end module m_get_obs
