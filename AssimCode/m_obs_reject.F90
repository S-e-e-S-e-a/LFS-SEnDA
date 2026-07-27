module m_obs_reject
   implicit none

   real, allocatable :: D0(:)
   logical, allocatable :: lobs0(:) 
   logical, allocatable :: lpass0(:)   
   real :: threshold

contains

subroutine obs_reject1(myid,obs0,nrobs0,psi_local,mean_ssh,depths,levels,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2,nrobs_local,nrobs)
   use mod_dimensions
   use mod_states
   use mod_measurement
   use m_Generate_element_Sij
   implicit none

   integer,           intent(in) :: myid,nrobs0,prep_rx,ry1,ry2
   type(sub_states),  intent(in) :: psi_local(nx_end-nx_begin+2*prep_rx,ny_end-ny_begin+ry1+ry2) 
   type(measurement), intent(in) :: obs0(nrobs0)  
   real, intent(in):: mean_ssh(nx,ny), depths(nx,ny), levels(nz)
   integer,           intent(in) :: nx_begin,nx_end,ny_begin,ny_end
   integer,           intent(out) :: nrobs_local,nrobs     
   
   integer :: i, j, m, num, num1
   real    :: Spsi_val

   if (allocated(D0)) deallocate(D0)
   if (allocated(lobs0)) deallocate(lobs0)
   if (allocated(lpass0)) deallocate(lpass0)
   
   allocate(D0(nrobs0))
   allocate(lobs0(nrobs0))
   allocate(lpass0(nrobs0))

   num  = 0
   num1 = 0

   do m = 1, nrobs0
      i = obs0(m)%ipiv
      j = obs0(m)%jpiv
      
      if(j>ny_begin-ry1+2 .and. j<ny_end+ry2-1 .and. &
        ((i>nx_begin-prep_rx+2 .and. i<nx_end+prep_rx-1).or. &
        (i+nx>nx_begin-prep_rx+2 .and. i+nx<nx_end+prep_rx-1).or. &
        (i-nx>nx_begin-prep_rx+2 .and. i-nx<nx_end+prep_rx-1))) then
         
         lobs0(m) = .true.
         Spsi_val = Generate_element_Sij(obs0(m),psi_local,mean_ssh,depths,levels,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
         D0(m)    = obs0(m)%d - Spsi_val
      else
         lobs0(m) = .false.
         D0(m)    = 99999. 
      endif

      select case (obs0(m)%id)
      case ('TEM', 'SST')
         threshold = 6.0
      case ('SAL')
         threshold = 5.0
      case ('SLA', 'SSH')
         threshold = 2.0
      case default
         threshold = 0.0
      end select

      ! 3. 记录通过状态 [cite: 77]
      if (abs(D0(m)) < threshold) then
         lpass0(m) = .true.
         num1 = num1 + 1
         if (lobs0(m)) num = num + 1
      else
         lpass0(m) = .false.
      endif
   enddo

   nrobs_local = num
   nrobs = num1
   print*, myid, ' Optimized obs_reject1 end, local_j:', num
end subroutine obs_reject1

subroutine obs_reject2(myid,obs0,nrobs0,obs_local,nrobs_local,obs,nrobs)
   use mod_measurement
   implicit none

   integer,  intent(in) :: myid,nrobs_local,nrobs0,nrobs        
   type(measurement), intent(out) :: obs_local(nrobs_local)  
   type(measurement), intent(out) :: obs(nrobs)
   type(measurement), intent(in) :: obs0(nrobs0)  
   
   integer :: m, num, num1

   num  = 0
   num1 = 0
   
   do m = 1, nrobs0
      if (lpass0(m)) then
         num1 = num1 + 1
         obs(num1) = obs0(m)
         
         if (lobs0(m)) then
            num = num + 1
            obs_local(num) = obs0(m)
         endif
      endif
   enddo
end subroutine obs_reject2

end module m_obs_reject
