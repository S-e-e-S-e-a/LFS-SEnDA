module m_Generate_element_Sij 
contains
! DESCRIPTION
! Function which generates the S matrix element by element. The function operates 
! for different observation types: temp, sal, u, v and SSH/SLA data. 
! The four first data types are found as progonostic parameters in the model
! While the SSH/SLA data has to be generated from the prognositc data. This
! Is done by a separte function in the case of this type of measurments.
! For the given location of the observation we combine the four surrounding 
! vertical model data using the biilinear coeffisients given for each obs(j).
! This is done for all types of observations. 
! In the  case of horizontal data we do not need to do more with the data
! While in case of vertical measurment profiles (data for depth .ne. 0) 
! the data has to be interporlated to a denser gridd in the vertical and for
! This we use a sline interpolation.

real function Generate_element_Sij(obs,mem,mean_ssh,depths,levels,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
   use mod_states
   use mod_measurement
   use mod_dimensions
!functions
   use m_modsubstate_point

   implicit none
   integer, intent(in) :: nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2
   type(measurement),    intent(in) :: obs
   type(sub_states),     intent(in) :: mem(nx_end-nx_begin+2*prep_rx,ny_end-ny_begin+ry1+ry2)
   real,                 intent(in) :: mean_ssh(nx,ny)
! Add by fuww
   real, intent(in):: depths(nx,ny), levels(nz)
   real deep(1)

   integer                          :: i,j,ip1,jp1
   integer                          :: k, klow, kup  !Add by fuww
   integer                          :: ii,jj, imin,imax,jmin,jmax, cnt
   real                             :: x0,y0,p1,p2,p3,p4

   real model_profile(2),  mod_value
   real obs_depths
   real, parameter :: undef=1.0E+35     ! land points have value huge()
   i   = obs%ipiv
   j   = obs%jpiv
      obs_depths= abs(obs%depths)  ! caveat by fuww depths must be positive  
      deep(1)=obs_depths
   if (deep(1) .lt. 5.0) then !zyz
     if(obs%ns .lt. 2) then ! point data : zero support
         ip1 = min(i+1,nx)
         jp1 = min(j+1,ny)
         p1=modsubstate_point(obs%id,mem,mean_ssh,i,j,1,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
         p2=modsubstate_point(obs%id,mem,mean_ssh,ip1,j,1,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
         p3=modsubstate_point(obs%id,mem,mean_ssh,ip1,jp1,1,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
         p4=modsubstate_point(obs%id,mem,mean_ssh,i,jp1,1,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
       
         Generate_element_Sij = p1*obs%a1 + p2*obs%a2 + p3*obs%a3 + p4*obs%a4
     else     ! data support assumed a square of 2ns*2ns grid cells
         imin = max(1,i-obs%ns)
         imax = min(nx,i+obs%ns)
         jmin = max(1,j-obs%ns)
         jmax = min(ny,j+obs%ns)
         cnt = 0                          ! counter
         Generate_element_Sij = 0.
         do jj= jmin, jmax
         do ii= imin, imax
            if (depths(ii,jj) > 20.0 .and.abs((depths(ii,jj)-undef)/undef)>0.01) then    ! Caveat by fuww
               mod_value=modsubstate_point(obs%id,mem,mean_ssh,ii,jj,1,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
               cnt=cnt+1
               Generate_element_Sij = Generate_element_Sij + mod_value
            endif
         enddo
         enddo
         if (cnt == 0) then
            print*, ' observation on land ', i,j, obs%d
            stop 'm_Generate_element_Sij: report bug to LB (laurentb@nersc.no)'
         endif
         Generate_element_Sij = Generate_element_Sij / float(cnt)
     endif
 else       ! in-situ data (in depth)
      model_profile(:)=0.0
      obs_depths= abs(obs%depths)  ! caveat by fuww depths must be positive  
      do k=1,nz  
      if(obs_depths>levels(k)-1.0.and.obs_depths<levels(k)+1.0)then
      klow = k
      exit
      endif      
      enddo
      ip1 = min(i+1,nx)
      jp1 = min(j+1,ny)
     model_profile(1) = modsubstate_point(obs%id,mem,mean_ssh,i,j,klow,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)*obs%a1 &
                       & +modsubstate_point(obs%id,mem,mean_ssh,ip1,j,klow,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)*obs%a2 &
                       & +modsubstate_point(obs%id,mem,mean_ssh,ip1,jp1,klow,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)*obs%a3&
                       & +modsubstate_point(obs%id,mem,mean_ssh,i,jp1,klow,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)*obs%a4
    Generate_element_Sij = model_profile(1)
endif

end function Generate_element_Sij
end module m_Generate_element_Sij
