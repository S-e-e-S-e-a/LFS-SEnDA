module m_prep_4_EnOI

! This subroutine uses the observation and ensembles from the model
! to prepare the input to the EnOI analysis scheme.
! The output from this routine is used directly in the global analysis
! while the output has to be run through a "filter" to be used in the
! local analysis scheme.

contains
subroutine prep_4_EnOI(myid,obs_local,nrobs_local,A_local,psi_local,mean_ssh,nrens,depths,levels,D,S_local,alpha,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
   use mod_dimensions
   use mod_states
   use mod_measurement
! Functions and subroutines
   use m_Generate_element_Sij

   implicit none

! Input variables
   integer,           intent(in) :: myid,nrobs_local,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2         ! Number of measurements
   integer,           intent(in) :: nrens         ! Size of ensemble
   type(sub_states),  intent(in) :: A_local(nx_end-nx_begin+2*prep_rx,ny_end-ny_begin+ry1+ry2,nrens)      ! Ensemble of model states 
!wwq
   type(sub_states),  intent(in) :: psi_local(nx_end-nx_begin+2*prep_rx,ny_end-ny_begin+ry1+ry2)           ! Model state 
   type(measurement), intent(in) :: obs_local(nrobs_local)    ! measurements
   real*4, intent(in):: mean_ssh(nx,ny)
   real*4, intent(in):: alpha
! Add by fuww
   real*4, intent(in):: depths(nx,ny), levels(nz)
! Local variables 
   real*4, intent(inout) :: D(nrobs_local)
   real*4, intent(inout) :: S_local(nrens,nrobs_local)
   real*4  Spsi(nrobs_local)
   real*4  meanS(nrobs_local)    ! Automatic array
   integer i,j,m,iens
!lzt

   meanS=0.0
   do m =1, nrobs_local  ! LB loops inverted
        do iens =1, nrens
            S_local(iens,m)  =  Generate_element_Sij(obs_local(m),A_local(:,:,iens),mean_ssh,depths,levels,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
            meanS(m)=meanS(m)+S_local(iens,m)
        enddo
        Spsi(m) = Generate_element_Sij(obs_local(m),psi_local, mean_ssh,depths,levels,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
   enddo
   meanS=(1.0/float(nrens))*meanS
   do j=1,nrens
      S_local(j,:)=S_local(j,:)-meanS(:)
   enddo 
   do m=1,nrobs_local
      D(m)=obs_local(m)%d 
   enddo

! Compute innovation D'=D-Hpsi
   D=D-Spsi
   print*,myid,'m_prep_4_EnOI:End Calculate innovations D=D-s'

end subroutine prep_4_EnOI
end module m_prep_4_EnOI
