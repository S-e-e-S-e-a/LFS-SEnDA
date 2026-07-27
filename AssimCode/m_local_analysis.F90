module m_local_analysis
contains 
subroutine local_analysis(myid,nrens,obs,nrobs,obs_local,nrobs_local,A_local,psi_local,radius_sla,radius_sal,radius_tem,radius_sst,modlon,modlat,depths,D,S_local,ndim,alpha,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2,assim_mask,anal_method,anal_solver,levels,vert_local)
   use mod_dimensions
   use mod_states
   use mod_measurement
   use m_getD ! add getDvec getpsi and m_getD4S special getD for S vector
   use m_active_obs
   use m_spherdist

   implicit none
   
   integer,                      intent(in)    :: myid,nrens,nrobs,nrobs_local,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2  ! Size of ensemble
   type(sub_states),             intent(in)    :: A_local(nx_end-nx_begin+2*prep_rx,ny_end-ny_begin+ry1+ry2,nrens)        ! static Ensemble matrix
   type(sub_states),             intent(inout) :: psi_local(nx_end-nx_begin+2*prep_rx,ny_end-ny_begin+ry1+ry2)             ! state vector
   type(measurement),            intent(in)    :: obs(nrobs),obs_local(nrobs_local)      ! measurements
   real*4, dimension(nx,ny),       intent(in)    :: modlon,modlat,depths,assim_mask
   real*4, dimension(nrobs_local),       intent(in)    :: D
   real*4, dimension(nrens,nrobs_local), intent(in)    :: S_local
   real*4,                         intent(in)    ::   radius_sla,radius_sal,radius_tem,radius_sst
   integer, intent(in) :: ndim            ! Number of elements per sample in analysis call
   logical, intent(in) :: vert_local
   real*4 :: alpha   
   type(sub_states)       subA(nrens)           ! local+static Ensemble matrix
   type(sub_states)       subpsi                ! local state vector
   logical                lobs(nrobs_local)           ! Active measurements in local analysis
   real*4                 wlobs(nrobs_local)           ! weight for active measurements in local analysis
   real*4                 levels(nz)
   real*4  R_trace, Lh
   integer nobs                                 ! Number of measurements in local analysis
   integer m,i,j,p,ii,jj,k,m1,m2,mm,f    !lzt 20120313
   real, parameter :: undef=1.0E+35     ! land points have value huge()

!   logical test
!lzt
   real meanD
!lzt 20120313
   integer, allocatable, dimension(:)  :: dep_obs_ind
   character(len=3), allocatable, dimension(:)  :: type_obs
   real*4, allocatable, dimension(:)   :: subD
   real*4, allocatable, dimension(:,:) :: subS,subE,R !geir changement
   character(len=3) :: anal_solver
   character(len=4) :: anal_method
!$OMP PARALLEL DO PRIVATE(j,i,m,mm,subA,lobs,wlobs,nobs,subD,subE,subS,R,type_obs,dep_obs_ind,subpsi,subA) SHARED(ndim,nrens,nrobs,D,E,S) SCHEDULE(STATIC,1)
   do j=ny_begin+1,ny_end
     do i=nx_begin+1,nx_end
         if (depths(i,j) > 0.0.and.abs((depths(i,j)-undef)/undef)>0.01) then 
            subpsi=psi_local(i-nx_begin+prep_rx,j-ny_begin+ry1)      ! The state at grid point i,j
            do m=1,nrens
              subA(m)=A_local(i-nx_begin+prep_rx,j-ny_begin+ry1,m)    ! The ensemble at grid point i,j
            enddo
            call active_obs(i,j,obs_local,nrobs_local,lobs,wlobs,nobs,radius_sla,radius_sal,radius_tem,radius_sst,modlon,modlat)
            if (nobs > 2) then
                 print*,'myid,i,j,nobs',myid,i,j,nobs
               allocate(subD(nobs))
               allocate(subS(nobs,nrens))
               call getDvecnew(D,subD,nrobs_local,lobs,wlobs,nobs) ! the innovations to use 
               call getD(S_local,subS,nrobs_local,nrens,lobs,nobs) ! the HA' to use
               allocate (R(nobs,nobs)) !geir changement
               R=0.
               allocate (dep_obs_ind(nobs))
               allocate (type_obs(nobs))
               !allocate (lon_obs(nobs))
               !allocate (lat_obs(nobs))
               !allocate (nz_obs(nobs))
               !lon_obs=0.
               !lat_obs=0.
               !nz_obs=0.
               dep_obs_ind=1
               type_obs=''
               mm = 0
               do m=1, nrobs_local
                 if (lobs(m)) then
                  mm = mm + 1
                  R(mm,mm) = obs_local(m)%var  !/ max(wlobs(m), 1e-5) ! R loc
                  type_obs(mm) = obs_local(m)%id
               !   lon_obs(mm) = obs_local(m)%lon
               !   lat_obs(mm) = obs_local(m)%lat
                  if (vert_local) then
                   dep_obs_ind(mm)=minloc( abs(levels - obs_local(m)%depths), dim=1)
                  endif
                 endif
               enddo
              if ( anal_method == 'EnOI' ) then
                 call analysis2(subA,subpsi, subD, R, subS, ndim, nz, nrens,nobs, alpha, assim_mask(i,j), modlon(i,j), modlat(i,j), anal_solver, vert_local, dep_obs_ind, type_obs)
              else if ( anal_method == 'ETKF' ) then
                 call analysis2_letkf(subA,subpsi, subD, R, subS, ndim, nz, nrens,nobs, alpha, assim_mask(i,j), modlon(i,j), modlat(i,j), anal_solver, vert_local, dep_obs_ind, type_obs)
              else
                 print *, '================================================'
                 print *, ' ERROR: anal_method "', trim(anal_method), '" is unavailable!'
                 print *, ' Please check your namelist or input settings.'
                 print *, '================================================'
              end if
              psi_local(i-nx_begin+prep_rx,j-ny_begin+ry1)=subpsi   ! The analysis for grid point i,j
              deallocate(subD, subS)
              deallocate(dep_obs_ind,type_obs)
              !deallocate(lon_obs, lat_obs, nz_obs)
              deallocate(R) !geir changement 
            endif              !nobs>2
         endif
     enddo
   enddo
   print *, ' Leaving m_local_analysis'

end subroutine local_analysis
end module m_local_analysis
