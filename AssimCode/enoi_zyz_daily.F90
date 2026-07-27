Module EnOI
  use mod_states
  use mod_measurement
  implicit none
  type(sub_states), allocatable    :: A_local(:,:,:)      ! static Ensemble matrix
  type(sub_states), allocatable    :: psi_local(:,:)

  integer nrens                                 ! Size of static ensemble
  real, dimension(nx,ny) :: modlon,modlat     ! is used in active_obs,to calculate the spherdist
  real, allocatable :: mean_ssh(:,:), assim_mask(:,:)
  real depths(nx,ny), levels(nz)
  integer :: nx_begin,nx_end,ny_begin,ny_end,nx_len,ny_len
  integer :: prep_rx,ry1,ry2

  type(measurement), allocatable :: obs0(:),obs(:),obs_local(:)    ! measurements
  type(measurement), allocatable :: obs01(:),obs02(:),obs03(:),obs04(:)    ! measurements
  integer nrobs0,nrobs,nrobs_local,nrobs01,nrobs02,nrobs03,nrobs04    ! Number of measurements

  real, allocatable, dimension(:,:) :: S_local   !S_plus,S
!  integer S_k
!  integer, allocatable, dimension(:) :: S_ks
  real, allocatable, dimension(:) :: D,D2
  real alpha
  real radius_sst,radius_tem,radius_sla,radius_sal
  character(len=9) rident
  character(len=5) rident_clim

  !Logical variables
  logical ex
  logical :: l_global_analysis=.false.
  logical :: l_local_analysis=.false.
  logical :: vert_local=.false.

  integer :: year, month, day
  character(len=3) :: obstype, anal_solver
  character(len=4) :: anal_method, assim_type
  character(len=256) :: filename, filename_clim, outfile
  character(len=256) :: filenames(4)
  
  !wwq
contains
subroutine read_A_psi(myid,numprocs,mympi_group_size,color,key,SplitWorld,ierr)
  use mod_dimensions
  use mod_states
!  use mpi
  implicit none
  include "/opt/hpc/software/mpi/hpcx/v2.7.4/intel-2017.5.239/include/mpif.h" !yuzp20211212

  integer, intent(in) :: myid,numprocs,mympi_group_size,color,key,SplitWorld,ierr
  integer type_substates_MPI,blocklens_global(1),offsets_global(1),oldtypes_global(1)
  type(sub_states4),allocatable :: substates4(:),substates4_clim(:)
  type(sub_states), allocatable :: psi_scatter(:)
  type(sub_states) :: substates
  integer reclA
  integer prep_ry,group_size
  integer yyear,mmonth,dday,i,j,jj,y1,y2
  integer, parameter :: NDperY = 15
  integer, parameter :: NYEAR = 15
  integer, dimension(NDperY) :: spm, spd
  logical :: special_days(12, 31)
  integer :: i_yr, i_day, i_ens
  type(sub_states), allocatable :: daily_mean(:)
  integer:: m,n,h,nn,sub_i
  integer*8 pos_begin
  integer time1,time2,timeA,timeB,time_sum,time_read1,time_read2,time_read
  integer, parameter :: maxday(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
    blocklens_global(1)=local_ndim
    offsets_global(1)=0
    oldtypes_global(1)=MPI_REAL
    call MPI_TYPE_STRUCT(1,blocklens_global,offsets_global,oldtypes_global,type_substates_MPI,ierr)
    call MPI_TYPE_COMMIT(type_substates_MPI,ierr)

  !calculate A_local_count
    !nrens=120
    nrens=NDperY*NYEAR
    group_size=ny/group            !linrp 20200213

    prep_ry=15                 !linrp revised 20200217
    prep_rx=prep_ry*2

    ny_begin=color*group_size
    ny_end=(color+1)*group_size

    nx_begin=key*(nx/mympi_group_size)
    nx_end=(key+1)*(nx/mympi_group_size)

    if(ny_begin-prep_ry<0) then    
      ry1=ny_begin
    else 
      ry1=prep_ry
    endif
    if(ny_end+prep_ry>ny) then
      ry2=ny-ny_end
    else 
      ry2=prep_ry
    endif
  
    nx_len=nx_end-nx_begin+2*prep_rx
    ny_len=ny_end-ny_begin+ry1+ry2
   
    open(103,file='date_type.dat',form='formatted',status='old',action='read')
      read(103,*) year, month, day
    close(103)

    open(104, file='sample_date.dat', form='formatted', status='old', action='read')
    special_days = .false.
    do i = 1, NDperY
        read(104, *) spm(i), spd(i)
        special_days(spm(i), spd(i)) = .true.
    end do
    close(104)


    !!read and sendrecv static_ensemble
    allocate(A_local(nx_len,ny_len,nrens))
    allocate(psi_local(nx_len,ny_len))
    write(filename,'(a)') '../input/static_ensemble_daily_new/static_ensemble_25km'
!    write(filename_clim,'(a)') '../input/static_ensemble_daily_clim/static_ensemble_25km_clim'
    call system_clock(time1)
    j=0
    if(key==0) then
      allocate(substates4(ny_len*nx))
      allocate(substates4_clim(ny_len*nx))
      allocate(psi_scatter((nx/mympi_group_size+2*prep_rx)*(ny_end-ny_begin+ry1+ry2)*mympi_group_size))
    endif
    time_sum=0
    time_read=0 

    do yyear=2007,2016
    do mmonth=1,12
    do dday=1,1
      j=j+1
      call system_clock(timeA)
      if(key==0) then
        call system_clock(time_read1)
        print*,'myid mmonth',myid,mmonth,'dday',dday
        write (rident,'(a,I4.4,I2.2,I2.2)') '_',yyear,mmonth,dday
        open(10,file=trim(filename)//rident//'_sub.uf',status='old',form='unformatted',access='stream',action='read')
        pos_begin=ny_begin-ry1
        read(10,pos=10+(pos_begin*nx*local_ndim*4)) substates4
        print*,'read success'
        call system_clock(time_read2)
        time_read=time_read+time_read2-time_read1
        do m=0,mympi_group_size-1
        do h=0,ny_end-ny_begin+ry1+ry2-1
        do n=m*(nx/mympi_group_size)-prep_rx+1,(m+1)*(nx/mympi_group_size)+prep_rx
          if(n<=0) then
            nn=n+nx
          else if(n>nx) then
            nn=n-nx
          else
            nn=n
          endif
          call substates4to8(substates,substates4(h*nx+nn))
          psi_scatter(m*(ny_end-ny_begin+ry1+ry2)*(nx/mympi_group_size+2*prep_rx)+h*(nx/mympi_group_size+2*prep_rx)+n-m*(nx/mympi_group_size)+prep_rx)=substates
        enddo
        enddo
        enddo
        close(10)
      print*,'read end'
      endif
      call system_clock(timeB)
      time_sum=time_sum+timeB-timeA
      call MPI_Barrier(SplitWorld,ierr)
      call MPI_Scatter(psi_scatter,(ny_end-ny_begin+ry1+ry2)*(nx/mympi_group_size+2*prep_rx),type_substates_MPI,A_local(1,1,j),(ny_end-ny_begin+ry1+ry2)*(nx/mympi_group_size+2*prep_rx),type_substates_MPI,0,SplitWorld,ierr)
    enddo
    enddo
    enddo
    
!    do yyear=1999,1999+NYEAR-1
!    do mmonth=1,12
!    do dday=1,31
!    if (special_days(mmonth, dday)) then
!      j=j+1
!      call system_clock(timeA)
!      if(key==0) then
!        call system_clock(time_read1)
!        print*,'myid mmonth',myid,mmonth,'dday',dday
!        write (rident,'(a,I4.4,I2.2,I2.2)') '_',yyear,mmonth,dday
!        write (rident_clim,'(a,I2.2,I2.2)') '_',mmonth,dday
!        open(10,file=trim(filename)//rident//'_sub.uf',status='old',form='unformatted',access='stream',action='read')
!!        open(11,file=trim(filename_clim)//rident_clim//'_sub.uf',status='old',form='unformatted',access='stream',action='read')
!        pos_begin=ny_begin-ry1
!        read(10,pos=10+(pos_begin*nx*local_ndim*4)) substates4
!!        read(11,pos=10+(pos_begin*nx*local_ndim*4)) substates4_clim
!!        do i=1, ny_len*nx
!!        substates4_clim(i) = nan_states(substates4_clim(i),1e30)
!!        substates4(i)=subtract_states(substates4(i), substates4_clim(i))
!!        enddo
!        print*,'read success'
!        call system_clock(time_read2)
!        time_read=time_read+time_read2-time_read1
!        do m=0,mympi_group_size-1
!        do h=0,ny_end-ny_begin+ry1+ry2-1
!        do n=m*(nx/mympi_group_size)-prep_rx+1,(m+1)*(nx/mympi_group_size)+prep_rx
!          if(n<=0) then
!            nn=n+nx
!          else if(n>nx) then
!            nn=n-nx
!          else 
!            nn=n
!          endif
!          call substates4to8(substates,substates4(h*nx+nn))
!          psi_scatter(m*(ny_end-ny_begin+ry1+ry2)*(nx/mympi_group_size+2*prep_rx)+h*(nx/mympi_group_size+2*prep_rx)+n-m*(nx/mympi_group_size)+prep_rx)=substates
!        enddo
!        enddo
!        enddo
!        close(10)
!        close(11)
!      print*,'read end'
!      endif
!      call system_clock(timeB)
!      time_sum=time_sum+timeB-timeA
!      call MPI_Barrier(SplitWorld,ierr)
!      call MPI_Scatter(psi_scatter,(ny_end-ny_begin+ry1+ry2)*(nx/mympi_group_size+2*prep_rx),type_substates_MPI,A_local(1,1,j),(ny_end-ny_begin+ry1+ry2)*(nx/mympi_group_size+2*prep_rx),type_substates_MPI,0,SplitWorld,ierr)
!    endif
!    enddo
!    enddo
!    enddo

    if(key==0) then
!      write(filename,'(a,I4.4,a,I2.2,a,I2.2,a)') '../../output/analysed_state_SLA-',year,'-',month,'-',day,'_sub-dvg-0608mod-OC.uf'
!      write(filename,'(a,I4.4,a,I4.4,I2.2,I2.2,a)') '../../../',year,'/static_ensemble_25km_',year,month,day,'.uf'
      write(filename,'(a,I4.4,I2.2,I2.2,a)') '../input/forecast/dvg-static_ensemble_25km_',year,month,day,'_OC.uf'
      open(10,file=trim(filename),status='old',form='unformatted',access='stream',action='read')
      pos_begin=ny_begin-ry1
      read(10,pos=10+(pos_begin*nx*local_ndim*4)) substates4
      do m=0,mympi_group_size-1
      do h=0,ny_end-ny_begin+ry1+ry2-1
      do n=m*(nx/mympi_group_size)-prep_rx+1,(m+1)*(nx/mympi_group_size)+prep_rx
        if(n<=0) then
          nn=n+nx
        else if(n>nx) then
          nn=n-nx
        else 
          nn=n
        endif
        call substates4to8(substates,substates4(h*nx+nn))
        psi_scatter(m*(ny_end-ny_begin+ry1+ry2)*(nx/mympi_group_size+2*prep_rx)+h*(nx/mympi_group_size+2*prep_rx)+n-m*(nx/mympi_group_size)+prep_rx)=substates
      enddo
      enddo
      enddo
      close(10)
    endif
    call MPI_Barrier(SplitWorld,ierr)
    call MPI_Scatter(psi_scatter,(ny_end-ny_begin+ry1+ry2)*(nx/mympi_group_size+2*prep_rx),type_substates_MPI,psi_local,(ny_end-ny_begin+ry1+ry2)*(nx/mympi_group_size+2*prep_rx),type_substates_MPI,0,SplitWorld,ierr) 
    
! clim
!    allocate(daily_mean(NDperY))
!     do i_day = 1, NDperY
!         do j = 1, ny_len
!             do i = 1, nx_len
!              substates= assign_states8(0.0d0)
!              do i_yr = 1, NYEAR
!                  i_ens = (i_yr - 1) * NDperY + i_day
!                  substates = add_states8(substates, A_local(i, j, i_ens))
!              end do
!              substates = multi_states8(substates, 1.0d0 / dble(NYEAR))
!              do i_yr = 1, NYEAR
!                  i_ens = (i_yr - 1) * NDperY + i_day
!                  A_local(i, j, i_ens) = subtract_states8(A_local(i, j, i_ens), substates)
!              end do
!          end do
!      end do
!  end do
!  deallocate(daily_mean)

    call system_clock(time2)
    print*,'myid',myid,'time2-time1 is',time2-time1,'time_sum',time_sum,'time_read',time_read
    if(key==0) then
      deallocate(substates4)
      deallocate(psi_scatter)   
    endif
end subroutine read_A_psi

subroutine read_obs_grid(myid,numprocs,ierr)
  use m_get_obs
  use m_obs_reject
  use m_get_mod_grid
  use m_get_mean_ssh
  use m_get_mask
  implicit none
  integer, intent(in) :: myid,numprocs,ierr
  logical exobs
  integer time1,time2,time3,time4,time5
  real sla_diff
  character(len=3) :: sla_ssh
  integer :: nrobs_each(4)
  logical :: exobs_each(4)
  integer :: total_obs, current_idx, i
  
  total_obs = 0
  nrobs_each = 0
  open(103,file='date_type.dat',form='formatted',status='old',action='read')
    read(103,*) year, month, day                                              
  close(103)

  open (10,file='assimilation.in')
  read (10,*) radius_sla
  read (10,*) radius_sal
  read (10,*) radius_tem
  read (10,*) radius_sst
  read (10,'(l1)') vert_local
  read (10,'(l1)') l_local_analysis
  read (10,'(l1)') l_global_analysis
  read (10,*) anal_method
  read (10,*) alpha
  read (10,*) anal_solver
  read (10,*) assim_type
  read (10,*) sla_ssh
  close(10)
  
  open (10,file='sla_diff')
  read (10,*) sla_diff
  close(10)

  write(filenames(1),'(a,I4.4,a,I2.2,a,I2.2,a)')    &
  '../input/obs_sla_uf_025/obs_SLA-',year,'-',month,'-',day,'.uf'
  write(filenames(2),'(a,I4.4,a,I2.2,a,I2.2,a)')    &
  '../input/obs_sst_uf_025/obs_SST-',year,'-',month,'-',day,'.uf'
  write(filenames(3),'(a,I4.4,a,I2.2,a,I2.2,a)')    &
  '../input/obs_ts_uf_025/obs_TEM-',year,'-',month,'-',day,'.uf'
  write(filenames(4),'(a,I4.4,a,I2.2,a,I2.2,a)')    &
  '../input/obs_ts_uf_025/obs_SAL-',year,'-',month,'-',day,'.uf'
  do i = 1, 4
     if (assim_type(i:i) == '1') then
        nrobs_each(i) = get_nrobs_d(filenames(i), exobs_each(i))
        if (exobs_each(i)) then
           total_obs = total_obs + nrobs_each(i)
        else
           print *, 'Warning: Flag set for obs type ', i, ' but file not found.'
        end if
     end if
  end do

  if (total_obs == 0) then
     print *, 'No observations selected or found on this day.'
     return
  end if

  nrobs0 = total_obs
  if (allocated(obs0)) deallocate(obs0)
  allocate(obs0(nrobs0))
  
  current_idx = 1
  do i = 1, 4
     if (assim_type(i:i) == '1' .and. nrobs_each(i) > 0) then
        call get_obs_d(obs0(current_idx : current_idx + nrobs_each(i) - 1), &
                       nrobs_each(i), filenames(i))
        current_idx = current_idx + nrobs_each(i)
     end if
  end do
  write(outfile,'(a,a,a,I4.4,a,I2.2,a,I2.2,a)') &
       '../output/analysed_state_', trim(assim_type), '-', year, '-', month, '-', day, '_sub-OC.uf'

  call system_clock(time3)
  call get_mod_grid(modlon,modlat,depths,levels)
  allocate (mean_ssh(nx,ny))
  allocate (assim_mask(nx,ny))
  mean_ssh(:,:)=0.0
  call get_mean_ssh(mean_ssh)

  mean_ssh=mean_ssh+sla_diff  !zyz auto diff with spatial distribution
  call get_mask(assim_mask)
  call system_clock(time4)
  call obs_reject1(myid,obs0,nrobs0,psi_local,mean_ssh,depths,levels,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2,nrobs_local,nrobs)
  print*,'obs_reject1 success',myid,nrobs,nrobs_local
  allocate(obs(nrobs))
  allocate(obs_local(nrobs_local))
  call obs_reject2(myid,obs0,nrobs0,obs_local,nrobs_local,obs,nrobs)
!  print*,'obs_reject2 success',nrobs,nrobs_local
  call system_clock(time5)
  deallocate(obs)
  deallocate(obs0)
!  print*,myid,'get_nrobs_d',time2-time1,time3-time2,time4-time3,time5-time4
end subroutine read_obs_grid

subroutine perform_enoi_data(myid)
   use m_prep_4_EnOI

   implicit none
   integer, intent(in) :: myid
   integer i,time1,time2,time3,time4,time5
   logical exobs
   logical, allocatable :: local_obs(:) 

   print*,'EnOI: Start calculations of input to the analysis'
   !allocate(local_obs(nrobs_local))
   !call count_S_k(myid,obs_local,nrobs_local,S_k,local_obs,nx_begin,ny_begin,nx_local,ny_local)
   allocate(D(nrobs_local))
   allocate(S_local(nrens,nrobs_local))
   call system_clock(time4)
   call prep_4_EnOI(myid,obs_local,nrobs_local,A_local,psi_local,mean_ssh,nrens,depths,levels,D,S_local,alpha,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2)
   !deallocate(local_obs)
   call system_clock(time5)
end subroutine perform_enoi_data

subroutine perform_enoi_analysis(myid)
   use m_local_analysis
   use mod_states
   implicit none
   integer, intent(in) :: myid
   integer i,j 

   if(l_local_analysis) then
      print *, 'EnOI: Computing local analysis'
      call local_analysis(myid,nrens,obs,nrobs,obs_local,nrobs_local,A_local,psi_local,radius_sla,radius_sal,radius_tem,radius_sst,modlon,modlat,depths,D,S_local,local_ndim,alpha,nx_begin,nx_end,ny_begin,ny_end,prep_rx,ry1,ry2,assim_mask,anal_method,anal_solver,levels,vert_local)
   endif
   print *,myid,'EnOI: Analysis done'
end subroutine perform_enoi_analysis

subroutine write_result(myid,numprocs)
    use mod_states
    implicit none
    integer , intent(in) :: myid,numprocs
    type(sub_states4),allocatable :: substates4_100(:)
    type(sub_states4) substates4
    type(sub_states)  substates
    integer*8 pos_begin
    integer i,j

  open(103,file='date_type.dat',form='formatted',status='old',action='read')
    read(103,*) year, month, day                                              
  close(103)
!    print*,myid,'write result'
    open(10,file=outfile,form='unformatted',access='stream',action='write')
    if(myid==0) then
      write(rident,'(a,I4.4,I2.2,I2.2)') '_',year,month,day
      write(10,pos=1) rident
    endif
    
    do j=ny_begin+1,ny_end
    do i=nx_begin+1,nx_end
      substates=psi_local(i-nx_begin+prep_rx,j-ny_begin+ry1) 
      call substates8to4(substates4,substates)
      pos_begin=j-1
      write(10,pos=10+(pos_begin*nx+i-1)*local_ndim*4) substates4
    enddo
    enddo
    close(10)
end subroutine write_result

subroutine perform_enoi_test
   use m_prep_4_EnOI
   implicit none
  deallocate(A_local)
  deallocate(mean_ssh)
  deallocate(assim_mask)
  deallocate(psi_local)
  deallocate(S_local)
  deallocate(D)
  print*,'EnOI: Finished'
end subroutine perform_enoi_test
end module enoi
