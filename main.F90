program do_enoi
  use mpi
  use mod_dimensions
  use enoi

    implicit none
    integer :: request,ierr, myid, numprocs,SplitWorld,status(MPI_STATUS_SIZE)
    integer :: mympi_group_size,color,key
    integer time1,time2,time3,time4,time5,time6,time7,time8,time9
    real :: t23, t35, t56, t68, t89
    real :: max_t23, max_t35, max_t56, max_t68, max_t89   
    print*,'begin'
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    mympi_group_size=numprocs/group
    color=myid/mympi_group_size
    key=mod(myid,mympi_group_size)
    call MPI_Comm_split(MPI_COMM_WORLD,color,key,SplitWorld,ierr)
    
    print*,myid,color,key
    call system_clock(time2)
    call read_A_psi(myid,numprocs,mympi_group_size,color,key,SplitWorld,ierr) 
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call system_clock(time3)
    call read_obs_grid(myid,numprocs,ierr)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call system_clock(time5)
    call perform_enoi_data(myid)
    call system_clock(time6)
    call perform_enoi_analysis(myid)
    print*,myid,'perform_enoi_analysis finished'
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
   
    print*,myid,'write_result begin'
    call system_clock(time8)
    call write_result(myid,numprocs)  
    !call perform_enoi_test
    call system_clock(time9)
    call MPI_Barrier(MPI_COMM_WORLD,ierr) 
!    print*,'myid is',myid,(time3-time2)/10000,(time5-time3)/10000,(time6-time5)/10000,(time8-time6)/10000,(time9-time8)/10000
    t23 = real(time3 - time2) / 10000.0
    t35 = real(time5 - time3) / 10000.0
    t56 = real(time6 - time5) / 10000.0
    t68 = real(time8 - time6) / 10000.0
    t89 = real(time9 - time8) / 10000.0

    call MPI_REDUCE(t23, max_t23, 1, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(t35, max_t35, 1, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(t56, max_t56, 1, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(t68, max_t68, 1, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(t89, max_t89, 1, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    if (myid == 0) then
        print*,'--- Maximum Execution Time Across All Procs (s) ---'
        print*,'Read_A_psi:  ', max_t23
        print*,'Read_obs:    ', max_t35
        print*,'EnOI_Prep:   ', max_t56
        print*,'EnOI_Analysis:', max_t68
        print*,'Write_Result: ', max_t89
    endif
    call MPI_FINALIZE(ierr)
end
