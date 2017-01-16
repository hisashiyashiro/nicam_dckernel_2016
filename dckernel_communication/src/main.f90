!-------------------------------------------------------------------------------
!
!+  Program dynamics kernel driver (setup operator)
!
!-------------------------------------------------------------------------------
program dckernel_communication
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_process, only: &
     PRC_MPIstart,    &
     PRC_LOCAL_setup, &
     PRC_MPIfinish
  use mod_adm, only: &
     ADM_setup
  use mod_comm, only: &
     COMM_setup
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  integer :: comm_world
  integer :: myrank
  logical :: ismaster
  !=============================================================================

  !---< MPI start >---
  call PRC_MPIstart( comm_world ) ! [OUT]

  !---< STDIO setup >---
  call IO_setup( 'NICAM-DC',         & ! [IN]
                 'communication.cnf' ) ! [IN]

  !---< Local process management setup >---
  call PRC_LOCAL_setup( comm_world, & ! [IN]
                        myrank,     & ! [OUT]
                        ismaster    ) ! [OUT]

  if( ismaster ) write(*,*) '##### [KERNEL] dckernel_communication #####'

  !---< Logfile setup >---
  call IO_LOG_setup( myrank,  & ! [IN]
                     ismaster ) ! [IN]

  !---< profiler module setup >---
  call PROF_setup

  !#############################################################################
  call PROF_setprefx('MAIN')
  call PROF_rapstart('Kernel_ALL',0)

  if( ismaster ) write(*,*) '##### start kernel #####'

  !---< admin module setup >---
  call ADM_setup

  !---< comm module setup >---
  call COMM_setup

  if( ismaster ) write(*,*) '##### finish kernel #####'

  call PROF_rapend('Kernel_ALL',0)
  !#############################################################################

  call PROF_rapreport

  !--- finalize all process
  call PRC_MPIfinish

  stop
end program dckernel_communication
!-------------------------------------------------------------------------------
