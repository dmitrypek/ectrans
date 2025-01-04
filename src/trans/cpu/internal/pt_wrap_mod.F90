module progress_thread

use mpi
use iso_c_binding
implicit none
!include 'mpif.h'

interface

      ! integer(C_INT) 
   subroutine start_mpi_helper() bind(C,name='start_mpi_helper')
     import
   end subroutine start_mpi_helper
   subroutine stop_mpi_helper() bind(C,name='stop_mpi_helper')
     import
   end subroutine stop_mpi_helper

   subroutine pt_reqset_register(count,array_of_requests,flags,gid,status) bind(C,name='pt_reqset_register_f')

      import 
      integer(C_INT) :: count,flags,gid,status
      integer :: array_of_requests(*)

    end subroutine pt_reqset_register

    subroutine pt_reqset_start(gid,status) bind(C,name='pt_reqset_start')
      import 
      integer(C_INT) :: gid,status
      
    end subroutine pt_reqset_start

    subroutine pt_reqset_wait(gid, status) bind(C,name='pt_reqset_wait')
      import
      integer(C_INT) :: gid,status  !,mpistatuses(*,MPI_STATUS_SIZE)
    end subroutine pt_reqset_wait

    subroutine pt_reqset_test(gid, flag,status) bind(C,name='pt_reqset_test')
      import
      integer(C_INT) :: gid,status,flag  !,mpistatuses(*,MPI_STATUS_SIZE)
    end subroutine pt_reqset_test

   end interface
 end module progress_thread
 
