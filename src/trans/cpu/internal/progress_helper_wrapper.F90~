module progress_thread

    use mpi
    use iso_c_binding
    implicit none

interface

    subroutine start_mpi_helper() bind(C, name='start_mpi_helper')
        import
    end subroutine start_mpi_helper
    subroutine stop_mpi_helper() bind(C, name='stop_mpi_helper')
        import
    end subroutine stop_mpi_helper

    subroutine pt_reqset_register(count, array_of_requests, flags, gid, ierr) bind(C, name='pt_reqset_register_f')

        import 
        integer(C_INT), VALUE, INTENT(IN) :: count, flags
        integer(C_INT), INTENT(OUT) :: gid
        integer(C_INT), OPTIONAL, INTENT(OUT) :: ierr
        integer(C_INT), INTENT(IN) :: array_of_requests(*)

    end subroutine pt_reqset_register

    subroutine pt_reqset_start(gid, ierr) bind(C, name='pt_reqset_start')
        import 
        integer(C_INT), VALUE, INTENT(IN) :: gid
        integer(C_INT), OPTIONAL, INTENT(OUT) :: ierr
      
    end subroutine pt_reqset_start

    subroutine pt_reqset_wait_f(gid, ierr) bind(C, name='pt_reqset_wait_f')
        import
        integer(C_INT), INTENT(IN) :: gid
!        integer(C_INT), OPTIONAL, INTENT(OUT) :: mpistatuses(:)        
        integer(C_INT), OPTIONAL, INTENT(OUT) :: ierr
    end subroutine pt_reqset_wait_f

!    subroutine pt_reqset_test_f(gid, flag, mpistatuses, ierr) bind(C, name='pt_reqset_test_f')
    subroutine pt_reqset_test_f(gid, flag, ierr) bind(C, name='pt_reqset_test_f')
        use mpi_f08
        integer(C_INT), INTENT(IN) :: gid
        integer(C_INT), INTENT(OUT) :: flag
        integer(C_INT), INTENT(OUT) :: ierr
!        type(MPI_STATUS), INTENT(OUT) :: mpistatuses(*)
    end subroutine pt_reqset_test_f

end interface

contains

    subroutine pt_reqset_wait(gid, mpistatuses, ierr)
        integer(C_INT), INTENT(IN) :: gid
        integer(C_INT), OPTIONAL, INTENT(OUT) :: ierr
        integer(C_INT), OPTIONAL, INTENT(OUT) :: mpistatuses(:)
        integer(C_INT) :: c_err

!        call pt_reqset_wait_f(gid, mpistatuses, c_err)
        call pt_reqset_wait_f(gid, c_err)
        if(present(ierr)) ierr = c_err
    end subroutine pt_reqset_wait

    subroutine pt_reqset_test(gid, flag, mpistatuses, ierr)

      use mpi_f08

      integer(C_INT), INTENT(IN) :: gid
        integer(C_INT), INTENT(OUT) :: flag
        integer(C_INT), OPTIONAL, INTENT(OUT) :: ierr
        type(MPI_STATUS), OPTIONAL, INTENT(OUT) :: mpistatuses(*)
        integer(C_INT) :: c_err
!        integer(C_INT), POINTER :: mpistats(:)

!        if(.not. present(mpistatuses)) then
!           nullify(mpistats)
!!        else
!           mpistats => mpistatuses
!        endif

        call pt_reqset_test_f(gid, flag, c_err)
        if(present(ierr)) ierr = c_err
    end subroutine pt_reqset_test

end module progress_thread
