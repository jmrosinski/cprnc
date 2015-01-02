module utils
  use netcdf
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! handle_err:
!   utility function to print an error msg and stop on netcdf error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine handle_err (ret)
    integer, intent(in) :: ret

    if (ret /= NF90_NOERR) then
      write(6,*) nf90_strerror (ret)
      stop 999
    end if
    return
  end subroutine handle_err
end module utils
