program cprnc
  use netcdf
  use utils
  use stats

  implicit none

  character(len=256) :: arg
  character(len=NF90_MAX_NAME) :: dimname
  character(len=NF90_MAX_NAME) :: name
  character(len=NF90_MAX_NAME) :: name2    ! for fields on file 2 not on file 1
  integer :: mode(2)
  integer :: ncid(2)
  integer :: ndims(2)
  integer :: ndims_var(2)                  ! number of dimensions for a variable
  integer :: nvars(2)
  integer :: natts(2)
  integer :: unlimdimid(2)
  integer :: nipdimid(2)
  integer :: ntime(2)
  integer :: varid(2)
  integer :: dimids(NF90_MAX_DIMS,2)
  integer :: nip(2)
  integer(8) :: totdiffs

  integer :: t                ! time index
  integer :: f                ! file index
  integer :: n                ! field index
  integer :: d                ! dimension index
  integer :: nfiles           ! 1 or 2

  nfiles = command_argument_count ()
  if (nfiles < 1 .or. nfiles > 2) then
    write(6,*)'Usage: cprnc file1.nc [file2.nc]'
    stop 999
  end if

  nipdimid(:) = -1
  do f=1,nfiles
    call get_command_argument (f, arg)
    call handle_err (nf90_open (arg, mode(f), ncid(f)))
    call handle_err (nf90_inquire (ncid(f), ndimensions=ndims(f), nvariables=nvars(f), &
                                   nattributes=natts(f), unlimiteddimid=unlimdimid(f)))
    write(6,100) f, trim(arg), nvars(f)
100 format ('File ',i1,' =',a,' nvars=',i3)

    call handle_err (nf90_inquire_dimension (ncid(f), unlimdimid(f), len=ntime(f)))
! Discover which dimension is named "nip"
    do d=1,ndims(f)
      call handle_err (nf90_inquire_dimension (ncid(f), d, name=dimname, len=nip(f)))
      if (trim (dimname) == 'nip') then
        nipdimid(f) = d
        exit
      end if
    end do

    if (nipdimid(f) == -1) then
      write(6,*) 'No nip dimension found for file=', f
      stop 999
    end if
  end do

! For now barf if these conditions are not all met. Not too hard to fix
! if any are ever not met

  if (nfiles == 2) then
    if (unlimdimid(1) /= unlimdimid(2)) then
      write(6,*) 'unlimdimid=', unlimdimid(:), ' do not match'
      stop 999
    end if
    if (nipdimid(1) /= nipdimid(2)) then
      write(6,*) 'nipdimid=', nipdimid(:), ' do not match'
      stop 999
    end if
    if (nip(1) /= nip(2)) then
      write(6,*) 'nip=', nip(:), ' do not match'
      stop 999
    end if
  end if

  do t=1,ntime(1)
    call print_heading (nfiles, t)
    totdiffs = 0
    do n=1,nvars(1)
      call handle_err (nf90_inquire_variable (ncid(1), n, name=name, ndims=ndims_var(1), &
                       dimids=dimids(:,1)))
! Only compare variables whose 1st dimension is nip and last dimension is time
! For now ASSUME that if var1 meets 1st and last criteria, so will var2
      if (dimids(1,1) == nipdimid(1) .and. dimids(ndims_var(1),1) == unlimdimid(1)) then
        if (nf90_inq_varid (ncid(2), name, varid(2)) == NF90_NOERR) then
          call printstats2 (ncid, name, n, t, nip(1), varid(2), totdiffs)
        else
          call printstats1 (ncid(1), name, n, t, nip(1))
        end if
      end if
    end do
!
! Now handle fields that are on the 2nd file but not the 1st (those also on the 1st were handled above)
!
    if (nfiles == 2) then
      do n=1,nvars(2)
        call handle_err (nf90_inquire_variable (ncid(2), n, name=name2, ndims=ndims_var(2), &
                         dimids=dimids(:,2)))
        if (nf90_inq_varid (ncid(1), name2, varid(1)) /= NF90_NOERR) then
          call printstats1 (ncid(2), name2, n, t, nip(2))
        end if
      end do
      write(6,101) t, totdiffs
101   format ('Total diffs time slice ',i3,' are ',i10)
    end if
  end do
  stop
end program cprnc
