module stats
  use netcdf
  use utils
  implicit none

  private
  public :: print_heading, printstats1, printstats2

contains

  subroutine print_heading (nfiles, t)
    integer, intent(in) :: nfiles
    integer, intent(in) :: t

    write(6,98)
    write(6,99) t
98  format ('-----------------------------------------------------')
99  format ('Fields for time slice=', i3,/)
    if (nfiles == 1) then
      write(6,100)
      write(6,101)
      write(6,102)
      write(6,105)
    else
      write(6,200)
      write(6,201)
      write(6,202)
      write(6,203)
      write(6,204)
      write(6,205)
    end if
100 format ('                 max_info        min_info')
101 format ('name  size      (     ipn1, k1) (     ipn1, k1)')
102 format ('                 case1_values    case1_values')
105 format ('Additional info',/,/)

200 format ('                 max_info        min_info        absolute_diff_info     relative_diff_info')
201 format ('name  size      (     ipn1, k1) (     ipn1, k1)  diffmax (     ipn,  k) rdiffmax  (     ipn,  k)')
202 format ('      num_diffs  case1_values    case1_values             case1_values             case1_values')
203 format ('                 case2_values    case2_values             case2_values             case2_values')
204 format ('                (     ipn2, k2) (     ipn2, k2)')
205 format ('Additional_info',/,/)
  end subroutine print_heading

  subroutine printstats1 (ncid, name, varid, t, nip)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer, intent(in) :: varid
    integer, intent(in) :: t
    integer, intent(in) :: nip

    integer :: ndims_var
    integer :: dimids(NF90_MAX_DIMS)
    integer :: start(3)
    integer :: count(3)
    integer :: nz
    integer :: ipnmax, kmax
    integer :: ipnmin, kmin
    real, allocatable :: values(:,:)
    real :: avg

    call handle_err (nf90_inquire_variable (ncid, varid, ndims=ndims_var, dimids=dimids))

    if (ndims_var == 2) then
      nz = 1
      allocate (values(nip,1))
      start(1:2) = (/1,t/)
      count(1:2) = (/nip,1/)
      call handle_err (nf90_get_var (ncid, varid, values(:,1), start(1:2), count(1:2)))
    else if (ndims_var == 3) then
      call handle_err (nf90_inquire_dimension (ncid, dimids(2), len=nz))
      allocate (values(nip,nz))
      start(1:3) = (/1,1,t/)
      count(1:3) = (/nip,nz,1/)
      call handle_err (nf90_get_var (ncid, varid, values(:,:), start(1:3), count(1:3)))
    else
      write(6,*)'printstats1: ndims_var=', ndims_var
      stop 999
    end if

    call get_stats (values, ipnmax, kmax, ipnmin, kmin, avg)

    write(6,100) name, nz*nip, ipnmax, kmax, ipnmin, kmin
    write(6,101) values(ipnmax,kmax), values(ipnmin,kmin)
    write(6,102) avg

100 format (a4,1x,i9,' (',i9,',',i3,')',' (',i9,',',i3,')')
101 format (1p,15x,       e14.7,2x,            e14.7)
102 format ('avg abs field values: ', 1pe14.7,/)

    deallocate (values)
  end subroutine printstats1

  subroutine printstats2 (ncid, name, varid, t, nip, varid2, totdiffs)
    integer, intent(in) :: ncid(2)
    character(len=*), intent(in) :: name
    integer, intent(in) :: varid
    integer, intent(in) :: t
    integer, intent(in) :: nip
    integer, intent(in) :: varid2
    integer(8), intent(inout) :: totdiffs

    integer :: ndims_var(2)
    integer :: dimids(NF90_MAX_DIMS,2)
    integer :: start(3)
    integer :: count(3)
    integer :: nz(2)
    integer :: ipnmax(2), kmax(2)
    integer :: ipnmin(2), kmin(2)
    integer :: ndiff
    integer :: ipn_diffmax, k_diffmax
    integer :: ipn_rdiffmax, k_rdiffmax
    integer :: n
    real :: diffmax, rdiffmax
    real :: rms
    real :: rdbar
    real :: digbar
    real :: digworst
    real, allocatable :: values(:,:,:)
    real :: avg(2)

    call handle_err (nf90_inquire_variable (ncid(1), varid,  ndims=ndims_var(1), dimids=dimids(:,1)))
    call handle_err (nf90_inquire_variable (ncid(2), varid2, ndims=ndims_var(2), dimids=dimids(:,2)))

    if (ndims_var(1) /= ndims_var(2)) then
      write(6,*) 'printstats2 name=', name,' ndims_var=',ndims_var(:),' do not match: skipping'
      return
    end if

    do n=1,max(ndims_var(1), ndims_var(2))
      if (dimids(n,1) /= dimids(n,2)) then
        write(6,*) 'printstats2 name=', trim(name),' n=',n,' dimids=',dimids(n,1),dimids(n,2),' do not match: skipping'
        return
      end if
    end do

    if (ndims_var(1) == 2) then
! Field is (nip,time): Use nz=1
      nz(1) = 1
      allocate (values(nip,1,2))
      start(1:2) = (/1,t/)
      count(1:2) = (/nip,1/)
      call handle_err (nf90_get_var (ncid(1), varid,  values(:,1,1), start(1:2), count(1:2)))
      call handle_err (nf90_get_var (ncid(2), varid2, values(:,1,2), start(1:2), count(1:2)))
    else if (ndims_var(1) == 3) then
! Field is (nip,nz,time): Discover nz
      call handle_err (nf90_inquire_variable (ncid(1), varid, dimids=dimids(:,1)))
      call handle_err (nf90_inquire_dimension (ncid(1), dimids(2,1), len=nz(1)))

      call handle_err (nf90_inquire_variable (ncid(2), varid2, dimids=dimids(:,2)))
      call handle_err (nf90_inquire_dimension (ncid(2), dimids(2,1), len=nz(2)))

      if (nz(1) /= nz(2)) then
        write(6,*) 'printstats2 name=', name,' nz=', nz(:), ' do not match: skipping'
        return
      end if

      allocate (values(nip,nz(1),2))
      start(1:3) = (/1,1,t/)
      count(1:3) = (/nip,nz(1),1/)
      call handle_err (nf90_get_var (ncid(1), varid, values(:,:,1), start(1:3), count(1:3)))
      call handle_err (nf90_get_var (ncid(2), varid2, values(:,:,2), start(1:3), count(1:3)))
    else
      write(6,*)'printstats2 ERROR: ndims_var=', ndims_var(:)
      stop 999
    end if

    call get_stats (values(:,:,1), ipnmax(1), kmax(1), ipnmin(1), kmin(1), avg(1))
    call get_stats (values(:,:,2), ipnmax(2), kmax(2), ipnmin(2), kmin(2), avg(2))
    call get_diffstats (values, ndiff, diffmax, ipn_diffmax, k_diffmax, rdiffmax, &
                        ipn_rdiffmax, k_rdiffmax, rms, rdbar)
    if (ndiff > 0) then
      digbar = log10 (1./rdbar)
      digworst = log10 (1./rdiffmax)

      write(6,100) name, nz(1)*nip, ipnmax(1), kmax(1), ipnmin(1), kmin(1), diffmax, &
                   ipn_diffmax, k_diffmax, rdiffmax, ipn_rdiffmax, k_rdiffmax
      write(6,101) ndiff, values(ipnmax(1),  kmax(1),  1), values(ipnmin(1),   kmin(1),   1), &
                          values(ipn_diffmax,k_diffmax,1), values(ipn_rdiffmax,k_rdiffmax,1) 
      write(6,102) values(ipnmax(2),  kmax(2),  2), values(ipnmin(2),   kmin(2),   2), &
                   values(ipn_diffmax,k_diffmax,2), values(ipn_rdiffmax,k_rdiffmax,2)
      write(6,103)                  ipnmax(1), kmax(1), ipnmin(1), kmin(1)
      write(6,104) avg(1), name, rms, rdbar
      write(6,105) avg(2), digbar, digworst
      totdiffs = totdiffs + ndiff
    else
      write(6,200) name, nz(1)*nip, ipnmax(1), kmax(1), ipnmin(1), kmin(1), diffmax
      write(6,201) ndiff, values(ipnmax(1),  kmax(1),  1), values(ipnmin(1),   kmin(1),   1)
      write(6,202)        values(ipnmax(2),  kmax(2),  2), values(ipnmin(2),   kmin(2),   2)
      write(6,203)                  ipnmax(1), kmax(1), ipnmin(1), kmin(1)
      write(6,204) avg(1), name, rms
      write(6,205) avg(2)
    end if

100 format (a4,1x,i9,1x,' (',i9,',',i3,')',' (',i9,',',i3,') ',1p,e8.1,' (',i9,',',i3,') ',e8.1,' (',i9,',',i3,') ')
101 format (5x,   i9,2x,  1p,e14.7,          2x,e14.7,                  11x,e14.7,               11x,e14.7)
102 format (16x,          1p,e14.7,          2x,e14.7,                  11x,e14.7,               11x,e14.7)
103 format (15x,        ' (',i9,',',i3,')',' (',i9,',',i3,') ')
104 format ('avg abs field values:',1p,e14.7,'    rms ',a4,':',e10.3,' avg rel diff(all points):',e14.7)
105 format (21x,                    1p,e14.7,23x,                    ' avg decimal digits(ndif):',0p,f4.1,' worst:',f4.1,/)

200 format (a4,1x,i9,1x,' (',i9,',',i3,')',' (',i9,',',i3,') ',1p,e8.1)
201 format (5x,   i9,2x,  1p,e14.7,            2x,e14.7)
202 format (16x,          1p,e14.7,            2x,e14.7)
203 format (15x,        ' (',i9,',',i3,')',' (',i9,',',i3,') ')
204 format ('avg abs field values:',1p,e14.7,'    rms ',a4,':',e10.3)
205 format (21x,                    1p,e14.7,/)

    deallocate (values)
  end subroutine printstats2

  subroutine get_stats (values, ipnmax, kmax, ipnmin, kmin, avg)
    real, intent(in) :: values(:,:)
    integer, intent(out) :: ipnmax, kmax
    integer, intent(out) :: ipnmin, kmin
    real, intent(out) :: avg

    integer :: nz
    integer :: nip
    integer :: k, ipn
    real :: valmax, valmin

    nip = size(values,1)
    nz = size(values,2)
    valmax = -1.e38
    valmin = +1.e38
    avg = 0.

    do k=1,nz
      do ipn=1,nip
        avg = avg + abs (values(ipn,k))
        if (values(ipn,k) > valmax) then
          valmax = values(ipn,k)
          ipnmax = ipn
          kmax = k
        end if
        if (values(ipn,k) < valmin) then
          valmin = values(ipn,k)
          ipnmin = ipn
          kmin = k
        end if
      end do
    end do
    avg = avg / (nz*nip)
  end subroutine get_stats

  subroutine get_diffstats (values, ndiff, diffmax, ipn_diffmax, k_diffmax, rdiffmax, &
                            ipn_rdiffmax, k_rdiffmax, rms, rdbar)
    real, intent(in) :: values(:,:,:)    ! field values case 1 and case 2
    integer, intent(out) :: ndiff        ! number of differences
    integer, intent(out) :: ipn_diffmax  ! ipn loc of biggest difference
    integer, intent(out) :: k_diffmax    ! k loc of biggest difference
    integer, intent(out) :: ipn_rdiffmax ! ipn loc of biggest relative difference
    integer, intent(out) :: k_rdiffmax   ! k loc of biggest relative difference
    real, intent(out) :: diffmax         ! max abs difference
    real, intent(out) :: rdiffmax        ! max relative difference
    real, intent(out) :: rms             ! RMS difference for the entire field
    real, intent(out) :: rdbar           ! mean relative difference (across points with non-zero diffs)

    integer :: nz         ! number of vertical levels
    integer :: nip        ! number of horizontal points
    integer :: k, ipn     ! loop indices
    integer :: npoints    ! nz*nip
    real :: diff          ! difference
    real :: rdiff         ! relative difference
    real :: denom         ! denominator of expression

    nip = size(values,1)
    nz = size(values,2)

    diffmax      = 0.
    rdiffmax     = 0.
    rms          = 0.
    rdbar        = 0.
    ipn_diffmax  = 1
    k_diffmax    = 1
    ipn_rdiffmax = 1
    k_rdiffmax   = 1
    ndiff        = 0

    do k=1,nz
      do ipn=1,nip
        diff = abs (values(ipn,k,1) - values(ipn,k,2))
        rms = rms + diff*diff
        if (diff > 0.) then
          ndiff = ndiff + 1
          denom = 2.*(max(abs(values(ipn,k,1)), abs(values(ipn,k,2))))
          rdiff = diff / denom
        else
          rdiff = 0.
        end if
        rdbar = rdbar + rdiff
        if (diff > diffmax) then
          diffmax = diff
          ipn_diffmax = ipn
          k_diffmax = k
        end if

        if (rdiff > rdiffmax) then
          rdiffmax = rdiff
          ipn_rdiffmax = ipn
          k_rdiffmax = k
        end if
      end do
    end do
    npoints = nz * nip
    rms = sqrt (rms / npoints)
    if (ndiff > 0) then
      rdbar = rdbar / npoints
    else
      rdbar = 0.
    end if
  end subroutine get_diffstats
end module stats
