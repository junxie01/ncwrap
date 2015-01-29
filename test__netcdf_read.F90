!! ----------------------------------------------------------------------------------------------------------------------------- !!
!> 
!! An example program to read netcdf files by calling m_netcdf module subroutines
!!
!! This program tries to read netcdf files created by a program test__netcdf_write.F90 .
!! Please execute this after test__netcdf_write.F90.
!<
!! ----
program test__netcdf_read
 
  use m_netcdf
  implicit none

  integer, parameter :: STDERR = 0
  integer, parameter :: PREC = 8
  real(PREC), allocatable :: lon(:), lat(:), tim(:)
  real(PREC), allocatable :: dat1(:,:), dat2(:,:)
  integer :: nlon, nlat, nt
  character(16) :: xname, yname, zname, xunit, yunit, zunit, title
  character(16) :: zname2(10), tname, zunit2(10), tunit
  integer :: i, j, k
  integer :: ncid
  integer :: ndim, nvar
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !!
  !! 1. An easy way to read single value netcdf file by calling a wrapper routine.
  !!
  !! --
  call netcdf__read_grd( 'test_1.nc', title, nlon, nlat, lon, lat, dat1, xname, yname, zname, xunit, yunit, zunit )
  write(STDERR,*) "file = test_1.nc"
  write(STDERR,*) "title = ", trim(title)
  write(STDERR,*) "(nlon, nlat) = (", nlon, ",", nlat, ")"
  write(STDERR,*) "export data to fort.101"
  do j=1, nlat
     do i=1, nlon
        write(101,*) lon(i), lat(j), dat1(i,j)
     end do
     write(101,*)
  end do

  deallocate( lon, lat, dat1 )
  write(STDERR,*)

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !!
  !! 2. Read multiple, time-independent dataset
  !!
  !! --

  !! inquire header information of netcdf files
  call netcdf__inq( 'test_2.nc', title, nlon, nlat, nt, ndim, nvar, xname, yname, tname, zname2, xunit, yunit, tunit, zunit2 )

  write(STDERR,*) "file = test_2.nc"
  write(STDERR,*) "dimension = ", ndim
  write(STDERR,*) "title = ", trim(title)
  write(STDERR,*) "(nlon, nlat) = (", nlon, ",", nlat, ")"
  do i=1, nvar
     write(STDERR,*) "Variable ", i, trim(zname2(i))
  end do
  allocate( lon(nlon), lat(nlat), dat1(nlon,nlat), dat2(nlon,nlat) )

  !! read data
  call netcdf__open('test_2.nc', nlon, nlat, lon, lat, ncid )
  call netcdf__inq( 'test_2.nc', title, nlon, nlat, nt, ndim, nvar, xname, yname, tname, zname2, xunit, yunit, tunit, zunit2 )
  write(STDERR,*) "file = test_2.nc"
  write(STDERR,*) "dimension = ", ndim
  write(STDERR,*) "title = ", trim(title)
  write(STDERR,*) "(nlon, nlat) = (", nlon, ",", nlat, ")"

  call netcdf__read_data( ncid, trim(zname2(1)), dat1 )
  call netcdf__read_data( ncid, trim(zname2(2)), dat2 )
  call netcdf__close( ncid )
  write(STDERR,*) "export data to fort.102"
  do j=1, nlat
     do i=1, nlon
        write(102,'(4F15.5)') lon(i), lat(j), dat1(i,j), dat2(i,j)
     end do
     write(102,*)
  end do
  
  deallocate( lon, lat, dat1, dat2 )
  write(STDERR,*)

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !!
  !! 3. Read time-dependent dataset
  !!
  !! --
  call netcdf__inq( 'test_3.nc', title, nlon, nlat, nt, ndim, nvar, xname, yname, tname, zname2, xunit, yunit, tunit, zunit2 )
  allocate( dat1(nlon,nlat), dat2(nlon,nlat) )
  call netcdf__open('test_3.nc', nlon, nlat, nt, lon, lat, tim, ncid )
  do k=1, nt
     call netcdf__read_data( ncid, zname2(1), k, dat1 )
     call netcdf__read_data( ncid, zname2(1), k, dat2 )
     do j=1, nlat
        do i=1, nlon
           write(200+k,'(4F15.5)') lon(i), lat(j), dat1(i,j), dat2(i,j)
        end do
        write(200+k,*)
     end do
  end do
  call netcdf__close(ncid)

end program test__netcdf_read
!! ----------------------------------------------------------------------------------------------------------------------------- !!
