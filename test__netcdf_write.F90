!! ----------------------------------------------------------------------------------------------------------------------------- !!
!> 
!! An example program to write netcdf files by calling m_netcdf module subroutines
!!
!<
!! --
program test__netcdf_write

  use m_netcdf
  
  implicit none

  integer, parameter :: STDERR = 0
  integer, parameter :: PREC = 8

  !! parameters for dummy data
  integer, parameter :: NLON =  360
  integer, parameter :: NLAT =  180
  real,    parameter :: DT   =  0.5
  integer, parameter :: NT   =   20
  !! 
  real(PREC) :: lon(NLON), lat(NLAT)
  real(PREC) :: dat1(NLON,NLAT,NT)
  real(PREC) :: dat2(NLON,NLAT,NT)
  integer :: ncid, varid1, varid2
  real(PREC) :: t
  integer :: k
  !! ----

  call prepare_dummy_data( nlon, nlat, nt, lon, lat, dat1, dat2 )

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !!
  !! 1. An easy way to generate single value netcdf file by calling a wrapper routine. 
  !!    arguments: filename, title, nx, ny, x(1:nx), y(1:ny), z(1:nx,1:ny)
  !!               x-label, y-label, z-label, x-unit, y-unit, z-unit
  !! ---- 
  call netcdf__write_grd( 'test_1.nc', 'test', NLON, NLAT, lon, lat, dat1(:,:,1), &
                          'longitude', 'latitude', 'amplitude', 'degree-east', 'degree-north', 'm' )


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !!
  !! 2. 2D time-independent, multiple dataset with some attributes
  !!
  !! ---- 

  !! first defines the file. with x&y-axes and units. It returns file id as ncid variable
  !! arguments: filename, title, nx, ny, x(1:nx), y(1:ny), xlabel, ylabel, tlabel, xunit, yunit, zunit, ID
  call netcdf__create( 'test_2.nc', 'test', NLON, NLAT, lon, lat, &
                        'longitude', 'latitude', 'degree-east', 'degree-north', ncid )

  !! define data fields with names and units. 
  !! Data are distinguished by variable ids (varid1, varid2)
  call netcdf__def_data( ncid, 'amplitude1', 'm', varid1)
  call netcdf__def_data( ncid, 'amplitude2', 'm', varid2)

  !! one can add any kinds of header information as "attribute" of netcdf file
  call netcdf__add_attribute( ncid, "ascii   attribute", "12345" ) 
  call netcdf__add_attribute( ncid, "real    attribute", 123.45  )
  call netcdf__add_attribute( ncid, "integer attribute", 12345   )

  !! add data to netcdf file
  call netcdf__add_data( ncid, varid1, dat1(:,:,1) )
  call netcdf__add_data( ncid, varid2, dat2(:,:,1) )
  
  !! close the file
  call netcdf__close( ncid )


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !!
  !! 3. 2D time-dependent, multiple dataset 
  !!
  !! ----

  !! first defines the file. with x,y and time-axes and units. It returns file id as ncid variable
  !! arguments: filename, title, nx, ny, x(1:nx), y(1:ny), xlabel, ylabel, tlabel, xunit, yunit, zunit, ID
  !! number of time grid is not necessary 
  call netcdf__create( 'test_3.nc', 'test', NLON, NLAT, lon(1:NLON), lat(1:NLAT), &
                        'longitude', 'latitude', 'time', 'degree-east', 'degree-north', 'sec', ncid )
  
  !! define data fields with names and units.
  !! Data are distinguished by variable ids (varid1, varid2)
  call netcdf__def_data( ncid, 'amplitude1','m', varid1)
  call netcdf__def_data( ncid, 'amplitude2','m', varid2 )

  !! one can add any kinds of header information as "attribute" of netcdf file
  call netcdf__add_attribute( ncid, "ascii   attribute", "12345" ) 
  call netcdf__add_attribute( ncid, "real    attribute", 123.45  )
  call netcdf__add_attribute( ncid, "integer attribute", 12345   )

  !! time loop
  t = 0
  do k=1, NT
     !! add space snapshot of the data at time t with time-index k
     call netcdf__add_data( ncid, varid1, k, t, dat1(:,:,k) )
     call netcdf__add_data( ncid, varid2, k, t, dat2(:,:,k) )
     t = t + dt
  end do

  !! close the file
  call netcdf__close( ncid )
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  

  
contains 

  !!
  subroutine prepare_dummy_data( nlon, nlat, nt, lon, lat, dat1, dat2 )
    integer, intent(in)  :: nlon, nlat, nt
    real(PREC),    intent(out) :: lon(nlon), lat(nlat)
    real(PREC),    intent(out) :: dat1(nlon,nlat,nt), dat2(nlon,nlat,nt)
    !!
    
    real(PREC),    parameter :: LON0 =    0.0
    real(PREC),    parameter :: LON1 =  359.0
    real(PREC),    parameter :: LAT0 =  -90.0
    real(PREC),    parameter :: LAT1 =   90.0
    real(PREC), parameter :: PI = atan(1.0)*4
    
    real    :: theta
    real    :: phi
    real    :: dlon, dlat
    integer :: i, j, k 
    !! ----
    
    !! axis data preparation
    dlon = ( LON1-LON0 ) / (nlon-1)
    do i=1, nlon
       lon(i) = (i-1) * dlon + LON0
    end do
    dlat = ( LAT1-LAT0 ) / (nlat-1)
    do i=1, nlat
       lat(i) = (i-1) * dlat + LAT0
    end do
    
    !! data
    !! monochromatic time-dependent oscillation of spherical harmonics Y_10^5
    !! real and imaginary parts are stored to dat1 and dat2
    do k=1, nt
       do j=1, nlat
          do i=1, nlon
             theta = ( 90 - lat(j) ) * PI / 180  !! colatitude
             phi   = lon(i) * PI / 180
             
             dat1(i,j,k) = cos( real(k-1)/real(NT)*(2*PI) ) &
                  * ( -3.0 / 256.0 ) * sqrt( 1001.0 / PI ) &
                  * cos( 5 * phi ) &
                  * sin( theta ) **5 * ( 323 * cos(theta)**5 - 170 * cos(theta)**3 + 15 * cos( theta ) )
             
             dat2(i,j,k) = cos( real(k-1)/real(NT)*(2*PI) ) &
                  * ( -3.0 / 256.0 ) * sqrt( 1001.0 / PI ) &
                  * sin( 5 * phi ) &
                  * sin( theta ) **5 * ( 323 * cos(theta)**5 - 170 * cos(theta)**3 + 15 * cos( theta ) )
             
          end do
       end do
    end do
  
  end subroutine prepare_dummy_data
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end program test__netcdf_write
!! ----------------------------------------------------------------------------------------------------------------------------- !!
