!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Wrapper module for reading/writing netcdf data
!!
!! This module is activated only if _NETCDF is defined in preprocessor by setting -D_NETCDF as an compiler option
!! Requires netcdf library from unidata (version>=3.6)
!!
!! ----
module m_ncwrap

#ifdef _NETCDF
  use netcdf
#endif

  implicit none
  private
  save

  !! -- Constant
  integer,     parameter :: DP = selected_real_kind(13) !< Double Precision
  integer,     parameter :: SP = selected_real_kind(5)  !< Single Precision
  integer,     parameter :: STDERR  = 0 !< Standard Error
  
  !! Write
  public :: ncwrap__create
  public :: ncwrap__def_data
  public :: ncwrap__add_data
  public :: ncwrap__add_attribute

  !! Read
  public :: ncwrap__open
  public :: ncwrap__read_data

  !! Common
  public :: ncwrap__close
  public :: ncwrap__inq

  !! Read/Write GMT grd file
  public :: ncwrap__read_grd
  public :: ncwrap__write_grd

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !!
  !! interfaces
  !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Add an attribute to netcdf file
  !!
  !! @par arguments 
  !! - ncid: netcdf file ID
  !! - name: attribute name
  !! - var : attribute value (single, double, int, char)
  !<
  !! --
  interface ncwrap__add_attribute

     module procedure ncwrap__add_attribute_s, ncwrap__add_attribute_d, ncwrap__add_attribute_i, ncwrap__add_attribute_a

  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  interface ncwrap__open
     
     module procedure ncwrap__open_2d_s,  ncwrap__open_2d_d, &
                      ncwrap__open_2dt_s, ncwrap__open_2dt_d

  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  interface ncwrap__read_data
     module procedure ncwrap__read_data_2d_s, ncwrap__read_data_2d_d, &
                      ncwrap__read_data_2dt_s, ncwrap__read_data_2dt_d
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  interface ncwrap__create
     module procedure ncwrap__create_2d_s, ncwrap__create_2d_d, &
                      ncwrap__create_2dt_s, ncwrap__create_2dt_d
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  interface ncwrap__add_data
     module procedure ncwrap__add_data_2d_s, ncwrap__add_data_2d_d, &
                      ncwrap__add_data_2dt_s, ncwrap__add_data_2dt_d
  end interface
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  interface ncwrap__write_grd
     module procedure ncwrap__write_grd_s, ncwrap__write_grd_d
  end interface
  interface ncwrap__read_grd
     module procedure ncwrap__read_grd_s, ncwrap__read_grd_d
  end interface
  

  !! internal memory for mamipulating several netcdf files at once
  
  integer, parameter :: MAX_FILE = 100   !< maximum number of files operated
  integer            :: ncidtbl(MAX_FILE)
  integer            :: dimid_x(MAX_FILE)
  integer            :: dimid_y(MAX_FILE)
  integer            :: dimid_t(MAX_FILE)
  integer            :: varid_x(MAX_FILE)
  integer            :: varid_y(MAX_FILE)
  integer            :: varid_t(MAX_FILE)
  integer            :: vartype(MAX_FILE) !< real (4) or double (8)
  integer            :: ndim(MAX_FILE)    !< xy (2) or xyt(3)
  integer            :: nx0(MAX_FILE)
  integer            :: ny0(MAX_FILE)
  integer            :: nt0(MAX_FILE)
  logical            :: is_open(MAX_FILE) = .false.


contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Open time-dependent 2d(+time) netcdf file for read. It also read axis arrays. Double Precision
  !<
  !! --
  subroutine ncwrap__open_2dt_d( fname, nx, ny, nt, x, y, t, ncid )

    character(*),          intent(in)  :: fname
    integer,               intent(out) :: nx
    integer,               intent(out) :: ny
    integer,               intent(out) :: nt
    real(DP), allocatable, intent(out) :: x(:)
    real(DP), allocatable, intent(out) :: y(:)
    real(DP), allocatable, intent(out) :: t(:)
    integer,               intent(out) :: ncid
    !! --

    integer :: id
    character(256) :: xname, yname, tname
    !! ----
    
#ifdef _NETCDF

    !! obtain internal file id
    call getid( id )
    call nc_chk(  nf90_open( fname, NF90_NOWRITE, ncid ) )
    is_open(id) = .true.
    ncidtbl(id) = ncid

    !! size & name
    call nc_chk( nf90_inquire_dimension( ncid, 1, xname, nx0(id) ) )
    call nc_chk( nf90_inquire_dimension( ncid, 2, yname, ny0(id) ) )
    call nc_chk( nf90_inquire_dimension( ncid, 3, tname, nt0(id) ) )
    nx = nx0(id)
    ny = ny0(id)
    nt = nt0(id)

    !! dimension ids (assume fixed value)
    dimid_x(id) = 1
    dimid_y(id) = 1
    dimid_t(id) = 1

    !! variable ids
    call nc_chk(  nf90_inq_varid        ( ncid, xname, varid_x(id) ) )
    call nc_chk(  nf90_inq_varid        ( ncid, yname, varid_y(id) ) )
    call nc_chk(  nf90_inq_varid        ( ncid, tname, varid_t(id) ) )

    if( .not. allocated( x ) ) allocate( x(nx) )
    if( .not. allocated( y ) ) allocate( y(ny) )
    if( .not. allocated( t ) ) allocate( t(nt) )

    !! read axis data
    call nc_chk(  nf90_get_var          ( ncid, varid_x(id), x(:) )    )
    call nc_chk(  nf90_get_var          ( ncid, varid_y(id), y(:) )    )
    call nc_chk(  nf90_get_var          ( ncid, varid_t(id), t(:) )    )

#endif

  end subroutine ncwrap__open_2dt_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Open time-dependent 2d(+time) netcdf file for read. It also read axis arrays. Single Precision
  !<
  !! --
  subroutine ncwrap__open_2dt_s( fname, nx, ny, nt, x, y, t, ncid )

    character(*),          intent(in)  :: fname
    integer,               intent(out) :: nx
    integer,               intent(out) :: ny
    integer,               intent(out) :: nt
    real(SP), allocatable, intent(out) :: x(:)
    real(SP), allocatable, intent(out) :: y(:)
    real(SP), allocatable, intent(out) :: t(:)
    integer,               intent(out) :: ncid
    !! --

    integer :: id
    character(256) :: xname, yname, tname
    !! ----
    
#ifdef _NETCDF

    !! obtain internal file id
    call getid( id )
    call nc_chk(  nf90_open( fname, NF90_NOWRITE, ncid ) )
    is_open(id) = .true.
    ncidtbl(id) = ncid

    !! size & name
    call nc_chk( nf90_inquire_dimension( ncid, 1, xname, nx0(id) ) )
    call nc_chk( nf90_inquire_dimension( ncid, 2, yname, ny0(id) ) )
    call nc_chk( nf90_inquire_dimension( ncid, 3, tname, nt0(id) ) )
    nx = nx0(id)
    ny = ny0(id)
    nt = nt0(id)

    !! dimension ids (assume fixed value)
    dimid_x(id) = 1
    dimid_y(id) = 1
    dimid_t(id) = 1

    !! variable ids
    call nc_chk(  nf90_inq_varid        ( ncid, xname, varid_x(id) ) )
    call nc_chk(  nf90_inq_varid        ( ncid, yname, varid_y(id) ) )
    call nc_chk(  nf90_inq_varid        ( ncid, tname, varid_t(id) ) )

    if( .not. allocated( x ) ) allocate( x(nx) )
    if( .not. allocated( y ) ) allocate( y(ny) )
    if( .not. allocated( t ) ) allocate( t(nt) )

    !! read axis data
    call nc_chk(  nf90_get_var          ( ncid, varid_x(id), x(:) )    )
    call nc_chk(  nf90_get_var          ( ncid, varid_y(id), y(:) )    )
    call nc_chk(  nf90_get_var          ( ncid, varid_t(id), t(:) )    )

#endif

  end subroutine ncwrap__open_2dt_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Open time-dependent 2d netcdf file for read. It also read axis arrays. Double Precision
  !<
  !! --
  subroutine ncwrap__open_2d_d( fname, nx, ny, x, y, ncid )

    character(*),          intent(in)  :: fname
    integer,               intent(out) :: nx
    integer,               intent(out) :: ny
    real(DP), allocatable, intent(out) :: x(:)
    real(DP), allocatable, intent(out) :: y(:)
    integer,               intent(out) :: ncid
    !! --

    integer :: id
    character(256) :: xname, yname
    !! ----
    
#ifdef _NETCDF

    !! obtain internal file id
    call getid( id )
    call nc_chk(  nf90_open( fname, NF90_NOWRITE, ncid ) )
    is_open(id) = .true.
    ncidtbl(id) = ncid

    !! size & name
    call nc_chk( nf90_inquire_dimension( ncid, 1, xname, nx0(id) ) )
    call nc_chk( nf90_inquire_dimension( ncid, 2, yname, ny0(id) ) )
    nx = nx0(id)
    ny = ny0(id)

    !! dimension ids (assume fixed value)
    dimid_x(id) = 1
    dimid_y(id) = 1

    !! variable ids
    call nc_chk(  nf90_inq_varid        ( ncid, xname, varid_x(id) ) )
    call nc_chk(  nf90_inq_varid        ( ncid, yname, varid_y(id) ) )

    if( .not. allocated( x ) ) allocate( x(nx) )
    if( .not. allocated( y ) ) allocate( y(ny) )

    !! read axis data
    call nc_chk(  nf90_get_var          ( ncid, varid_x(id), x(:) )    )
    call nc_chk(  nf90_get_var          ( ncid, varid_y(id), y(:) )    )

#endif

  end subroutine ncwrap__open_2d_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Open time-dependent 2d netcdf file for read. It also read axis arrays. Single Precision
  !<
  !! --
  subroutine ncwrap__open_2d_s( fname, nx, ny, x, y, ncid )

    character(*),          intent(in)  :: fname
    integer,               intent(out) :: nx
    integer,               intent(out) :: ny
    real(SP), allocatable, intent(out) :: x(:)
    real(SP), allocatable, intent(out) :: y(:)
    integer,               intent(out) :: ncid
    !! --

    integer :: id
    character(256) :: xname, yname
    !! ----
    
#ifdef _NETCDF

    !! obtain internal file id
    call getid( id )
    call nc_chk(  nf90_open( fname, NF90_NOWRITE, ncid ) )
    is_open(id) = .true.
    ncidtbl(id) = ncid

    !! size & name
    call nc_chk( nf90_inquire_dimension( ncid, 1, xname, nx0(id) ) )
    call nc_chk( nf90_inquire_dimension( ncid, 2, yname, ny0(id) ) )
    nx = nx0(id)
    ny = ny0(id)

    !! dimension ids (assume fixed value)
    dimid_x(id) = 1
    dimid_y(id) = 1

    !! variable ids
    call nc_chk(  nf90_inq_varid        ( ncid, xname, varid_x(id) ) )
    call nc_chk(  nf90_inq_varid        ( ncid, yname, varid_y(id) ) )

    if( .not. allocated( x ) ) allocate( x(nx) )
    if( .not. allocated( y ) ) allocate( y(ny) )

    !! read axis data
    call nc_chk(  nf90_get_var          ( ncid, varid_x(id), x(:) )    )
    call nc_chk(  nf90_get_var          ( ncid, varid_y(id), y(:) )    )

#endif

  end subroutine ncwrap__open_2d_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine ncwrap__read_data_2d_d( ncid, varnm, z )
    
    integer,      intent(in)  :: ncid
    character(*), intent(in)  :: varnm  !< variable name
    real(DP),     intent(out) :: z(:,:) 
    !! 
    integer :: id
    logical :: found
    integer :: zid
    !! ----
#ifdef _NETCDF
    
    call idsearch( ncid, id, found )
    if( .not. found ) then
       write(STDERR,*) "ncwrap__read_data: no such ncid"
       return
    end if

    call nc_chk( nf90_inq_varid( ncid, varnm, zid ) )
    call nc_chk( nf90_get_var( ncid, zid, z ) )

#endif

  end subroutine ncwrap__read_data_2d_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine ncwrap__read_data_2d_s( ncid, varnm, z )
    
    integer,      intent(in)  :: ncid
    character(*), intent(in)  :: varnm  !< variable name
    real(SP),     intent(out) :: z(:,:) 
    !! 
    integer :: id
    integer :: zid
    logical :: found
    !! ----

#ifdef _NETCDF
    
    call idsearch( ncid, id, found )
    if( .not. found ) then
       write(STDERR,*) "ncwrap__read_data: no such ncid"
       return
    end if

    call nc_chk( nf90_inq_varid( ncid, varnm, zid ) )
    call nc_chk( nf90_get_var( ncid, zid, z ) )

#endif

  end subroutine ncwrap__read_data_2d_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine ncwrap__read_data_2dt_d( ncid, varnm, it, z )
    
    integer,      intent(in)  :: ncid
    character(*), intent(in)  :: varnm !< variable name
    integer,      intent(in)  :: it    !< time snapshot number
    real(DP),     intent(out) :: z(:,:)

    integer :: count(3)
    integer :: start(3)
    logical :: found 
    integer :: id
    integer :: zid
    !! ----

#ifdef _NETCDF

    call idsearch( ncid, id, found )
    if( .not. found ) then
       write(STDERR,*) "ncwrap__read_data: no such ncid"
       return
    end if

    count = (/ nx0(id), ny0(id), 1/)
    start = (/ 1, 1, it /)

    call nc_chk( nf90_inq_varid( ncid, varnm, zid ) )
    call nc_chk( nf90_get_var( ncid, zid, z, start=start, count=count ))
    
#endif

  end subroutine ncwrap__read_data_2dt_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
 
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine ncwrap__read_data_2dt_s( ncid, varnm, it, z )
    
    integer,      intent(in)  :: ncid
    character(*), intent(in)  :: varnm !< variable name
    integer,      intent(in)  :: it    !< time snapshot number
    real(SP),     intent(out) :: z(:,:)

    integer :: count(3)
    integer :: start(3)
    logical :: found 
    integer :: id
    integer :: zid
    !! ----

#ifdef _NETCDF

    call idsearch( ncid, id, found )
    if( .not. found ) then
       write(STDERR,*) "ncwrap__read_data: no such ncid"
       return
    end if

    count = (/ nx0(id), ny0(id), 1/)
    start = (/ 1, 1, it /)

    call nc_chk( nf90_inq_varid( ncid, varnm, zid ) )
    call nc_chk( nf90_get_var( ncid, zid, z, start=start, count=count ))
    
#endif

  end subroutine ncwrap__read_data_2dt_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Inquire grid size of 2D array NETCDF file
  !<
  !! --
  subroutine ncwrap__inq( fn_in, title, nx, ny, nt, ndim, nvar, xname, yname, tname, zname, xunit, yunit, tunit, zunit )

    !! -- Argumetns
    character(*), intent(in)  :: fn_in
    character(*), intent(out) :: title
    integer,      intent(out) :: nx
    integer,      intent(out) :: ny
    integer,      intent(out) :: nt       !< number of time grid. zero for 2D data
    integer,      intent(out) :: nvar     !< number of dependent variables
    integer,      intent(out) :: ndim     !< dimensions. 2 for 2D data, 3 for time-dependent 2D data
    character(*), intent(out) :: xname    !<  name of x axis
    character(*), intent(out) :: yname    !<  name of y axis
    character(*), intent(out) :: tname    !<  name of time axis. void for 2D data
    character(*), intent(out) :: zname(:)
    character(*), intent(out) :: xunit
    character(*), intent(out) :: yunit
    character(*), intent(out) :: tunit
    character(*), intent(out) :: zunit(:)
    !! --
    integer :: ncid
    integer :: i
    !! ----

#ifdef _NETCDF
    
    call nc_chk   (  nf90_open( trim(fn_in), NF90_NOWRITE, ncid )     )

    call nc_chk   (  nf90_inquire( ncid, ndim, nvar )                 )
    
    !! assume first ndim variables are x, y, (t) axes in variable field. The rest is dependent variables
    nvar = nvar - ndim

    call nc_chk   (  nf90_inquire_dimension( ncid, 1, xname, nx )     )
    call nc_chk   (  nf90_inquire_dimension( ncid, 2, yname, ny )     )

    call nc_chk   (  nf90_get_att( ncid, NF90_GLOBAL, 'title', title ))
    call nc_chk   (  nf90_get_att( ncid, 1, 'units', xunit )          )
    call nc_chk   (  nf90_get_att( ncid, 2, 'units', yunit )          )
    
    if( ndim == 3 ) then                                              
       call nc_chk(  nf90_inquire_dimension( ncid, 3, tname, nt )     )
       call nc_chk(  nf90_get_att( ncid, 3, 'units', tunit ) )
    else
       nt = 0
       tname = ''
    end if
    
    !! unit
    do i=1, nvar
       call nc_chk( nf90_inquire_variable ( ncid, ndim+i, zname(i) ) )
       call nc_chk( nf90_get_att          ( ncid, ndim+i, 'units', zunit(i) ) )
    end do
    call nc_chk   (  nf90_close            ( ncid )                   )
    
#endif

  end subroutine ncwrap__inq
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Create new netcdf file for time-dependent 2D array
  !!
  !! This routine create the netcdf-fomatted file fn_out, and defines axis information (name&unit)
  !! It returns netcdf file ID number ncid.
  !!
  !! @par Units
  !! default units for latitude, longitude are degree-north/degree-east for better compatibility with GMT commands. 
  !<
  !! ----
  subroutine ncwrap__create_2dt_d( fname, title, nx, ny, x, y, xname, yname, tname, xunit, yunit, tunit, ncid )

    character(*), intent(in)  :: fname  !< Filename
    character(*), intent(in)  :: title  !< Title of the file
    integer,      intent(in)  :: nx     !< x-axis size
    integer,      intent(in)  :: ny     !< y-axis size
    real(DP),     intent(in)  :: x(nx)  !< dependent variable
    real(DP),     intent(in)  :: y(ny)  !< dependent variable
    character(*), intent(in)  :: xname, yname, tname !< axis long names
    character(*), intent(in)  :: xunit, yunit, tunit !< axis units
    integer,      intent(out) :: ncid
    !!
    integer :: id
    !! ----

#ifdef _NETCDF

    !! internal file id
    call getid( id )

    !! create file
    !! CLOBBER: create new or replace old file
    call nc_chk( nf90_create( trim(fname), NF90_CLOBBER, ncidtbl(id) ))
    ncid = ncidtbl(id)
    nx0(id) = nx
    ny0(id) = ny
    is_open(id) = .true.
    vartype(id) = DP
    ndim(id)    = 3
    
    !!
    !! dimensions
    !!
    call nc_chk( nf90_def_dim( ncid, trim(xname), nx,             dimid_x(id) ) )
    call nc_chk( nf90_def_dim( ncid, trim(yname), ny,             dimid_y(id) ) )
    call nc_chk( nf90_def_dim( ncid, trim(tname), NF90_UNLIMITED, dimid_t(id) ) )

    !!
    !! variables
    !!
    call nc_chk( nf90_def_var( ncid, trim(xname), NF90_DOUBLE, dimid_x(id), varid_x(id) ) )
    call nc_chk( nf90_def_var( ncid, trim(yname), NF90_DOUBLE, dimid_y(id), varid_y(id) ) )
    call nc_chk( nf90_def_var( ncid, trim(tname), NF90_DOUBLE, dimid_t(id), varid_t(id) ) )

    !!
    !! default global attributes
    !!
    call nc_chk( nf90_put_att( ncid, NF90_GLOBAL, 'title', trim(title) ) )

    !!
    !! default variable attributes
    !!
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'long_name',    trim(xname) ) )
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'units',        trim(xunit) ) )
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'actual_range', (/ minval(x), maxval(x) /) ) )

    call nc_chk( nf90_put_att( ncid, varid_y(id), 'long_name',    trim(yname) ) )
    call nc_chk( nf90_put_att( ncid, varid_y(id), 'units',        trim(yunit) ) )
    call nc_chk( nf90_put_att( ncid, varid_y(id), 'actual_range', (/ minval(y), maxval(y) /) ) )

    call nc_chk( nf90_put_att( ncid, varid_t(id), 'long_name',    trim(tname) ) )
    call nc_chk( nf90_put_att( ncid, varid_t(id), 'units',        trim(tunit) ) )

    call nc_chk( nf90_enddef( ncid ) )

    !!
    !! add independent variables
    !!
    call nc_chk( nf90_put_var( ncid, varid_x(id), x ) )
    call nc_chk( nf90_put_var( ncid, varid_y(id), y ) )

#endif

  end subroutine ncwrap__create_2dt_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Create new netcdf file for time-dependent 2D array
  !!
  !! This routine create the netcdf-fomatted file fn_out, and defines axis information (name&unit)
  !! It returns netcdf file ID number ncid.
  !!
  !! @par Units
  !! default units for latitude, longitude are degree-north/degree-east for better compatibility with GMT commands. 
  !<
  !! ----
  subroutine ncwrap__create_2dt_s( fname, title, nx, ny, x, y, xname, yname, tname, xunit, yunit, tunit, ncid )

    character(*), intent(in)  :: fname  !< Filename
    character(*), intent(in)  :: title  !< Title of the file
    integer,      intent(in)  :: nx     !< x-axis size
    integer,      intent(in)  :: ny     !< y-axis size
    real(SP),     intent(in)  :: x(nx)  !< dependent variable
    real(SP),     intent(in)  :: y(ny)  !< dependent variable
    character(*), intent(in)  :: xname, yname, tname !< axis long names
    character(*), intent(in)  :: xunit, yunit, tunit !< axis units
    integer,      intent(out) :: ncid
    !!
    integer :: id
    !! ----

#ifdef _NETCDF

    !! internal file id
    call getid( id )

    !! create file
    !! CLOBBER: create new or replace old file
    call nc_chk( nf90_create( trim(fname), NF90_CLOBBER, ncidtbl(id) ))
    ncid = ncidtbl(id)
    nx0(id) = nx
    ny0(id) = ny
    is_open(id) = .true.
    vartype(id) = SP
    ndim(id)    = 3
    
    !!
    !! dimensions
    !!
    call nc_chk( nf90_def_dim( ncid, trim(xname), nx,             dimid_x(id) ) )
    call nc_chk( nf90_def_dim( ncid, trim(yname), ny,             dimid_y(id) ) )
    call nc_chk( nf90_def_dim( ncid, trim(tname), NF90_UNLIMITED, dimid_t(id) ) )

    !!
    !! variables
    !!
    call nc_chk( nf90_def_var( ncid, trim(xname), NF90_REAL, dimid_x(id), varid_x(id) ) )
    call nc_chk( nf90_def_var( ncid, trim(yname), NF90_REAL, dimid_y(id), varid_y(id) ) )
    call nc_chk( nf90_def_var( ncid, trim(tname), NF90_REAL, dimid_t(id), varid_t(id) ) )

    !!
    !! default global attributes
    !!
    call nc_chk( nf90_put_att( ncid, NF90_GLOBAL, 'title', trim(title) ) )

    !!
    !! default variable attributes
    !!
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'long_name',    trim(xname) ) )
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'units',        trim(xunit) ) )
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'actual_range', (/ minval(x), maxval(x) /) ) )

    call nc_chk( nf90_put_att( ncid, varid_y(id), 'long_name',    trim(yname) ) )
    call nc_chk( nf90_put_att( ncid, varid_y(id), 'units',        trim(yunit) ) )
    call nc_chk( nf90_put_att( ncid, varid_y(id), 'actual_range', (/ minval(y), maxval(y) /) ) )

    call nc_chk( nf90_put_att( ncid, varid_t(id), 'long_name',    trim(tname) ) )
    call nc_chk( nf90_put_att( ncid, varid_t(id), 'units',        trim(tunit) ) )

    call nc_chk( nf90_enddef( ncid ) )

    !!
    !! add independent variables
    !!
    call nc_chk( nf90_put_var( ncid, varid_x(id), x ) )
    call nc_chk( nf90_put_var( ncid, varid_y(id), y ) )

#endif

  end subroutine ncwrap__create_2dt_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Create new netcdf file for 2D array (GMT grd file)
  !!
  !! Assume unlimited number of 2D array whose spatial dimensions are (nx,ny)
  !! This routine create the netcdf-fomatted file fn_out, and defines axis information (name&unit)
  !! It returns netcdf file ID number ncid.
  !!
  !! @par Units
  !! default units for latitude, longitude are degree-north/degree-east for better compatibility with GMT. 
  !<
  !! --
  subroutine ncwrap__create_2d_d( fname, title, nx, ny, x, y, xname, yname, xunit, yunit, ncid )
    
    character(*), intent(in)  :: fname        !< filename
    character(*), intent(in)  :: title        !< title of the file
    integer,      intent(in)  :: nx           !< x-axis size
    integer,      intent(in)  :: ny           !< y-axis size
    real(DP),     intent(in)  :: x(nx)        !< dependent variable
    real(DP),     intent(in)  :: y(ny)        !< dependent variable
    character(*), intent(in)  :: xname, yname !< axis long names
    character(*), intent(in)  :: xunit, yunit !< axis units
    integer,      intent(out) :: ncid         !< netcdf file id
    integer :: id

#ifdef _NETCDF
    
    !! internal file id
    call getid( id )

    !! create file
    !! CLOBBER: create new or replace old file
    call nc_chk( nf90_create( trim(fname), NF90_CLOBBER, ncidtbl(id) ))

    !! remember ncid and size
    ncid     = ncidtbl(id)
    nx0(id) = nx
    ny0(id) = ny
    is_open(id) = .true.
    vartype(id) = DP
    ndim(id)    = 2

    !!
    !! dimensions
    !!
    call nc_chk( nf90_def_dim( ncid, trim(xname), nx,             dimid_x(id) ) )
    call nc_chk( nf90_def_dim( ncid, trim(yname), ny,             dimid_y(id) ) )

    !!
    !! variables
    !!
    call nc_chk( nf90_def_var( ncid, trim(xname), NF90_DOUBLE, dimid_x(id), varid_x(id) ) )
    call nc_chk( nf90_def_var( ncid, trim(yname), NF90_DOUBLE, dimid_y(id), varid_y(id) ) )

    !!
    !! default global attributes
    !!
    call nc_chk( nf90_put_att( ncid, NF90_GLOBAL, 'title', trim(title) ) )

    !!
    !! default variable attributes
    !!
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'long_name',    trim(xname) ) )
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'units',        trim(xunit) ) )
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'actual_range', (/ minval(x), maxval(x) /) ) )

    call nc_chk( nf90_put_att( ncid, varid_y(id), 'long_name',    trim(yname) ) )
    call nc_chk( nf90_put_att( ncid, varid_y(id), 'units',        trim(yunit) ) )
    call nc_chk( nf90_put_att( ncid, varid_y(id), 'actual_range', (/ minval(y), maxval(y) /) ) )

    call nc_chk( nf90_enddef( ncid ) )

    !!
    !! add independent variables
    !!
    call nc_chk( nf90_put_var( ncid, varid_x(id), x ) )
    call nc_chk( nf90_put_var( ncid, varid_y(id), y ) )

#endif
    
  end subroutine ncwrap__create_2d_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Create new netcdf file for 2D array (GMT grd file)
  !!
  !! Assume unlimited number of 2D array whose spatial dimensions are (nx,ny)
  !! This routine create the netcdf-fomatted file fn_out, and defines axis information (name&unit)
  !! It returns netcdf file ID number ncid.
  !!
  !! @par Units
  !! default units for latitude, longitude are degree-north/degree-east for better compatibility with GMT. 
  !<
  !! --
  subroutine ncwrap__create_2d_s ( fname, title, nx, ny, x, y, xname, yname, xunit, yunit, ncid )
    
    character(*), intent(in)  :: fname        !< filename
    character(*), intent(in)  :: title        !< title of the file
    integer,      intent(in)  :: nx           !< x-axis size
    integer,      intent(in)  :: ny           !< y-axis size
    real(SP),     intent(in)  :: x(nx)        !< dependent variable
    real(SP),     intent(in)  :: y(ny)        !< dependent variable
    character(*), intent(in)  :: xname, yname !< axis long names
    character(*), intent(in)  :: xunit, yunit !< axis units
    integer,      intent(out) :: ncid         !< netcdf file id
    integer :: id

#ifdef _NETCDF
    
    !! internal file id
    call getid( id )

    !! create file
    !! CLOBBER: create new or replace old file
    call nc_chk( nf90_create( trim(fname), NF90_CLOBBER, ncidtbl(id) ))

    !! remember ncid and size
    ncid     = ncidtbl(id)
    nx0(id) = nx
    ny0(id) = ny
    is_open(id) = .true.
    vartype(id) = SP
    ndim(id)    = 2

    !!
    !! dimensions
    !!
    call nc_chk( nf90_def_dim( ncid, trim(xname), nx,             dimid_x(id) ) )
    call nc_chk( nf90_def_dim( ncid, trim(yname), ny,             dimid_y(id) ) )

    !!
    !! variables
    !!
    call nc_chk( nf90_def_var( ncid, trim(xname), NF90_REAL, dimid_x(id), varid_x(id) ) )
    call nc_chk( nf90_def_var( ncid, trim(yname), NF90_REAL, dimid_y(id), varid_y(id) ) )

    !!
    !! default global attributes
    !!
    call nc_chk( nf90_put_att( ncid, NF90_GLOBAL, 'title', trim(title) ) )

    !!
    !! default variable attributes
    !!
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'long_name',    trim(xname) ) )
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'units',        trim(xunit) ) )
    call nc_chk( nf90_put_att( ncid, varid_x(id), 'actual_range', (/ minval(x), maxval(x) /) ) )

    call nc_chk( nf90_put_att( ncid, varid_y(id), 'long_name',    trim(yname) ) )
    call nc_chk( nf90_put_att( ncid, varid_y(id), 'units',        trim(yunit) ) )
    call nc_chk( nf90_put_att( ncid, varid_y(id), 'actual_range', (/ minval(y), maxval(y) /) ) )

    call nc_chk( nf90_enddef( ncid ) )

    !!
    !! add independent variables
    !!
    call nc_chk( nf90_put_var( ncid, varid_x(id), x ) )
    call nc_chk( nf90_put_var( ncid, varid_y(id), y ) )

#endif
    
  end subroutine ncwrap__create_2d_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Define new data field. Returns variable id. 
  !<
  !! --
  subroutine ncwrap__def_data( ncid, zname, zunit, zid )
    
    integer,      intent(in)  :: ncid  !< netcdf file id
    character(*), intent(in)  :: zname !< variable longname
    character(*), intent(in)  :: zunit !< unit
    integer,      intent(out) :: zid   !< return variable id

    !!--
    integer :: id
    logical :: found
    integer :: dimids_3(3), dimids_2(2)

    !! ----
    
#ifdef _NETCDF

    call idsearch( ncid, id, found )
    if( .not. found ) then
       write(STDERR,*) "ncwrap__def_dat_2dt: no such ncid"
       return
    end if
    

    !!
    !! enter the define mode
    !!
    call nc_chk( nf90_redef  ( ncid ) )

    !!
    !! define new variable
    !!
    if( ndim(id) == 3 ) then
       dimids_3( 1:3 ) = (/ dimid_x(id), dimid_y(id), dimid_t(id) /)
       if( vartype(id) == SP ) then
          call nc_chk( nf90_def_var( ncid, trim(zname), NF90_REAL, dimids_3, zid ) )
       else
          call nc_chk( nf90_def_var( ncid, trim(zname), NF90_DOUBLE, dimids_3, zid ) )
       end if
    else if ( ndim(id) == 2 ) then
       dimids_2( 1:2 ) = (/ dimid_x(id), dimid_y(id) /)
       if( vartype(id) == SP ) then
          call nc_chk( nf90_def_var( ncid, trim(zname), NF90_REAL, dimids_2, zid ) )
       else
          call nc_chk( nf90_def_var( ncid, trim(zname), NF90_DOUBLE, dimids_2, zid ) )
       end if
    end if
    

    !!
    !! add default variable attributes
    !!
    call nc_chk( nf90_put_att( ncid, zid, 'long_name', trim(zname) ) )
    call nc_chk( nf90_put_att( ncid, zid, 'units',     trim(zunit) ) )

    !!
    !! close the define mode
    !!
    call nc_chk( nf90_enddef ( ncid ) )

#endif
    
  end subroutine ncwrap__def_data
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Add data snapshot at time t attributed by varid to netcdf datafile.
  !<
  !! --
  subroutine ncwrap__add_data_2dt_d( ncid, varid, it, t, z )
    
    integer,  intent(in) :: ncid
    integer,  intent(in) :: varid
    integer,  intent(in) :: it
    real(DP), intent(in) :: t
    real(DP), intent(in) :: z(:,:)
    !! --
    integer :: count(3)
    integer :: start(3)
    integer :: id
    logical :: found
    !! ----

#ifdef _NETCDF
    
    !! seek file id
    call idsearch( ncid, id, found )
    if( .not. found ) then
       write(STDERR,*) "ncwrap__add_data: no such ncid"
       return
    end if

    
    count = (/ nx0(id), ny0(id), 1/)
    start = (/ 1, 1, it /)

    call nc_chk( nf90_put_var( ncid, varid, z, start, count ) )
    call nc_chk( nf90_put_var( ncid, varid_t(id), t, start=(/ it /) ) )
    call nc_chk( nf90_sync   ( ncid ) )

#endif
    
  end subroutine ncwrap__add_data_2dt_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Add data snapshot at time t attributed by varid to netcdf datafile.
  !<
  !! --
  subroutine ncwrap__add_data_2dt_s( ncid, varid, it, t, z )
    
    integer,  intent(in) :: ncid
    integer,  intent(in) :: varid
    integer,  intent(in) :: it
    real(SP), intent(in) :: t
    real(SP), intent(in) :: z(:,:)
    !! --
    integer :: count(3)
    integer :: start(3)
    integer :: id
    logical :: found
    !! ----

#ifdef _NETCDF
    
    !! seek file id
    call idsearch( ncid, id, found )
    if( .not. found ) then
       write(STDERR,*) "ncwrap__add_data: no such ncid"
       return
    end if

    
    count = (/ nx0(id), ny0(id), 1/)
    start = (/ 1, 1, it /)

    call nc_chk( nf90_put_var( ncid, varid, z, start, count ) )
    call nc_chk( nf90_put_var( ncid, varid_t(id), t, start=(/ it /) ) )
    call nc_chk( nf90_sync   ( ncid ) )

#endif
    
  end subroutine ncwrap__add_data_2dt_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Add data attributed by varid to netcdf datafile.
  !<
  !! --
  subroutine ncwrap__add_data_2d_d( ncid, varid, z )
    
    integer,  intent(in) :: ncid
    integer,  intent(in) :: varid
    real(DP), intent(in) :: z(:,:)
    !! --
    integer :: id
    logical :: found
    !! ----

#ifdef _NETCDF
    
    !! seek file id
    call idsearch( ncid, id, found )
    if( .not. found ) then
       write(STDERR,*) "ncwrap__add_data_2d: no such ncid"
       return
    end if
    
    call nc_chk( nf90_redef  ( ncid ) )
    call nc_chk( nf90_put_att( ncid, varid, 'actual_range', (/ minval(z), maxval(z) /) ) )
    call nc_chk( nf90_enddef ( ncid ) )

    call nc_chk( nf90_put_var( ncid, varid, z ) )
    call nc_chk( nf90_sync   ( ncid ) )

#endif
    
  end subroutine ncwrap__add_data_2d_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Add data attributed by varid to netcdf datafile.
  !<
  !! --
  subroutine ncwrap__add_data_2d_s( ncid, varid, z )
    
    integer,  intent(in) :: ncid
    integer,  intent(in) :: varid
    real(SP), intent(in) :: z(:,:)
    !! --
    integer :: id
    logical :: found
    !! ----

#ifdef _NETCDF
    
    !! seek file id
    call idsearch( ncid, id, found )
    if( .not. found ) then
       write(STDERR,*) "ncwrap__add_data_2d: no such ncid"
       return
    end if
    
    call nc_chk( nf90_redef  ( ncid ) )
    call nc_chk( nf90_put_att( ncid, varid, 'actual_range', (/ minval(z), maxval(z) /) ) )
    call nc_chk( nf90_enddef ( ncid ) )

    call nc_chk( nf90_put_var( ncid, varid, z ) )
    call nc_chk( nf90_sync   ( ncid ) )

#endif
    
  end subroutine ncwrap__add_data_2d_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Add real-valued attribute to netcdf file
  !<
  !! --
  subroutine ncwrap__add_attribute_s( ncid, name, var )

    integer,      intent(in) :: ncid
    character(*), intent(in) :: name
    real(SP),     intent(in) :: var

    !! ----

#ifdef _NETCDF
    call nc_chk( nf90_redef( ncid ) )
    call nc_chk( nf90_put_att( ncid, NF90_GLOBAL, trim(name), var ) )
    call nc_chk( nf90_enddef( ncid ) )
#endif

  end subroutine ncwrap__add_attribute_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Add real-valued attribute to netcdf file
  !<
  !! --
  subroutine ncwrap__add_attribute_d( ncid, name, var )

    integer,      intent(in) :: ncid
    character(*), intent(in) :: name
    real(DP),     intent(in) :: var

    !! ----

#ifdef _NETCDF
    call nc_chk( nf90_redef( ncid ) )
    call nc_chk( nf90_put_att( ncid, NF90_GLOBAL, trim(name), var ) )
    call nc_chk( nf90_enddef( ncid ) )
#endif

  end subroutine ncwrap__add_attribute_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Add integer-valued attribute to netcdf file
  !<
  !! --
  subroutine ncwrap__add_attribute_i( ncid, name, var )

    integer,      intent(in) :: ncid
    character(*), intent(in) :: name
    integer,      intent(in) :: var

    !! ----
#ifdef _NETCDF
    call nc_chk( nf90_redef( ncid ) )
    call nc_chk( nf90_put_att( ncid, NF90_GLOBAL, trim(name), var ) )
    call nc_chk( nf90_enddef( ncid ) )
#endif
  end subroutine ncwrap__add_attribute_i
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Add character attribute to netcdf file
  !<
  !! --
  subroutine ncwrap__add_attribute_a( ncid, name, var )

    integer,      intent(in) :: ncid
    character(*), intent(in) :: name
    character(*), intent(in) :: var

    !! ---- 
#ifdef _NETCDF   
    call nc_chk( nf90_redef( ncid ) )
    call nc_chk( nf90_put_att( ncid, NF90_GLOBAL, trim(name), trim(var) ) )
    call nc_chk( nf90_enddef( ncid ) )
#endif
  end subroutine ncwrap__add_attribute_a
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Close netcdf file
  !<
  !! --
  subroutine ncwrap__close( ncid )

    integer, intent(in) :: ncid
    integer :: id
    logical :: found
    !! ----

#ifdef _NETCDF
    call idsearch( ncid, id, found )
    if( .not. found ) then
       write(STDERR,*) "ncwrap__close: no such ncid"
    end if
    
    call nc_chk( nf90_close( ncid ) )
    is_open( id ) = .false.
#endif    
  end subroutine ncwrap__close
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! A set of subroutine call for creating 2D time-independent netcdf file (GMT grd file)
  !<
  !! --
  subroutine ncwrap__write_grd_s( fname, title, nx, ny, x, y, z, xname, yname, zname, xunit, yunit, zunit )
    
    character(*), intent(in) :: fname
    character(*), intent(in) :: title  !< Title of the file
    integer,      intent(in) :: nx     !< x-axis size
    integer,      intent(in) :: ny     !< y-axis size
    real(SP),     intent(in) :: x(nx)  !< dependent variable
    real(SP),     intent(in) :: y(ny)  !< dependent variable
    real(SP),     intent(in) :: z(nx,ny)
    character(*), intent(in) :: xname, yname, zname !< axis long names
    character(*), intent(in) :: xunit, yunit, zunit !< axis units
    integer :: ncid
    integer :: zid

#ifdef _NETCDF
    
    call ncwrap__create  ( fname, title, nx, ny, x, y, xname, yname, xunit, yunit, ncid )
    call ncwrap__def_data( ncid, zname, zunit, zid )
    call ncwrap__add_data( ncid, zid, z )
    call ncwrap__close( ncid )

#endif

  end subroutine ncwrap__write_grd_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! A set of subroutine call for creating 2D time-independent netcdf file (GMT grd file)
  !<
  !! --
  subroutine ncwrap__write_grd_d( fname, title, nx, ny, x, y, z, xname, yname, zname, xunit, yunit, zunit )
    
    character(*), intent(in) :: fname
    character(*), intent(in) :: title  !< Title of the file
    integer,      intent(in) :: nx     !< x-axis size
    integer,      intent(in) :: ny     !< y-axis size
    real(DP),     intent(in) :: x(nx)  !< dependent variable
    real(DP),     intent(in) :: y(ny)  !< dependent variable
    real(DP),     intent(in) :: z(nx,ny)
    character(*), intent(in) :: xname, yname, zname !< axis long names
    character(*), intent(in) :: xunit, yunit, zunit !< axis units
    integer :: ncid
    integer :: zid

#ifdef _NETCDF
    
    call ncwrap__create  ( fname, title, nx, ny, x, y, xname, yname, xunit, yunit, ncid )
    call ncwrap__def_data( ncid, zname, zunit, zid )
    call ncwrap__add_data( ncid, zid, z )
    call ncwrap__close( ncid )

#endif
    
  end subroutine ncwrap__write_grd_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! A set of subroutine call for reading 2d time-independent netcdf file (GMT grd file)
  !<
  !! --
  subroutine ncwrap__read_grd_s( fname, title, nx, ny, x, y, z, xname, yname, zname, xunit, yunit, zunit )
    character(*),          intent(in)  :: fname
    character(*),          intent(out) :: title  !< Title of the file
    integer,               intent(out) :: nx     !< x-axis size
    integer,               intent(out) :: ny     !< y-axis size
    real(SP), allocatable, intent(out) :: x(:)  !< dependent variable
    real(SP), allocatable, intent(out) :: y(:)  !< dependent variable
    real(SP), allocatable, intent(out) :: z(:,:)
    character(*),          intent(out) :: xname, yname, zname !< axis long names
    character(*),          intent(out) :: xunit, yunit, zunit !< axis units
    !! --
    integer :: ncid
    integer :: nt
    character(80) :: tname, tunit ! dummy 
    integer :: nvar, ndim
    character(80) :: zname2(10), zunit2(10)

#ifdef _NETCDF

    call ncwrap__inq( fname, title, nx, ny, nt, ndim, nvar, xname, yname, tname, zname2, xunit, yunit, tunit, zunit2 )
    ! assume first record
    zname = zname2(1)
    zunit = zunit2(1)
    allocate( x(nx), y(ny), z(nx,ny) )

    call ncwrap__open( fname, nx, ny, x, y, ncid )
    call ncwrap__read_data( ncid, zname2(1), z )
    call ncwrap__close( ncid )
#endif

  end subroutine ncwrap__read_grd_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! A set of subroutine call for reading 2d time-independent netcdf file (GMT grd file)
  !<
  !! --
  subroutine ncwrap__read_grd_d( fname, title, nx, ny, x, y, z, xname, yname, zname, xunit, yunit, zunit )
    character(*),          intent(in)  :: fname
    character(*),          intent(out) :: title  !< Title of the file
    integer,               intent(out) :: nx     !< x-axis size
    integer,               intent(out) :: ny     !< y-axis size
    real(DP), allocatable, intent(out) :: x(:)  !< dependent variable
    real(DP), allocatable, intent(out) :: y(:)  !< dependent variable
    real(DP), allocatable, intent(out) :: z(:,:)
    character(*),          intent(out) :: xname, yname, zname !< axis long names
    character(*),          intent(out) :: xunit, yunit, zunit !< axis units
    !! --
    integer :: ncid
    integer :: nt
    character(80) :: tname, tunit ! dummy 
    integer :: nvar, ndim
    character(80) :: zname2(10), zunit2(10)

#ifdef _NETCDF
    call ncwrap__inq( fname, title, nx, ny, nt, ndim, nvar, xname, yname, tname, zname2, xunit, yunit, tunit, zunit2 )
    ! assume first record
    zname(:) = zname2(1)(:)
    zunit(:) = zunit2(1)(:)
    allocate( x(nx), y(ny), z(nx,ny) )

    call ncwrap__open( fname, nx, ny, x, y, ncid )
    call ncwrap__read_data( ncid, zname2(1), z )
    call ncwrap__close( ncid )
#endif

  end subroutine ncwrap__read_grd_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!
    
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! An internal subroutine to find the file id from given netcdf ncid. 
  !<
  !! --
  subroutine idsearch( ncid, id, found )

    integer, intent(in)  :: ncid
    integer, intent(out) :: id
    logical, intent(out) :: found

    !! --
    integer :: i

    !! ----

#ifdef _NETCDF
    found = .false.
    do i=1, MAX_FILE
       if( ncidtbl(i) == ncid ) then
          if( is_open(i) ) then
             found = .true.
             id = i
             exit
          end if
       end if
    end do
#endif
    
  end subroutine idsearch
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine getid( id )
    
    integer, intent(out) :: id
    integer :: i
#ifdef _NETCDF
    do i=1, MAX_FILE
       if( .not. is_open(i) ) then
          id = i
          exit
       end if
    end do
#endif
  end subroutine getid
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! An internal subroutine to check error in netcdf function calls
  !<
  !! --
  subroutine nc_chk( ierr )

    integer, intent(in) :: ierr
    !! ----

#ifdef _NETCDF
    if( ierr /= NF90_NOERR )  write(STDERR,*) NF90_STRERROR( ierr )
#endif

  end subroutine nc_chk
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  
end module m_ncwrap
!! ----------------------------------------------------------------------------------------------------------------------------- !!
