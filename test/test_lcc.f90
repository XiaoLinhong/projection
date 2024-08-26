module test_ppm
  use projection, only: fp, proj_type, create_proj
  use testdrive, only : error_type, unittest_type, new_unittest, check

  implicit none
  private

  public :: collect_ppm

contains

  !> Collect all exported unit tests
  subroutine collect_ppm(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [new_unittest("lcc", test_lcc)]
  end subroutine collect_ppm

  !> Check substitution of a single line
  subroutine test_lcc(error)
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
  
    real(fp)  :: lat1 = 30.    !! 中心纬度
    real(fp)  :: lon1 = 115.    !! 中心经度
    real(fp)  :: truelat1 = 20. !! First true latitude
    real(fp)  :: truelat2 = 40. !! Second true latitude

    real(fp)  :: lat, lon
    real(fp)  :: x, y

    type(proj_type) :: p ! 投影信息
    p = create_proj(1, lon1, lat1, truelat1, truelat2)

    ! import numpy as np
    ! import cartopy.crs as ccrs
    ! RADIUS = 6370000
    ! globe = ccrs.Globe(ellipse='WGS84', semimajor_axis=RADIUS, semiminor_axis=RADIUS)
    ! param=dict(
    !     central_latitude=30.5, # 同namelist.wps中的ref_lat
    !     central_longitude=115.0, # 同namelist.wps中的ref_lon
    !     standard_parallels=(20.0, 40.0), # 同namelist.wps中的truelat1，truelat2
    ! )
    ! lambert = ccrs.LambertConformal(globe=globe, **param)
    ! latlon = ccrs.PlateCarree()
    ! lambert.transform_points(latlon, np.array([120]), np.array([35])).astype('float32')

    ! Convert lat/lon to x/y
    lat = 35.0
    lon = 120.0
    call p%ll_to_ij(lon, lat, x, y)
    print *, "Cartesian coordinates: ", x, y

    call check(error, real(abs(x-449911.12)) < 1.0, "x should be near 449911.12")
    call check(error, real(abs(y-557935.4)) < 1.0, "y should be near 557935.4")

    ! Convert x/y back to lat/lon
    call p%ij_to_ll(x, y, lon, lat)
    print *, "Geographical coordinates: ", lon, lat

    call check(error, real(lat) == 35.0, "lat must be 35")
    call check(error, real(lon) == 120.0, "lon must be 35")
  
  end subroutine test_lcc
end module test_ppm

program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite
  use test_ppm, only : collect_ppm
  implicit none
  integer :: stat

  stat = 0
  call run_testsuite(collect_ppm, error_unit, stat)

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

end program tester
