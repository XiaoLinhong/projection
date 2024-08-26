program main
    !! This is our program
    use projection, only: fp, proj_type, create_proj
    implicit none

    real(fp)  :: lat1 = 30.    !! 中心纬度
    real(fp)  :: lon1 = 115.    !! 中心经度
    real(fp)  :: truelat1 = 20. !! First true latitude
    real(fp)  :: truelat2 = 40. !! Second true latitude

    real(fp) :: lat, lon
    real(fp) :: x, y

    type(proj_type) :: p ! 投影信息
    p = create_proj(1, lon1, lat1, truelat1, truelat2)

    ! Convert lat/lon to x/y
    lat = 35.0
    lon = 120.0
    call p%ll_to_ij(lon, lat, x, y)
    print *, "Cartesian coordinates: ", x, y

    ! Convert x/y back to lat/lon
    call p%ij_to_ll(x, y, lon, lat)
    print *, "Geographical coordinates: ", lon, lat

end program main
