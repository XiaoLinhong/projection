module mod_proj
   !! projected coordinate reference system
   !!
   !! Module that defines constants, data structures, and
   !! subroutines used to convert grid indices to lat/lon
   !! and vice versa.
   !!
   !! SUPPORTED PROJECTIONS
   !! ---------------------
   !! Lambert Conformal (code = 1)
   !!
   !!  Data Definitions:
   !!       1. Any arguments that are a latitude value are expressed in
   !!          degrees north with a valid range of -90 -> 90
   !!       2. Any arguments that are a longitude value are expressed in
   !!          degrees east with a valid range of -180 -> 180.
   !!       3. Distances are in meters and are always positive.

   implicit none
   private
   public fp, proj_type, create_proj

   integer, parameter :: fp = 8

   real(fp), parameter :: PI = 3.141592653589793_fp
   real(fp), parameter :: DEG_PER_RAD = 180./PI
   real(fp), parameter :: RAD_PER_DEG = PI/180.
   real(fp), parameter :: EARTH_RADIUS_M = 6370000. !! same as MM5 system

   integer, parameter :: LCC = 1

   type globe ! 地球参数，基面
      real(fp)  :: a = EARTH_RADIUS_M ! 长轴半径
      real(fp)  :: flattening = 0. ! flattenging of the ellipsoid (a-b)/a
      real(fp)  :: eccentricity = 0. ! eccentricity of the ellipsoid sqrt(2f-f**2)
   end type globe

   type proj_type ! 投影数据结构
      integer :: code = 0 ! 投影编号
      ! degree
      real(fp)  :: lat1 = 30.    !! 中心纬度
      real(fp)  :: lon1 = 115.    !! 中心经度
      real(fp)  :: truelat1 = 20. !! First true latitude
      real(fp)  :: truelat2 = 40. !! Second true latitude

      ! Radians
      real(fp)  :: phi1 = 0. !! First true latitude
      real(fp)  :: phi2 = 0. !! Second true latitude
      real(fp)  :: phi = 0.  !! Latitude of false origin
      real(fp)  :: lambda = 0. !! Longitude of false origin
      ! 目前没用到
      real(fp)  :: xshift = 0.
      real(fp)  :: yshift = 0.
      ! datum
      type(globe) :: ellipsoid !! 基面

   contains
      procedure :: ij_to_ll ! lambert to latlon
      procedure :: ll_to_ij ! latlon to lambert
   end type proj_type

contains

   type(proj_type) function create_proj(code, lon1, lat1, truelat1, truelat2, flattening, xshift, yshift) result(p)
      !! 支持Lambert投影
      implicit none

      integer, intent(in) :: code !! 1 for lcc, 0 for latlon
      real(fp), intent(in), optional  :: lon1     !! 中心经度，默认是本初子午线(lcc), 左下角纬度（latlon）
      real(fp), intent(in), optional  :: lat1     !! 中心纬度(lcc), 左下角纬度（latlon）
      real(fp), intent(in), optional  :: truelat1 !! First true latitude
      real(fp), intent(in), optional  :: truelat2 !! Second true latitude
      real(fp), intent(in), optional  :: flattening !! eccentricity of the ellipsoid
      real(fp), intent(in), optional  :: xshift     !! Easting at false origin
      real(fp), intent(in), optional  :: yshift     !! Northing at false origin

      p%code = code
      ! lcc
      if (present(lat1)) p%lat1 = lat1
      if (present(lon1)) p%lon1 = lon1
      if (present(xshift)) p%xshift = xshift
      if (present(yshift)) p%yshift = yshift
      if (present(truelat1)) p%truelat1 = truelat1
      if (present(truelat2)) p%truelat2 = truelat2
      if (present(flattening)) p%ellipsoid%flattening = flattening

      p%ellipsoid%eccentricity = sqrt(p%ellipsoid%flattening*(2-p%ellipsoid%flattening))

      p%phi1 = p%truelat1*RAD_PER_DEG
      p%phi2 = p%truelat2*RAD_PER_DEG
      p%phi = p%lat1*RAD_PER_DEG
      p%lambda = p%lon1*RAD_PER_DEG

   end function create_proj

   real(fp) function calculate_m(eccentricity, x) result(m)
      !! IOGP Publication 373-7-2
      implicit none
      real(fp), intent(in) :: eccentricity
      real(fp), intent(in) :: x
      m = cos(x) / (1 - eccentricity**2 * sin(x)**2)**0.5
   end function calculate_m

   real(fp) function calculate_t(eccentricity, x) result(t)
      !! IOGP Publication 373-7-2
      implicit none
      real(fp), intent(in) :: eccentricity
      real(fp), intent(in) :: x
      t = tan(PI/4. - x/2.) / ( (1-eccentricity*sin(x))/(1+eccentricity*sin(x)))**(eccentricity/2)
   end function calculate_t

   subroutine ij_to_ll(this, i, j, lon, lat)
      ! Subroutine to compute the geographical latitude and longitude values
      ! to the cartesian x/y on a Lambert Conformal projection.
      ! IOGP Publication 373-7-2
      implicit none
      ! Input Args
      class(proj_type), intent(in)  :: this !! Projection info structure
      real(fp), intent(in) :: i        !! Cartesian X coordinate, unit: m
      real(fp), intent(in) :: j        !! Cartesian Y coordinate, unit: m
      ! Output Args
      real(fp), intent(out)  :: lon      !! Longitude (-180->180 E)
      real(fp), intent(out)  :: lat      !! Latitude (-90->90 N)
      ! Output Args

      ! Locals
      integer :: k
      real(fp)   :: phi ! Latitude
      real(fp)   :: lamda ! Longitude
      real(fp)   :: m1, m2
      real(fp)   :: t1, t2, t3
      real(fp)   :: n, F, r
      real(fp)   :: r_dash, t_dash, theta_dash, eccentricity

      eccentricity = this%ellipsoid%eccentricity
      m1 = calculate_m(eccentricity, this%phi1)
      m2 = calculate_m(eccentricity, this%phi2)

      t1 = calculate_t(eccentricity, this%phi1)
      t2 = calculate_t(eccentricity, this%phi2)
      t3 = calculate_t(eccentricity, this%phi)

      n = (log(m1)-log(m2))/(log(t1)-log(t2))
      F = m1/(n * sign(1., t1)*abs(t1)**n)
      r = this%ellipsoid%a * F * (sign(1., t3)*abs(t3)**n)

      ! Iterative solution for the latitude
      r_dash = sign(1., n) * sqrt((i-this%xshift)**2 + (r - (j-this%yshift))**2)
      t_dash = sign(1., r_dash/(this%ellipsoid%a*F)) * abs(r_dash/(this%ellipsoid%a*F))**(1/n)
      phi = PI/2 - 2*ATAN(t_dash)  ! Initial value
      do k = 1, 20
         phi = PI/2 - 2*ATAN(t_dash * ((1-eccentricity*sin(phi))/(1 + eccentricity*sin(phi)))**(eccentricity/2))
      end do

      ! Longitude
      theta_dash = ATAN2(i-this%xshift, r-(j-this%yshift))
      lamda = theta_dash / n + this%lambda

      ! rad to deg
      lat = phi*DEG_PER_RAD
      lon = lamda*DEG_PER_RAD
      if (lon > +180.) lon = lon - 360.
      if (lon < -180.) lon = lon + 360.

   end subroutine ij_to_ll

   subroutine ll_to_ij(this, lon, lat, i, j)
      !! Subroutine to compute the geographical latitude and longitude values
      !! to the cartesian x/y on a Lambert Conformal projection.
      !! IOGP Publication 373-7-2
      implicit none
      ! Input Args
      class(proj_type), intent(in) :: this !! Projection info structure
      real(fp), intent(in)  :: lon      !! Longitude (-180->180 E)
      real(fp), intent(in)  :: lat      !! Latitude (-90->90 N)
      ! Output Args
      real(fp), intent(out) :: i        !! Cartesian X coordinate, unit: m
      real(fp), intent(out) :: j        !! Cartesian Y coordinate, unit: m

      ! Locals
      real(fp)   :: phi ! Latitude
      real(fp)   :: lamda ! Longitude
      real(fp)   :: m1, m2
      real(fp)   :: t, t1, t2, t_f
      real(fp)   :: n, F, r, r_f, lon1
      real(fp)   :: theta, eccentricity

      ! 处理维度, 0 - 360;
      lon1 = lon
      if (lon1 < 0.) lon1 = lon1 + 360.

      phi  = lat  * RAD_PER_DEG
      lamda = lon1 * RAD_PER_DEG

      eccentricity = this%ellipsoid%eccentricity
      m1 = calculate_m(eccentricity, this%phi1)
      m2 = calculate_m(eccentricity, this%phi2)

      t = calculate_t(eccentricity, phi)
      t1 = calculate_t(eccentricity, this%phi1)
      t2 = calculate_t(eccentricity, this%phi2)
      t_f = calculate_t(eccentricity, this%phi)

      n = (log(m1)-log(m2))/(log(t1)-log(t2))
      F = m1/(n*t1**n)
      r = this%ellipsoid%a * F * (sign(1., t)*abs(t)**n)
      r_f = this%ellipsoid%a * F * (sign(1., t_f)*abs(t_f)**n)

      theta = n*(lamda - this%lambda)

      i = this%xshift + r*sin(theta)  ! Easting (x-axis)
      j = this%yshift + r_f - r*cos(theta)  ! Northing (y-axis)

   end subroutine ll_to_ij

end module mod_proj
