!*******************************************************************************

!  Module:      HEALPixModule

!  Description: HEALPix routines for generating nested pixels on the 2-Sphere,
!  modified to take HEALpix pixel array from Particles_Data

module HEALPixModule

!===============================================================================

contains

  subroutine mk_pix2xy()
    !=======================================================================
    !     constructs the array giving x and y in the face from pixel number
    !     for the nested (quad-cube like) ordering of pixels
    !
    !     the bits corresponding to x and y are interleaved in the pixel number
    !     one breaks up the pixel number by even and odd bits
    !=======================================================================
		!get 
    use   Particles_rayData, only: x2pix, y2pix, pix2x, pix2y, x2pix1, y2pix1
    implicit none
    INTEGER*4 ::  kpix, jpix, ix, iy, ip, id

    !cc cf block data      data      pix2x(1023) /0/
    !-----------------------------------------------------------------------
    !      print *, 'initiate pix2xy'
    do kpix=0,1023          ! pixel number
       jpix = kpix
       IX = 0
       IY = 0
       IP = 1               ! bit position (in x and y)
!        do while (jpix/=0) ! go through all the bits
       do
          if (jpix == 0) exit ! go through all the bits
          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in ix
          jpix = jpix/2
          IX = ID*IP+IX

          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in iy
          jpix = jpix/2
          IY = ID*IP+IY

          IP = 2*IP         ! next bit (in x and y)
       enddo
       pix2x(kpix) = IX     ! in 0,31
       pix2y(kpix) = IY     ! in 0,31
    enddo

    return
  end subroutine mk_pix2xy

	!only difference is that the loop starts at i=0
  subroutine mk_xy2pix1()
    !=======================================================================
    !     sets the array giving the number of the pixel lying in (x,y)
    !     x and y are in {1,128}
    !     the pixel number is in {0,128**2-1}
    !
    !     if  i-1 = sum_p=0  b_p * 2^p
    !     then ix = sum_p=0  b_p * 4^p
    !          iy = 2*ix
    !     ix + iy in {0, 128**2 -1}
    !=======================================================================
    use  Particles_rayData, only: x2pix, y2pix, pix2x, pix2y, x2pix1, y2pix1
    implicit none
    INTEGER*4:: k,ip,i,j,id
    !=======================================================================
    do i = 0,127           !for converting x,y into
       j  = i           !pixel numbers
       k  = 0
       ip = 1

       do
          if (j==0) then
             x2pix1(i) = k
             y2pix1(i) = 2*k
             exit
          else
             id = MODULO(J,2)
             j  = j/2
             k  = ip*id+k
             ip = ip*4
          endif
       enddo

    enddo

    return
  end subroutine mk_xy2pix1
  !=======================================================================

  subroutine mk_xy2pix()
		use  Particles_rayData, only: x2pix, y2pix, pix2x, pix2y, x2pix1, y2pix1
		implicit none
    !=======================================================================
    !     sets the array giving the number of the pixel lying in (x,y)
    !     x and y are in {1,128}
    !     the pixel number is in {0,128**2-1}
    !
    !     if  i-1 = sum_p=0  b_p * 2^p
    !     then ix = sum_p=0  b_p * 4^p
    !          iy = 2*ix
    !     ix + iy in {0, 128**2 -1}
    !=======================================================================
    INTEGER*4:: k,ip,i,j,id
    !=======================================================================

    do i = 1,128           !for converting x,y into
       j  = i-1            !pixel numbers
       k  = 0
       ip = 1

       do
          if (j==0) then
             x2pix(i) = k
             y2pix(i) = 2*k
             exit
          else
             id = MODULO(J,2)
             j  = j/2
             k  = ip*id+k
             ip = ip*4
          endif
       enddo

    enddo

    return
  end subroutine mk_xy2pix

! gives normal vector for a pixel number

!=======================================================================
!     pix2vec_nest
!
!     renders vector (x,y,z) coordinates of the nominal pixel center
!     for the pixel number ipix (NESTED scheme)
!     given the map resolution parameter nside
!     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
!     in the order N,W,S,E
!=======================================================================
  subroutine pix2vec_nest  (nside, ipix, vector, vertex)
    use  Particles_rayData, only: x2pix, y2pix, pix2x, pix2y, x2pix1, y2pix1
    implicit none

    INTEGER*4, INTENT(IN) :: nside
    INTEGER*8, INTENT(IN) :: ipix
    REAL*8,     INTENT(OUT), dimension(1:) :: vector
    REAL*8,     INTENT(OUT), dimension(1:,1:), optional :: vertex

    INTEGER*4 :: npix, npface, ipf
    INTEGER*8 :: ip_low, ip_trunc, ip_med, ip_hi
    INTEGER*8 :: face_num, ix, iy, kshift, scale, i, ismax
    INTEGER*8 :: jrt, jr, nr, jpt, jp, nl4
    REAL*8     :: z, fn, fact1, fact2, sth, phi

    ! coordinate of the lowest corner of each face
    INTEGER*4, dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER*4, dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2

    real*8 :: phi_nv, phi_wv, phi_sv, phi_ev, phi_up, phi_dn, sin_phi, cos_phi
    real*8 :: z_nv, z_sv, sth_nv, sth_sv
    real*8 :: hdelta_phi
    integer*8 :: iphi_mod, iphi_rat
    logical :: do_vertex
    integer*8 :: diff_phi,ns_max4 = 8192

    !-----------------------------------------------------------------------

    if (nside > ns_max4) call fatal_error("nside out of range")

    npix = 12*nside*nside       ! total number of points

    if (ipix <0 .or. ipix>npix-1) then 
			print*,ipix,npix
			call fatal_error("ipix out of range")
		endif

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    npface = nside * int(nside, 4)
    nl4    = 4*nside

    !     finds the face, and the number in the face
    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    do_vertex = .false.
    if (present(vertex)) then
       if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
          do_vertex = .true.
       else
          call fatal_error(" pix2vec_ring : vertex array has wrong size ")
       endif
    endif
    fn = real(nside, 8)
    fact1 = 1.0d0/(3.0d0*fn*fn)
    fact2 = 2.0d0/(3.0d0*fn)

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    if (nside <= ns_max4) then
       ip_low = iand(ipf,1023)       ! content of the last 10 bits
       ip_trunc =    ipf/1024        ! truncation of the last 10 bits
       ! fix for gfortran 9+, based on https://stackoverflow.com/questions/59756525/iand-with-different-kind-parameters-using-new-gfortran-version
       ip_med = iand(int(ip_trunc, kind(1023)),1023)  ! content of the next 10 bits
       !ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
       ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

       ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
       iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    else
       ix = 0
       iy = 0
       scale = 1
       ismax = 4
       do i=0, ismax
          ip_low = iand(ipf,1023)
          ix = ix + scale * pix2x(ip_low)
          iy = iy + scale * pix2y(ip_low)
          scale = scale * 32
          ipf   = ipf/1024
       enddo
       ix = ix + scale * pix2x(ipf)
       iy = iy + scale * pix2y(ipf)
    endif

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}


    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0d0 - nr*fact1*nr
       kshift = 0
       if (do_vertex) then
          z_nv = 1.0d0 - (nr-1)*fact1*(nr-1)
          z_sv = 1.0d0 - (nr+1)*fact1*(nr+1)
       endif

    else if (jr <= 3*nside) then ! equatorial region
       nr = nside
       z  = (2*nside-jr)*fact2
       ! fix for gfortran 9+, based on https://stackoverflow.com/questions/59756525/iand-with-different-kind-parameters-using-new-gfortran-version
       !kshift = iand(jr - nside, 1)
       kshift = iand(int(jr - nside, kind(1)), 1)
       if (do_vertex) then
          z_nv = (2*nside-jr+1)*fact2
          z_sv = (2*nside-jr-1)*fact2
          if (jr == nside) then ! northern transition
             z_nv =  1.0d0 - (nside-1) * fact1 * (nside-1)
          elseif (jr == 3*nside) then  ! southern transition
             z_sv = -1.0d0 + (nside-1) * fact1 * (nside-1)
          endif
       endif

    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1.0d0 + nr*fact1*nr
       kshift = 0
       if (do_vertex) then
          z_nv = - 1.0d0 + (nr+1)*fact1*(nr+1)
          z_sv = - 1.0d0 + (nr-1)*fact1*(nr-1)
       endif
    endif

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    phi = (jp - (kshift+1)*0.5d0) * (1.57079632679490 / nr)

    ! pixel center
    sth = SQRT((1.0d0-z)*(1.0d0+z))
    cos_phi = cos(phi)
    sin_phi = sin(phi)
    vector(1) = sth * cos_phi
    vector(2) = sth * sin_phi
    vector(3) = z

    if (do_vertex) then
       phi_nv = phi
       phi_sv = phi
       diff_phi = 0 ! phi_nv = phi_sv = phi

       phi_up = 0.0d0
       iphi_mod = MODULO(jp-1, nr) ! in {0,1,... nr-1}
       iphi_rat = (jp-1) / nr      ! in {0,1,2,3}
       if (nr > 1) phi_up = 1.57079632679490 * (iphi_rat +  iphi_mod   /real(nr-1,8))
       phi_dn             = 1.57079632679490 * (iphi_rat + (iphi_mod+1)/real(nr+1,8))
       if (jr < nside) then            ! North polar cap
          phi_nv = phi_up
          phi_sv = phi_dn
          diff_phi = 3 ! both phi_nv and phi_sv different from phi
       else if (jr > 3*nside) then     ! South polar cap
          phi_nv = phi_dn
          phi_sv = phi_up
          diff_phi = 3 ! both phi_nv and phi_sv different from phi
       else if (jr == nside) then      ! North transition
          phi_nv = phi_up
          diff_phi = 1
       else if (jr == 3*nside) then    ! South transition
          phi_sv = phi_up
          diff_phi = 2
       endif

       hdelta_phi = 3.14159265358979 / (4.0d0*nr)

       ! west vertex
       phi_wv      = phi - hdelta_phi
       vertex(1,2) = sth * COS(phi_wv)
       vertex(2,2) = sth * SIN(phi_wv)
       vertex(3,2) = z

       ! east vertex
       phi_ev      = phi + hdelta_phi
       vertex(1,4) = sth * COS(phi_ev)
       vertex(2,4) = sth * SIN(phi_ev)
       vertex(3,4) = z

       ! north and south vertices
       sth_nv = SQRT((1.0d0-z_nv)*(1.0d0+z_nv))
       sth_sv = SQRT((1.0d0-z_sv)*(1.0d0+z_sv))
       if (diff_phi == 0) then
          vertex(1,1) = sth_nv * cos_phi
          vertex(2,1) = sth_nv * sin_phi
          vertex(1,3) = sth_sv * cos_phi
          vertex(2,3) = sth_sv * sin_phi
       else
          vertex(1,1) = sth_nv * COS(phi_nv)
          vertex(2,1) = sth_nv * SIN(phi_nv)
          vertex(1,3) = sth_sv * COS(phi_sv)
          vertex(2,3) = sth_sv * SIN(phi_sv)
       endif
       vertex(3,1) = z_nv
       vertex(3,3) = z_sv
    endif

    return

  end subroutine pix2vec_nest

! gives pixel number for a normal vector

!=======================================================================
!     vec2pix_nest
!
!     renders the pixel number ipix (NESTED scheme) for a pixel which contains
!     a point on a sphere at coordinate vector (=x,y,z), given the map
!     resolution parameter nside
!
! 2009-03-10: calculations done directly at nside rather than ns_max
!=======================================================================

  subroutine vec2pix_nest  (nside, vector, ipix)
    use Particles_rayData, only: x2pix, y2pix, pix2x, pix2y, x2pix1, y2pix1
    implicit none
    INTEGER*4, INTENT(IN)                 :: nside
    REAL*8,     INTENT(IN), dimension(1:) :: vector
    INTEGER*4, INTENT(OUT)                :: ipix

    integer*4 :: ipf, scale, scale_factor,ns_max4 = 8192

    REAL*8    ::  z, za, tt, tp, tmp, dnorm, phi
    INTEGER*4 ::  jp, jm, ifp, ifm, face_num, &
         &     ix, iy, ix_low, iy_low, ntt, i, ismax, ipix4

    !-----------------------------------------------------------------------

    if (nside<1 .or. nside>ns_max4) call fatal_error("nside out of range")

    dnorm = SQRT(vector(1)**2+vector(2)**2+vector(3)**2)
    z = vector(3) / dnorm
    phi = 0.0d0
    if (vector(1) /= 0.0d0 .or. vector(2) /= 0.0d0) &
         &     phi = ATAN2(vector(2),vector(1)) ! phi in ]-pi,pi]

    za = ABS(z)
    if (phi < 0.0)    phi = phi + 6.28318530717959 ! phi in [0,2pi[
    tt = phi / 1.57079632679490 ! in [0,4[
    if (x2pix1(127) <= 0) call mk_xy2pix1()

!		changed twothird to 0.666666666, CB
    if (za <= 0.666666666666666) then ! equatorial region

       !        (the index of edge lines increase when the longitude=phi goes up)
       jp = INT(nside*(0.5d0 + tt - z*0.75d0)) !  ascending edge line index
       jm = INT(nside*(0.5d0 + tt + z*0.75d0)) ! descending edge line index

       !        finds the face
       ifp = jp / nside  ! in {0,4}
       ifm = jm / nside
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = iand(ifp,3) + 4
       else if (ifp < ifm) then     ! (half-)faces 0 to 3
          face_num = iand(ifp,3)
       else                            ! (half-)faces 8 to 11
          face_num = iand(ifm,3) + 8
       endif

       ix =         iand(jm, nside-1)
       iy = nside - iand(jp, nside-1) - 1

    else ! polar region, za > 2/3

       ntt = INT(tt)
       if (ntt >= 4) ntt = 3
       tp = tt - ntt
       tmp = SQRT( 3.0d0*(1.0d0 - za) )  ! in ]0,1]

       !        (the index of edge lines increase when distance from the closest pole goes up)
       jp = INT( nside * tp          * tmp ) ! line going toward the pole as phi increases
       jm = INT( nside * (1.0d0 - tp) * tmp ) ! that one goes away of the closest pole
       jp = MIN(nside-1, jp) ! for points too close to the boundary
       jm = MIN(nside-1, jm)

       !        finds the face and pixel's (x,y)
       if (z >= 0) then
          face_num = ntt  ! in {0,3}
          ix = nside - jm - 1
          iy = nside - jp - 1
       else
          face_num = ntt + 8 ! in {8,11}
          ix =  jp
          iy =  jm
       endif

    endif

    if (nside <= ns_max4) then 
       ix_low = iand(ix, 127)
       iy_low = iand(iy, 127)
       ipf =     x2pix1(ix_low) + y2pix1(iy_low) &
            & + (x2pix1(ix/128) + y2pix1(iy/128)) * 16384
    else
       scale = 1
       scale_factor = 16384 ! 128*128
       ipf = 0
       ismax = 1 ! for nside in [2^14, 2^20]
       if (nside >  1048576 ) ismax = 3
       do i=0, ismax
          ix_low = iand(ix, 127) ! last 7 bits
          iy_low = iand(iy, 127) ! last 7 bits
          ipf = ipf + (x2pix1(ix_low)+y2pix1(iy_low)) * scale
          scale = scale * scale_factor
          ix  =     ix / 128 ! truncate out last 7 bits
          iy  =     iy / 128
       enddo
       ipf =  ipf + (x2pix1(ix)+y2pix1(iy)) * scale
    endif
    ipix = ipf + face_num* int(nside,4) * nside    ! in {0, 12*nside**2 - 1}

    return
end subroutine vec2pix_nest

function nside2npix(nside) result(npix_result)
  !=======================================================================
  ! given nside, returns npix such that npix = 12*nside^2
  !  nside should be a power of 2 smaller than ns_max
  !  if not, -1 is returned
  ! EH, Feb-2000
  ! 2009-03-04: returns i8b result, faster
  !=======================================================================
  implicit none
  INTEGER*8             :: npix_result
  INTEGER*4, INTENT(IN) :: nside
  INTEGER*4							:: ns_max4 = 8192

  INTEGER*8 :: npix
  CHARACTER(LEN=*), PARAMETER :: code = "nside2npix"
  !=======================================================================

  npix = (12*nside)*nside
  if (nside < 1 .or. nside > ns_max4 ) then !.or. iand(nside-1,nside) /= 0) then
     print*,code//": Nside=",nside," is not a power of 2."
     npix = -1
  endif
  npix_result = npix

  return
end function nside2npix

!=======================================================================
!     pix2ang_nest
!
!     renders theta and phi coordinates of the nominal pixel center
!     for the pixel number ipix (NESTED scheme)
!     given the map resolution parameter nside
!=======================================================================

subroutine pix2ang_nest  (nside, ipix, theta, phi)
    use  Particles_rayData, only: x2pix, y2pix, pix2x, pix2y, x2pix1, y2pix1
    implicit none
    INTEGER*4, INTENT(IN)  :: nside
    INTEGER*4, INTENT(IN)  :: ipix
    REAL*8,     INTENT(OUT) :: theta, phi

    INTEGER*4 :: npix, npface, ipf,ns_max4=8192
    INTEGER*4 :: ip_low, ip_trunc, ip_med, ip_hi, &
         &     jrt, jr, nr, jpt, jp, kshift, nl4, scale, i, ismax
    INTEGER*4 :: ix, iy, face_num
    REAL*8     :: z, fn, fact1, fact2

    ! coordinate of the lowest corner of each face
    INTEGER*4, dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER*4, dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
    !-----------------------------------------------------------------------

    if (nside > ns_max4) call fatal_error("nside out of range")

    npix = 12*nside*nside       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error ("ipix out of range")

    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    npface = nside * int(nside, 4)
    nl4    = 4*nside

    !     finds the face, and the number in the face
    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    fn = real(nside, 8)
    fact1 = 1.0d0/(3.0d0*fn*fn)
    fact2 = 2.0d0/(3.0d0*fn)

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    if (nside <= ns_max4) then
       ip_low = iand(ipf,1023)       ! content of the last 10 bits
       ip_trunc =    ipf/1024        ! truncation of the last 10 bits
       ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
       ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

       ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
       iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    else
       ix = 0
       iy = 0
       scale = 1
       ismax = 4
       do i=0, ismax
          ip_low = iand(ipf,1023)
          ix = ix + scale * pix2x(ip_low)
          iy = iy + scale * pix2y(ip_low)
          scale = scale * 32
          ipf   = ipf/1024
       enddo
       ix = ix + scale * pix2x(ipf)
       iy = iy + scale * pix2y(ipf)
    endif

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0d0 - nr * fact1 * nr
       kshift = 0

    else if (jr <= 3*nside) then ! equatorial region
       nr = nside
       z  = (2*nside-jr)*fact2
       kshift = iand(jr - nside, 1)

    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1.0d0 + nr * fact1 * nr
       kshift = 0
    endif

    theta = ACOS(z)

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4

    phi = (jp - (kshift+1)*0.5d0) * (1.57079632679490 / nr)

    return

  end subroutine pix2ang_nest

!	HEALPix error handlers
  subroutine fatal_error (msg)
    character(len=*), intent(in), optional :: msg

    if (present(msg)) then
       print *,'Fatal error: ', trim(msg)
    else
       print *,'Fatal error'
    endif
    call exit_with_status(1)
  end subroutine fatal_error


!  subroutine fatal_error_msg (msg)
!    character(len=*), intent(in) :: msg
!       print *,'Fatal error: ', trim(msg)
!    call exit_with_status(1)
!  end subroutine fatal_error_msg

  ! ===========================================================
  subroutine exit_with_status (code, msg)
    ! ===========================================================
    integer, intent(in) :: code
    character (len=*), intent(in), optional :: msg
    ! ===========================================================

    if (present(msg)) print *,trim(msg)
    print *,'program exits with exit code ', code

    call exit (code)

  end subroutine exit_with_status


end module HEALPixModule
