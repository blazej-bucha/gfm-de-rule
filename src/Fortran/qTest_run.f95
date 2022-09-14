program qTest_run
!
! ==========================================================================
!
! This program computes the gravitation potential of the Bennu asteroid in
! quadruple precision. The shape of the asteroid is given by a surface 
! spherical harmonic coefficients (in metres; the structure of the file is 
! explained, for instance, in "./qanm_bnm.f95")
!
! ../../data/Bennu_Shape_SHCs_to15.txt
!
! and the mass of the asteroid is assumed to be constant. Three evaluation
! points (one located inside and two located outside the body) are defined in
!
! ../../data/Bennu_Computing_points.txt
!
! The file contains three columns representing the spherical latitude, 
! the spherical longitude (both in degrees) and the spherical radius 
! (in metres) of the evaluation points.
!
! In this example, the same number of digits in input data files is used for
! double precision computations. This, however, can be tailored by
! the user.
!
! The loop over the evaluation points is paralellized using OpenMP. By default,
! this program employs automatically all available threads (the value can
! be changed, however; see the variable "nmaxthreads").
!
!
! NOTE: In case of parallel computation using OpenMP, the user may need to set
!       the value of the "OMP_STACKSIZE" environment variable to a sufficiently
!       large value as well as to set the OS limits. For instance, the 
!       following settings may be sufficient (Linux bash commands):
!
!       ulimit -s unlimited
!       export OMP_STACKSIZE=2g
!
!
! REFERENCES: Fukushima, T., 2017. Precise and fast computation of the
!                 gravitational field of a general finite body and its
!                 application to the gravitational study of asteroid Eros. The
!                 Astronomical Journal 154(145):15pp,
!                 doi: 10.3847/1538-3881/aa88b8
!
! ==========================================================================


    ! Import modules
    ! ===============================================
    use vartypes  ! Module defining floating point double and quadruple
                  ! precision numbers
    use constants ! Module defining constants (pi, etc.)
    !$ use omp_lib, only: omp_get_max_threads ! Module to get the maximum
                                              ! number of available threads. By
                                              ! default, the maximum number is
                                              ! used for parallel parts of the
                                              ! synthesis and 
                                              ! analysis. However, the value 
                                              ! can be changed in the
                                              ! code below (look for
                                              ! "omp_set_num_threads").
    ! ===============================================


    ! Declaration of variables
    ! ===============================================
    implicit none

    integer :: nmax = 15 ! Maximum harmonic degree of the surface spherical
                         ! harmonic expansion of the body's surface (must
                         ! correspond to the imported "cs" matrix, see below)

    integer :: npoints = 3 ! Number of evaluation points in the data file
                           ! containing the evaluation points

    real(qp) :: rb = 0.0_qp ! Radius (in metres) of a ball, the gravitational
                            ! potential of which will be subtracted from the
                            ! output potential.
                            ! If "rb = 0.0_qp", the gravitational field of the
                            ! entire body (as defined by its spherical harmonic
                            ! coefficients) is computed.

    real(qp) :: rho   = 1260.0_qp ! Constant mass density of the body

    real(qp) :: gc    = 6.67384e-11_qp ! Newton's gravitational constant

    real(qp) :: delta = 1e-32_qp ! Relative error tolerance (see
                                 ! Fukushima 2017). For instance,
                                 ! "delta = 1e-16_qp" ensures about 15- or
                                 ! 16-digit accuracy of the output 
                                 ! gravitational potential.

    real(qp), allocatable :: lat(:) ! Spherical latitudes of the evaluation
                                    ! points (input values in degrees)
    real(qp), allocatable :: lon(:) ! Spherical longitudes of the evaluation
                                    ! points (input values in degrees)
    real(qp), allocatable :: r(:)   ! Spherical radii of the evaluation points
                                    ! (input values in metres)

    real(qp) :: lati, loni, ri ! The "i"th spherical latitude, longitude and
                               ! spherical radius

    real(qp), allocatable :: v(:) ! Output gravitational potential at
                                  ! evaluation points in m**2 * s**-2

    real(qp) :: vi1, vi2 ! Partial contributions to the gravitational potential
                         ! of the "i"th evaluation point

    real(qp), allocatable :: cs(:,:) ! Spherical harmonic coefficients defining
                                     ! the gravitating body (for the structure
                                     ! of the matrix, see, e.g., the 
                                     ! "dqde1.f95" subroutine). The size of the 
                                     ! matrix must correspond to the specified 
                                     ! "nmax" value

    real(qp), allocatable :: anm(:,:), bnm(:,:) ! Auxiliary coefficients to
                                                ! compute Legendre functions

    real(qp) :: r0 ! Reference radius in metres. This value will be set to the
                   ! body's mean sphere later as an effective choice

    integer :: fid ! Unit number for file reading

    integer :: nmaxthreads ! Number of threads used for the parallel loop over
                           ! evaluation points

    integer :: n, m ! Harmonic degree and order, respectively
    integer :: i    ! Loop variable for evaluation points

    character(len = 256) :: path ! Relative path to the directory, where the
                                 ! input data files are stored.
    ! ===============================================


    ! Main computational part of the program
    ! ===============================================

    ! Relative path to the directory, where the input data files are stored.
    ! -----------------------------------------------
    path = '../data'
    ! -----------------------------------------------


    ! Get and then set the maximum number of available threads for parallel
    ! parts of the package.
    ! -----------------------------------------------
    !$ nmaxthreads = omp_get_max_threads()
    !$ call omp_set_num_threads(nmaxthreads)
    ! -----------------------------------------------


    ! Import spherical harmonic coefficients defining the gravitating body
    ! -----------------------------------------------
    print *, 'Reading spherical harmonic coefficients...'

    open(newunit = fid, file = &
         trim(adjustl(path))//'/Bennu_Shape_SHCs_to15.txt', status = 'old')

    allocate(cs(nmax + 1, nmax + 1))

    ! Loop over harmonic degrees (rows) in the input data file with spherical
    ! harmonic coefficients
    do n = 1, (nmax + 1)
        read(fid, *) cs(n, :)
    end do

    close(fid)
    ! -----------------------------------------------


    ! Import evaluation points
    ! -----------------------------------------------
    print *, 'Reading evaluation points...'

    open(newunit = fid, file = &
         trim(adjustl(path))//'/Bennu_Computing_points.txt', &
         status = 'old')

    allocate(lat(npoints), lon(npoints), r(npoints))

    ! Loop over evaluation points in the input data file with evaluation points
    do i = 1, npoints
        read(fid, *) lat(i), lon(i), r(i)
    end do

    close(fid)
    ! -----------------------------------------------


    ! Transform latitudes and longitudes of the evaluation points to radians
    ! -----------------------------------------------
    lat = lat * pi_qp / 180.0_qp
    lon = lon * pi_qp / 180.0_qp
    ! -----------------------------------------------


    ! Set the reference radius to the mean sphere
    ! -----------------------------------------------
    r0 = cs(1,1)
    ! -----------------------------------------------


    ! Precompute "anm" and "bnm" coefficients for Legendre recurrences
    ! -----------------------------------------------
    print *, 'Computing "anm" and "bnm" coefficients...'

    allocate(anm(nmax + 1, nmax + 1), bnm(nmax + 1, nmax + 1))

    anm = 0.0_qp
    bnm = 0.0_qp

    do n = 0, nmax ! Loop over harmonic degrees
        do m = 0, (n - 1) ! Loop over harmonic orders
            anm(n + 1, m + 1) = sqrt((2.0_qp * real(n, qp) - 1.0_qp) &
                                     * (2.0_qp * real(n, qp) + 1.0_qp) &
                                     / ((real(n, qp) - real(m, qp)) &
                                     * (real(n, qp) + real(m, qp))))

            bnm(n + 1, m + 1) = sqrt(((2.0_qp * real(n, qp) + 1.0_qp) &
                                     * (real(n, qp) + real(m, qp) - 1.0_qp) &
                                     * (real(n, qp) - real(m, qp) - 1.0_qp)) &
                                     / ((real(n, qp) - real(m, qp)) &
                                     * (real(n, qp) + real(m, qp)) &
                                     * (2.0_qp * real(n, qp) - 3.0_qp)))
        end do ! End of the loop over harmonic orders
    end do ! End of the loop over harmonic degrees
    ! -----------------------------------------------


    ! Computation of the gravitational potential
    ! -----------------------------------------------
    print *, 'Computing the gravitational potential...'

    allocate(v(npoints))

    ! Parallel loop over the evaluation points
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& SHARED(npoints, lat, lon, r, delta, cs, nmax, rb, r0) &
    !$OMP& SHARED(anm, bnm, gc, rho, v) &
    !$OMP& PRIVATE(i, lati, loni, ri, vi1, vi2)
    do i = 1, npoints

        lati = lat(i)
        loni = lon(i)
        ri = r(i)

        call qqde1(-pi_qp / 2.0_qp - lati, 0.0_qp, delta, lati, loni, ri, &
                   cs, nmax, rb, r0, anm, bnm, &
                   vi1)

        call qqde1(0.0_qp, pi_qp / 2.0_qp - lati, delta, lati, loni, ri, &
                   cs, nmax, rb, r0, anm, bnm, &
                   vi2)

        v(i) = vi1 + vi2

    end do ! End of the loop over evaluation points
    !$OMP END PARALLEL DO

    v = gc * r0**2 * rho * v
    ! -----------------------------------------------
    ! ===============================================


    ! Print outputs
    ! ===============================================
    print *, ''
    print *, ''
    print *, 'Validation'
    print *, '========================================'
    print *, 'Latitudes of the evaluation points (in degrees)'
    print *, (lat * 180.0_qp / pi_qp)
    print *, ''
    print *, 'Longitudes of the evaluation points (in degrees)'
    print *, (lon * 180.0_qp / pi_qp)
    print *, ''
    print *, 'Radii of the evaluation points (in metres)'
    print *, r
    print *, ''
    print *, 'Computed gravitational potential (in m**2 * s**-2)'
    print *, v
    print *, ''
    print *, 'Reference values of the gravitational potential (in m**2 * s**-2)'
    print *, '  3.01907463537325614667197010432098283E-0002&
              &   1.73060538290977726180610558316772352E-0002&
              &   5.23751595802099559263452974360152990E-0003'
    print *, '========================================'

    print *, ''
    print *, ''
    print *, 'End of the program'
    ! ===============================================

end program qTest_run
