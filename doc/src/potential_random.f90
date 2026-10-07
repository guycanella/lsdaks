!> Module for random (disorder) potentials
!!
!! Implements disordered potentials modeling Anderson localization.
!! Two distributions are available:
!! - Uniform distribution: V(i) ~ U[-W/2, W/2]
!! - Gaussian distribution: V(i) ~ N(0, σ²)
!!
!! Physical interpretation:
!! - Models disorder in condensed matter systems
!! - For large W or σ: Anderson localization (exponential decay of wavefunctions)
!! - For small W or σ: Extended states (weak disorder)
!! - Average: ⟨V(i)⟩ = 0 (typically)
!!
!! Reference: P.W. Anderson, "Absence of Diffusion in Certain Random Lattices" (1958)
module potential_random
    use lsda_constants, only: dp, PI
    use lsda_errors, only: ERROR_SUCCESS, ERROR_NEGATIVE_VALUE
    implicit none
    private
    
    public :: potential_random_uniform
    public :: potential_random_gaussian
    public :: box_muller_pair

contains

    !> Box-Muller transform of one pair of uniform deviates
    !!
    !! Maps a pair `(u1, u2)` drawn from the uniform distribution on [0, 1) -
    !! which is exactly the range documented for `random_number` - onto a pair
    !! of independent standard normal deviates.
    !!
    !! The radial factor is `sqrt(-2*log(1 - u1))` and NOT the textbook
    !! `sqrt(-2*log(u1))`: because the source interval is half-open and closed
    !! on the left, `u1` may legitimately come back as exactly 0, and
    !! `log(0) = -Inf` makes the radius +Inf and poisons the whole potential
    !! array with non-finite values (which then surfaced far from here, as
    !! ERROR_NOT_A_NUMBER inside validate_hamiltonian_inputs). Using
    !! `1 - u1 ∈ (0, 1]` keeps the argument of the logarithm strictly positive
    !! for every value `random_number` can produce, while sampling the same
    !! distribution: if u1 is uniform on [0,1), so is 1 - u1 on (0,1].
    !!
    !! This is a separate pure routine precisely so that the `u1 = 0` corner
    !! can be exercised deterministically by a test instead of being left to
    !! the (astronomically unlikely) chance of the generator emitting it.
    !!
    !! @param[in]  u1  First uniform deviate in [0, 1)
    !! @param[in]  u2  Second uniform deviate in [0, 1)
    !! @param[out] z1  First standard normal deviate, N(0, 1)
    !! @param[out] z2  Second standard normal deviate, N(0, 1)
    pure subroutine box_muller_pair(u1, u2, z1, z2)
        real(dp), intent(in) :: u1, u2
        real(dp), intent(out) :: z1, z2

        real(dp) :: radius

        radius = sqrt(-2.0_dp * log(1.0_dp - u1))
        z1 = radius * cos(2.0_dp * PI * u2)
        z2 = radius * sin(2.0_dp * PI * u2)
    end subroutine box_muller_pair

    !> Random potential with uniform distribution: V(i) ~ U[-W, W]
    !!
    !! Creates a disordered potential with uniform distribution.
    !! Each site has an independent random value drawn from [-W, W].
    !! Matches C++ random1 implementation: V*(2*rand - 1).
    !!
    !! Physical properties:
    !! - Mean: ⟨V(i)⟩ = 0
    !! - Variance: σ² = W²/3
    !! - Models box disorder (flat probability distribution)
    !!
    !! @param[in]  W      Disorder strength (half-width of distribution), W > 0
    !! @param[in]  L      Number of lattice sites
    !! @param[in]  seed   Random seed for reproducibility (use system time if < 0)
    !! @param[out] V      Potential array V(i) for i = 1..L
    !! @param[out] ierr   Error flag (ERROR_SUCCESS or ERROR_NEGATIVE_VALUE)
    !!
    !! @note For W → 0: No disorder (all sites have V ≈ 0)
    !! @note For W >> t (hopping): Strong disorder, Anderson localization
    !! @note Range is [-W, +W] with total width = 2W (matches C++ random1)
    subroutine potential_random_uniform(W, L, seed, V, ierr)
        real(dp), intent(in) :: W
        integer, intent(in) :: L, seed
        real(dp), dimension(L), intent(out) :: V
        integer, intent(out) :: ierr
        integer :: i, seed_size
        integer, allocatable :: seed_array(:)
        real(dp) :: rand_val

        ierr = ERROR_SUCCESS

        if (W < 0.0_dp) then
            ierr = ERROR_NEGATIVE_VALUE
            V = 0.0_dp
            return
        end if

        ! Special case: W = 0 means no disorder. `abs(W) < tiny(W)` is written
        ! instead of `W == 0` for readability and to avoid the equality-comparison
        ! warning, NOT because it is a looser test: for every normal number the
        ! two conditions are identical, and they differ only for subnormal W
        ! (|W| < ~2.2e-308), which the branch then also treats as no disorder.
        ! That is deliberate - a subnormal disorder width is zero for every
        ! physical purpose - but it is the whole of the difference.
        if (abs(W) < tiny(W)) then
            V = 0.0_dp
            return
        end if

        if (seed >= 0) then
            call random_seed(size=seed_size)
            allocate(seed_array(seed_size))
            seed_array = seed
            call random_seed(put=seed_array)
            deallocate(seed_array)
        else
            call random_seed()
        end if

        do i = 1, L
            call random_number(rand_val)
            ! C++ formula: V*(2.0*rand - 1.0) gives [-V, +V]
            ! Matches C++ random1 implementation exactly
            V(i) = W * (2.0_dp * rand_val - 1.0_dp)  ! Maps [0,1] → [-W, +W]
        end do
    end subroutine potential_random_uniform

    !> Random potential with Gaussian distribution: V(i) ~ N(0, σ²)
    !!
    !! Creates a disordered potential with Gaussian (normal) distribution.
    !! Each site has an independent random value drawn from N(0, σ²).
    !!
    !! Physical properties:
    !! - Mean: ⟨V(i)⟩ = 0
    !! - Standard deviation: σ
    !! - Models thermal/quantum fluctuations
    !! - More realistic than uniform disorder for many physical systems
    !!
    !! @param[in]  sigma  Standard deviation of disorder, σ > 0
    !! @param[in]  L      Number of lattice sites
    !! @param[in]  seed   Random seed for reproducibility (use system time if < 0)
    !! @param[out] V      Potential array V(i) for i = 1..L
    !! @param[out] ierr   Error flag (ERROR_SUCCESS or ERROR_NEGATIVE_VALUE)
    !!
    !! @note Uses the Box-Muller transform (see `box_muller_pair`) to generate
    !!       Gaussian random numbers; the transform is written so that it stays
    !!       finite for every value `random_number` can return, including an
    !!       exact 0
    !! @note For σ → 0: No disorder (all sites have V ≈ 0)
    !! @note For σ >> t (hopping): Strong disorder, Anderson localization
    subroutine potential_random_gaussian(sigma, L, seed, V, ierr)
        real(dp), intent(in) :: sigma
        integer, intent(in) :: L, seed
        real(dp), dimension(L), intent(out) :: V
        integer, intent(out) :: ierr
        integer :: i, seed_size
        integer, allocatable :: seed_array(:)
        real(dp) :: u1, u2, z1, z2

        ierr = ERROR_SUCCESS

        if (sigma < 0.0_dp) then
            ierr = ERROR_NEGATIVE_VALUE
            V = 0.0_dp
            return
        end if

        ! Special case: sigma = 0 means no disorder. Same threshold as for W
        ! above, and the same caveat: for normal numbers this is exactly
        ! `sigma == 0`; only subnormal sigma is additionally caught.
        if (abs(sigma) < tiny(sigma)) then
            V = 0.0_dp
            return
        end if

        if (seed >= 0) then
            call random_seed(size=seed_size)
            allocate(seed_array(seed_size))
            seed_array = seed
            call random_seed(put=seed_array)
            deallocate(seed_array)
        else
            call random_seed()
        end if

        ! Generate Gaussian random values using Box-Muller transform
        do i = 1, L, 2
            call random_number(u1)
            call random_number(u2)

            call box_muller_pair(u1, u2, z1, z2)

            V(i) = sigma * z1
            if (i + 1 <= L) then
                V(i + 1) = sigma * z2
            end if
        end do
    end subroutine potential_random_gaussian

end module potential_random