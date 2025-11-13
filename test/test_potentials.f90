!> Unit tests for potential modules
program test_potentials
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_potential_tests())

contains

    function get_potential_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("uniform_constant_value", test_uniform_constant_value), &
            test("harmonic_symmetry", test_harmonic_symmetry), &
            test("harmonic_center_minimum", test_harmonic_center_minimum), &
            test("impurity_single_position", test_impurity_single_position), &
            test("impurity_single_bounds", test_impurity_single_bounds), &
            test("impurity_multiple_positions", test_impurity_multiple_positions), &
            test("impurity_multiple_overlap", test_impurity_multiple_overlap), &
            test("impurity_random_concentration", test_impurity_random_concentration), &
            test("random_uniform_mean", test_random_uniform_mean), &
            test("random_gaussian_mean", test_random_gaussian_mean), &
            test("barrier_single_width", test_barrier_single_width), &
            test("barrier_single_bounds", test_barrier_single_bounds), &
            test("barrier_double_well_separation", test_barrier_double_well_separation), &
            test("barrier_double_no_overlap", test_barrier_double_no_overlap), &
            test("factory_uniform", test_factory_uniform), &
            test("factory_harmonic", test_factory_harmonic), &
            test("factory_invalid_type", test_factory_invalid_type) &
        ])
    end function get_potential_tests

    !> Test uniform potential returns constant value
    !!
    !! Physics: A uniform potential V(i) = V₀ represents a global energy shift.
    !! This does not affect the physics of the system, only the absolute energy scale.
    !! All sites should have the same potential value.
    subroutine test_uniform_constant_value()
        use fortuno_serial, only: check => serial_check
        use potential_uniform, only: apply_potential_uniform
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10
        real(dp) :: V(L)
        real(dp), parameter :: V0 = 2.5_dp
        integer :: i, ierr

        call apply_potential_uniform(V0, L, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")

        do i = 1, L
            call check(abs(V(i) - V0) < TOL, "All sites should have value V0")
        end do
    end subroutine test_uniform_constant_value

    !> Test harmonic potential has parity symmetry
    !!
    !! Physics: The harmonic trap V(i) = 0.5·k·(i-center)² creates a parabolic
    !! confining potential that models optical traps in cold atom systems.
    !! Due to the (i-center)² dependence, the potential must have parity symmetry:
    !! V(center+d) = V(center-d) for any displacement d.
    subroutine test_harmonic_symmetry()
        use fortuno_serial, only: check => serial_check
        use potential_harmonic, only: apply_potential_harmonic
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 11  ! Odd for symmetric center
        real(dp) :: V(L), k
        integer :: i, center, ierr

        k = 0.1_dp
        center = (L + 1) / 2

        call apply_potential_harmonic(k, L, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")

        ! Check parity symmetry: V(center+d) = V(center-d)
        do i = 1, center - 1
            call check(abs(V(i) - V(L + 1 - i)) < TOL, &
                       "Harmonic potential should have parity symmetry")
        end do
    end subroutine test_harmonic_symmetry

    !> Test harmonic potential has minimum at center
    !!
    !! Physics: The confining nature of the harmonic trap means particles
    !! are attracted to the center (lowest energy). The potential minimum
    !! must be at the center position i_center = (L+1)/2.
    subroutine test_harmonic_center_minimum()
        use fortuno_serial, only: check => serial_check
        use potential_harmonic, only: apply_potential_harmonic
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 21
        real(dp) :: V(L), k
        integer :: center, ierr

        k = 0.2_dp
        center = (L + 1) / 2

        call apply_potential_harmonic(k, L, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")

        ! Center should be the minimum
        call check(V(center) < V(1), "Center should have lower energy than edge")
        call check(V(center) < V(L), "Center should have lower energy than edge")
        call check(abs(V(center)) < TOL, "Center should have V ≈ 0")
    end subroutine test_harmonic_center_minimum

    !> Test single impurity is at correct position
    !!
    !! Physics: A point impurity creates a localized perturbation at a single site.
    !! For V_imp > 0 (repulsive), particles are scattered. For V_imp < 0 (attractive),
    !! bound states can form. The impurity should only affect the specified site.
    subroutine test_impurity_single_position()
        use fortuno_serial, only: check => serial_check
        use potential_impurity, only: potential_impurity_single
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10, i_imp = 5
        real(dp) :: V(L)
        real(dp), parameter :: V_imp = 3.0_dp
        integer :: ierr, i

        call potential_impurity_single(V_imp, i_imp, L, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")
        call check(abs(V(i_imp) - V_imp) < TOL, "Impurity site should have V_imp")

        do i = 1, L
            if (i /= i_imp) then
                call check(abs(V(i)) < TOL, "Non-impurity sites should be zero")
            end if
        end do
    end subroutine test_impurity_single_position

    !> Test single impurity bounds checking
    !!
    !! Physics: The impurity position must be within the physical lattice [1, L].
    !! Positions outside this range are unphysical and should be rejected.
    subroutine test_impurity_single_bounds()
        use fortuno_serial, only: check => serial_check
        use potential_impurity, only: potential_impurity_single
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_OUT_OF_BOUNDS

        integer, parameter :: L = 10
        real(dp) :: V(L)
        integer :: ierr

        ! Test i_imp = 0 (below bounds)
        call potential_impurity_single(1.0_dp, 0, L, V, ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "i_imp = 0 should fail")

        ! Test i_imp = L+1 (above bounds)
        call potential_impurity_single(1.0_dp, L + 1, L, V, ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "i_imp > L should fail")
    end subroutine test_impurity_single_bounds

    !> Test multiple impurities are placed correctly
    !!
    !! Physics: Multiple impurities model disorder or multiple scattering centers.
    !! Each impurity independently perturbs the electronic structure, and their
    !! effects can interfere (constructive or destructive) in transport properties.
    subroutine test_impurity_multiple_positions()
        use fortuno_serial, only: check => serial_check
        use potential_impurity, only: potential_impurity_multiple
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 20, N_imp = 3
        real(dp) :: V(L), V_imp_array(N_imp)
        integer :: imp_positions(N_imp)
        integer :: ierr

        V_imp_array = [2.0_dp, -1.5_dp, 3.0_dp]
        imp_positions = [5, 10, 15]

        call potential_impurity_multiple(V_imp_array, imp_positions, L, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")
        call check(abs(V(5) - 2.0_dp) < TOL, "First impurity at position 5")
        call check(abs(V(10) + 1.5_dp) < TOL, "Second impurity at position 10")
        call check(abs(V(15) - 3.0_dp) < TOL, "Third impurity at position 15")
    end subroutine test_impurity_multiple_positions

    !> Test overlapping impurities add their amplitudes
    !!
    !! Physics: When two impurities occupy the same site, their potentials
    !! superpose (linear addition). This is a consequence of the linearity
    !! of the one-particle Hamiltonian.
    subroutine test_impurity_multiple_overlap()
        use fortuno_serial, only: check => serial_check
        use potential_impurity, only: potential_impurity_multiple
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10, N_imp = 2
        real(dp) :: V(L), V_imp_array(N_imp)
        integer :: imp_positions(N_imp)
        integer :: ierr

        ! Two impurities at the same position
        V_imp_array = [2.0_dp, 3.0_dp]
        imp_positions = [5, 5]

        call potential_impurity_multiple(V_imp_array, imp_positions, L, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")
        call check(abs(V(5) - 5.0_dp) < TOL, "Overlapping impurities should add")
    end subroutine test_impurity_multiple_overlap

    !> Test random impurity concentration is correct
    !!
    !! Physics: Random impurities with a specified concentration model dilute
    !! magnetic impurities or defects in condensed matter. The concentration
    !! parameter determines the disorder strength (number of scattering centers).
    subroutine test_impurity_random_concentration()
        use fortuno_serial, only: check => serial_check
        use potential_impurity, only: potential_impurity_random
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 100, seed = 42
        real(dp) :: V(L), concentration
        integer, allocatable :: imp_positions(:)
        integer :: ierr, N_imp, count_nonzero

        concentration = 10.0_dp  ! 10% of sites

        call potential_impurity_random(2.0_dp, concentration, L, seed, V, imp_positions, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")

        N_imp = size(imp_positions)
        call check(N_imp == 10, "Should have 10 impurities for 10% concentration")

        ! Count non-zero sites
        count_nonzero = count(abs(V) > TOL)
        call check(count_nonzero == N_imp, "Number of non-zero sites should match N_imp")

        deallocate(imp_positions)
    end subroutine test_impurity_random_concentration

    !> Test random uniform potential has zero mean
    !!
    !! Physics: Random disorder with V(i) ~ U[-W/2, W/2] models Anderson localization.
    !! The uniform distribution ensures ⟨V⟩ = 0 (no systematic bias).
    !! For large W/t, wavefunctions become exponentially localized (Anderson insulator).
    subroutine test_random_uniform_mean()
        use fortuno_serial, only: check => serial_check
        use potential_random, only: potential_random_uniform
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10000, seed = 123
        real(dp) :: V(L), W, mean_V
        integer :: ierr

        W = 2.0_dp

        call potential_random_uniform(W, L, seed, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")

        ! Check mean is close to zero (within statistical fluctuations)
        mean_V = sum(V) / real(L, dp)
        call check(abs(mean_V) < 0.01_dp, "Mean should be close to zero")

        ! Check all values are within [-W/2, W/2]
        call check(all(V >= -W/2.0_dp .and. V <= W/2.0_dp), &
                   "All values should be in [-W/2, W/2]")
    end subroutine test_random_uniform_mean

    !> Test random Gaussian potential has zero mean
    !!
    !! Physics: Gaussian disorder V(i) ~ N(0, σ²) is more realistic than uniform
    !! disorder for many systems (thermal fluctuations, quantum fluctuations).
    !! The central limit theorem ensures ⟨V⟩ = 0 for large L.
    subroutine test_random_gaussian_mean()
        use fortuno_serial, only: check => serial_check
        use potential_random, only: potential_random_gaussian
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10000, seed = 456
        real(dp) :: V(L), sigma, mean_V, std_V
        integer :: ierr

        sigma = 1.0_dp

        call potential_random_gaussian(sigma, L, seed, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")

        ! Check mean is close to zero (within statistical fluctuations: σ/√N)
        mean_V = sum(V) / real(L, dp)
        call check(abs(mean_V) < 0.02_dp, "Mean should be close to zero")

        ! Check standard deviation is close to sigma
        std_V = sqrt(sum((V - mean_V)**2) / real(L, dp))
        call check(abs(std_V - sigma) < 0.05_dp, "Std dev should be close to sigma")
    end subroutine test_random_gaussian_mean

    !> Test single barrier has correct width
    !!
    !! Physics: A rectangular barrier V(i) = V_bar for i ∈ [i_start, i_end]
    !! models quantum tunneling. The barrier width w = i_end - i_start + 1
    !! determines the tunneling probability T ~ exp(-2κw) where κ² ~ V_bar - E.
    subroutine test_barrier_single_width()
        use fortuno_serial, only: check => serial_check
        use potential_barrier, only: potential_barrier_single
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 20, i_start = 8, i_end = 12
        real(dp) :: V(L), V_bar
        integer :: ierr, i, width

        V_bar = 5.0_dp
        width = i_end - i_start + 1

        call potential_barrier_single(V_bar, i_start, i_end, L, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")
        call check(width == 5, "Barrier width should be 5")

        ! Check barrier region
        do i = i_start, i_end
            call check(abs(V(i) - V_bar) < TOL, "Barrier sites should have V_bar")
        end do

        ! Check outside barrier
        do i = 1, i_start - 1
            call check(abs(V(i)) < TOL, "Sites before barrier should be zero")
        end do
        do i = i_end + 1, L
            call check(abs(V(i)) < TOL, "Sites after barrier should be zero")
        end do
    end subroutine test_barrier_single_width

    !> Test single barrier bounds checking
    !!
    !! Physics: The barrier must be within the physical lattice [1, L].
    !! Invalid bounds are unphysical and should be rejected.
    subroutine test_barrier_single_bounds()
        use fortuno_serial, only: check => serial_check
        use potential_barrier, only: potential_barrier_single
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_OUT_OF_BOUNDS

        integer, parameter :: L = 10
        real(dp) :: V(L)
        integer :: ierr

        ! Test i_start = 0 (invalid)
        call potential_barrier_single(1.0_dp, 0, 5, L, V, ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "i_start = 0 should fail")

        ! Test i_end > L (invalid)
        call potential_barrier_single(1.0_dp, 5, L + 1, L, V, ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "i_end > L should fail")

        ! Test i_end < i_start (invalid)
        call potential_barrier_single(1.0_dp, 8, 5, L, V, ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "i_end < i_start should fail")
    end subroutine test_barrier_single_bounds

    !> Test double barrier quantum well separation
    !!
    !! Physics: A double barrier creates a quantum well between the two barriers.
    !! The well width d = i2_start - i1_end - 1 determines the energy levels:
    !! E_n ~ n²/d² (particle in a box). Resonant tunneling occurs when the
    !! incident energy matches a well energy level.
    subroutine test_barrier_double_well_separation()
        use fortuno_serial, only: check => serial_check
        use potential_barrier, only: potential_barrier_double
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 30
        integer, parameter :: i1_start = 8, i1_end = 10
        integer, parameter :: i2_start = 18, i2_end = 20
        real(dp) :: V(L), V_bar
        integer :: ierr, well_width, i

        V_bar = 4.0_dp
        well_width = i2_start - i1_end - 1

        call potential_barrier_double(V_bar, i1_start, i1_end, i2_start, i2_end, L, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")
        call check(well_width == 7, "Well width should be 7")

        ! Check first barrier
        do i = i1_start, i1_end
            call check(abs(V(i) - V_bar) < TOL, "First barrier should have V_bar")
        end do

        ! Check well region (should be zero)
        do i = i1_end + 1, i2_start - 1
            call check(abs(V(i)) < TOL, "Well region should be zero")
        end do

        ! Check second barrier
        do i = i2_start, i2_end
            call check(abs(V(i) - V_bar) < TOL, "Second barrier should have V_bar")
        end do
    end subroutine test_barrier_double_well_separation

    !> Test double barrier rejects overlapping barriers
    !!
    !! Physics: Overlapping barriers are ill-defined for a double barrier
    !! quantum well configuration. This would not create the characteristic
    !! Fabry-Pérot resonances expected from a true double barrier structure.
    subroutine test_barrier_double_no_overlap()
        use fortuno_serial, only: check => serial_check
        use potential_barrier, only: potential_barrier_double
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_OUT_OF_BOUNDS

        integer, parameter :: L = 30
        real(dp) :: V(L)
        integer :: ierr

        ! i2_start <= i1_end (overlapping)
        call potential_barrier_double(1.0_dp, 5, 10, 10, 15, L, V, ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "Overlapping barriers should fail")

        call potential_barrier_double(1.0_dp, 5, 10, 8, 15, L, V, ierr)
        call check(ierr == ERROR_OUT_OF_BOUNDS, "Overlapping barriers should fail")
    end subroutine test_barrier_double_no_overlap

    !> Test factory creates uniform potential
    !!
    !! Physics: The factory pattern provides a unified interface for creating
    !! any potential type from a string identifier and parameter array.
    subroutine test_factory_uniform()
        use fortuno_serial, only: check => serial_check
        use potential_factory, only: create_potential
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 10
        real(dp) :: V(L), params(1)
        integer :: ierr, i

        params(1) = 3.0_dp

        call create_potential("uniform", params, L, -1, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")
        do i = 1, L
            call check(abs(V(i) - params(1)) < TOL, "Should match uniform potential")
        end do
    end subroutine test_factory_uniform

    !> Test factory creates harmonic potential
    !!
    !! Physics: Factory-created harmonic trap should be identical to direct call.
    subroutine test_factory_harmonic()
        use fortuno_serial, only: check => serial_check
        use potential_factory, only: create_potential
        use potential_harmonic, only: apply_potential_harmonic
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 11
        real(dp) :: V_factory(L), V_direct(L), params(1)
        integer :: ierr, ierr2, i

        params(1) = 0.2_dp

        call create_potential("harmonic", params, L, -1, V_factory, ierr)
        call apply_potential_harmonic(params(1), L, V_direct, ierr2)

        call check(ierr == ERROR_SUCCESS, "Should succeed")
        do i = 1, L
            call check(abs(V_factory(i) - V_direct(i)) < TOL, &
                       "Factory and direct call should match")
        end do
    end subroutine test_factory_harmonic

    !> Test factory rejects invalid potential type
    !!
    !! Physics: Invalid potential types should be caught at runtime to prevent
    !! silent failures in simulation setup.
    subroutine test_factory_invalid_type()
        use fortuno_serial, only: check => serial_check
        use potential_factory, only: create_potential
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_INVALID_INPUT

        integer, parameter :: L = 10
        real(dp) :: V(L), params(1)
        integer :: ierr

        params(1) = 1.0_dp

        call create_potential("invalid_type", params, L, -1, V, ierr)

        call check(ierr == ERROR_INVALID_INPUT, "Invalid type should fail")
    end subroutine test_factory_invalid_type

end program test_potentials
