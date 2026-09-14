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
            test("quasiperiodic_golden_ratio", test_quasiperiodic_golden_ratio), &
            test("quasiperiodic_phase_shift", test_quasiperiodic_phase_shift), &
            test("quasiperiodic_critical_point", test_quasiperiodic_critical_point), &
            test("quasiperiodic_localization", test_quasiperiodic_localization), &
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
    !! Physics: Random disorder with V(i) ~ U[-W, W] models Anderson localization.
    !! The uniform distribution ensures ⟨V⟩ = 0 (no systematic bias).
    !! For large W/t, wavefunctions become exponentially localized (Anderson insulator).
    !!
    !! Regression guard for Bug #2 (see CLAUDE.md): the disorder amplitude must
    !! follow the C++ reference, V = W*(2*rand - 1), which spans the FULL width
    !! [-W, W]. The old Fortran formula V = W*(rand - 0.5) only spanned
    !! [-W/2, W/2], i.e. half the intended disorder strength. Checking the upper
    !! bound alone cannot detect that regression (values in [-W/2, W/2] also
    !! satisfy it), so this test additionally asserts that the sample actually
    !! reaches close to both ends of [-W, W].
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
        ! Tolerance is statistical, not exact: for V ~ U[-W, W] the standard
        ! error of the mean is W/sqrt(3*L) ~= 0.01155 for W = 2, L = 10000.
        ! A 0.01 bound would sit below one standard error and could fail on a
        ! different compiler's random_number sequence; 0.06 is ~5 standard
        ! errors. Detecting a halved amplitude (Bug #2) is the job of the
        ! extreme-value checks below, not of this bound.
        call check(abs(mean_V) < 0.06_dp, "Mean should be close to zero")

        ! Check all values are within [-W, W] (C++ formula V = W*(2*rand - 1))
        call check(all(V >= -W .and. V <= W), &
                   "All values should be in [-W, W]")

        ! Check the sample spans the full width: with L = 10000 draws the
        ! extremes must land within 1% of the bounds. This fails if the
        ! amplitude is halved (max would only reach W/2).
        call check(maxval(V) > 0.99_dp * W, &
                   "Maximum should approach +W (full disorder amplitude)")
        call check(minval(V) < -0.99_dp * W, &
                   "Minimum should approach -W (full disorder amplitude)")
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

    !> Test double barrier quantum well geometry (barrier - well - barrier)
    !!
    !! Physics: A double barrier creates a quantum well between the two barriers.
    !! The well width determines the quasi-bound level spacing (E_n ~ n²/L_well²)
    !! and resonant tunneling occurs when the incident energy matches such a level.
    !!
    !! Geometry (mirrors the C++ double_barrier, lsda_potential.cc:166-194):
    !!   L = 30 (even)  =>  x0 = 30/2 + 0.5 = 15.5
    !!   L_well = 7     =>  x_1 = 12.0, x1 = 19.0
    !!   L_bar  = 3     =>  x_2 =  9.0, x2 = 22.0
    !! Because the barrier branches are tested first and use strict inequalities
    !! widened by SMALL = 1e-10, sites sitting exactly on the nominal well edges
    !! (i = 12 and i = 19) are classified as BARRIER, not well. Hence:
    !!   barrier sites: i = 9..12 and i = 19..22  -> 8 sites with V = V_bar
    !!   well sites:    i = 13..18                -> 6 sites with V = V_well
    !!   zero sites:    the remaining 16 sites
    subroutine test_barrier_double_well_separation()
        use fortuno_serial, only: check => serial_check
        use potential_barrier, only: potential_barrier_double
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 30
        real(dp), parameter :: V_bar = 4.0_dp, L_bar = 3.0_dp
        real(dp), parameter :: V_well = -3.0_dp, L_well = 7.0_dp
        real(dp) :: V(L)
        integer :: ierr, i, n_bar, n_well, n_zero

        call potential_barrier_double(V_bar, L_bar, V_well, L_well, L, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")

        n_bar = 0
        n_well = 0
        n_zero = 0
        do i = 1, L
            if (abs(V(i) - V_bar) < TOL) then
                n_bar = n_bar + 1
            else if (abs(V(i) - V_well) < TOL) then
                n_well = n_well + 1
            else if (abs(V(i)) < TOL) then
                n_zero = n_zero + 1
            end if
        end do

        call check(n_bar == 8, "Should have 8 barrier sites (i=9..12 and i=19..22)")
        call check(n_well == 6, "Should have 6 well sites (i=13..18)")
        call check(n_zero == L - 14, "Remaining sites should be zero")

        ! Explicit region boundaries
        call check(abs(V(8)) < TOL, "Site 8 is outside the left barrier")
        call check(abs(V(9) - V_bar) < TOL, "Site 9 starts the left barrier")
        call check(abs(V(12) - V_bar) < TOL, "Site 12 is barrier (well edge falls on it)")
        call check(abs(V(13) - V_well) < TOL, "Site 13 starts the well")
        call check(abs(V(18) - V_well) < TOL, "Site 18 ends the well")
        call check(abs(V(19) - V_bar) < TOL, "Site 19 is barrier (well edge falls on it)")
        call check(abs(V(22) - V_bar) < TOL, "Site 22 ends the right barrier")
        call check(abs(V(23)) < TOL, "Site 23 is outside the right barrier")
    end subroutine test_barrier_double_well_separation

    !> Test that the three double-barrier regions never overlap
    !!
    !! Physics: barrier, well and field-free regions must partition the lattice;
    !! a site cannot be simultaneously barrier and well, otherwise the
    !! Fabry-Pérot resonance structure would be ill-defined.
    !!
    !! Geometry chosen so that no boundary lands on an integer site:
    !!   L = 30 => x0 = 15.5; L_well = 4 => x_1 = 13.5, x1 = 17.5;
    !!   L_bar = 5 => x_2 = 8.5, x2 = 22.5
    !!   barriers: i = 9..13 and i = 18..22 (5 + 5 sites)
    !!   well:     i = 14..17 (4 sites)
    !!   zero:     16 sites
    !!
    !! The second part checks the oversized case L_well + 2*L_bar > L: the C++
    !! original performs no validation and neither does the Fortran port, so the
    !! call must still succeed, simply clipping the regions to the lattice.
    subroutine test_barrier_double_no_overlap()
        use fortuno_serial, only: check => serial_check
        use potential_barrier, only: potential_barrier_double
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 30
        real(dp), parameter :: V_bar = 1.0_dp, V_well = -2.0_dp
        real(dp) :: V(L)
        real(dp) :: V_small(10)
        integer :: ierr, i, n_bar, n_well, n_zero

        call potential_barrier_double(V_bar, 5.0_dp, V_well, 4.0_dp, L, V, ierr)
        call check(ierr == ERROR_SUCCESS, "Non-overlapping geometry should succeed")

        n_bar = 0
        n_well = 0
        n_zero = 0
        do i = 1, L
            if (abs(V(i) - V_bar) < TOL) then
                n_bar = n_bar + 1
            else if (abs(V(i) - V_well) < TOL) then
                n_well = n_well + 1
            else if (abs(V(i)) < TOL) then
                n_zero = n_zero + 1
            end if
        end do

        ! Every site belongs to exactly one region
        call check(n_bar + n_well + n_zero == L, "Regions must partition the lattice")
        call check(n_bar == 10, "Should have 10 barrier sites (i=9..13 and i=18..22)")
        call check(n_well == 4, "Should have 4 well sites (i=14..17)")
        call check(n_zero == 16, "Should have 16 field-free sites")

        ! Oversized geometry: L_well + 2*L_bar = 22 > L = 10
        ! L = 10 => x0 = 5.5; x_1 = 2.5, x1 = 8.5; x_2 = -5.5, x2 = 16.5
        ! => barriers at i = 1,2 and i = 9,10; well at i = 3..8; no field-free site
        call potential_barrier_double(V_bar, 8.0_dp, V_well, 6.0_dp, 10, V_small, ierr)
        call check(ierr == ERROR_SUCCESS, "Oversized geometry is not rejected (matches C++)")
        call check(abs(V_small(1) - V_bar) < TOL, "Site 1 is clipped left barrier")
        call check(abs(V_small(2) - V_bar) < TOL, "Site 2 is clipped left barrier")
        call check(abs(V_small(3) - V_well) < TOL, "Site 3 starts the well")
        call check(abs(V_small(8) - V_well) < TOL, "Site 8 ends the well")
        call check(abs(V_small(9) - V_bar) < TOL, "Site 9 is clipped right barrier")
        call check(abs(V_small(10) - V_bar) < TOL, "Site 10 is clipped right barrier")
    end subroutine test_barrier_double_no_overlap

    !> Test quasiperiodic potential with golden ratio
    !!
    !! Physics: The Aubry-André-Harper (AAH) model with β = golden ratio exhibits
    !! maximum incommensurability, preventing periodic repetition. For λ < 2, all
    !! states are extended (delocalized). This test verifies that the potential
    !! is continuous and bounded by [-λ, λ], characteristic of the cosine modulation.
    subroutine test_quasiperiodic_golden_ratio()
        use fortuno_serial, only: check => serial_check
        use potential_quasiperiodic, only: apply_potential_quasiperiodic
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 100
        real(dp) :: V(L)
        real(dp), parameter :: lambda = 1.5_dp
        real(dp), parameter :: beta = 0.5_dp * (sqrt(5.0_dp) - 1.0_dp)  ! Golden ratio
        real(dp), parameter :: phi = 0.0_dp
        integer :: ierr
        real(dp) :: v_min, v_max

        call apply_potential_quasiperiodic(lambda, beta, phi, L, V, ierr)

        call check(ierr == ERROR_SUCCESS, "Should succeed")

        v_min = minval(V)
        v_max = maxval(V)
        call check(v_min >= -lambda - TOL, "V_min should be >= -lambda")
        call check(v_max <= lambda + TOL, "V_max should be <= lambda")

        call check(abs(v_max - v_min) > 0.1_dp, "Potential should vary")
    end subroutine test_quasiperiodic_golden_ratio

    !> Test quasiperiodic potential with phase shift
    !!
    !! Physics: The phase φ in V(i) = λcos(2πβi + φ) shifts the potential pattern
    !! along the lattice without changing the physics (gauge freedom). A phase shift
    !! of π inverts the potential: V_new(i) = -V_old(i). This test verifies that
    !! φ = 0 and φ = π produce opposite potentials, as expected.
    subroutine test_quasiperiodic_phase_shift()
        use fortuno_serial, only: check => serial_check
        use potential_quasiperiodic, only: apply_potential_quasiperiodic
        use lsda_constants, only: dp, PI
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 50
        real(dp) :: V_phi0(L), V_phi_pi(L)
        real(dp), parameter :: lambda = 2.0_dp
        real(dp), parameter :: beta = 0.5_dp * (sqrt(5.0_dp) - 1.0_dp)
        integer :: ierr, i

        ! Calculate with phi = 0
        call apply_potential_quasiperiodic(lambda, beta, 0.0_dp, L, V_phi0, ierr)
        call check(ierr == ERROR_SUCCESS, "phi=0 should succeed")

        ! Calculate with phi = π
        call apply_potential_quasiperiodic(lambda, beta, PI, L, V_phi_pi, ierr)
        call check(ierr == ERROR_SUCCESS, "phi=π should succeed")

        ! Check that V(φ=π) ≈ -V(φ=0)
        do i = 1, L
            call check(abs(V_phi_pi(i) + V_phi0(i)) < TOL, &
                       "Phase shift of π should invert potential")
        end do
    end subroutine test_quasiperiodic_phase_shift

    !> Test quasiperiodic potential at critical point
    !!
    !! Physics: The AAH model undergoes a localization transition at λ_c = 2
    !! (for β = golden ratio). At this critical point, the system exhibits:
    !! - Critical wavefunctions (multifractal, neither extended nor exponentially localized)
    !! - Subdiffusive transport with anomalous diffusion exponent
    !! - Self-similar fractal structure in the density of states
    !! This test verifies that λ = 2 produces a potential with correct amplitude.
    subroutine test_quasiperiodic_critical_point()
        use fortuno_serial, only: check => serial_check
        use potential_quasiperiodic, only: apply_potential_quasiperiodic
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 200
        real(dp) :: V(L)
        real(dp), parameter :: lambda_c = 2.0_dp
        real(dp), parameter :: beta = 0.5_dp * (sqrt(5.0_dp) - 1.0_dp)
        real(dp), parameter :: phi = 0.0_dp
        integer :: ierr
        real(dp) :: v_min, v_max

        call apply_potential_quasiperiodic(lambda_c, beta, phi, L, V, ierr)
        call check(ierr == ERROR_SUCCESS, "Should succeed at critical point")

        v_min = minval(V)
        v_max = maxval(V)

        call check(v_min < -1.8_dp, "Min should be close to -2")
        call check(v_max > 1.8_dp, "Max should be close to 2")
        call check(abs(v_max + v_min) < 0.1_dp, "Potential should be symmetric around 0")
    end subroutine test_quasiperiodic_critical_point

    !> Test quasiperiodic potential in localized regime
    !!
    !! Physics: For λ > 2 (with β = golden ratio), the AAH model is in the localized
    !! phase where all eigenstates are exponentially localized. The localization length
    !! ξ decreases as λ increases beyond 2. In this regime:
    !! - Transport is suppressed (DC conductivity → 0)
    !! - Density of states shows pure point spectrum
    !! - All wavefunctions decay exponentially: |ψ(x)| ~ exp(-|x|/ξ)
    !! This test verifies that strong potentials (λ >> 2) are correctly calculated.
    subroutine test_quasiperiodic_localization()
        use fortuno_serial, only: check => serial_check
        use potential_quasiperiodic, only: apply_potential_quasiperiodic
        use lsda_constants, only: dp
        use lsda_errors, only: ERROR_SUCCESS

        integer, parameter :: L = 150
        real(dp) :: V(L)
        real(dp), parameter :: lambda = 5.0_dp
        real(dp), parameter :: beta = 0.5_dp * (sqrt(5.0_dp) - 1.0_dp)
        real(dp), parameter :: phi = 0.0_dp
        integer :: ierr
        real(dp) :: v_min, v_max, v_range

        call apply_potential_quasiperiodic(lambda, beta, phi, L, V, ierr)
        call check(ierr == ERROR_SUCCESS, "Should succeed in localized regime")

        v_min = minval(V)
        v_max = maxval(V)
        v_range = v_max - v_min

        call check(v_range > 8.0_dp, "Potential range should be large for strong λ")
        call check(v_min >= -lambda - TOL, "V_min should be >= -lambda")
        call check(v_max <= lambda + TOL, "V_max should be <= lambda")

        ! Verify approximate bounds are reached (within 90%)
        call check(v_min < -0.9_dp * lambda, "Should reach close to -lambda")
        call check(v_max > 0.9_dp * lambda, "Should reach close to +lambda")
    end subroutine test_quasiperiodic_localization

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
