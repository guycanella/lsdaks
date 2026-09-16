!> Unit tests for density_calculator module
program test_density_calculator
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp, PI
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_density_tests())

contains

    function get_density_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("single_electron_density", test_single_electron_density), &
            test("half_filling_unpolarized", test_half_filling_unpolarized), &
            test("particle_number_conservation", test_particle_number_conservation), &
            test("density_positivity", test_density_positivity), &
            test("physical_bounds", test_physical_bounds), &
            test("density_from_harmonic_trap", test_density_from_harmonic_trap), &
            test("occupations_closed_shell", test_occupations_closed_shell), &
            test("occupations_open_shell", test_occupations_open_shell), &
            test("occupations_consecutive_window", test_occupations_consecutive_window), &
            test("occupations_spectrum_edges", test_occupations_spectrum_edges), &
            test("occupations_invalid_input", test_occupations_invalid_input), &
            test("open_shell_density_is_uniform", test_open_shell_density_is_uniform), &
            test("open_shell_density_complex", test_open_shell_density_complex), &
            test("empty_channel_density_is_zero", test_empty_channel_density_is_zero) &
        ])
    end function get_density_tests

    !> Closed Fermi shell: occupations are plain integers
    !!
    !! With L = 8 and periodic BC the U = 0 spectrum is -2cos(2 pi m / 8), i.e.
    !! -2, -sqrt(2), -sqrt(2), 0, 0, sqrt(2), sqrt(2), 2. N = 5 puts the Fermi
    !! level in the middle of the (4, 5) degenerate pair, but N = 3 closes the
    !! (2, 3) pair exactly: there the partial fill must come out as 1.0 and the
    !! occupations must be indistinguishable from integer filling.
    subroutine test_occupations_closed_shell()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_occupations
        use lsda_constants, only: DEG_TOL
        use lsda_errors, only: ERROR_SUCCESS
        integer, parameter :: n_levels = 8
        real(dp) :: eigvals(n_levels), occ(n_levels)
        integer :: ierr

        eigvals = [-2.0_dp, -sqrt(2.0_dp), -sqrt(2.0_dp), 0.0_dp, &
                    0.0_dp, sqrt(2.0_dp), sqrt(2.0_dp), 2.0_dp]

        call compute_occupations(eigvals, 3, DEG_TOL, occ, ierr)

        call check(ierr == ERROR_SUCCESS, "Closed shell: computation should succeed")
        call check(all(abs(occ(1:3) - 1.0_dp) < TOL), &
                   "Closed shell: the three lowest levels must be fully occupied")
        call check(all(abs(occ(4:n_levels)) < TOL), &
                   "Closed shell: levels above the Fermi level must be empty")
        call check(abs(sum(occ) - 3.0_dp) < TOL, &
                   "Closed shell: sum(occ) must equal N")
    end subroutine test_occupations_closed_shell

    !> Open (partially filled) degenerate Fermi shell
    !!
    !! Same L = 8 periodic spectrum, N = 4: levels 4 and 5 are degenerate and
    !! together hold a single electron, so each must carry weight 1/2. This is
    !! the case that used to be filled with integer occupation, making the
    !! density depend on the arbitrary LAPACK basis inside the shell.
    !!
    !! (N = 5 in the same spectrum is a CLOSED shell: the same window then holds
    !! two electrons in two states, partial_fill = 1.)
    subroutine test_occupations_open_shell()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_occupations
        use lsda_constants, only: DEG_TOL
        use lsda_errors, only: ERROR_SUCCESS
        integer, parameter :: n_levels = 8
        real(dp) :: eigvals(n_levels), occ(n_levels)
        integer :: ierr

        eigvals = [-2.0_dp, -sqrt(2.0_dp), -sqrt(2.0_dp), 0.0_dp, &
                    0.0_dp, sqrt(2.0_dp), sqrt(2.0_dp), 2.0_dp]

        call compute_occupations(eigvals, 4, DEG_TOL, occ, ierr)

        call check(ierr == ERROR_SUCCESS, "Open shell: computation should succeed")
        call check(all(abs(occ(1:3) - 1.0_dp) < TOL), &
                   "Open shell: levels below the shell must be fully occupied")
        call check(abs(occ(4) - 0.5_dp) < TOL .and. abs(occ(5) - 0.5_dp) < TOL, &
                   "Open shell: the two degenerate levels must share the electron equally")
        call check(all(abs(occ(6:n_levels)) < TOL), &
                   "Open shell: levels above the shell must be empty")
        call check(abs(sum(occ) - 4.0_dp) < TOL, &
                   "Open shell: sum(occ) must equal N with no renormalisation")

        ! The closed-shell sibling: N = 5 fills the same window completely.
        call compute_occupations(eigvals, 5, DEG_TOL, occ, ierr)
        call check(ierr == ERROR_SUCCESS, "Open shell: N = 5 should succeed")
        call check(all(abs(occ(1:5) - 1.0_dp) < TOL) .and. all(abs(occ(6:n_levels)) < TOL), &
                   "Open shell: N = 5 closes the shell, partial_fill = 1")
    end subroutine test_occupations_open_shell

    !> The shell window walks over CONSECUTIVE neighbours, as the C++ does
    !!
    !! The C++ `update_degen` widens the window while |e(j) - e(j-1)| < tol, not
    !! while |e(j) - e(N)| < tol. The two criteria disagree on a cascade of
    !! levels each within tol of the previous one but far from e(N): here four
    !! levels spaced by tol/2 span 1.5*tol, so the consecutive criterion takes
    !! all four while the distance-to-e(N) criterion would take only two.
    subroutine test_occupations_consecutive_window()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_occupations
        use lsda_constants, only: DEG_TOL
        use lsda_errors, only: ERROR_SUCCESS
        integer, parameter :: n_levels = 6
        real(dp) :: eigvals(n_levels), occ(n_levels)
        real(dp) :: step
        integer :: ierr

        step = 0.5_dp * DEG_TOL
        eigvals = [-1.0_dp, 0.0_dp, step, 2.0_dp * step, 3.0_dp * step, 1.0_dp]

        ! N = 3 sits inside the cascade (levels 2..5).
        call compute_occupations(eigvals, 3, DEG_TOL, occ, ierr)

        call check(ierr == ERROR_SUCCESS, "Cascade: computation should succeed")
        call check(abs(occ(1) - 1.0_dp) < TOL, &
                   "Cascade: the isolated bottom level stays fully occupied")
        call check(all(abs(occ(2:5) - 0.5_dp) < TOL), &
                   "Cascade: all four chained levels share the 2 remaining electrons (2/4 each)")
        call check(abs(occ(6)) < TOL, "Cascade: the isolated top level stays empty")
        call check(abs(sum(occ) - 3.0_dp) < TOL, "Cascade: sum(occ) must equal N")
    end subroutine test_occupations_consecutive_window

    !> The window search must not read outside the spectrum
    !!
    !! The C++ reads erg[level_min - 1] and erg[level_max + 1] unguarded and
    !! relies on asking for Ne + 5 eigenvalues to make that harmless. Here the
    !! whole spectrum is degenerate and N = n_levels, so both ends of the window
    !! hit the array bounds; the routine must clamp instead of running off.
    subroutine test_occupations_spectrum_edges()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_occupations
        use lsda_constants, only: DEG_TOL
        use lsda_errors, only: ERROR_SUCCESS
        integer, parameter :: n_levels = 4
        real(dp) :: eigvals(n_levels), occ(n_levels)
        integer :: ierr

        eigvals = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]

        call compute_occupations(eigvals, n_levels, DEG_TOL, occ, ierr)

        call check(ierr == ERROR_SUCCESS, "Edges: full band should succeed")
        call check(all(abs(occ - 1.0_dp) < TOL), &
                   "Edges: a completely filled degenerate band has occ = 1 everywhere")

        ! N = 1 in the same fully degenerate band: one electron over four states.
        call compute_occupations(eigvals, 1, DEG_TOL, occ, ierr)
        call check(ierr == ERROR_SUCCESS, "Edges: N = 1 should succeed")
        call check(all(abs(occ - 0.25_dp) < TOL), &
                   "Edges: one electron in a 4-fold degenerate band gives occ = 1/4")
    end subroutine test_occupations_spectrum_edges

    !> Input validation of compute_occupations
    !!
    !! n_elec = 0 is a legal, fully polarised channel and must return zeros
    !! WITHOUT touching eigvals(0); only n_elec < 0 (and n_elec > n_levels) are
    !! rejected.
    subroutine test_occupations_invalid_input()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_occupations
        use lsda_constants, only: DEG_TOL
        use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT, ERROR_SIZE_MISMATCH
        integer, parameter :: n_levels = 4
        real(dp) :: eigvals(n_levels), occ(n_levels), occ_short(n_levels - 1)
        integer :: ierr

        eigvals = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]

        call compute_occupations(eigvals, 0, DEG_TOL, occ, ierr)
        call check(ierr == ERROR_SUCCESS, "Validation: n_elec = 0 must be accepted")
        call check(all(abs(occ) < TOL), "Validation: n_elec = 0 gives all-zero occupations")

        call compute_occupations(eigvals, -1, DEG_TOL, occ, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Validation: n_elec < 0 must be rejected")

        call compute_occupations(eigvals, n_levels + 1, DEG_TOL, occ, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Validation: n_elec > n_levels must be rejected")

        call compute_occupations(eigvals, 2, -1.0_dp, occ, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "Validation: negative deg_tol must be rejected")

        call compute_occupations(eigvals, 2, DEG_TOL, occ_short, ierr)
        call check(ierr == ERROR_SIZE_MISMATCH, &
                   "Validation: occ must have the same length as eigvals")
    end subroutine test_occupations_invalid_input

    !> REGRESSION (T7): open degenerate shell must give a uniform density
    !!
    !! L = 8, N_sigma = 4, periodic BC, U = 0, zero potential. The spectrum is
    !! -2, -sqrt(2), -sqrt(2), 0, 0, sqrt(2), sqrt(2), 2: level 4 is degenerate
    !! with level 5, so the Fermi shell is open and holds one electron. With the
    !! translational symmetry unbroken the density must be exactly N/L = 0.5 at
    !! every site.
    !!
    !! Without fractional occupation this test FAILS: LAPACK returns the real
    !! cos/sin combinations of the k = +-pi/2 pair, occupying only the cos one
    !! with weight 1 produces a period-2 density wave of amplitude 0.25, and the
    !! second assertion below (which checks integer filling is NOT uniform) is
    !! there to prove the test has teeth rather than passing for free.
    subroutine test_open_shell_density_is_uniform()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_density_spin, compute_occupations
        use hamiltonian_builder, only: build_hamiltonian
        use boundary_conditions, only: BC_PERIODIC
        use lapack_wrapper, only: diagonalize_symmetric_real
        use lsda_constants, only: DEG_TOL
        use lsda_errors, only: ERROR_SUCCESS
        integer, parameter :: L = 8, n_elec = 4
        real(dp) :: V_zero(L), H(L, L), eigvals(L), eigvecs(L, L)
        real(dp) :: occ(L), density(L), density_int(L)
        integer :: ierr

        V_zero = 0.0_dp

        call build_hamiltonian(L, V_zero, V_zero, BC_PERIODIC, 0.0_dp, H, ierr)
        call check(ierr == ERROR_SUCCESS, "Open shell density: Hamiltonian build should succeed")

        call diagonalize_symmetric_real(H, L, eigvals, eigvecs, ierr)
        call check(ierr == ERROR_SUCCESS, "Open shell density: diagonalization should succeed")

        ! Precondition: level 4 and level 5 really are degenerate.
        call check(abs(eigvals(5) - eigvals(4)) < DEG_TOL, &
                   "Open shell density: levels 4 and 5 must be degenerate (precondition)")

        call compute_occupations(eigvals, n_elec, DEG_TOL, occ, ierr)
        call check(ierr == ERROR_SUCCESS, "Open shell density: occupations should succeed")

        call compute_density_spin(eigvecs, L, occ, density, ierr)
        call check(ierr == ERROR_SUCCESS, "Open shell density: density should succeed")

        call check(maxval(abs(density - 0.5_dp)) < 1.0e-12_dp, &
                   "Open shell density: n(i) must be exactly N/L = 0.5 at every site")
        call check(abs(sum(density) - real(n_elec, dp)) < 1.0e-12_dp, &
                   "Open shell density: particle number must be conserved")

        ! Teeth: integer occupation of the same eigenvectors is NOT uniform.
        call compute_density_spin(eigvecs, L, n_elec, density_int, ierr)
        call check(ierr == ERROR_SUCCESS, "Open shell density: integer filling should succeed")
        call check(maxval(abs(density_int - 0.5_dp)) > 1.0e-3_dp, &
                   "Open shell density: integer filling must break translational symmetry")
    end subroutine test_open_shell_density_is_uniform

    !> REGRESSION (T7): same open-shell uniformity on the COMPLEX path
    !!
    !! The twisted-BC branch of the SCF cycle needs fractional occupation too
    !! (the C++ has the same logic in `update_dens_twist`). Here the degenerate
    !! pair is built explicitly as the two plane waves k = +-pi/2 of an L = 8
    !! ring, mixed by an arbitrary unitary rotation inside the subspace to mimic
    !! what a diagonalizer may return.
    subroutine test_open_shell_density_complex()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_density_spin
        use lsda_errors, only: ERROR_SUCCESS
        integer, parameter :: L = 8, n_elec = 4
        complex(dp) :: eigvecs(L, L), plus(L), minus(L)
        real(dp) :: occ(L), density(L), k
        integer :: i, ierr

        eigvecs = (0.0_dp, 0.0_dp)

        ! Levels 1..3: plane waves k = 0, and the k = +-pi/4 pair. Only their
        ! modulus matters here (all are uniform), so any normalised plane wave
        ! does; what matters is that levels 4 and 5 are an arbitrary rotation of
        ! the degenerate k = +-pi/2 pair.
        do i = 1, L
            eigvecs(i, 1) = cmplx(1.0_dp / sqrt(real(L, dp)), 0.0_dp, kind=dp)
        end do
        do i = 1, L
            k = 2.0_dp * PI / real(L, dp)
            eigvecs(i, 2) = exp(cmplx(0.0_dp, k * real(i, dp), kind=dp)) / sqrt(real(L, dp))
            eigvecs(i, 3) = exp(cmplx(0.0_dp, -k * real(i, dp), kind=dp)) / sqrt(real(L, dp))
        end do

        k = PI / 2.0_dp
        do i = 1, L
            plus(i) = exp(cmplx(0.0_dp, k * real(i, dp), kind=dp)) / sqrt(real(L, dp))
            minus(i) = exp(cmplx(0.0_dp, -k * real(i, dp), kind=dp)) / sqrt(real(L, dp))
        end do

        ! Arbitrary unitary rotation inside the degenerate subspace: the real
        ! cos/sin combinations, which individually are NOT uniform.
        eigvecs(:, 4) = (plus + minus) / sqrt(2.0_dp)
        eigvecs(:, 5) = (plus - minus) / cmplx(0.0_dp, sqrt(2.0_dp), kind=dp)

        occ = 0.0_dp
        occ(1:3) = 1.0_dp
        occ(4:5) = 0.5_dp

        call compute_density_spin(eigvecs, L, occ, density, ierr)
        call check(ierr == ERROR_SUCCESS, "Complex open shell: density should succeed")
        call check(maxval(abs(density - 0.5_dp)) < 1.0e-12_dp, &
                   "Complex open shell: n(i) must be exactly 0.5 at every site")
        call check(abs(sum(density) - real(n_elec, dp)) < 1.0e-12_dp, &
                   "Complex open shell: particle number must be conserved")

        ! Teeth: filling only the cos combination with weight 1 is not uniform.
        occ(4) = 1.0_dp
        occ(5) = 0.0_dp
        call compute_density_spin(eigvecs, L, occ, density, ierr)
        call check(ierr == ERROR_SUCCESS, "Complex open shell: integer filling should succeed")
        call check(maxval(abs(density - 0.5_dp)) > 1.0e-3_dp, &
                   "Complex open shell: integer filling must break translational symmetry")
    end subroutine test_open_shell_density_complex

    !> REGRESSION (T7): a fully polarised channel (N_sigma = 0) is legal
    !!
    !! validate_kohn_sham_cycle_inputs already accepts N_sigma = 0, but the
    !! density routines used to reject n_elec <= 0 with ERROR_INVALID_INPUT, so
    !! a fully polarised run died downstream of a successful validation. Both
    !! the integer and the occupation-array entry points must now return a zero
    !! density instead.
    subroutine test_empty_channel_density_is_zero()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_density_spin
        use lsda_errors, only: ERROR_SUCCESS
        integer, parameter :: L = 6
        real(dp) :: eigvecs(L, L), density(L), occ(L)
        complex(dp) :: eigvecs_c(L, L)
        integer :: i, ierr

        eigvecs = 0.0_dp
        do i = 1, L
            eigvecs(i, i) = 1.0_dp
        end do
        eigvecs_c = cmplx(eigvecs, 0.0_dp, kind=dp)
        occ = 0.0_dp

        call compute_density_spin(eigvecs, L, 0, density, ierr)
        call check(ierr == ERROR_SUCCESS, "Empty channel: n_elec = 0 must be accepted (real)")
        call check(all(abs(density) < TOL), "Empty channel: density must be zero (real)")

        call compute_density_spin(eigvecs_c, L, 0, density, ierr)
        call check(ierr == ERROR_SUCCESS, "Empty channel: n_elec = 0 must be accepted (complex)")
        call check(all(abs(density) < TOL), "Empty channel: density must be zero (complex)")

        call compute_density_spin(eigvecs, L, occ, density, ierr)
        call check(ierr == ERROR_SUCCESS, "Empty channel: zero occ must be accepted (real)")
        call check(all(abs(density) < TOL), "Empty channel: zero occ gives zero density (real)")

        call compute_density_spin(eigvecs_c, L, occ, density, ierr)
        call check(ierr == ERROR_SUCCESS, "Empty channel: zero occ must be accepted (complex)")
        call check(all(abs(density) < TOL), "Empty channel: zero occ gives zero density (complex)")
    end subroutine test_empty_channel_density_is_zero

    !> Test density for single electron system
    !!
    !! Physics: For a single electron in a 1D box (open BC), the ground state
    !! is ψ₁(i) = √(2/(L+1)) sin(πi/(L+1)). The density should be n(i) = |ψ₁(i)|²,
    !! which is maximal at the center and zero at the edges.
    !! This tests that compute_density_spin correctly calculates |ψ|² at each site.
    subroutine test_single_electron_density()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_density_spin
        integer, parameter :: L = 5
        real(dp) :: eigvecs(L, L), density(L)
        real(dp) :: norm_factor
        integer :: i, ierr

        ! Create ground state wavefunction for particle in a box (OBC)
        ! ψ₁(i) = √(2/(L+1)) sin(πi/(L+1))
        norm_factor = sqrt(2.0_dp / real(L + 1, dp))
        do i = 1, L
            eigvecs(i, 1) = norm_factor * sin(PI * real(i, dp) / real(L + 1, dp))
        end do

        ! Fill remaining columns with dummy data (won't be used)
        eigvecs(:, 2:L) = 0.0_dp

        ! Compute density for 1 electron
        call compute_density_spin(eigvecs, L, 1, density, ierr)

        ! Check success
        call check(ierr == 0, "Single electron: computation should succeed")

        ! Check density is |ψ₁(i)|²
        do i = 1, L
            call check(abs(density(i) - eigvecs(i, 1)**2) < TOL, &
                       "Single electron: n(i) should equal |ψ₁(i)|²")
        end do

        ! Check normalization: Σ n(i) = 1
        call check(abs(sum(density) - 1.0_dp) < TOL, &
                   "Single electron: total density should equal 1")

        ! Check density is maximal near center (i=3 for L=5)
        call check(density(3) > density(1) .and. density(3) > density(5), &
                   "Single electron: density should be maximal at center")
    end subroutine test_single_electron_density

    !> Test half-filling with unpolarized system
    !!
    !! Physics: At half-filling (N = L) with periodic BC and U=0 (non-interacting),
    !! the eigenstates are plane waves ψₖ(i) = (1/√L) exp(ikxᵢ) with k = 2πn/L.
    !! For unpolarized case (N_up = N_down = L/2), the density should be uniform:
    !! n(i) = N/L = 1 at all sites due to translational symmetry.
    subroutine test_half_filling_unpolarized()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_density_spin, compute_total_density
        use hamiltonian_builder, only: build_hamiltonian_free
        use boundary_conditions, only: BC_PERIODIC
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 10
        real(dp) :: H(L, L), eigenvals(L), eigvecs(L, L)
        real(dp) :: density_up(L), density_dw(L), density_total(L)
        integer :: n_up, n_dw, ierr, i

        ! Half-filling: N_up = N_down = L/2
        n_up = L / 2
        n_dw = L / 2

        ! Build free Hamiltonian with periodic BC
        call build_hamiltonian_free(L, BC_PERIODIC, 0.0_dp, H, ierr)
        call check(ierr == 0, "Half-filling: Hamiltonian construction should succeed")

        ! Diagonalize to get eigenvectors
        call diagonalize_symmetric_real(H, L, eigenvals, eigvecs, ierr)
        call check(ierr == 0, "Half-filling: Diagonalization should succeed")

        ! Compute spin-resolved densities (same eigenvectors for unpolarized)
        call compute_density_spin(eigvecs, L, n_up, density_up, ierr)
        call check(ierr == 0, "Half-filling: density_up computation should succeed")

        call compute_density_spin(eigvecs, L, n_dw, density_dw, ierr)
        call check(ierr == 0, "Half-filling: density_dw computation should succeed")

        ! Compute total density
        call compute_total_density(density_up, density_dw, density_total, ierr)
        call check(ierr == 0, "Half-filling: total density computation should succeed")

        ! Check uniform density: n(i) = 1 for all i
        do i = 1, L
            call check(abs(density_total(i) - 1.0_dp) < TOL, &
                       "Half-filling: density should be uniform n(i) = 1")
        end do

        ! Check spin symmetry: n_up = n_down = 0.5
        do i = 1, L
            call check(abs(density_up(i) - 0.5_dp) < TOL, &
                       "Half-filling: n_up(i) should equal 0.5")
            call check(abs(density_dw(i) - 0.5_dp) < TOL, &
                       "Half-filling: n_down(i) should equal 0.5")
        end do
    end subroutine test_half_filling_unpolarized

    !> Test particle number conservation
    !!
    !! Physics: The total number of particles must be conserved:
    !! Σᵢ n_σ(i) = N_σ for each spin, and Σᵢ n(i) = N_total.
    !! This is a fundamental requirement from the normalization of wavefunctions.
    subroutine test_particle_number_conservation()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_density_spin, verify_particle_number
        use hamiltonian_builder, only: build_hamiltonian_free
        use boundary_conditions, only: BC_OPEN
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 8
        integer, parameter :: n_elec = 3
        real(dp) :: H(L, L), eigenvals(L), eigvecs(L, L), density(L)
        logical :: is_conserved
        integer :: ierr

        ! Build free Hamiltonian with open BC
        call build_hamiltonian_free(L, BC_OPEN, 0.0_dp, H, ierr)
        call check(ierr == 0, "Conservation: Hamiltonian construction should succeed")

        ! Diagonalize
        call diagonalize_symmetric_real(H, L, eigenvals, eigvecs, ierr)
        call check(ierr == 0, "Conservation: Diagonalization should succeed")

        ! Compute density for n_elec electrons
        call compute_density_spin(eigvecs, L, n_elec, density, ierr)
        call check(ierr == 0, "Conservation: density computation should succeed")

        ! Verify particle number conservation
        call verify_particle_number(density, L, n_elec, is_conserved, ierr)
        call check(ierr == 0, "Conservation: verification should succeed")
        call check(is_conserved, "Conservation: Σ n(i) should equal N")

        ! Explicit check: sum(density) ≈ n_elec
        call check(abs(sum(density) - real(n_elec, dp)) < TOL, &
                   "Conservation: explicit sum check should pass")
    end subroutine test_particle_number_conservation

    !> Test that density is always non-negative
    !!
    !! Physics: Electron density n(i) = Σⱼ |ψⱼ(i)|² is always non-negative
    !! since it's a sum of squared magnitudes. Negative densities are unphysical.
    subroutine test_density_positivity()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_density_spin
        use hamiltonian_builder, only: build_hamiltonian_free
        use boundary_conditions, only: BC_PERIODIC
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 12
        integer, parameter :: n_elec = 6
        real(dp) :: H(L, L), eigenvals(L), eigvecs(L, L), density(L)
        integer :: ierr, i

        ! Build and diagonalize free Hamiltonian
        call build_hamiltonian_free(L, BC_PERIODIC, 0.0_dp, H, ierr)
        call diagonalize_symmetric_real(H, L, eigenvals, eigvecs, ierr)

        ! Compute density
        call compute_density_spin(eigvecs, L, n_elec, density, ierr)
        call check(ierr == 0, "Positivity: computation should succeed")

        ! Check all densities are non-negative
        do i = 1, L
            call check(density(i) >= -TOL, &
                       "Positivity: n(i) should be non-negative")
        end do
    end subroutine test_density_positivity

    !> Test physical bounds on density
    !!
    !! Physics: For fermions with spin-1/2, each site can hold at most 2 electrons
    !! (one spin-up, one spin-down). Therefore, the physical bounds are:
    !! - 0 ≤ n_σ(i) ≤ 1 for each spin
    !! - 0 ≤ n(i) = n_up(i) + n_down(i) ≤ 2 for total density
    subroutine test_physical_bounds()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_density_spin, check_density_bounds
        use hamiltonian_builder, only: build_hamiltonian_free
        use boundary_conditions, only: BC_OPEN
        use lapack_wrapper, only: diagonalize_symmetric_real
        integer, parameter :: L = 10
        integer, parameter :: n_up = 3, n_dw = 3
        real(dp) :: H(L, L), eigenvals(L), eigvecs(L, L)
        real(dp) :: density_up(L), density_dw(L)
        logical :: all_valid
        integer :: ierr, i

        ! Build and diagonalize
        call build_hamiltonian_free(L, BC_OPEN, 0.0_dp, H, ierr)
        call diagonalize_symmetric_real(H, L, eigenvals, eigvecs, ierr)

        ! Compute densities
        call compute_density_spin(eigvecs, L, n_up, density_up, ierr)
        call compute_density_spin(eigvecs, L, n_dw, density_dw, ierr)

        ! Check bounds using utility function
        call check_density_bounds(density_up, density_dw, L, all_valid, ierr)
        call check(ierr == 0, "Bounds: check should succeed")
        call check(all_valid, "Bounds: all densities should be physical")

        ! Explicit checks
        do i = 1, L
            call check(density_up(i) >= -TOL, "Bounds: n_up(i) ≥ 0")
            call check(density_up(i) <= 1.0_dp + TOL, "Bounds: n_up(i) ≤ 1")
            call check(density_dw(i) >= -TOL, "Bounds: n_down(i) ≥ 0")
            call check(density_dw(i) <= 1.0_dp + TOL, "Bounds: n_down(i) ≤ 1")
            call check(density_up(i) + density_dw(i) <= 2.0_dp + TOL, &
                       "Bounds: n_total(i) ≤ 2")
        end do
    end subroutine test_physical_bounds

    !> Test density from harmonic trap potential
    !!
    !! Physics: In a harmonic trap V(i) = k·(i - i₀)², the particles are pulled
    !! toward the trap centre, so the density at the centre is larger than at
    !! either edge. The envelope is Thomas-Fermi-like, but it is NOT monotonic:
    !! a finite number of occupied Kohn-Sham levels produces Friedel
    !! oscillations (n_up + n_dw = 10 particles gives visible shell structure),
    !! so local dips such as n(centre+2) > n(centre) are physically expected and
    !! must not be asserted against. What is exact here is particle number
    !! conservation, sum_i n(i) = N, which is checked to 1e-10.
    subroutine test_density_from_harmonic_trap()
        use fortuno_serial, only: check => serial_check
        use density_calculator, only: compute_density_spin, compute_total_density
        use hamiltonian_builder, only: build_hamiltonian
        use boundary_conditions, only: BC_OPEN
        use lapack_wrapper, only: diagonalize_symmetric_real
        use potential_harmonic, only: apply_potential_harmonic
        integer, parameter :: L = 20
        integer, parameter :: n_up = 5, n_dw = 5
        real(dp), parameter :: spring_const = 0.1_dp
        ! Windows for the confinement assertions, spelled out for L = 20 to keep
        ! them free of truncating integer division:
        integer, parameter :: MID_LO = 6, MID_HI = 15     ! middle half, sites 6..15
        integer, parameter :: OUT_LEFT = 3                ! outer sixth: sites 1..3
        integer, parameter :: OUT_RIGHT = 18              ! outer sixth: sites 18..20
        real(dp) :: V_ext(L), V_xc(L), H(L, L)
        real(dp) :: eigenvals(L), eigvecs(L, L)
        real(dp) :: density_up(L), density_dw(L), density_total(L)
        real(dp) :: center_density, edge_density
        real(dp) :: n_total, n_middle_half, n_outer_sixth
        integer :: center, ierr

        ! Create harmonic potential (automatically centered at (L+1)/2)
        V_xc = 0.0_dp  ! No XC potential (U=0)
        call apply_potential_harmonic(spring_const, L, V_ext, ierr)
        call check(ierr == 0, "Harmonic: potential creation should succeed")

        ! Build Hamiltonian with harmonic trap
        call build_hamiltonian(L, V_ext, V_xc, BC_OPEN, 0.0_dp, H, ierr)
        call check(ierr == 0, "Harmonic: Hamiltonian construction should succeed")

        ! Diagonalize
        call diagonalize_symmetric_real(H, L, eigenvals, eigvecs, ierr)
        call check(ierr == 0, "Harmonic: Diagonalization should succeed")

        ! Compute densities
        call compute_density_spin(eigvecs, L, n_up, density_up, ierr)
        call compute_density_spin(eigvecs, L, n_dw, density_dw, ierr)
        call compute_total_density(density_up, density_dw, density_total, ierr)

        ! Check shell structure: density at center > density at edges
        center = L / 2
        center_density = density_total(center)
        edge_density = (density_total(1) + density_total(L)) / 2.0_dp

        call check(center_density > edge_density, &
                   "Harmonic: density at center should exceed mean edge density")

        ! Confinement: the center must beat each edge individually, not just
        ! their average (an average can be dominated by a single large edge).
        call check(density_total(center) > density_total(1), &
                   "Harmonic: density at center should exceed left edge density")
        call check(density_total(center) > density_total(L), &
                   "Harmonic: density at center should exceed right edge density")

        ! Confinement with actual content: comparing the centre with the (nearly
        ! empty) edges is almost vacuous, so measure how much of the particle
        ! number the trap actually holds near the middle.
        !   - middle half  = sites 6..15
        !   - outer sixth  = sites 1..3 and 18..20
        n_total = real(n_up + n_dw, dp)
        n_middle_half = sum(density_total(MID_LO:MID_HI))
        n_outer_sixth = sum(density_total(1:OUT_LEFT)) + sum(density_total(OUT_RIGHT:L))

        call check(n_middle_half > 0.95_dp * n_total, &
                   "Harmonic: middle half of the trap must hold more than 95% of the particles")
        call check(n_outer_sixth < 1.0e-3_dp * n_total, &
                   "Harmonic: outer sixth of the trap must hold less than 0.1% of the particles")

        ! Check particle conservation: exact, independent of Friedel oscillations
        call check(abs(sum(density_total) - n_total) < 1.0e-10_dp, &
                   "Harmonic: total particle number should be conserved to 1e-10")
    end subroutine test_density_from_harmonic_trap

end program test_density_calculator
