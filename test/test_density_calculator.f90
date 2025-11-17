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
            test("density_from_harmonic_trap", test_density_from_harmonic_trap) &
        ])
    end function get_density_tests

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
    !! Physics: In a harmonic trap V(i) = 0.5·k·(i - i₀)², the density should
    !! exhibit "shell structure" - higher density at the center and decreasing
    !! toward the edges. This is analogous to cold atoms in optical traps.
    !! For non-interacting fermions, we expect a Thomas-Fermi-like profile.
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
        real(dp) :: V_ext(L), V_xc(L), H(L, L)
        real(dp) :: eigenvals(L), eigvecs(L, L)
        real(dp) :: density_up(L), density_dw(L), density_total(L)
        real(dp) :: center_density, edge_density
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
                   "Harmonic: density at center should exceed edge density")

        ! Check qualitative Gaussian-like profile (peak at center)
        call check(density_total(center) > density_total(center - 2), &
                   "Harmonic: density should decrease away from center (left)")
        call check(density_total(center) > density_total(center + 2), &
                   "Harmonic: density should decrease away from center (right)")

        ! Check particle conservation
        call check(abs(sum(density_total) - real(n_up + n_dw, dp)) < TOL, &
                   "Harmonic: total particle number should be conserved")
    end subroutine test_density_from_harmonic_trap

end program test_density_calculator
