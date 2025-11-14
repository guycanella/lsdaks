!> Unit tests for hamiltonian_builder module
program test_hamiltonian_builder
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp
    implicit none

    real(dp), parameter :: TOL = 1.0e-9_dp

    call execute_serial_cmd_app(get_hamiltonian_tests())

contains

    function get_hamiltonian_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("validate_inputs_valid", test_validate_inputs_valid), &
            test("validate_inputs_invalid_L", test_validate_inputs_invalid_L), &
            test("validate_inputs_size_mismatch_vext", test_validate_inputs_size_mismatch_vext), &
            test("validate_inputs_size_mismatch_vxc", test_validate_inputs_size_mismatch_vxc), &
            test("validate_inputs_nan_in_vext", test_validate_inputs_nan_in_vext), &
            test("validate_inputs_inf_in_vxc", test_validate_inputs_inf_in_vxc), &
            test("compute_veff_simple", test_compute_veff_simple), &
            test("compute_veff_size_mismatch", test_compute_veff_size_mismatch), &
            test("build_free_hamiltonian_structure", test_build_free_hamiltonian_structure), &
            test("build_free_hamiltonian_open_bc", test_build_free_hamiltonian_open_bc), &
            test("build_free_hamiltonian_periodic_bc", test_build_free_hamiltonian_periodic_bc), &
            test("build_hamiltonian_diagonal", test_build_hamiltonian_diagonal), &
            test("build_hamiltonian_offdiagonal", test_build_hamiltonian_offdiagonal), &
            test("build_hamiltonian_open_bc", test_build_hamiltonian_open_bc), &
            test("build_hamiltonian_periodic_bc", test_build_hamiltonian_periodic_bc), &
            test("build_hamiltonian_symmetry", test_build_hamiltonian_symmetry), &
            test("build_hamiltonian_complex_real_parts", test_build_hamiltonian_complex_real_parts), &
            test("build_hamiltonian_size_mismatch", test_build_hamiltonian_size_mismatch) &
        ])
    end function get_hamiltonian_tests

    !> Test validation with correct inputs
    !!
    !! Physics: Before constructing the Hamiltonian, we must verify that all
    !! input arrays are properly sized and contain valid numerical values.
    !! This prevents silent errors that could propagate through calculations.
    subroutine test_validate_inputs_valid()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: validate_hamiltonian_inputs
        integer, parameter :: L = 10
        real(dp) :: V_ext(L), V_xc(L)
        integer :: ierr

        V_ext = 0.0_dp
        V_xc = 0.0_dp

        call validate_hamiltonian_inputs(L, V_ext, V_xc, ierr)
        call check(ierr == 0, "Valid inputs should pass validation")
    end subroutine test_validate_inputs_valid

    !> Test validation failure for invalid system size
    !!
    !! Physics: A lattice must have at least one site. L ≤ 0 is unphysical
    !! and indicates a programming error that should be caught early.
    subroutine test_validate_inputs_invalid_L()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: validate_hamiltonian_inputs
        use lsda_errors, only: ERROR_INVALID_INPUT
        real(dp) :: V_ext(1), V_xc(1)
        integer :: ierr

        V_ext = 0.0_dp
        V_xc = 0.0_dp

        call validate_hamiltonian_inputs(0, V_ext, V_xc, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "L=0 should fail")

        call validate_hamiltonian_inputs(-1, V_ext, V_xc, ierr)
        call check(ierr == ERROR_INVALID_INPUT, "L<0 should fail")
    end subroutine test_validate_inputs_invalid_L

    !> Test validation failure for V_ext size mismatch
    !!
    !! Physics: The external potential must be defined at every lattice site.
    !! A mismatch between array size and system size indicates inconsistent
    !! problem setup.
    subroutine test_validate_inputs_size_mismatch_vext()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: validate_hamiltonian_inputs
        use lsda_errors, only: ERROR_SIZE_MISMATCH
        integer, parameter :: L = 10
        real(dp) :: V_ext(L-1), V_xc(L)
        integer :: ierr

        V_ext = 0.0_dp
        V_xc = 0.0_dp

        call validate_hamiltonian_inputs(L, V_ext, V_xc, ierr)
        call check(ierr == ERROR_SIZE_MISMATCH, "V_ext size mismatch should fail")
    end subroutine test_validate_inputs_size_mismatch_vext

    !> Test validation failure for V_xc size mismatch
    !!
    !! Physics: The exchange-correlation potential from DFT must also be
    !! defined at every site. Size mismatch indicates an error in the
    !! self-consistent cycle or table interpolation.
    subroutine test_validate_inputs_size_mismatch_vxc()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: validate_hamiltonian_inputs
        use lsda_errors, only: ERROR_SIZE_MISMATCH
        integer, parameter :: L = 10
        real(dp) :: V_ext(L), V_xc(L+1)
        integer :: ierr

        V_ext = 0.0_dp
        V_xc = 0.0_dp

        call validate_hamiltonian_inputs(L, V_ext, V_xc, ierr)
        call check(ierr == ERROR_SIZE_MISMATCH, "V_xc size mismatch should fail")
    end subroutine test_validate_inputs_size_mismatch_vxc

    !> Test validation failure for NaN in V_ext
    !!
    !! Physics: NaN (Not a Number) values in the external potential indicate
    !! numerical errors in potential construction (e.g., 0/0 in quasiperiodic
    !! potential). These must be caught before attempting diagonalization.
    subroutine test_validate_inputs_nan_in_vext()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: validate_hamiltonian_inputs
        use lsda_errors, only: ERROR_NOT_A_NUMBER
        integer, parameter :: L = 10
        real(dp) :: V_ext(L), V_xc(L), zero, nan
        integer :: ierr

        zero = 0.0_dp
        nan = zero / zero

        V_ext = 0.0_dp
        V_ext(5) = nan
        V_xc = 0.0_dp

        call validate_hamiltonian_inputs(L, V_ext, V_xc, ierr)
        call check(ierr == ERROR_NOT_A_NUMBER, "NaN in V_ext should fail")
    end subroutine test_validate_inputs_nan_in_vext

    !> Test validation failure for Inf in V_xc
    !!
    !! Physics: Infinite values in V_xc indicate overflow in functional
    !! evaluation or spline extrapolation. Such potentials lead to unphysical
    !! eigenstates and must be rejected.
    subroutine test_validate_inputs_inf_in_vxc()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: validate_hamiltonian_inputs
        use lsda_errors, only: ERROR_NOT_A_NUMBER
        integer, parameter :: L = 10
        real(dp) :: V_ext(L), V_xc(L), zero, inf
        integer :: ierr

        zero = 0.0_dp
        inf = 1.0_dp / zero

        V_ext = 0.0_dp
        V_xc = 0.0_dp
        V_xc(3) = inf

        call validate_hamiltonian_inputs(L, V_ext, V_xc, ierr)
        call check(ierr == ERROR_NOT_A_NUMBER, "Inf in V_xc should fail")
    end subroutine test_validate_inputs_inf_in_vxc

    !> Test effective potential computation
    !!
    !! Physics: In Kohn-Sham DFT, the effective single-particle potential
    !! is V_eff = V_ext + V_Hartree + V_xc. For the 1D Hubbard model,
    !! V_Hartree is absorbed into V_xc, so V_eff = V_ext + V_xc is the
    !! total potential seen by electrons.
    subroutine test_compute_veff_simple()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: compute_effective_potential
        integer, parameter :: L = 5
        real(dp) :: V_ext(L), V_xc(L), V_eff(L)
        integer :: ierr

        V_ext = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
        V_xc = [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp]

        call compute_effective_potential(V_ext, V_xc, V_eff, ierr)

        call check(ierr == 0, "Should succeed")
        call check(abs(V_eff(1) - 1.5_dp) < TOL, "V_eff(1) = V_ext(1) + V_xc(1)")
        call check(abs(V_eff(3) - 3.5_dp) < TOL, "V_eff(3) = V_ext(3) + V_xc(3)")
        call check(abs(V_eff(5) - 5.5_dp) < TOL, "V_eff(5) = V_ext(5) + V_xc(5)")
    end subroutine test_compute_veff_simple

    !> Test effective potential size mismatch
    !!
    !! Physics: All potential arrays must have the same length (number of sites).
    !! Mismatched sizes indicate inconsistent data that could lead to
    !! accessing uninitialized memory.
    subroutine test_compute_veff_size_mismatch()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: compute_effective_potential
        use lsda_errors, only: ERROR_SIZE_MISMATCH
        real(dp) :: V_ext(5), V_xc(6), V_eff(5)
        integer :: ierr

        V_ext = 0.0_dp
        V_xc = 0.0_dp

        call compute_effective_potential(V_ext, V_xc, V_eff, ierr)
        call check(ierr == ERROR_SIZE_MISMATCH, "Size mismatch should fail")
    end subroutine test_compute_veff_size_mismatch

    !> Test free Hamiltonian has tridiagonal structure
    !!
    !! Physics: The tight-binding Hamiltonian with only nearest-neighbor
    !! hopping (and no potentials) produces a tridiagonal matrix. Elements
    !! beyond the main diagonal and first off-diagonals must be zero,
    !! except for boundary terms in periodic BC.
    subroutine test_build_free_hamiltonian_structure()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: build_hamiltonian_free
        use boundary_conditions, only: BC_OPEN
        integer, parameter :: L = 5
        real(dp) :: H(L,L)
        integer :: ierr, i, j

        call build_hamiltonian_free(L, BC_OPEN, H=H, ierr=ierr)

        call check(ierr == 0, "Should succeed")

        ! Check that only tridiagonal elements are non-zero (for open BC)
        do i = 1, L
            do j = 1, L
                if (abs(i - j) > 1) then
                    call check(abs(H(i,j)) < TOL, "Non-tridiagonal elements should be zero")
                end if
            end do
        end do
    end subroutine test_build_free_hamiltonian_structure

    !> Test free Hamiltonian with open BC
    !!
    !! Physics: Open boundary conditions for free particles create a
    !! finite chain with hard-wall boundaries. The Hamiltonian is purely
    !! tridiagonal: H(i,i±1) = -t, H(i,i) = 0, with no connections at edges.
    subroutine test_build_free_hamiltonian_open_bc()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: build_hamiltonian_free
        use boundary_conditions, only: BC_OPEN
        integer, parameter :: L = 5
        real(dp) :: H(L,L)
        integer :: ierr

        call build_hamiltonian_free(L, BC_OPEN, H=H, ierr=ierr)

        call check(ierr == 0, "Should succeed")

        ! Diagonal should be zero (no potential)
        call check(abs(H(1,1)) < TOL, "Diagonal should be zero")
        call check(abs(H(3,3)) < TOL, "Diagonal should be zero")

        ! Off-diagonal should be -1 (hopping t=1)
        call check(abs(H(1,2) + 1.0_dp) < TOL, "H(1,2) should be -1")
        call check(abs(H(2,1) + 1.0_dp) < TOL, "H(2,1) should be -1")

        ! No boundary connections
        call check(abs(H(1,L)) < TOL, "H(1,L) should be zero for open BC")
        call check(abs(H(L,1)) < TOL, "H(L,1) should be zero for open BC")
    end subroutine test_build_free_hamiltonian_open_bc

    !> Test free Hamiltonian with periodic BC
    !!
    !! Physics: Periodic boundary conditions connect the last site back to
    !! the first, creating a ring. The Hamiltonian is circulant with
    !! H(1,L) = H(L,1) = -t, enabling Bloch wave solutions with momentum k.
    subroutine test_build_free_hamiltonian_periodic_bc()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: build_hamiltonian_free
        use boundary_conditions, only: BC_PERIODIC
        integer, parameter :: L = 5
        real(dp) :: H(L,L)
        integer :: ierr

        call build_hamiltonian_free(L, BC_PERIODIC, H=H, ierr=ierr)

        call check(ierr == 0, "Should succeed")

        ! Check periodic connections
        call check(abs(H(1,L) + 1.0_dp) < TOL, "H(1,L) should be -1 for periodic BC")
        call check(abs(H(L,1) + 1.0_dp) < TOL, "H(L,1) should be -1 for periodic BC")
    end subroutine test_build_free_hamiltonian_periodic_bc

    !> Test Hamiltonian diagonal elements contain potentials
    !!
    !! Physics: The diagonal elements of the Kohn-Sham Hamiltonian represent
    !! on-site energies, which are the sum of external potential (traps,
    !! barriers) and exchange-correlation potential from many-body effects.
    subroutine test_build_hamiltonian_diagonal()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: build_hamiltonian
        use boundary_conditions, only: BC_OPEN
        integer, parameter :: L = 5
        real(dp) :: V_ext(L), V_xc(L), H(L,L)
        integer :: ierr

        V_ext = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
        V_xc = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp]

        call build_hamiltonian(L, V_ext, V_xc, BC_OPEN, H=H, ierr=ierr)

        call check(ierr == 0, "Should succeed")

        ! Check diagonal = V_ext + V_xc
        call check(abs(H(1,1) - 1.1_dp) < TOL, "H(1,1) = V_ext(1) + V_xc(1)")
        call check(abs(H(2,2) - 2.2_dp) < TOL, "H(2,2) = V_ext(2) + V_xc(2)")
        call check(abs(H(5,5) - 5.5_dp) < TOL, "H(5,5) = V_ext(5) + V_xc(5)")
    end subroutine test_build_hamiltonian_diagonal

    !> Test Hamiltonian off-diagonal hopping elements
    !!
    !! Physics: Off-diagonal elements represent nearest-neighbor hopping
    !! with amplitude -t (t=1 in our units). This kinetic energy term
    !! allows electrons to delocalize and form bands. For hermitian H,
    !! we must have H(i,j) = H(j,i).
    subroutine test_build_hamiltonian_offdiagonal()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: build_hamiltonian
        use boundary_conditions, only: BC_OPEN
        integer, parameter :: L = 5
        real(dp) :: V_ext(L), V_xc(L), H(L,L)
        integer :: ierr

        V_ext = 0.0_dp
        V_xc = 0.0_dp

        call build_hamiltonian(L, V_ext, V_xc, BC_OPEN, H=H, ierr=ierr)

        call check(ierr == 0, "Should succeed")

        ! Check hopping elements
        call check(abs(H(1,2) + 1.0_dp) < TOL, "H(1,2) should be -t = -1")
        call check(abs(H(2,3) + 1.0_dp) < TOL, "H(2,3) should be -t = -1")
        call check(abs(H(4,5) + 1.0_dp) < TOL, "H(4,5) should be -t = -1")
    end subroutine test_build_hamiltonian_offdiagonal

    !> Test Hamiltonian with open BC has no edge connections
    !!
    !! Physics: Open boundaries model a finite system with no current flow
    !! at the edges. This creates localized edge states and standing waves,
    !! crucial for studying surface physics and quantum confinement.
    subroutine test_build_hamiltonian_open_bc()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: build_hamiltonian
        use boundary_conditions, only: BC_OPEN
        integer, parameter :: L = 6
        real(dp) :: V_ext(L), V_xc(L), H(L,L)
        integer :: ierr

        V_ext = 1.0_dp
        V_xc = 0.5_dp

        call build_hamiltonian(L, V_ext, V_xc, BC_OPEN, H=H, ierr=ierr)

        call check(ierr == 0, "Should succeed")
        call check(abs(H(1,L)) < TOL, "No connection at edges for open BC")
        call check(abs(H(L,1)) < TOL, "No connection at edges for open BC")
    end subroutine test_build_hamiltonian_open_bc

    !> Test Hamiltonian with periodic BC has edge connections
    !!
    !! Physics: Periodic boundaries eliminate edge effects and restore
    !! translational symmetry. Essential for bulk properties, Bethe Ansatz
    !! calculations, and when finite-size effects must be minimized.
    subroutine test_build_hamiltonian_periodic_bc()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: build_hamiltonian
        use boundary_conditions, only: BC_PERIODIC
        integer, parameter :: L = 6
        real(dp) :: V_ext(L), V_xc(L), H(L,L)
        integer :: ierr

        V_ext = 1.0_dp
        V_xc = 0.5_dp

        call build_hamiltonian(L, V_ext, V_xc, BC_PERIODIC, H=H, ierr=ierr)

        call check(ierr == 0, "Should succeed")
        call check(abs(H(1,L) + 1.0_dp) < TOL, "Periodic connection at edges")
        call check(abs(H(L,1) + 1.0_dp) < TOL, "Periodic connection at edges")
    end subroutine test_build_hamiltonian_periodic_bc

    !> Test Hamiltonian is Hermitian (symmetric for real H)
    !!
    !! Physics: The Hamiltonian must be Hermitian (H† = H) to ensure real
    !! eigenvalues (energies) and unitary time evolution. For real matrices,
    !! Hermiticity reduces to symmetry: H(i,j) = H(j,i). This is fundamental
    !! to quantum mechanics.
    subroutine test_build_hamiltonian_symmetry()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: build_hamiltonian
        use boundary_conditions, only: BC_PERIODIC
        integer, parameter :: L = 6
        real(dp) :: V_ext(L), V_xc(L), H(L,L)
        integer :: ierr, i, j

        V_ext = [1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp, 3.5_dp]
        V_xc = 0.5_dp

        call build_hamiltonian(L, V_ext, V_xc, BC_PERIODIC, H=H, ierr=ierr)

        call check(ierr == 0, "Should succeed")

        ! Check symmetry H(i,j) = H(j,i)
        do i = 1, L
            do j = 1, L
                call check(abs(H(i,j) - H(j,i)) < TOL, &
                           "Hamiltonian should be symmetric (Hermitian)")
            end do
        end do
    end subroutine test_build_hamiltonian_symmetry

    !> Test complex Hamiltonian has correct real parts
    !!
    !! Physics: Even for twisted BC with complex hopping, the on-site
    !! potentials remain real (no local phase). Only the hopping terms
    !! acquire complex phases from the Aharonov-Bohm effect. The diagonal
    !! must be purely real.
    subroutine test_build_hamiltonian_complex_real_parts()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: build_hamiltonian_complex
        use boundary_conditions, only: BC_TWISTED
        use lsda_constants, only: PI
        integer, parameter :: L = 5
        real(dp) :: V_ext(L), V_xc(L), theta
        complex(dp) :: H(L,L)
        integer :: ierr, i

        V_ext = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
        V_xc = 0.5_dp
        theta = PI / 4.0_dp

        call build_hamiltonian_complex(L, V_ext, V_xc, BC_TWISTED, theta, H, ierr)

        call check(ierr == 0, "Should succeed")

        ! Diagonal should be real (on-site potentials)
        do i = 1, L
            call check(abs(aimag(H(i,i))) < TOL, "Diagonal should be real")
            call check(abs(real(H(i,i)) - (V_ext(i) + V_xc(i))) < TOL, &
                       "Diagonal should equal V_ext + V_xc")
        end do
    end subroutine test_build_hamiltonian_complex_real_parts

    !> Test Hamiltonian construction with size mismatch
    !!
    !! Physics: The Hamiltonian matrix must be L×L to represent L sites.
    !! A mismatch indicates memory allocation errors or incorrect problem
    !! setup. This must be caught before attempting expensive diagonalization.
    subroutine test_build_hamiltonian_size_mismatch()
        use fortuno_serial, only: check => serial_check
        use hamiltonian_builder, only: build_hamiltonian
        use boundary_conditions, only: BC_OPEN
        use lsda_errors, only: ERROR_SIZE_MISMATCH
        integer, parameter :: L = 5
        real(dp) :: V_ext(L), V_xc(L), H(L+1,L+1)
        integer :: ierr

        V_ext = 0.0_dp
        V_xc = 0.0_dp

        call build_hamiltonian(L, V_ext, V_xc, BC_OPEN, H=H, ierr=ierr)
        call check(ierr == ERROR_SIZE_MISMATCH, "H size mismatch should fail")
    end subroutine test_build_hamiltonian_size_mismatch

end program test_hamiltonian_builder
