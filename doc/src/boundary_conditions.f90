!> @brief Boundary conditions module for 1D tight-binding Hamiltonian
!!
!! This module implements different boundary conditions for the 1D Hubbard model:
!! - Open Boundary Conditions (OBC): Confined chain with no hopping at edges
!! - Periodic Boundary Conditions (PBC): Ring topology with momentum conservation
!! - Twisted Boundary Conditions (TBC): Generalized PBC with phase twist (magnetic flux)
!!
!! The tight-binding Hamiltonian in matrix form is:
!! \[ H_{ij} = -t \delta_{i,j\pm 1} + V_i \delta_{ij} \]
!!
!! Boundary conditions modify the edge terms (i=1, i=L):
!! - OBC: H(1,L) = H(L,1) = 0 (no connection)
!! - PBC: H(1,L) = H(L,1) = -t (periodic ring)
!! - TBC: H(1,L) = -t e^{iθ}, H(L,1) = -t e^{-iθ} (twisted ring)
!!
!! Physical applications:
!! - OBC: Edge states, surface physics, finite-size effects
!! - PBC: Bulk properties, momentum conservation, Bethe Ansatz calculations
!! - TBC: Persistent currents, Aharonov-Bohm effect, flux threading
!!
!! @note All functions assume hopping parameter t = 1 (energy units)
!! @note TBC requires complex Hamiltonian matrix
!!
module boundary_conditions
    use lsda_constants, only: dp, PI, TWOPI
    use lsda_errors, only: ERROR_INVALID_INPUT, ERROR_OUT_OF_BOUNDS, ERROR_SIZE_MISMATCH
    implicit none
    private

    enum, bind(c)
        enumerator :: BC_OPEN = 1
        enumerator :: BC_PERIODIC = 2
        enumerator :: BC_TWISTED = 3
    end enum

    public :: BC_OPEN, BC_PERIODIC, BC_TWISTED
    public :: apply_boundary_conditions
    public :: apply_boundary_conditions_complex
    public :: validate_bc_parameters
    public :: get_free_particle_eigenvalues

contains
    !> @brief Validate boundary condition parameters
    !!
    !! Checks if BC type and associated parameters are physically valid:
    !! - BC type must be BC_OPEN, BC_PERIODIC, or BC_TWISTED
    !! - System size L must be > 1
    !! - For TBC, theta must be present and in range [0, 2π)
    !!
    !! @param[in] bc_type Boundary condition type (BC_OPEN/BC_PERIODIC/BC_TWISTED)
    !! @param[in] theta Optional twist angle in radians (required for BC_TWISTED)
    !! @param[in] L Number of lattice sites
    !! @param[out] ierr Error code (0=success)
    !!
    !! @note Twist angle θ ∈ [0, 2π) corresponds to magnetic flux Φ/Φ₀ = θ/(2π)
    subroutine validate_bc_parameters(bc_type, theta, L, ierr)
        integer, intent(in) :: bc_type
        real(dp), intent(in), optional :: theta
        integer, intent(in) :: L
        integer, intent(out) :: ierr
        integer, parameter :: allowed(3) = [BC_OPEN, BC_PERIODIC, BC_TWISTED]

        ierr = 0

        if (L <= 1) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (.not. any(bc_type == allowed)) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (bc_type == BC_OPEN) then
            return
        end if

        if (bc_type == BC_TWISTED) then
            if (.not. present(theta)) then
                ierr = ERROR_INVALID_INPUT
                return
            end if

            if (theta < 0.0_dp .or. theta >= TWOPI) then
                ierr = ERROR_OUT_OF_BOUNDS
                return
            end if
        end if
    end subroutine validate_bc_parameters

    !> @brief Apply boundary conditions to real Hamiltonian matrix
    !!
    !! Modifies edge elements of tight-binding Hamiltonian matrix H(L×L)
    !! according to specified boundary conditions. Assumes hopping t = 1.
    !!
    !! Modifications:
    !! - BC_OPEN: No modification (H already tridiagonal)
    !! - BC_PERIODIC: Sets H(1,L) = H(L,1) = -1
    !! - BC_TWISTED: Not supported for real H (returns error)
    !!
    !! @param[inout] H Hamiltonian matrix (L×L, real)
    !! @param[in] L Number of lattice sites
    !! @param[in] bc_type Boundary condition type
    !! @param[in] theta Optional twist angle (ignored for OBC/PBC)
    !! @param[out] ierr Error code (0=success)
    !!
    !! @note For TBC, use apply_boundary_conditions_complex instead
    !! @note Matrix H must be allocated with dimensions (L,L) before calling
    subroutine apply_boundary_conditions(H, L, bc_type, theta, ierr)
        real(dp), intent(inout) :: H(:,:)
        integer, intent(in) :: L
        integer, intent(in) :: bc_type
        real(dp), intent(in), optional :: theta
        integer, intent(out) :: ierr

        ierr = 0
        call validate_bc_parameters(bc_type, theta, L, ierr)

        if (ierr /= 0) return

        if (size(H, 1) /= L .or. size(H, 2) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        select case (bc_type)
            case (BC_OPEN)
                return
            case (BC_PERIODIC)
                H(1,L) = -1.0_dp ! Assuming hopping parameter t = 1
                H(L,1) = -1.0_dp ! Assuming hopping parameter t = 1
            case (BC_TWISTED)
                ierr = ERROR_INVALID_INPUT
                return
            case default
            ierr = ERROR_INVALID_INPUT
        end select
    end subroutine apply_boundary_conditions

    !> @brief Apply boundary conditions to complex Hamiltonian matrix
    !!
    !! Modifies edge elements of tight-binding Hamiltonian matrix H(L×L)
    !! for all boundary condition types, including twisted BC.
    !! Assumes hopping t = 1.
    !!
    !! Modifications:
    !! - BC_OPEN: No modification
    !! - BC_PERIODIC: Sets H(1,L) = H(L,1) = -1
    !! - BC_TWISTED: Sets H(1,L) = -e^{iθ}, H(L,1) = -e^{-iθ}
    !!
    !! @param[inout] H Hamiltonian matrix (L×L, complex)
    !! @param[in] L Number of lattice sites
    !! @param[in] bc_type Boundary condition type
    !! @param[in] theta Optional twist angle in radians (required for BC_TWISTED)
    !! @param[out] ierr Error code (0=success)
    !!
    !! @note Twist angle θ represents Aharonov-Bohm phase from magnetic flux
    !! @note For θ=0, TBC reduces to PBC
    !! @note For θ=π, antiperiodic BC (useful for odd electron numbers)
    subroutine apply_boundary_conditions_complex(H, L, bc_type, theta, ierr)
        complex(dp), intent(inout) :: H(:,:)
        integer, intent(in) :: L
        integer, intent(in) :: bc_type
        real(dp), intent(in), optional :: theta
        complex(dp), allocatable :: H_complex(:,:)
        integer, intent(out) :: ierr
        real(dp) :: theta_val

        complex(dp), parameter :: IMAG_UNIT = (0.0_dp, 1.0_dp)

        ierr = 0
        call validate_bc_parameters(bc_type, theta, L, ierr)

        if (ierr /= 0) return

        if (size(H, 1) /= L .or. size(H, 2) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        select case (bc_type)
            case (BC_OPEN)
                return
            case (BC_PERIODIC)
                H(1, L) = cmplx(-1.0_dp, 0.0_dp, kind=dp) ! Assuming hopping parameter t = 1
                H(L, 1) = cmplx(-1.0_dp, 0.0_dp, kind=dp) ! Assuming hopping parameter t = 1
            case (BC_TWISTED)
                theta_val = theta

                H(1,L) = -1.0_dp * exp(IMAG_UNIT * theta_val)
                H(L,1) = -1.0_dp * exp(-IMAG_UNIT * theta_val)

            case default
            ierr = ERROR_INVALID_INPUT
        end select
    end subroutine apply_boundary_conditions_complex

    !> @brief Calculate analytical eigenvalues for free particles (U=0, V=0)
    !!
    !! Returns exact eigenvalues for non-interacting tight-binding model
    !! with specified boundary conditions. Assumes hopping t = 1.
    !!
    !! Eigenvalue formulas:
    !! - OBC: E_n = -2cos(nπ/(L+1)), n = 1,...,L
    !! - PBC: E_k = -2cos(2πk/L), k = 0,...,L-1
    !! - TBC: E_k(θ) = -2cos((2πk + θ)/L), k = 0,...,L-1
    !!
    !! @param[in] L Number of lattice sites
    !! @param[in] bc_type Boundary condition type
    !! @param[in] theta Optional twist angle (required for BC_TWISTED)
    !! @param[out] eigenvalues Array of L eigenvalues (ordered)
    !! @param[out] ierr Error code (0=success)
    !!
    !! @note Useful for validating numerical diagonalization
    !! @note For PBC/TBC, momentum k is a good quantum number
    !! @note For OBC, eigenvalues correspond to standing wave solutions
    subroutine get_free_particle_eigenvalues(L, bc_type, theta, eigenvalues, ierr)
        integer, intent(in) :: L
        integer, intent(in) :: bc_type
        real(dp), intent(in), optional :: theta
        real(dp), intent(out) :: eigenvalues(:)
        integer, intent(out) :: ierr
        integer :: n, k
        real(dp) :: theta_val

        ierr = 0
        call validate_bc_parameters(bc_type, theta, L, ierr)

        if (ierr /= 0) return

        if (size(eigenvalues) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        select case (bc_type)
            case (BC_OPEN)
                do n = 1, L
                    eigenvalues(n) = -2.0_dp * cos(n * PI / real(L + 1, dp)) ! Assuming hopping parameter t = 1
                end do
            case (BC_PERIODIC)
                do k = 0, L-1
                    eigenvalues(k+1) = -2.0_dp * cos(TWOPI * k / real(L, dp)) ! Assuming hopping parameter t = 1
                end do
            case (BC_TWISTED)
                theta_val = theta
                do k = 0, L-1
                    eigenvalues(k+1) = -2.0_dp * cos((TWOPI * k + theta_val) / real(L, dp)) ! Assuming hopping parameter t = 1
                end do
            case default
                ierr = ERROR_INVALID_INPUT
                return
        end select
    end subroutine get_free_particle_eigenvalues
end module boundary_conditions