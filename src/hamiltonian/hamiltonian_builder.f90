module hamiltonian_builder
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_INVALID_INPUT, ERROR_SIZE_MISMATCH, ERROR_NOT_A_NUMBER
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use boundary_conditions, only: apply_boundary_conditions, &
                                    apply_boundary_conditions_complex
    implicit none
    private

    public :: validate_hamiltonian_inputs
    public :: build_hamiltonian
    public :: build_hamiltonian_complex
    public :: build_hamiltonian_free
    public :: compute_effective_potential

contains
    
    !> @brief Validate Hamiltonian construction inputs
    !!
    !! Checks that system size and potential arrays are valid:
    !! - L > 0
    !! - Array sizes match L
    !! - No NaN or Inf values in potentials
    !!
    !! @param[in] L Number of lattice sites
    !! @param[in] V_ext External potential array
    !! @param[in] V_xc Exchange-correlation potential array
    !! @param[out] ierr Error code (0=success)
    subroutine validate_hamiltonian_inputs(L, V_ext, V_xc, ierr)
        ! No hopping parameter as input - assuming t = 1
        integer, intent(in) :: L
        real(dp), intent(in) :: V_ext(:)
        real(dp), intent(in) :: V_xc(:)
        integer, intent(out) :: ierr

        ierr = 0

        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(V_ext) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (size(V_xc) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        if (.not. ALL(ieee_is_finite(V_ext))) then
            ierr = ERROR_NOT_A_NUMBER
            return
        end if

        if (.not. ALL(ieee_is_finite(V_xc))) then
            ierr = ERROR_NOT_A_NUMBER
            return
        end if
    end subroutine validate_hamiltonian_inputs

    !> @brief Build tight-binding Hamiltonian with external and XC potentials (real)
    !!
    !! Constructs the Kohn-Sham Hamiltonian matrix:
    !! H(i,j) = (V_ext(i) + V_xc(i)) δ_ij - t δ_{i,j±1}
    !!
    !! Assumes hopping parameter t = 1 (energy units).
    !!
    !! @param[in] L Number of lattice sites
    !! @param[in] V_ext External potential (length L)
    !! @param[in] V_xc Exchange-correlation potential (length L)
    !! @param[in] bc_type Boundary condition type
    !! @param[in] theta Optional twist angle (for validation only; TBC not supported)
    !! @param[out] H Hamiltonian matrix (L×L, real)
    !! @param[out] ierr Error code (0=success)
    !!
    !! @note For twisted BC, use build_hamiltonian_complex instead
    subroutine build_hamiltonian(L, V_ext, V_xc, bc_type, theta, H, ierr)
        integer, intent(in) :: L
        real(dp), intent(in) :: V_ext(:), V_xc(:)
        integer, intent(in) :: bc_type
        real(dp), intent(in), optional :: theta
        real(dp), intent(out) :: H(:,:)
        integer, intent(out) :: ierr
        integer :: i

        call validate_hamiltonian_inputs(L, V_ext, V_xc, ierr)
        if (ierr /= 0) return

        if (size(H, 1) /= L .or. size(H, 2) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        H = 0.0_dp

        do i = 1, L
            H(i,i) = V_ext(i) + V_xc(i)
        end do

        do i = 1, L - 1
            H(i,i+1) = -1.0_dp ! Assuming hopping parameter t = 1
            H(i+1,i) = -1.0_dp ! Assuming hopping parameter t = 1
        end do

        call apply_boundary_conditions(H, L, bc_type, theta, ierr)
    end subroutine build_hamiltonian

    !> @brief Build tight-binding Hamiltonian with external and XC potentials (complex)
    !!
    !! Complex version supporting all boundary conditions including twisted BC.
    !! 
    !! @param[in] L Number of lattice sites
    !! @param[in] V_ext External potential (length L)
    !! @param[in] V_xc Exchange-correlation potential (length L)
    !! @param[in] bc_type Boundary condition type
    !! @param[in] theta Optional twist angle (required for BC_TWISTED)
    !! @param[out] H Hamiltonian matrix (L×L, complex)
    !! @param[out] ierr Error code (0=success)
    subroutine build_hamiltonian_complex(L, V_ext, V_xc, bc_type, theta, H, ierr)
        integer, intent(in) :: L
        real(dp), intent(in) :: V_ext(:), V_xc(:)
        integer, intent(in) :: bc_type
        real(dp), intent(in) :: theta
        complex(dp), intent(out) :: H(:,:)
        integer, intent(out) :: ierr
        integer :: i

        call validate_hamiltonian_inputs(L, V_ext, V_xc, ierr)
        if (ierr /= 0) return

        if (size(H, 1) /= L .or. size(H, 2) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        H = cmplx(0.0_dp, 0.0_dp, kind=dp)

        do i = 1, L
            H(i,i) = cmplx(V_ext(i) + V_xc(i), 0.0_dp, kind=dp)
        end do

        do i = 1, L - 1
            H(i,i+1) = cmplx(-1.0_dp, 0.0_dp, kind=dp) ! Assuming hopping parameter t = 1
            H(i+1,i) = cmplx(-1.0_dp, 0.0_dp, kind=dp) ! Assuming hopping parameter t = 1
        end do

        call apply_boundary_conditions_complex(H, L, bc_type, theta, ierr)
    end subroutine build_hamiltonian_complex

    !> @brief Build free-particle Hamiltonian (U=0, V=0)
    !!
    !! Constructs tight-binding Hamiltonian with no potentials.
    !! Useful for testing and validation against analytical results.
    !!
    !! @param[in] L Number of lattice sites
    !! @param[in] bc_type Boundary condition type
    !! @param[in] theta Optional twist angle
    !! @param[out] H Hamiltonian matrix (L×L, real)
    !! @param[out] ierr Error code (0=success)
    subroutine build_hamiltonian_free(L, bc_type, theta, H, ierr)
        integer, intent(in) :: L
        integer, intent(in) :: bc_type
        real(dp), intent(in), optional :: theta
        real(dp), intent(out) :: H(:,:)
        integer, intent(out) :: ierr
        integer :: i

        ierr = 0
        if (L <= 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        if (size(H, 1) /= L .or. size(H, 2) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        H = 0.0_dp

        do i = 1, L - 1
            H(i,i+1) = -1.0_dp ! Assuming hopping parameter t = 1
            H(i+1,i) = -1.0_dp ! Assuming hopping parameter t = 1
        end do

        call apply_boundary_conditions(H, L, bc_type, theta, ierr)
    end subroutine build_hamiltonian_free

    !> @brief Compute effective potential V_eff = V_ext + V_xc
    !!
    !! Simple helper function to combine external and XC potentials.
    !!
    !! @param[in] V_ext External potential
    !! @param[in] V_xc Exchange-correlation potential
    !! @param[out] V_eff Effective potential
    !! @param[out] ierr Error code (0=success)
    subroutine compute_effective_potential(V_ext, V_xc, V_eff, ierr)
        real(dp), intent(in) :: V_ext(:)
        real(dp), intent(in) :: V_xc(:)
        real(dp), intent(out) :: V_eff(:)
        integer, intent(out) :: ierr
        integer :: L

        ierr = 0

        L = size(V_ext)

        if (size(V_xc) /= L .or. size(V_eff) /= L) then
            ierr = ERROR_SIZE_MISMATCH
            return
        end if

        V_eff = V_ext + V_xc
    end subroutine compute_effective_potential
end module hamiltonian_builder