!> Adaptive mixing scheme for SCF convergence
!!
!! Implements the adaptive mixing algorithm from the original C++ code.
!! The mixing parameter dynamically adjusts based on SCF convergence behavior:
!!
!! - If energy oscillates within a band for CountSCmax iterations → UpMix (more conservative)
!! - If energy keeps increasing/decreasing for too long → DwMix (more aggressive)
!!
!! This replicates the logic from lsda_stop.cc (Convergencia class).
module adaptive_mixing
    use lsda_constants, only: dp, INITIAL_MIX, ITER_MAX
    use lsda_errors, only: ERROR_SUCCESS
    implicit none
    private

    !> Adaptive mixing control parameters
    type, public :: adaptive_mix_t
        integer :: iter = 0                      !< Current iteration
        integer :: count_sc = 0                  !< Count within energy band
        integer :: count_bot = 0                 !< Count when hitting bottom
        integer :: count_top = 0                 !< Count when hitting top
        integer :: count_sc_max = 10             !< Max iterations in band before UpMix
        real(dp) :: mix = INITIAL_MIX            !< Current mixing parameter (C++ convention)
        real(dp) :: energy_top = 0.0_dp          !< Upper energy bound
        real(dp) :: energy_bot = 0.0_dp          !< Lower energy bound
        real(dp) :: energy_new = 0.0_dp          !< Current energy
        real(dp) :: energy_old = 0.0_dp          !< Previous energy
        real(dp) :: tol = 1.0e-8_dp              !< Convergence tolerance
        logical :: converged = .false.           !< Convergence flag
    end type adaptive_mix_t

    public :: adaptive_mix_init
    public :: adaptive_mix_update
    public :: adaptive_mix_get_alpha
    public :: adaptive_mix_reset

contains

    !> Initialize adaptive mixing
    !!
    !! @param[out] mix_ctrl Adaptive mixing control structure
    !! @param[in]  tol      Convergence tolerance (optional, default 1e-8)
    subroutine adaptive_mix_init(mix_ctrl, tol)
        type(adaptive_mix_t), intent(out) :: mix_ctrl
        real(dp), intent(in), optional :: tol

        mix_ctrl%iter = 0
        mix_ctrl%count_sc = 0
        mix_ctrl%count_bot = 0
        mix_ctrl%count_top = 0
        mix_ctrl%count_sc_max = 10
        mix_ctrl%mix = INITIAL_MIX
        mix_ctrl%energy_top = 0.0_dp
        mix_ctrl%energy_bot = 0.0_dp
        mix_ctrl%energy_new = 0.0_dp
        mix_ctrl%energy_old = 0.0_dp
        mix_ctrl%converged = .false.

        if (present(tol)) then
            mix_ctrl%tol = tol
        else
            mix_ctrl%tol = 1.0e-8_dp
        end if
    end subroutine adaptive_mix_init

    !> Update adaptive mixing based on new energy
    !!
    !! This implements the C++ Convergencia::Update() logic.
    !!
    !! @param[inout] mix_ctrl Adaptive mixing control structure
    !! @param[in]    energy   New SCF energy
    subroutine adaptive_mix_update(mix_ctrl, energy)
        type(adaptive_mix_t), intent(inout) :: mix_ctrl
        real(dp), intent(in) :: energy

        real(dp) :: error, band_error

        ! Increment iteration
        mix_ctrl%iter = mix_ctrl%iter + 1

        ! Update energies
        mix_ctrl%energy_old = mix_ctrl%energy_new
        mix_ctrl%energy_new = energy

        ! First iteration: initialize bounds
        if (mix_ctrl%iter == 1) then
            mix_ctrl%energy_top = energy
            mix_ctrl%energy_bot = energy
            mix_ctrl%energy_old = energy
            mix_ctrl%converged = .false.
            return
        end if

        ! Calculate relative error
        if (abs(mix_ctrl%energy_new) > 1.0e-15_dp) then
            error = abs((mix_ctrl%energy_new - mix_ctrl%energy_old) / mix_ctrl%energy_new)
        else
            error = abs(mix_ctrl%energy_new - mix_ctrl%energy_old)
        end if

        ! Calculate band error
        if (abs(mix_ctrl%energy_bot) > 1.0e-15_dp) then
            band_error = abs((mix_ctrl%energy_top - mix_ctrl%energy_bot) / mix_ctrl%energy_bot)
        else
            band_error = abs(mix_ctrl%energy_top - mix_ctrl%energy_bot)
        end if

        ! Check if energy is within current band [Bot, Top]
        if (mix_ctrl%energy_bot <= mix_ctrl%energy_new .and. &
            mix_ctrl%energy_new <= mix_ctrl%energy_top) then

            ! Energy within band
            mix_ctrl%count_sc = mix_ctrl%count_sc + 1
            mix_ctrl%count_top = 0
            mix_ctrl%count_bot = 0

            ! Check convergence
            if (mix_ctrl%count_sc >= mix_ctrl%count_sc_max .and. &
                error < mix_ctrl%tol .and. band_error < mix_ctrl%tol) then
                mix_ctrl%converged = .true.
                return
            end if

            ! Energy in band but not converged → UpMix (more conservative)
            if (mix_ctrl%count_sc >= mix_ctrl%count_sc_max .and. &
                (error >= mix_ctrl%tol .or. band_error >= mix_ctrl%tol)) then
                call up_mix(mix_ctrl)
                call reset_counts(mix_ctrl)
            end if

        else if (mix_ctrl%energy_new > mix_ctrl%energy_top) then
            ! Energy increased above top
            mix_ctrl%energy_top = mix_ctrl%energy_new
            mix_ctrl%count_top = mix_ctrl%count_top + 1
            mix_ctrl%count_bot = 0
            mix_ctrl%count_sc = 0

        else if (mix_ctrl%energy_new < mix_ctrl%energy_bot) then
            ! Energy decreased below bottom
            mix_ctrl%energy_bot = mix_ctrl%energy_new
            mix_ctrl%count_bot = mix_ctrl%count_bot + 1
            mix_ctrl%count_top = 0
            mix_ctrl%count_sc = 0
        end if

        ! If energy only increases or decreases for too long → DwMix (more aggressive)
        ! IMPORTANT: C++ checks Mix > 0.35 to prevent Mix from becoming too small
        if ((mix_ctrl%count_bot > mix_ctrl%count_sc_max * 5 .or. &
             mix_ctrl%count_top > mix_ctrl%count_sc_max * 5) .and. &
             mix_ctrl%mix > 0.35_dp) then
            call dw_mix(mix_ctrl)

            ! Clamp Mix to prevent it from going too negative
            ! (C++ doesn't clamp explicitly, but the check above prevents problems)
            if (mix_ctrl%mix < 0.0_dp) then
                mix_ctrl%mix = 0.0_dp
            end if

            call reset_counts(mix_ctrl)
        end if

        ! Max iterations reached
        if (mix_ctrl%iter >= ITER_MAX) then
            mix_ctrl%converged = .false.
        end if

    end subroutine adaptive_mix_update

    !> Get alpha (Fortran convention) from mix (C++ convention)
    !!
    !! Fortran: n_new = (1-α)*n_old + α*n_calc  (α = weight of NEW)
    !! C++:     v_new = Mix*v_old + (1-Mix)*v_calc  (Mix = weight of OLD)
    !!
    !! Therefore: α = 1 - Mix
    !!
    !! IMPORTANT: Clamps alpha to (0, 1] to prevent linear_mixing errors
    !!
    !! @param[in] mix_ctrl Adaptive mixing control structure
    !! @return alpha Fortran mixing parameter (clamped to valid range)
    function adaptive_mix_get_alpha(mix_ctrl) result(alpha)
        type(adaptive_mix_t), intent(in) :: mix_ctrl
        real(dp) :: alpha

        alpha = 1.0_dp - mix_ctrl%mix

        ! Clamp to valid range (0, 1] to prevent linear_mixing from failing
        if (alpha <= 0.0_dp) then
            alpha = 1.0e-10_dp  ! Very small but positive
        else if (alpha > 1.0_dp) then
            alpha = 1.0_dp
        end if
    end function adaptive_mix_get_alpha

    !> Reset adaptive mixing counters (keep Mix value)
    !!
    !! @param[inout] mix_ctrl Adaptive mixing control structure
    subroutine adaptive_mix_reset(mix_ctrl)
        type(adaptive_mix_t), intent(inout) :: mix_ctrl

        call reset_counts(mix_ctrl)

        ! Reset energy bounds to current energy
        mix_ctrl%energy_top = mix_ctrl%energy_new
        mix_ctrl%energy_bot = mix_ctrl%energy_new
        mix_ctrl%energy_old = mix_ctrl%energy_new
    end subroutine adaptive_mix_reset

    !> Increase mixing parameter (more conservative)
    !!
    !! C++ implementation: NewMix = Mix + (1.0 - Mix)/1.5
    !! Capped at 0.999999999
    !!
    !! @param[inout] mix_ctrl Adaptive mixing control structure
    subroutine up_mix(mix_ctrl)
        type(adaptive_mix_t), intent(inout) :: mix_ctrl
        real(dp) :: new_mix

        new_mix = mix_ctrl%mix + (1.0_dp - mix_ctrl%mix) / 1.5_dp

        if (new_mix < 0.999999999_dp) then
            mix_ctrl%mix = new_mix
        end if
    end subroutine up_mix

    !> Decrease mixing parameter (more aggressive)
    !!
    !! C++ implementation: Mix = Mix - (1.0 - Mix)*1.9
    !! NO CLAMP - Let it go negative and wrap to 0 if needed
    !!
    !! @param[inout] mix_ctrl Adaptive mixing control structure
    subroutine dw_mix(mix_ctrl)
        type(adaptive_mix_t), intent(inout) :: mix_ctrl

        mix_ctrl%mix = mix_ctrl%mix - (1.0_dp - mix_ctrl%mix) * 1.9_dp
    end subroutine dw_mix

    !> Reset counters only (internal helper)
    !!
    !! @param[inout] mix_ctrl Adaptive mixing control structure
    subroutine reset_counts(mix_ctrl)
        type(adaptive_mix_t), intent(inout) :: mix_ctrl

        mix_ctrl%count_sc = 0
        mix_ctrl%count_bot = 0
        mix_ctrl%count_top = 0
    end subroutine reset_counts

end module adaptive_mixing
