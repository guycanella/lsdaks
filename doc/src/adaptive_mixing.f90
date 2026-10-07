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
    use lsda_constants, only: dp, INITIAL_MIX, MIX_ALPHA_MIN
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
        real(dp) :: alpha_min = MIX_ALPHA_MIN    !< Floor for alpha = 1 - mix
        logical :: converged = .false.           !< Convergence flag
    end type adaptive_mix_t

    !> Cap applied by up_mix to the raw C++ formula (lsda_stop.cc:266-272).
    !! The effective cap is the tighter of this value and 1 - alpha_min.
    real(dp), parameter, public :: MIX_FORMULA_CAP = 0.999999999_dp

    public :: adaptive_mix_init
    public :: adaptive_mix_update
    public :: adaptive_mix_get_alpha
    public :: adaptive_mix_reset

contains

    !> Initialize adaptive mixing
    !!
    !! @param[out] mix_ctrl      Adaptive mixing control structure
    !! @param[in]  tol           Convergence tolerance (optional, default 1e-8)
    !! @param[in]  alpha_min     Floor for alpha = 1 - mix (optional, default
    !!                           MIX_ALPHA_MIN). Only 0 < alpha_min <= 1 is
    !!                           accepted; anything else falls back to the
    !!                           default.
    !! @param[in]  initial_alpha Starting mixing weight in the Fortran convention
    !!                           (optional). When present the controller starts at
    !!                           mix = 1 - initial_alpha instead of the hard-coded
    !!                           INITIAL_MIX, so that a user-supplied mixing_alpha
    !!                           is honoured even with the adaptive controller on.
    !!                           Values outside (0, 1] are ignored.
    subroutine adaptive_mix_init(mix_ctrl, tol, alpha_min, initial_alpha)
        type(adaptive_mix_t), intent(out) :: mix_ctrl
        real(dp), intent(in), optional :: tol
        real(dp), intent(in), optional :: alpha_min
        real(dp), intent(in), optional :: initial_alpha

        mix_ctrl%iter = 0
        mix_ctrl%count_sc = 0
        mix_ctrl%count_bot = 0
        mix_ctrl%count_top = 0
        mix_ctrl%count_sc_max = 10
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

        if (present(alpha_min)) then
            ! A non-positive floor would reopen the alpha -> 0 pathology, and
            ! linear_mixing rejects alpha <= 0 outright. A floor ABOVE 1 is just
            ! as invalid: clamp_mix would pin mix at 0 and adaptive_mix_get_alpha
            ! would return alpha_min itself (the `alpha > 1` branch is an
            ! else-if and is never re-evaluated after the floor is applied), so
            ! the controller would hand out a mixing weight greater than 1.
            if (alpha_min > 0.0_dp .and. alpha_min <= 1.0_dp) then
                mix_ctrl%alpha_min = alpha_min
            else
                mix_ctrl%alpha_min = MIX_ALPHA_MIN
            end if
        else
            mix_ctrl%alpha_min = MIX_ALPHA_MIN
        end if

        mix_ctrl%mix = INITIAL_MIX
        if (present(initial_alpha)) then
            if (initial_alpha > 0.0_dp .and. initial_alpha <= 1.0_dp) then
                mix_ctrl%mix = 1.0_dp - initial_alpha
            end if
        end if

        call clamp_mix(mix_ctrl)
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
                ! C++ Convergencia::Reset() (lsda_stop.cc:260-264) also collapses
                ! the band onto the current energy (Top = Bot = Old = New).
                ! Zeroing the counters alone would let energy_top only grow and
                ! energy_bot only shrink for the whole run, making the band error
                ! monotonically non-decreasing and the convergence test below
                ! permanently unreachable.
                call adaptive_mix_reset(mix_ctrl)
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

            ! Same rationale as after up_mix: the C++ Reset() collapses the band.
            call adaptive_mix_reset(mix_ctrl)
        end if

        ! NOTE: there is deliberately no iteration cap here. This module only
        ! decides HOW MUCH to mix; the number of SCF iterations belongs to the
        ! caller, which loops up to scf_params%max_iter. The previous code
        ! compared mix_ctrl%iter against the global constant ITER_MAX (10000)
        ! and, on reaching it, forced converged = .false. That was wrong twice
        ! over: it ignored the user's max_iter, and it did the opposite of the
        ! C++ original (lsda_stop.cc, which sets its Stop flag, i.e. terminates
        ! the loop). Since the SCF cycle no longer consults mix_ctrl%converged
        ! at all - convergence is decided by residual_V plus energy stability -
        ! the only observable effect left was to silently clear a flag on
        ! iteration 10000. Dropping it keeps the mixing state a pure function of
        ! the energy sequence.
    end subroutine adaptive_mix_update

    !> Get alpha (Fortran convention) from mix (C++ convention)
    !!
    !! Fortran: n_new = (1-α)*n_old + α*n_calc  (α = weight of NEW)
    !! C++:     v_new = Mix*v_old + (1-Mix)*v_calc  (Mix = weight of OLD)
    !!
    !! Therefore: α = 1 - Mix
    !!
    !! IMPORTANT: Clamps alpha to [alpha_min, 1] to prevent linear_mixing errors
    !! and, above all, to keep the SCF moving. Repeated up_mix drives mix towards
    !! 1, i.e. alpha towards 0; with an unbounded alpha -> 0 each step changes the
    !! effective potential by a vanishing amount, ||delta_n|| collapses in step
    !! with alpha and a density-based criterion declares a convergence that never
    !! happened. The floor alpha_min (default MIX_ALPHA_MIN) forbids that regime.
    !!
    !! @param[in] mix_ctrl Adaptive mixing control structure
    !! @return alpha Fortran mixing parameter (clamped to valid range)
    function adaptive_mix_get_alpha(mix_ctrl) result(alpha)
        type(adaptive_mix_t), intent(in) :: mix_ctrl
        real(dp) :: alpha

        alpha = 1.0_dp - mix_ctrl%mix

        ! Clamp to [alpha_min, 1]
        if (alpha < mix_ctrl%alpha_min) then
            alpha = mix_ctrl%alpha_min
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
    !! C++ implementation: NewMix = Mix + (1.0 - Mix)/1.5, capped at
    !! MIX_FORMULA_CAP = 0.999999999 (lsda_stop.cc:266-272).
    !!
    !! On top of that formula cap, `mix` itself is clamped to 1 - alpha_min (see
    !! clamp_mix): otherwise the state and the mixing weight actually used would
    !! disagree, since adaptive_mix_get_alpha never returns less than alpha_min.
    !!
    !! @param[inout] mix_ctrl Adaptive mixing control structure
    subroutine up_mix(mix_ctrl)
        type(adaptive_mix_t), intent(inout) :: mix_ctrl
        real(dp) :: new_mix

        new_mix = mix_ctrl%mix + (1.0_dp - mix_ctrl%mix) / 1.5_dp

        if (new_mix < MIX_FORMULA_CAP) then
            mix_ctrl%mix = new_mix
        end if

        call clamp_mix(mix_ctrl)
    end subroutine up_mix

    !> Clamp `mix` to [0, 1 - alpha_min]
    !!
    !! adaptive_mix_get_alpha already clamps alpha = 1 - mix from below at
    !! alpha_min, but without this the stored `mix` could keep climbing towards
    !! 1 (the formula cap 0.999999999) while the alpha actually applied stayed
    !! pinned at alpha_min. The controller would then be stuck: dw_mix computes
    !! Mix - (1 - Mix)*1.9, whose correction (1 - Mix) is ~1e-9 at the formula
    !! cap, so the "be more aggressive" branch could no longer move the mixing
    !! weight at all. Clamping the state keeps 1 - mix equal to the alpha in use
    !! and leaves dw_mix a finite step (alpha_min*1.9) to recover with.
    !!
    !! @param[inout] mix_ctrl Adaptive mixing control structure
    subroutine clamp_mix(mix_ctrl)
        type(adaptive_mix_t), intent(inout) :: mix_ctrl
        real(dp) :: mix_cap

        mix_cap = max(0.0_dp, 1.0_dp - mix_ctrl%alpha_min)

        if (mix_ctrl%mix > mix_cap) then
            mix_ctrl%mix = mix_cap
        else if (mix_ctrl%mix < 0.0_dp) then
            mix_ctrl%mix = 0.0_dp
        end if
    end subroutine clamp_mix

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
