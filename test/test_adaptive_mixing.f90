!> Unit tests for the adaptive mixing controller
!!
!! The adaptive mixing module reproduces the C++ Convergencia class: it watches
!! the SCF energy and retunes the mixing parameter.
!!
!!   - Energy bouncing inside the band [Bot, Top] for count_sc_max iterations
!!     without converging => UpMix (keep more of the old potential, be safer).
!!   - Energy drifting monotonically past the band for 5*count_sc_max
!!     iterations => DwMix (keep less of the old potential, be more aggressive).
!!
!! Conventions: `mix` is the C++ convention (weight of the OLD value) and
!! `alpha = 1 - mix` is the Fortran convention (weight of the NEW value).
!! The SCF driver divides by nothing but does require alpha > 0, so
!! adaptive_mix_get_alpha must never return a non-positive value.
program test_adaptive_mixing
    use fortuno_serial, only: execute_serial_cmd_app
    use lsda_constants, only: dp, INITIAL_MIX
    use adaptive_mixing
    implicit none

    !> Cap applied by up_mix: a candidate mix is only accepted if below this
    real(dp), parameter :: MIX_CAP = 0.999999999_dp

    !> Tolerance for comparing mixing parameters
    real(dp), parameter :: TOL = 1.0e-12_dp

    !> Band used by the in-band scenarios: [BAND_BOT, BAND_TOP]
    real(dp), parameter :: BAND_TOP = -0.9_dp
    real(dp), parameter :: BAND_BOT = -1.1_dp

    call execute_serial_cmd_app(get_adaptive_mixing_tests())

contains

    function get_adaptive_mixing_tests() result(tests)
        use fortuno_serial, only: test_list, test => serial_case_item
        type(test_list) :: tests

        tests = test_list([ &
            test("init_defaults", test_init_defaults), &
            test("up_mix_formula", test_up_mix_formula), &
            test("up_mix_respects_cap", test_up_mix_respects_cap), &
            test("converged_in_band", test_converged_in_band), &
            test("not_converged_before_count_sc_max", test_not_converged_before_count_sc_max), &
            test("dw_mix_on_monotonic_drift", test_dw_mix_on_monotonic_drift), &
            test("dw_mix_clamped_non_negative", test_dw_mix_clamped_non_negative), &
            test("get_alpha_always_positive", test_get_alpha_always_positive), &
            test("reset_keeps_mix", test_reset_keeps_mix) &
        ])
    end function get_adaptive_mixing_tests


    ! ------------------------------------------------------------------
    ! Scenario helpers
    ! ------------------------------------------------------------------

    !> Establish a non-degenerate energy band [BAND_BOT, BAND_TOP]
    !!
    !! After adaptive_mix_init the band is empty (Top = Bot). Three updates are
    !! needed to open it: one to seed it, one above Top, one below Bot. This
    !! leaves band_error ~ 0.18, i.e. far above any sensible tolerance, so the
    !! controller cannot declare convergence while inside this band.
    !!
    !! @param[inout] ctrl Adaptive mixing controller
    subroutine open_band(ctrl)
        type(adaptive_mix_t), intent(inout) :: ctrl

        call adaptive_mix_update(ctrl, -1.0_dp)     ! seeds Top = Bot = -1.0
        call adaptive_mix_update(ctrl, BAND_TOP)    ! pushes Top up
        call adaptive_mix_update(ctrl, BAND_BOT)    ! pushes Bot down
    end subroutine open_band

    !> Feed n strictly decreasing energies (a monotonic below-band drift)
    !!
    !! Each update lowers energy_bot, so the controller takes the below-band
    !! branch every time and count_bot grows; the 51st consecutive step trips the
    !! dw_mix guard. Because the sequence is monotonic, this behaves identically
    !! whether the post-dw_mix reset only zeroes the counters (current Fortran)
    !! or also collapses the band (C++ lsda_stop.cc:260-264).
    !!
    !! @param[inout] ctrl           Adaptive mixing controller
    !! @param[inout] energy         Running energy value, decreased in place
    !! @param[in]    n_steps        Number of updates to feed
    !! @param[inout] alpha_positive Set to .false. if alpha ever becomes non-positive
    subroutine drive_drift(ctrl, energy, n_steps, alpha_positive)
        type(adaptive_mix_t), intent(inout) :: ctrl
        real(dp), intent(inout) :: energy
        integer, intent(in) :: n_steps
        logical, intent(inout) :: alpha_positive
        integer :: i

        do i = 1, n_steps
            energy = energy - 0.001_dp
            call adaptive_mix_update(ctrl, energy)
            if (adaptive_mix_get_alpha(ctrl) <= 0.0_dp) alpha_positive = .false.
        end do
    end subroutine drive_drift

    !> Re-open the band around the current energy and trigger exactly one up_mix
    !!
    !! KNOWN DEVIATION FROM THE C++ REFERENCE (resolution belongs to T3):
    !! after UpMix/DwMix the C++ Convergencia::Reset() collapses the band
    !! (original/lsda_stop.cc:260-264 `Top = Bot = Old = New`), while
    !! src/convergence/adaptive_mixing.f90:128 and :159 only zero the counters.
    !! The state left behind by an up_mix therefore differs between the two
    !! semantics. To keep this test about the up_mix FORMULA only (which does
    !! match lsda_stop.cc:266-272), we normalise the state ourselves with
    !! adaptive_mix_reset and then re-open the band explicitly. Under either
    !! reset semantics this sequence performs exactly one up_mix.
    !!
    !! @param[inout] ctrl Adaptive mixing controller
    subroutine trigger_up_mix(ctrl)
        type(adaptive_mix_t), intent(inout) :: ctrl

        call adaptive_mix_reset(ctrl)               ! Top = Bot = current energy
        call adaptive_mix_update(ctrl, BAND_TOP)    ! pushes Top up
        call adaptive_mix_update(ctrl, BAND_BOT)    ! pushes Bot down
        call run_in_band_cycle(ctrl)                ! count_sc reaches count_sc_max
    end subroutine trigger_up_mix

    !> Feed count_sc_max in-band energies, which triggers exactly one up_mix
    !!
    !! The energies alternate between -1.00 and -1.05: both lie strictly inside
    !! [BAND_BOT, BAND_TOP], so count_sc keeps incrementing, while the relative
    !! step error (~5e-2) stays far above the tolerance, so the controller must
    !! choose up_mix instead of declaring convergence.
    !!
    !! @param[inout] ctrl Adaptive mixing controller
    subroutine run_in_band_cycle(ctrl)
        type(adaptive_mix_t), intent(inout) :: ctrl
        integer :: i

        do i = 1, ctrl%count_sc_max
            if (mod(i, 2) == 0) then
                call adaptive_mix_update(ctrl, -1.00_dp)
            else
                call adaptive_mix_update(ctrl, -1.05_dp)
            end if
        end do
    end subroutine run_in_band_cycle


    ! ------------------------------------------------------------------
    ! Tests
    ! ------------------------------------------------------------------

    !> Initialization must produce the documented starting state
    subroutine test_init_defaults()
        use fortuno_serial, only: check => serial_check
        type(adaptive_mix_t) :: ctrl

        call adaptive_mix_init(ctrl)

        call check(ctrl%iter == 0, "iter should start at 0")
        call check(ctrl%count_sc == 0, "count_sc should start at 0")
        call check(ctrl%count_sc_max == 10, "count_sc_max should default to 10")
        call check(abs(ctrl%mix - INITIAL_MIX) < TOL, "mix should start at INITIAL_MIX")
        call check(.not. ctrl%converged, "should not start converged")
        call check(abs(ctrl%tol - 1.0e-8_dp) < TOL, "tol should default to 1e-8")

        call adaptive_mix_init(ctrl, tol=1.0e-5_dp)
        call check(abs(ctrl%tol - 1.0e-5_dp) < TOL, "explicit tol should be honoured")
    end subroutine test_init_defaults


    !> up_mix must apply exactly mix + (1 - mix)/1.5
    subroutine test_up_mix_formula()
        use fortuno_serial, only: check => serial_check
        type(adaptive_mix_t) :: ctrl
        real(dp) :: expected

        call adaptive_mix_init(ctrl)
        call open_band(ctrl)

        call check(abs(ctrl%mix - INITIAL_MIX) < TOL, &
                   "Opening the band alone must not change mix")

        call run_in_band_cycle(ctrl)

        expected = INITIAL_MIX + (1.0_dp - INITIAL_MIX) / 1.5_dp
        call check(abs(ctrl%mix - expected) < TOL, &
                   "One in-band cycle should apply mix + (1-mix)/1.5")
        call check(ctrl%mix > INITIAL_MIX, "up_mix must increase mix (more conservative)")
        call check(.not. ctrl%converged, &
                   "Oscillating in a wide band must not count as converged")
        call check(ctrl%count_sc == 0, "Counters should be reset after up_mix")

        ! A second up_mix applies the same rule again, from the new mix.
        ! trigger_up_mix normalises the controller state first, so this assertion
        ! only depends on the formula, not on what reset does to the band.
        expected = ctrl%mix + (1.0_dp - ctrl%mix) / 1.5_dp
        call trigger_up_mix(ctrl)
        call check(abs(ctrl%mix - expected) < TOL, &
                   "Second up_mix should apply the formula from the updated mix")
    end subroutine test_up_mix_formula


    !> Repeated up_mix must never push mix to or past the 0.999999999 cap
    !!
    !! Each up_mix shrinks (1 - mix) by a factor of 3, so ~17 cycles from
    !! INITIAL_MIX = 0.95 would reach 1e-9 were the cap absent. 40 cycles is
    !! comfortably past that point.
    subroutine test_up_mix_respects_cap()
        use fortuno_serial, only: check => serial_check
        type(adaptive_mix_t) :: ctrl
        integer :: cycle_idx
        logical :: cap_respected, alpha_positive

        call adaptive_mix_init(ctrl)
        call open_band(ctrl)

        cap_respected = .true.
        alpha_positive = .true.

        do cycle_idx = 1, 40
            call trigger_up_mix(ctrl)
            if (ctrl%mix >= MIX_CAP) cap_respected = .false.
            if (adaptive_mix_get_alpha(ctrl) <= 0.0_dp) alpha_positive = .false.
        end do

        call check(cap_respected, "mix must always stay below the 0.999999999 cap")
        call check(alpha_positive, "alpha must stay strictly positive while up_mix saturates")
        call check(ctrl%mix < 1.0_dp, "mix must stay below 1")

        ! 40 forced up_mix calls do saturate the formula at the cap. This says
        ! something about up_mix alone (parity with lsda_stop.cc:266-272); it is
        ! NOT a statement that a real SCF run should ever get here. Whether the
        ! SCF can be driven into this corner depends on the reset semantics
        ! discussed in trigger_up_mix, which is T3 territory.
        call check(ctrl%mix > 0.999_dp, "40 forced up_mix calls should saturate the formula near the cap")
    end subroutine test_up_mix_respects_cap


    !> count_sc_max in-band iterations with a sub-tolerance error converge
    !!
    !! The band is opened with a width of 1e-12 relative, and the energies fed
    !! afterwards move by less than that, so both the step error and the band
    !! error stay below the default tol = 1e-8.
    subroutine test_converged_in_band()
        use fortuno_serial, only: check => serial_check
        type(adaptive_mix_t) :: ctrl
        integer :: i
        real(dp) :: e_base

        call adaptive_mix_init(ctrl)

        e_base = -1.0_dp
        call adaptive_mix_update(ctrl, e_base)                  ! Top = Bot = -1.0
        call adaptive_mix_update(ctrl, e_base - 1.0e-12_dp)     ! opens a 1e-12 wide band

        do i = 1, ctrl%count_sc_max
            call adaptive_mix_update(ctrl, e_base - 5.0e-13_dp)
        end do

        call check(ctrl%converged, &
                   "10 in-band iterations with relative error below tol should converge")
        call check(abs(ctrl%mix - INITIAL_MIX) < TOL, &
                   "Converging must not have triggered up_mix")
    end subroutine test_converged_in_band


    !> Convergence must not be declared before count_sc_max in-band iterations
    subroutine test_not_converged_before_count_sc_max()
        use fortuno_serial, only: check => serial_check
        type(adaptive_mix_t) :: ctrl
        integer :: i
        real(dp) :: e_base

        call adaptive_mix_init(ctrl)

        e_base = -1.0_dp
        call adaptive_mix_update(ctrl, e_base)
        call adaptive_mix_update(ctrl, e_base - 1.0e-12_dp)

        ! One iteration short of the threshold
        do i = 1, ctrl%count_sc_max - 1
            call adaptive_mix_update(ctrl, e_base - 5.0e-13_dp)
            call check(.not. ctrl%converged, &
                       "Must not converge before count_sc_max in-band iterations")
        end do

        call adaptive_mix_update(ctrl, e_base - 5.0e-13_dp)
        call check(ctrl%converged, "Must converge once the threshold is reached")
    end subroutine test_not_converged_before_count_sc_max


    !> A monotonic drift longer than 5*count_sc_max must trigger dw_mix
    !!
    !! The guard in adaptive_mix_update is count_bot > 5*count_sc_max, i.e. the
    !! 51st consecutive below-band energy. Counting the update that only seeds
    !! the band, that is the 52nd call.
    subroutine test_dw_mix_on_monotonic_drift()
        use fortuno_serial, only: check => serial_check
        type(adaptive_mix_t) :: ctrl
        integer :: i, drift_limit
        real(dp) :: expected

        call adaptive_mix_init(ctrl)
        drift_limit = 5 * ctrl%count_sc_max     ! 50

        call adaptive_mix_update(ctrl, -1.0_dp)  ! seeds Top = Bot

        ! 50 consecutive decreases: still one short of the trigger
        do i = 1, drift_limit
            call adaptive_mix_update(ctrl, -1.0_dp - 0.001_dp * real(i, dp))
        end do
        call check(abs(ctrl%mix - INITIAL_MIX) < TOL, &
                   "50 monotonic iterations must not yet trigger dw_mix")

        ! The 51st consecutive decrease trips the guard
        expected = INITIAL_MIX - (1.0_dp - INITIAL_MIX) * 1.9_dp
        call adaptive_mix_update(ctrl, -1.0_dp - 0.001_dp * real(drift_limit + 1, dp))

        call check(abs(ctrl%mix - expected) < TOL, &
                   "51 monotonic iterations should apply mix - (1-mix)*1.9")
        call check(ctrl%mix < INITIAL_MIX, "dw_mix must decrease mix (more aggressive)")
        call check(ctrl%count_bot == 0, "Counters should be reset after dw_mix")
        call check(.not. ctrl%converged, "A monotonic drift must not be reported as converged")
    end subroutine test_dw_mix_on_monotonic_drift


    !> Repeated dw_mix follows the formula and is then blocked by the mix > 0.35 guard
    !!
    !! From INITIAL_MIX = 0.95 the successive dw_mix values are 0.855 and 0.5795.
    !! The third application would give -0.21945.
    !!
    !! KNOWN DEVIATION FROM THE C++ REFERENCE (resolution belongs to T3):
    !! original/lsda_stop.cc:274-276 applies `Mix = Mix - (1.0-Mix)*1.9;` with no
    !! clamp, so the C++ really does reach Mix = -0.21945 and uses it as
    !! deliberate over-relaxation (weight 1.219 on the new potential,
    !! original/lsdaks.cc:633 and :641). src/convergence/adaptive_mixing.f90:154-156
    !! clamps mix to 0 instead, and adaptive_mixing.f90:186-191 clamps alpha to 1.
    !! This test therefore asserts only what both variants share: the two
    !! applications that stay inside the guard follow the formula exactly, and
    !! once mix has dropped to or below 0.35 the guard freezes it, whatever the
    !! third application produced.
    subroutine test_dw_mix_clamped_non_negative()
        use fortuno_serial, only: check => serial_check
        type(adaptive_mix_t) :: ctrl
        real(dp) :: energy, expected, mix_after_third
        logical :: alpha_positive

        call adaptive_mix_init(ctrl)
        energy = -1.0_dp
        call adaptive_mix_update(ctrl, energy)      ! seeds Top = Bot

        alpha_positive = .true.

        ! First dw_mix: 0.95 -> 0.855
        expected = INITIAL_MIX - (1.0_dp - INITIAL_MIX) * 1.9_dp
        call drive_drift(ctrl, energy, 51, alpha_positive)
        call check(abs(ctrl%mix - expected) < TOL, "First dw_mix applies mix - (1-mix)*1.9")

        ! Second dw_mix: 0.855 -> 0.5795 (0.855 is still above the 0.35 guard)
        expected = ctrl%mix - (1.0_dp - ctrl%mix) * 1.9_dp
        call drive_drift(ctrl, energy, 51, alpha_positive)
        call check(abs(ctrl%mix - expected) < TOL, "Second dw_mix applies the formula again")
        call check(ctrl%mix > 0.35_dp, "Precondition: 0.5795 is still above the mix > 0.35 guard")

        ! Third dw_mix fires and takes mix below the guard
        call drive_drift(ctrl, energy, 51, alpha_positive)
        mix_after_third = ctrl%mix
        call check(mix_after_third <= 0.35_dp, "Third dw_mix must take mix below the 0.35 guard")

        ! From here on the guard must block every further dw_mix
        call drive_drift(ctrl, energy, 400, alpha_positive)
        call check(abs(ctrl%mix - mix_after_third) < TOL, &
                   "mix must stay frozen once the mix > 0.35 guard blocks dw_mix")
        call check(alpha_positive, "alpha must stay strictly positive throughout dw_mix")
    end subroutine test_dw_mix_clamped_non_negative


    !> adaptive_mix_get_alpha must always return a strictly positive value
    !!
    !! The real requirement on alpha is alpha > 0: linear_mixing rejects a
    !! non-positive weight for the new potential, so mix = 1 (and any mix > 1)
    !! must not be handed through as alpha <= 0.
    !!
    !! KNOWN DEVIATION FROM THE C++ REFERENCE (resolution belongs to T3):
    !! the UPPER clamp in src/convergence/adaptive_mixing.f90:186-191
    !! (alpha > 1 -> alpha = 1) has no counterpart in the C++. There
    !! original/lsda_stop.cc:274-276 applies `Mix = Mix - (1.0-Mix)*1.9;` with no
    !! clamp, so Mix legitimately becomes negative (-0.21945 from Mix = 0.95) and
    !! original/lsdaks.cc:633 and :641 then use `v_eff = Mix*v_eff + (1-Mix)*v_novo`
    !! with weight 1.21945 on the new potential: deliberate over-relaxation.
    !! This test therefore does NOT assert alpha <= 1; if T3 restores the C++
    !! over-relaxation, only the lower bound below must keep holding.
    subroutine test_get_alpha_always_positive()
        use fortuno_serial, only: check => serial_check
        type(adaptive_mix_t) :: ctrl
        real(dp) :: alpha

        call adaptive_mix_init(ctrl)
        alpha = adaptive_mix_get_alpha(ctrl)
        call check(abs(alpha - (1.0_dp - INITIAL_MIX)) < TOL, &
                   "alpha should be 1 - mix in the ordinary case")

        ! mix = 1 would give alpha = 0, which linear_mixing rejects
        ctrl%mix = 1.0_dp
        alpha = adaptive_mix_get_alpha(ctrl)
        call check(alpha > 0.0_dp, "alpha must stay positive when mix = 1")

        ! mix > 1 would give a negative alpha
        ctrl%mix = 2.5_dp
        alpha = adaptive_mix_get_alpha(ctrl)
        call check(alpha > 0.0_dp, "alpha must stay positive when mix > 1")

        ! mix < 0 is the over-relaxation regime of the C++ (alpha = 1 - mix > 1).
        ! Only positivity is asserted here; see the deviation note above.
        ctrl%mix = -3.0_dp
        alpha = adaptive_mix_get_alpha(ctrl)
        call check(alpha > 0.0_dp, "alpha must stay positive when mix < 0")
    end subroutine test_get_alpha_always_positive


    !> adaptive_mix_reset clears the counters and band but keeps the tuned mix
    subroutine test_reset_keeps_mix()
        use fortuno_serial, only: check => serial_check
        type(adaptive_mix_t) :: ctrl
        real(dp) :: mix_before

        call adaptive_mix_init(ctrl)
        call open_band(ctrl)
        call run_in_band_cycle(ctrl)

        mix_before = ctrl%mix
        call check(mix_before > INITIAL_MIX, "Precondition: mix should have been tuned up")

        call adaptive_mix_reset(ctrl)

        call check(abs(ctrl%mix - mix_before) < TOL, "reset must preserve the tuned mix")
        call check(ctrl%count_sc == 0, "reset must clear count_sc")
        call check(ctrl%count_top == 0, "reset must clear count_top")
        call check(ctrl%count_bot == 0, "reset must clear count_bot")
        call check(abs(ctrl%energy_top - ctrl%energy_new) < TOL, &
                   "reset must collapse Top onto the current energy")
        call check(abs(ctrl%energy_bot - ctrl%energy_new) < TOL, &
                   "reset must collapse Bot onto the current energy")
    end subroutine test_reset_keeps_mix

end program test_adaptive_mixing
