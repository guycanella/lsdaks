!> Resolution of the random seed of the stochastic external potentials
!!
!! The disorder generators (`potential_random_uniform`,
!! `potential_random_gaussian`, `potential_impurity_random`) accept a seed and
!! treat any negative value as "draw one from the system clock", by calling
!! `random_seed()` with no arguments. That is the right *default* for a
!! disordered system - the physical observable is the average over realisations,
!! so every run must sample a new one - but it made each run irreproducible even
!! for the person who ran it: `random_seed()` does not report what it drew, so
!! the realisation could not be recovered, and two identical invocations gave
!! measurably different ground-state energies (E/site = -3.204508987 and
!! -3.220780358 for the same input file).
!!
!! This module closes that gap by drawing the seed EXPLICITLY, one level above
!! the generators: the caller resolves the requested seed into a non-negative
!! integer, hands that integer to the generator (which then takes its
!! deterministic `seed >= 0` branch) and records it in the output provenance.
!! The run stays random by default and becomes replayable a posteriori by
!! feeding the recorded value back as `pot_seed`.
module potential_seed
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private

    public :: resolve_random_seed

    !> Exclusive upper bound of the drawn seeds.
    !!
    !! Keeps the result comfortably inside the default integer range, so that
    !! the value can be printed, stored in an input file and read back without
    !! any risk of overflow.
    integer, parameter :: SEED_MODULUS = 1000000000

    !> Odd increment applied to successive draws inside one process.
    !!
    !! `system_clock` can return the same count for two calls that are close
    !! enough together, which would hand the same realisation to two different
    !! potentials of the same run. Advancing by a fixed odd stride guarantees
    !! distinct seeds without pretending to add entropy.
    integer, parameter :: DRAW_STRIDE = 7919

    !> Number of seeds already drawn in this process (see DRAW_STRIDE).
    integer, save :: n_drawn = 0

contains

    !> Turn a requested seed into the non-negative seed actually to be used
    !!
    !! @param[in]  requested_seed Seed as given by the user; any negative value
    !!                            (conventionally -1) means "draw one"
    !! @param[out] effective_seed Seed to hand to the potential generator and to
    !!                            record in the output provenance. Equal to
    !!                            `requested_seed` when that is non-negative,
    !!                            otherwise a clock-derived value in
    !!                            [0, SEED_MODULUS).
    subroutine resolve_random_seed(requested_seed, effective_seed)
        integer, intent(in) :: requested_seed
        integer, intent(out) :: effective_seed

        integer(int64) :: clock_count

        if (requested_seed >= 0) then
            effective_seed = requested_seed
            return
        end if

        call system_clock(count=clock_count)

        ! abs() before the modulus: a platform is free to return a negative
        ! count, and a negative seed would be interpreted by the generators as
        ! "draw one", i.e. exactly the irreproducibility this routine removes.
        effective_seed = int(mod(abs(clock_count), int(SEED_MODULUS, int64)))

        ! Distinct seeds for successive draws within the same process, even if
        ! the clock did not advance between them.
        n_drawn = n_drawn + 1
        effective_seed = modulo(effective_seed + n_drawn * DRAW_STRIDE, SEED_MODULUS)
    end subroutine resolve_random_seed

end module potential_seed
