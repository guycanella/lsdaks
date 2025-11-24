!> Module for potential factory pattern
!!
!! Provides a unified interface to create any type of potential dynamically.
!! This simplifies potential creation in simulation codes by allowing
!! string-based selection of potential types.
!!
!! Supported potentials:
!! - "uniform": Constant potential
!! - "harmonic": Harmonic trap
!! - "impurity_single": Single point impurity
!! - "impurity_multiple": Multiple impurities
!! - "impurity_random": Random impurities with concentration
!! - "random_uniform": Random disorder (uniform distribution)
!! - "random_gaussian": Random disorder (Gaussian distribution)
!! - "barrier_single": Single rectangular barrier
!! - "barrier_double": Double barrier (quantum well)
!! - "quasiperiodic": Aubry-André-Harper quasiperiodic potential
!!
!! Usage:
!!   call create_potential("harmonic", params, L, seed, V, ierr)
module potential_factory
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_SUCCESS, ERROR_INVALID_INPUT
    use potential_uniform, only: apply_potential_uniform
    use potential_harmonic, only: apply_potential_harmonic
    use potential_impurity, only: potential_impurity_single, &
                                   potential_impurity_multiple, &
                                   potential_impurity_random
    use potential_random, only: potential_random_uniform, &
                                potential_random_gaussian
    use potential_barrier, only: potential_barrier_single, &
                                 potential_barrier_double
    use potential_quasiperiodic, only: apply_potential_quasiperiodic
    implicit none
    private
    
    public :: create_potential
    public :: get_potential_info

contains

    !> Factory function to create any potential type
    !!
    !! Creates a potential based on a string identifier and parameter array.
    !! This provides a unified interface for potential creation.
    !!
    !! Parameter array format (params):
    !! - "uniform": params(1) = V0
    !! - "harmonic": params(1) = k
    !! - "impurity_single": params(1) = V_imp, params(2) = i_imp (as real)
    !! - "random_uniform": params(1) = W
    !! - "random_gaussian": params(1) = sigma
    !! - "barrier_single": params(1) = V_bar, params(2) = i_start, params(3) = i_end
    !! - "barrier_double": params(1) = V_bar, params(2) = L_bar, params(3) = V_well, params(4) = L_well
    !! - "quasiperiodic": params(1) = lambda, params(2) = beta, params(3) = phi
    !!
    !! @param[in]  potential_type  String identifying potential type
    !! @param[in]  params          Parameter array (size depends on potential type)
    !! @param[in]  L               Number of lattice sites
    !! @param[in]  seed            Random seed (for random potentials, use -1 for system time)
    !! @param[out] V               Potential array V(i) for i = 1..L
    !! @param[out] ierr            Error flag
    !!
    !! @note For multiple impurities, use potential_impurity_multiple directly
    !! @note For random impurities, use potential_impurity_random directly
    subroutine create_potential(potential_type, params, L, seed, V, ierr)
        character(len=*), intent(in) :: potential_type
        real(dp), intent(in) :: params(:)
        integer, intent(in) :: L, seed
        real(dp), dimension(L), intent(out) :: V
        integer, intent(out) :: ierr

        ! Initialize
        V = 0.0_dp
        ierr = ERROR_SUCCESS

        ! Select potential type
        select case (trim(adjustl(potential_type)))
        
        case ("uniform")
            ! Uniform potential: params(1) = V0
            if (size(params) < 1) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
            call apply_potential_uniform(params(1), L, V, ierr)

        case ("harmonic")
            ! Harmonic trap: params(1) = k
            if (size(params) < 1) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
            call apply_potential_harmonic(params(1), L, V, ierr)
            
        case ("impurity_single")
            ! Single impurity: params(1) = V_imp, params(2) = i_imp
            if (size(params) < 2) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
            call potential_impurity_single(params(1), int(params(2)), L, V, ierr)
            
        case ("random_uniform")
            ! Random uniform: params(1) = W
            if (size(params) < 1) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
            call potential_random_uniform(params(1), L, seed, V, ierr)
            
        case ("random_gaussian")
            ! Random Gaussian: params(1) = sigma
            if (size(params) < 1) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
            call potential_random_gaussian(params(1), L, seed, V, ierr)
            
        case ("barrier_single")
            ! Single barrier: params(1) = V_bar, params(2) = i_start, params(3) = i_end
            if (size(params) < 3) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
            call potential_barrier_single(params(1), int(params(2)), int(params(3)), L, V, ierr)
            
        case ("barrier_double")
            ! Double barrier: params(1) = V_bar, params(2) = L_bar, params(3) = V_well, params(4) = L_well
            ! Matches C++ double_barrier(Na, Vb, Lb, Vwell, Lwell, v_ext)
            if (size(params) < 4) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
            call potential_barrier_double(params(1), params(2), params(3), params(4), L, V, ierr)

        case ("quasiperiodic")
            ! Quasiperiodic: params(1) = lambda, params(2) = beta, params(3) = phi
            if (size(params) < 3) then
                ierr = ERROR_INVALID_INPUT
                return
            end if
            call apply_potential_quasiperiodic(params(1), params(2), params(3), L, V, ierr)

        case default
            ierr = ERROR_INVALID_INPUT
        end select
    end subroutine create_potential

    !> Get information about a potential type
    !!
    !! Returns a descriptive string explaining a potential type and its parameters.
    !!
    !! @param[in]  potential_type  String identifying potential type
    !! @return     info            Descriptive information string
    function get_potential_info(potential_type) result(info)
        character(len=*), intent(in) :: potential_type
        character(len=512) :: info

        select case (trim(adjustl(potential_type)))
        case ("uniform")
            info = "Uniform potential: V(i) = V0. Parameters: [V0]"
        case ("harmonic")
            info = "Harmonic trap: V(i) = k*(i-center)^2. Parameters: [k]"
        case ("impurity_single")
            info = "Single impurity: V(i) = V_imp at i_imp. Parameters: [V_imp, i_imp]"
        case ("impurity_multiple")
            info = "Multiple impurities: Use potential_impurity_multiple directly"
        case ("impurity_random")
            info = "Random impurities: Use potential_impurity_random directly"
        case ("random_uniform")
            info = "Random uniform disorder: V(i) ~ Uniform[-W, +W]. Parameters: [W]"
        case ("random_gaussian")
            info = "Random Gaussian: V(i) ~ N(0, sigma^2). Parameters: [sigma]"
        case ("barrier_single")
            info = "Single barrier: V(i) = V_bar in [i_start, i_end]. Parameters: [V_bar, i_start, i_end]"
        case ("barrier_double")
            info = "Double barrier: Two barriers with well. Parameters: [V_bar, L_bar, V_well, L_well]"
        case ("quasiperiodic")
            info = "Quasiperiodic AAH: V(i) = lambda*cos(2*pi*beta*i + phi). Parameters: [lambda, beta, phi]"
        case default
            info = "Unknown potential type"
        end select
    end function get_potential_info

end module potential_factory