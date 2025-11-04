module bethe_equations
    use lsda_constants, only: dp, PI, TWOPI, U_SMALL
    implicit none
    private

    public :: theta, Theta_capital
    public :: dtheta_dx, dTheta_capital_dx
    public :: initialize_quantum_numbers
    public :: compute_residual
    public :: compute_jacobian

contains

    !> charge-spin scattering function θ(x, U)
    !!
    !! θ(x, U) = 2·arctan(2x/U)
    !!
    !! @param[in] x   Difference in rapidities (k - Λ)
    !! @param[in] U   Hubbard interaction
    !! @return        Value of θ(x, U)
    !!
    !! @note For U → 0, θ → π·sign(x)
    function theta(x, U) result(res)
        real(dp), intent(in) :: x, U
        real(dp) :: res

        if (abs(U) < U_SMALL) then
            res = PI * sign(1.0_dp, x)
        else
            res = 2.0_dp * atan( (2.0_dp * x) / U )
        end if
    end function theta

    !> spin-spin scattering function Θ(x, U)
    !!
    !! Θ(x, U) = 2·arctan(x/U)
    !!
    !! @param[in] x   Difference in rapidities (Λ - Λ')
    !! @param[in] U   Hubbard interaction
    !! @return        Value of Θ(x, U)
    !!
    !! @note For U → 0, Θ → π·sign(x)
    function Theta_capital(x, U) result(res)
        real(dp), intent(in) :: x, U
        real(dp) :: res

        if (abs(U) < U_SMALL) then
            res = PI * sign(1.0_dp, x)
        else
            res = 2.0_dp * atan2( x, U )
        end if
    end function Theta_capital

    !> Derivative of the charge-spin scattering function θ with respect to x
    !!
    !! dθ/dx = 4U / (U² + 4x²)
    !!
    !! @param[in] x   Difference in rapidities
    !! @param[in] U   Hubbard interaction
    !! @return        Value of dθ/dx
    !!
    !! @note For U → 0, dθ/dx → 0 (except at x=0 where it is singular)
    function dtheta_dx(x, U) result(res)
        real(dp), intent(in) :: x, U
        real(dp) :: res, denom

        denom = U**2 + 4.0_dp * x**2

        if (abs(U) < U_SMALL) then
            res = 0.0_dp
        else
            res = (4.0_dp * U) / denom
        end if
    end function dtheta_dx

    !> Derivative of the function Θ with respect to x
    !!
    !! dΘ/dx = 2U / (U² + x²)
    !!
    !! @param[in] x   Difference in rapidities
    !! @param[in] U   Hubbard interaction
    !! @return        Value of dΘ/dx
    !!
    !! @note For U → 0, dΘ/dx → 0 (except at x=0 where it is singular)
    function dTheta_capital_dx(x, U) result(res)
        real(dp), intent(in) :: x, U
        real(dp) :: res, denom

        denom = U**2 + x**2

        if (abs(U) < U_SMALL) then
            res = 0.0_dp
        else
            res = (2.0_dp * U) / denom
        end if
    end function dTheta_capital_dx

    !> Initialize quantum numbers for the ground state
    !!
    !! Computes I_k and J_α using Fermi distribution:
    !!   I_k = k - (N_up + 1)/2,    k = 1, ..., N_up
    !!   J_α = α - (M + 1)/2,       α = 1, ..., M
    !!
    !! @param[in]  N_up   Number of spin-up electrons
    !! @param[in]  M      Number of spin rapidities (= N_down)
    !! @param[out] I      Array of charge quantum numbers [1:N_up]
    !! @param[out] J      Array of spin quantum numbers [1:M]
    !!
    !! @note Arrays I and J must be pre-allocated by the caller
    !!
    !! Example:
    !!   N_up = 5  →  I = [-2, -1, 0, 1, 2]
    !!   M = 3     →  J = [-1, 0, 1]
    subroutine initialize_quantum_numbers(Nup, M, I, J_capital)
        integer, intent(in) :: Nup, M
        real(dp), intent(out) :: I(:), J_capital(:)
        real(dp) :: offset_I, offset_J
        integer :: j, alpha

        offset_I = (Nup + 1) / 2.0_dp
        offset_J = (M + 1) / 2.0_dp

        do j = 1, Nup
            I(j) = real(j, dp) - offset_I
        end do

        do alpha = 1, M
            J_capital(alpha) = real(alpha, dp) - offset_J
        end do
    end subroutine initialize_quantum_numbers

    !> Computes the residual vector F(x) for the Lieb-Wu equations of the 1D Hubbard model.
    !!
    !! This function evaluates the nonlinear system F(x) = 0, where x = [k, Λ] contains
    !! the charge and spin rapidities. The residual measures how far the current guess
    !! is from satisfying the Bethe Ansatz equations.
    !!
    !! The system consists of two sets of equations:
    !!
    !! **Charge equations** (for j = 1, ..., N↑):
    !! \[ F_j^k = k_j - \frac{2\pi}{L} I_j - \frac{1}{L} \sum_{\alpha=1}^{M} \theta(k_j - \Lambda_\alpha, U) \]
    !!
    !! **Spin equations** (for α = 1, ..., M = N↓):
    !! \[ F_\alpha^\Lambda = \frac{2\pi}{L} J_\alpha - \sum_{j=1}^{N_\uparrow} \theta(\Lambda_\alpha - k_j, U) 
    !!    + \sum_{\beta \neq \alpha} \Theta(\Lambda_\alpha - \Lambda_\beta, U) \]
    !!
    !! @param[in] k         Charge rapidities (size N↑)
    !! @param[in] Lambda    Spin rapidities (size M = N↓)
    !! @param[in] I         Charge quantum numbers (size N↑, can be semi-integers)
    !! @param[in] J_capital Spin quantum numbers (size M, can be semi-integers)
    !! @param[in] L         Number of lattice sites
    !! @param[in] U         Hubbard interaction strength
    !! @return F            Residual vector (size N↑ + M)
    !!
    !! @note The residual is zero (F = 0) when the rapidities satisfy the Lieb-Wu equations.
    !! @note Quantum numbers I and J should be initialized using `initialize_quantum_numbers`.
    !!
    !! @see theta, Theta_capital, initialize_quantum_numbers
    !!
    function compute_residual(k, Lambda, I, J_capital, L, U) result(F)
        real(dp), intent(in) :: k(:), Lambda(:), I(:), J_capital(:), U
        integer, intent(in) :: L
        integer :: j, alpha, beta, Nup, M
        real(dp) :: summ, summ1, summ2
        real(dp) :: F(size(k) + size(Lambda))

        Nup = size(k)
        M = size(Lambda)
        
        ! Charge equations: F^k
        do j = 1, Nup
            summ = 0.0_dp
            do alpha = 1, M
                summ = summ + theta(k(j) - Lambda(alpha), U)
            end do

            F(j) = k(j) - (TWOPI * I(j))/L - summ/L
        end do

        ! Spin equations: F^Lambda
        do alpha = 1, M
            summ1 = 0.0_dp
            summ2 = 0.0_dp

            do j = 1, Nup
                summ1 = summ1 + theta(Lambda(alpha) - k(j), U)
            end do

            do beta = 1, M
                if (beta /= alpha) then
                    summ2 = summ2 + Theta_capital(Lambda(alpha) - Lambda(beta), U)
                end if
            end do

            F(Nup + alpha) = (TWOPI * J_capital(alpha))/L - summ1 + summ2
        end do
    end function compute_residual

    !> Computes the Jacobian matrix J = ∂F/∂x for the Lieb-Wu equations.
    !!
    !! The Jacobian is a square matrix of size (N↑ + M) × (N↑ + M) containing all
    !! partial derivatives of the residual F with respect to the rapidities x = [k, Λ].
    !! It has a 4-block structure:
    !!
    !! \[ J = \begin{bmatrix}
    !!   \frac{\partial F^k}{\partial k} & \frac{\partial F^k}{\partial \Lambda} \\
    !!   \frac{\partial F^\Lambda}{\partial k} & \frac{\partial F^\Lambda}{\partial \Lambda}
    !! \end{bmatrix} \]
    !!
    !! **Block A** (N↑ × N↑, diagonal):
    !! \[ J_{jj} = 1 - \frac{1}{L} \sum_{\alpha=1}^{M} \frac{4U}{U^2 + 4(k_j - \Lambda_\alpha)^2} \]
    !! \[ J_{ji} = 0 \quad \text{for } i \neq j \]
    !!
    !! **Block B** (N↑ × M):
    !! \[ J_{j,N_\uparrow+\beta} = \frac{1}{L} \cdot \frac{4U}{U^2 + 4(k_j - \Lambda_\beta)^2} \]
    !!
    !! **Block C** (M × N↑):
    !! \[ J_{N_\uparrow+\alpha,i} = \frac{4U}{U^2 + 4(\Lambda_\alpha - k_i)^2} \]
    !!
    !! **Block D** (M × M):
    !! \[ J_{N_\uparrow+\alpha,N_\uparrow+\gamma} = \begin{cases}
    !!   \sum_j \frac{4U}{U^2 + 4(\Lambda_\alpha - k_j)^2} + \sum_{\beta \neq \alpha} \frac{2U}{U^2 + (\Lambda_\alpha - \Lambda_\beta)^2} & \text{if } \gamma = \alpha \\
    !!   -\frac{2U}{U^2 + (\Lambda_\alpha - \Lambda_\gamma)^2} & \text{if } \gamma \neq \alpha
    !! \end{cases} \]
    !!
    !! @param[in] k       Charge rapidities (size N↑)
    !! @param[in] Lambda  Spin rapidities (size M = N↓)
    !! @param[in] L       Number of lattice sites
    !! @param[in] U       Hubbard interaction strength
    !! @return Jacobian   Jacobian matrix (size (N↑+M) × (N↑+M))
    !!
    !! @note This Jacobian is used in Newton-Raphson method to solve J·Δx = -F.
    !! @note The matrix is dense and fully populated (no sparsity exploitation).
    !!
    !! @see compute_residual, dtheta_dx, dTheta_dx
    !!
    function compute_jacobian(k, Lambda, L, U) result(Jacobian)
        real(dp), intent(in) :: k(:), Lambda(:), U
        integer, intent(in) :: L
        integer :: Nup, M, i, j, alpha, beta, gamma
        real(dp) :: Jacobian(size(k) + size(Lambda), size(k) + size(Lambda))
        real(dp) :: summ1, summ2

        Nup = size(k)
        M = size(Lambda)

        !! Jacobian
        !!
        !!            | Block A  Block B |
        !! Jacobian = |                  |
        !!            | Block C  Block D |
        
        !! Block A: dF^k_j/dk_i
        do j = 1, Nup
            do i = 1, Nup
                if (i == j) then
                    summ1 = 0.0_dp
                    do alpha = 1, M
                        summ1 = summ1 + (4.0_dp * U) / (U**2 + 4.0_dp * (k(j) - Lambda(alpha))**2)
                    end do
                    Jacobian(j, i) = 1.0_dp - (1.0_dp / L) * summ1
                else
                    Jacobian(j, i) = 0.0_dp
                end if
            end do
        end do

        !! Block B: dF^k_j/dLambda_beta
        do j = 1, Nup
            do beta = 1, M
                Jacobian(j, Nup + beta) = (1.0_dp / L) * (4.0_dp * U) / (U**2 + 4.0_dp * (k(j) - Lambda(beta))**2)
            end do
        end do

        !! Block C: dF^Lambda_alpha/dk_i
        do alpha = 1, M
            do i = 1, Nup
                Jacobian(Nup + alpha, i) = (4.0_dp * U) / (U**2 + 4.0_dp * (Lambda(alpha) - k(i))**2)
            end do
        end do

        !! Block D: dF^Lambda_alpha/dLambda_gamma
        do alpha = 1, M
            do gamma = 1, M
                if (alpha == gamma) then
                    summ1 = 0.0_dp
                    do j = 1, Nup
                        summ1 = summ1 + (4.0_dp * U) / (U**2 + 4.0_dp * (Lambda(alpha) - k(j))**2)
                    end do

                    summ2 = 0.0_dp
                    do beta = 1, M
                        if (beta /= alpha) then
                            summ2 = summ2 + (2.0_dp * U) / (U**2 + (Lambda(alpha) - Lambda(beta))**2)
                        end if
                    end do

                    Jacobian(Nup + alpha, Nup + gamma) = summ1 + summ2
                else
                    Jacobian(Nup + alpha, Nup + gamma) = - (2.0_dp * U) / (U**2 + (Lambda(alpha) - Lambda(gamma))**2)
                end if
            end do
        end do

    end function
end module bethe_equations