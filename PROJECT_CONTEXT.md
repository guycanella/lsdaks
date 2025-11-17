# PROJECT_CONTEXT.md

## 📚 Visão Geral

**LSDA-Hubbard-Fortran** é uma reimplementação moderna em Fortran do código LSDA (Local Spin Density Approximation) para cálculos de Bethe Ansatz-LSDA no modelo de Hubbard 1D. Este projeto migra um código legado em C++ para Fortran moderno (2008/2018), com foco em:

- **Arquitetura limpa e modular**
- **Testes unitários abrangentes**
- **Performance otimizada** (LAPACK, OpenMP)
- **Código bem documentado e mantível**

---

## 🎯 Motivação

### Por que Fortran?

1. **Adequação natural para física computacional**
   - Arrays multidimensionais nativos
   - Álgebra linear de alto desempenho (LAPACK/BLAS)
   - Tipos numéricos precisos e bem definidos

2. **Modernidade** (Fortran 2008/2018)
   - Módulos e namespaces
   - Tipos derivados (similar a structs/classes)
   - Alocação automática
   - Interoperabilidade com C

3. **Eliminação de problemas do código original**
   - Arrays 1-indexed manualmente (offset `+1`)
   - Gerenciamento manual de memória (`new[]`/`delete[]`)
   - Arquitetura monolítica
   - Falta de testes
   - Bethe Ansatz não implementado (apenas leitura de tabelas)

---

## 🏗️ Arquitetura do Projeto

### Estrutura de Diretórios

```
lsda-hubbard/
├── fpm.toml                    # Configuração do Fortran Package Manager
├── README.md                   # Documentação de uso
├── PROJECT_CONTEXT.md          # Este arquivo (contexto técnico)
├── LICENSE                     # Licença do projeto
│
├── src/                        # Código-fonte principal
│   ├── lsda_main.f90          # Programa principal
│   │
│   ├── types/                  # ✅ COMPLETO
│   │   ├── lsda_types.f90     # Tipos principais (SystemParams, State, etc)
│   │   ├── lsda_constants.f90 # Constantes físicas e numéricas
│   │   └── lsda_errors.f90    # ✅ COMPLETO - Sistema de erros centralizado
│   │
│   ├── io/                     # 🔜 TODO
│   │   ├── input_parser.f90   # Parse de argumentos e arquivos
│   │   ├── output_writer.f90  # Escrita de resultados
│   │   └── logger.f90         # Sistema de logging
│   │
│   ├── bethe_ansatz/          # ✅ COMPLETO (Fases 1 & 2)
│   │   ├── bethe_equations.f90      # ✅ COMPLETO - Equações de Lieb-Wu
│   │   ├── nonlinear_solvers.f90    # ✅ COMPLETO - Newton-Raphson
│   │   ├── continuation.f90         # ✅ COMPLETO - Sweep em U
│   │   ├── table_io.f90             # ✅ COMPLETO - I/O tabelas (ASCII/binário)
│   │   ├── bethe_tables.f90         # ✅ COMPLETO - Geração de tabelas XC
│   │   └── table_manager.f90        # 🔜 TODO - Cache híbrido (opcional)
│   │
│   ├── xc_functional/         # ✅ COMPLETO (Fase 3)
│   │   ├── spline2d.f90       # ✅ COMPLETO - Interpolação bicúbica 2D
│   │   └── xc_lsda.f90        # ✅ COMPLETO - Interface exc, Vxc_up, Vxc_dw
│   │
│   ├── potentials/            # ✅ COMPLETO (Fase 4)
│   │   ├── potential_uniform.f90      # ✅ COMPLETO - Potencial uniforme V(i) = V₀
│   │   ├── potential_harmonic.f90     # ✅ COMPLETO - Armadilha harmônica
│   │   ├── potential_impurity.f90     # ✅ COMPLETO - Impurezas (single/multiple/random)
│   │   ├── potential_random.f90       # ✅ COMPLETO - Desordem (uniform/Gaussian)
│   │   ├── potential_barrier.f90      # ✅ COMPLETO - Barreiras (single/double)
│   │   ├── potential_quasiperiodic.f90 # ✅ COMPLETO - Aubry-André-Harper (AAH)
│   │   └── potential_factory.f90      # ✅ COMPLETO - Factory pattern (7 tipos)
│   │
│   ├── hamiltonian/           # ✅ COMPLETO (Fase 5 - partial)
│   │   ├── hamiltonian_builder.f90 # ✅ COMPLETO - Tight-binding com Veff
│   │   ├── boundary_conditions.f90 # ✅ COMPLETO - Open, periodic, twisted
│   │   └── symmetry.f90            # 🔜 TODO - Exploração de simetria de paridade (opcional)
│   │
│   ├── diagonalization/       # ✅ COMPLETO (Fase 5)
│   │   ├── lapack_wrapper.f90      # ✅ COMPLETO - Wrappers DSYEVD/ZHEEVD
│   │   └── degeneracy_handler.f90  # ✅ COMPLETO - QR/Gram-Schmidt
│   │
│   ├── density/               # 🔄 EM PROGRESSO (Fase 6 - 20% completo)
│   │   └── density_calculator.f90  # ✅ COMPLETO - Densidade de autoestados KS
│   │
│   ├── convergence/           # 🔜 TODO (Fase 6 - 80% restante)
│   │   └── scf_mixer.f90      # Mixing linear, Broyden, Anderson
│   │
│   └── kohn_sham/             # 🔜 TODO (Fase 6)
│       └── ks_cycle.f90       # Loop SCF completo
│
├── app/                        # 🔄 EM PROGRESSO
│   ├── main.f90               # Ponto de entrada (placeholder)
│   └── convert_tables.f90     # ✅ COMPLETO - Utilitário conversão tabelas
│
├── test/                       # 🔄 EM PROGRESSO (164 testes, 100% passando)
│   ├── test_bethe_equations.f90       # ✅ COMPLETO - 17 testes
│   ├── test_nonlinear_solvers.f90     # ✅ COMPLETO - 9 testes
│   ├── test_continuation.f90          # ✅ COMPLETO - 5 testes
│   ├── test_table_io.f90              # ✅ COMPLETO - 10 testes
│   ├── test_bethe_tables.f90          # ✅ COMPLETO - 6 testes
│   ├── test_spline2d.f90              # ✅ COMPLETO - 5 testes
│   ├── test_xc_lsda.f90               # ✅ COMPLETO - 6 testes
│   ├── test_potentials.f90            # ✅ COMPLETO - 21 testes
│   ├── test_lsda_errors.f90           # ✅ COMPLETO - 13 testes
│   ├── test_boundary_conditions.f90   # ✅ COMPLETO - 17 testes
│   ├── test_hamiltonian_builder.f90   # ✅ COMPLETO - 18 testes
│   ├── test_lapack_wrapper.f90        # ✅ COMPLETO - 18 testes
│   ├── test_degeneracy_handler.f90    # ✅ COMPLETO - 13 testes
│   ├── test_density_calculator.f90    # ✅ COMPLETO - 6 testes
│   └── test_ks_cycle.f90              # 🔜 TODO
│
├── examples/                   # 🔜 TODO
│   ├── harmonic_trap.f90
│   ├── double_barrier.f90
│   └── half_filling.f90
│
└── data/                       # 🔜 TODO
    ├── potential_params/       # Parâmetros de potenciais
    ├── reference_results/      # Resultados de referência (validação)
    └── tables/                 # Diretório de cache
        └── lsda_hub_u4.00      # Tabelas
```

---

## 🧮 Física do Problema

### Modelo de Hubbard 1D

Hamiltoniano:
$$
H = -t \sum_{i,j,\sigma} (c_{i\sigma}^\dagger c_{j\sigma} + \text{h.c.}) + U \sum_i n_{i\uparrow} n_{i\downarrow} + \sum_{i,\sigma} V_i^{\text{ext}} n_{i\sigma}
$$

- **t = 1**: Hopping (unidade de energia)
- **U**: Interação on-site (Hubbard U)
- **$V_i^{\text{ext}}$**: Potencial externo (armadilha, impurezas, etc)

### Teoria do Funcional da Densidade (DFT)

Kohn-Sham equations:
$$
\left[-\nabla^2 + V_{\text{ext}}(r) + V_H(r) + V_{\text{xc}}(r)\right] \psi_i(r) = \varepsilon_i \psi_i(r)
$$

Para o Hubbard 1D:
$$
H_{\text{KS}} = H_0 + V_{\text{ext}} + U \cdot n_{-\sigma} + V_{\text{xc}}
$$

### Bethe Ansatz de Lieb-Wu

Solução exata para o estado fundamental via Ansatz de Bethe.

**Equações para rapidities de carga** ($j = 1, 2, \ldots, N_{\uparrow}$):
$$
F_j^k = k_j - \frac{2\pi}{L} I_j - \frac{1}{L} \sum_{\alpha=1}^{M} \theta(k_j - \Lambda_\alpha, U) = 0
$$

**Equações para rapidities de spin** ($\alpha = 1, 2, \ldots, M$ onde $M = N_{\downarrow}$):
$$
F_\alpha^\Lambda = \frac{2\pi}{L} J_\alpha - \sum_{j=1}^{N_{\uparrow}} \theta(\Lambda_\alpha - k_j, U) + \sum_{\substack{\beta=1 \\ \beta \neq \alpha}}^{M} \Theta(\Lambda_\alpha - \Lambda_\beta, U) = 0
$$

Onde:
- $\theta(x, U) = 2\arctan(2x/U)$ (espalhamento carga-spin)
- $\Theta(x, U) = 2\arctan(x/U)$ (espalhamento spin-spin)
- $\{I_j\}$: números quânticos de carga (inteiros/semi-inteiros)
- $\{J_\alpha\}$: números quânticos de spin (inteiros/semi-inteiros)

**Energia:**
$$
E = -2 \sum_j \cos(k_j)
$$

**Funcional XC:**
$$
\begin{align}
E_{\text{xc}} &= E_{\text{BA}}[n_\uparrow, n_\downarrow] - E_0[n_\uparrow, n_\downarrow] \\
V_{\text{xc}}^\uparrow &= \frac{\partial E_{\text{xc}}}{\partial n_\uparrow} \\
V_{\text{xc}}^\downarrow &= \frac{\partial E_{\text{xc}}}{\partial n_\downarrow}
\end{align}
$$

---

## 🔧 Decisões Técnicas

### 1. Build System: Fortran Package Manager (fpm)

**Por quê?**
- Gerenciamento automático de dependências
- Build simplificado: `fpm build`, `fpm test`, `fpm run`
- Integração com Fortuno (testes)
- Padrão moderno da comunidade Fortran

**Configuração (`fpm.toml`):**
```toml
name = "lsda-hubbard"
version = "0.1.0"
license = "MIT"
author = "Guilherme Canella"
maintainer = "guycanella@gmail.com"

[build]
auto-executables = true
auto-tests = true
auto-examples = true

[dependencies]

[dev-dependencies]
fortuno-serial = { git = "https://github.com/fortuno-repos/fortuno.git" }

[library]
source-dir = "src"

[[executable]]
name = "lsda"
source-dir = "app"
main = "main.f90"
```

### 2. Diagonalização: LAPACK

**Rotinas utilizadas:**
- `DSYEV`: Autovalores + autovetores (matriz simétrica, método QR)
- `DSYEVD`: Versão divide-and-conquer (mais rápida para N > 100)

**Por quê não Givens?**
- LAPACK é ~10-100x mais rápido
- Implementação otimizada (BLAS nível 3)
- Bem testada e mantida
- Padrão industrial

**Interface:**
```fortran
subroutine diagonalize(H, eigenvalues, eigenvectors, n)
    real(real64), intent(in) :: H(:,:)
    real(real64), intent(out) :: eigenvalues(:)
    real(real64), intent(out) :: eigenvectors(:,:)
    integer, intent(in) :: n
    
    ! Wrapper para DSYEVD
    call dsyevd('V', 'U', n, H_copy, n, eigenvalues, work, lwork, &
                iwork, liwork, info)
end subroutine
```

### 3. Splines 2D: Implementação Própria

**Método:** Bicúbica separável (spline 1D em cada direção)

**Estratégia:**
1. Para $(n_{\uparrow}, n_{\downarrow})$ → calcular $(n=n_\uparrow+n_\downarrow, m=n_\uparrow-n_\downarrow)$
2. Para cada $n_i$ fixo: interpolação spline 1D em $m$
3. Com valores interpolados: segunda interpolação spline 1D em $n$

**Vantagens:**
- Código limpo e moderno (sem dependências pesadas)
- Facilmente testável
- Totalmente controlado
- $C^2$ contínuo

**Estrutura:**
```fortran
type :: spline2d_t
    real(real64), allocatable :: n_grid(:)      ! Grid de densidades
    real(real64), allocatable :: m_grid(:,:)    ! Grid de magnetizações
    real(real64), allocatable :: f(:,:)         ! Valores da função
    real(real64), allocatable :: d2f_dm2(:,:)   ! Derivadas segundas em m
    real(real64), allocatable :: d2f_dn2(:)     ! Derivadas segundas em n
end type
```

### 4. Testes: Fortuno

**Por quê?**
- Framework moderno para Fortran
- Sintaxe clara e expressiva
- Integração nativa com fpm
- Suporte a fixtures e parametrização

**Exemplo:**
```fortran
module test_bethe_ansatz
    use fortuno_serial, only: test_item, check => serial_check
    use bethe_ansatz
    implicit none
    
contains
    
    function test_free_fermions() result(tests)
        type(test_item) :: tests
        
        real(real64) :: E_calculated, E_exact
        
        ! U=0 → Fermi gas livre
        call solve_lieb_wu(N=10, U=0.0_real64, E=E_calculated)
        E_exact = -4.0_real64 / pi * integrate_fermi_sea(N, L)
        
        call check(abs(E_calculated - E_exact) < 1e-10, &
                   "Free fermions energy")
    end function
    
end module
```

### 5. Convenções de Código

**Indexação:**
- Arrays 1-indexed para sítios físicos (convenção natural do problema)
- 0-indexed quando apropriado (ex: índices auxiliares)

**Precisão:**
```fortran
use, intrinsic :: iso_fortran_env, only: real64
integer, parameter :: dp = real64
real(dp) :: variable  ! Float64 (double precision)
```

**Nomes:**
- Módulos: `snake_case` (ex: `bethe_ansatz`)
- Tipos: `snake_case_t` (ex: `system_params_t`)
- Funções/subrotinas: `snake_case` (ex: `solve_newton`)
- Constantes: `UPPER_SNAKE_CASE` (ex: `MAX_ITER`)

**Documentação (FORD-compliant):**
```fortran
!> Resolve as equações de Lieb-Wu para o estado fundamental
!! do modelo de Hubbard 1D usando o método de Newton-Raphson.
!!
!! @param[in]  n_up    Número de elétrons spin-up
!! @param[in]  n_dn    Número de elétrons spin-down
!! @param[in]  L       Número de sítios
!! @param[in]  U       Interação de Hubbard
!! @param[out] k       Rapidities de carga
!! @param[out] Lambda  Rapidities de spin
!! @param[out] energy  Energia do estado fundamental
subroutine solve_lieb_wu(n_up, n_dn, L, U, k, Lambda, energy)
```

---

## 🧩 Solução das Equações de Lieb-Wu

### Estratégia Híbrida Newton-Broyden

#### 1. Formulação do Problema

**Variáveis:** $\mathbf{x} = [k_1, k_2, \ldots, k_{N_\uparrow}, \Lambda_1, \Lambda_2, \ldots, \Lambda_M]$ onde $M = N_\downarrow$

**Sistema não-linear:** $\mathbf{F}(\mathbf{x}) = 0$

$$
F_j^k = k_j - \frac{2\pi}{L} I_j - \frac{1}{L} \sum_{\alpha=1}^{M} \theta(k_j - \Lambda_\alpha, U)
$$

$$
F_\alpha^\Lambda = \frac{2\pi}{L} J_\alpha - \sum_{j=1}^{N_\uparrow} \theta(\Lambda_\alpha - k_j, U) + \sum_{\substack{\beta=1 \\ \beta \neq \alpha}}^{M} \Theta(\Lambda_\alpha - \Lambda_\beta, U)
$$

**Jacobiano analítico:**
$$
\mathbf{J} = \begin{bmatrix}
\frac{\partial F^k}{\partial k} & \frac{\partial F^k}{\partial \Lambda} \\[1em]
\frac{\partial F^\Lambda}{\partial k} & \frac{\partial F^\Lambda}{\partial \Lambda}
\end{bmatrix}
$$

Com derivadas:
$$
\frac{d\theta}{dx} = \frac{4U}{U^2 + 4x^2}, \quad \frac{d\Theta}{dx} = \frac{2U}{U^2 + x^2}
$$

#### 2. Escolha do Método

**Heurística:**
```fortran
if (N_up + N_dn < 100) then
    ! Sistema pequeno → Newton com Jacobiano analítico
    call solve_newton(x, F, J)
    
else if (N_up + N_dn < 500) then
    ! Sistema médio → Híbrido
    call solve_broyden(x, F, n_iter=5)  ! Warm-up
    call solve_newton(x, F, J)          ! Finaliza
    
else
    ! Sistema grande → Broyden puro (economiza memória)
    call solve_broyden(x, F)
end if
```

#### 3. Newton com Line Search

**Algoritmo:**
```
1. Resolver J·Δx = -F  (via DGESV/DGETRF)
2. Line search: encontrar α ∈ (0,1] tal que
   ||F(x + α·Δx)|| < (1-c·α)·||F(x)||  (Armijo condition)
3. Atualizar: x ← x + α·Δx
4. Repetir até convergência
```

**Implementação:**
```fortran
subroutine newton_with_linesearch(x, tol, max_iter)
    do iter = 1, max_iter
        call compute_residual(x, F)
        call compute_jacobian(x, J)
        
        call solve_linear_system(J, -F, dx)  ! LAPACK
        
        ! Line search (backtracking)
        alpha = 1.0_real64
        do ls_iter = 1, 20
            x_trial = x + alpha * dx
            call compute_residual(x_trial, F_trial)
            
            if (norm2(F_trial) < 0.9_real64 * norm2(F)) exit
            alpha = 0.5_real64 * alpha
        end do
        
        x = x_trial
        
        if (norm2(F) < tol) exit
    end do
end subroutine
```

#### 4. Broyden (Quasi-Newton)

**Ideia:** Aproximar $\mathbf{J}^{-1}$ iterativamente sem recalcular Jacobiano

**Atualização de Broyden:**
$$
\mathbf{B}_{k+1} = \mathbf{B}_k + \frac{(\Delta\mathbf{x} - \mathbf{B}_k \cdot \Delta\mathbf{F}) \otimes \Delta\mathbf{F}^T}{\Delta\mathbf{F}^T \cdot \Delta\mathbf{F}}
$$

Onde $\mathbf{B} \approx \mathbf{J}^{-1}$ (inversa aproximada).

**Vantagens:**
- Não precisa calcular $\mathbf{J}$ a cada iteração
- Custo $O(n^2)$ por iteração (vs $O(n^3)$ do Newton)

**Desvantagens:**
- Convergência superlinear (vs quadrática do Newton)
- Precisa de bom chute inicial

#### 5. Continuação em U (Preditor-Corretor)

**Objetivo:** Resolver para vários valores de U reutilizando soluções

**Algoritmo:**
```
1. Resolver para U₀ "fácil" (ex: U/t = 6 ou U → ∞)

2. Para cada passo i:
   a) Preditor: x_guess(Uᵢ₊₁) = x(Uᵢ) + ΔU · (dx/dU)
      - Estimar dx/dU via diferença finita ou implicitamente
   
   b) Corretor: Resolver F(x, Uᵢ₊₁) = 0 com chute x_guess
      - Usar Newton ou Broyden
   
   c) Adaptar ΔU baseado em número de iterações:
      - Se convergiu rápido (< 4 iter): ΔU ← 1.2·ΔU
      - Se demorou (> 8 iter): ΔU ← 0.5·ΔU

3. Checkpoints: salvar x(U) a cada 10 pontos
```

**Sweep bidirecional:**
```
Forward:  U = 0 → 2 → 4 → 6 → ... → 10
Backward: U = 10 → 8 → 6 → ... → 0  (refinamento)
```

#### 6. Normalização e Escalonamento

**Problema:** $k \in [-\pi,\pi]$, $\Lambda \sim O(U)$ → mal-condicionado para U grande

**Solução:**
```fortran
! Escalonar variáveis
k_scaled = k / pi              ! k ∈ [-1, 1]
Lambda_scaled = Lambda / U     ! Λ ∈ O(1)

! Resolver sistema escalonado
call solve_newton(x_scaled, F_scaled, J_scaled)

! De-escalonar resultado
k = k_scaled * pi
Lambda = Lambda_scaled * U
```

#### 7. Critérios de Convergência

**Duplo critério:**
```fortran
converged = (norm2(F) < 1e-10_real64) .and. &
            (norm2(dx) / max(norm2(x), 1.0_real64) < 1e-12_real64)
```

**Fail-safes:**
```fortran
! Divergência
if (norm2(F) > 1e6 .or. any(ieee_is_nan(F))) then
    error stop "Solver divergiu!"
end if

! Estagnação
if (iter > max_iter) then
    print *, "Warning: máximo de iterações atingido"
    exit
end if
```

#### 8. Números Quânticos (Estado Fundamental)

**Distribuição de Fermi:**
```fortran
! Carga (spin-up)
do j = 1, N_up
    I(j) = j - (N_up + 1) / 2  ! Centrado, consecutivo
end do

! Spin (spin-down)
M = N_dn
do alpha = 1, M
    J(alpha) = alpha - (M + 1) / 2  ! Centrado, consecutivo
end do
```

Para $N_\uparrow=5$: $I = [-2, -1, 0, 1, 2]$  
Para $N_\downarrow=3$: $J = [-1, 0, 1]$

#### 9. Validação

**Casos limite (testes unitários):**

```fortran
! 1. U=0 (Fermi gas livre)
E_exact = -4/π ∫₀^(πn/2) cos(k) dk

! 2. U→∞ (forte acoplamento)
E_exact ~ -4·J₀(πn)  ! Função de Bessel

! 3. Half-filling (n=1, m=0)
E_exact = -4·J₀(U)

! 4. Polarizado (n=m)
E_exact = -2·cos(πn/2) - U·n²/4
```

#### 10. Paralelização

Grid $(n, m, U)$ é **embaraçosamente paralelo**:

```fortran
!$omp parallel do schedule(dynamic) private(solution)
do i_U = 1, n_U_points
    do i_n = 1, n_density_points
        do i_m = 1, n_mag_points(i_n)
            call solve_lieb_wu(n(i_n), m(i_m), U(i_U), solution)
            
            exc_table(i_n, i_m, i_U) = solution%exc
            Vxc_up_table(i_n, i_m, i_U) = solution%Vxc_up
            Vxc_dn_table(i_n, i_m, i_U) = solution%Vxc_dn
        end do
    end do
end do
!$omp end parallel do
```

---

## 📅 Roadmap de Desenvolvimento

### Fase 0: Infraestrutura ✅ 100% COMPLETO

- [x] Criar estrutura fpm
- [x] Módulo de tipos (`lsda_types.f90`)
- [x] Módulo de constantes (`lsda_constants.f90`)
- [x] Configurar Fortuno (dependência instalada)
- [x] Programa principal placeholder (`app/main.f90`)

**Status:** ✅ 

---

### Fase 1: Bethe Ansatz ✅ 100% COMPLETO

#### ✅ Completo:
- [x] **`bethe_equations.f90`** (487 linhas, 100% testado):
  - [x] Funções θ e Θ (espalhamento carga-spin e spin-spin)
  - [x] Derivadas dθ/dx, dΘ/dx (analíticas, validadas numericamente)
  - [x] Derivadas dθ/dU, dΘ/dU (para continuation method)
  - [x] `initialize_quantum_numbers()` - Estado fundamental (distribuição de Fermi)
  - [x] `compute_residual()` - Vetor F(x) das equações de Lieb-Wu
  - [x] `compute_jacobian()` - Matriz Jacobiana analítica (4 blocos)
  - [x] `compute_dFdU()` - Derivada do resíduo para preditor-corretor
  - [x] `compute_energy()` - Energia do estado fundamental E = -2·Σcos(k_j)
  - [x] Tratamento especial para U=0 (Fermi gas livre)

- [x] **`nonlinear_solvers.f90`** (303 linhas, 100% testado):
  - [x] `solve_linear_system()` - Wrapper LAPACK DGESV com LU decomposition
  - [x] `line_search()` - Backtracking com condição de Armijo
  - [x] `solve_newton()` - Newton-Raphson com line search adaptativo
  - [x] Tratamento especial para U=0 (solução analítica do Fermi gas)
  - [x] Detecção de estagnação, divergência e convergência
  - [x] Robustez: NaN checking, singular matrix handling

- [x] **`continuation.f90`** (369 linhas, 100% testado):
  - [x] `estimate_dxdU()` - Estimativa de dx/dU via diferenças finitas
  - [x] `sweep_U_forward()` - Sweep forward (U_min → U_max) com preditor linear
  - [x] `sweep_U_backward()` - Sweep backward (U_max → U_min) para refinamento
  - [x] `sweep_U_bidirectional()` - Média de forward + backward (maior precisão)
  - [x] Predictor-corrector: típico speedup de 5-10x vs soluções independentes

- [x] **`test/test_bethe_equations.f90`** (446 linhas, 17 testes ✅):
  - [x] Funções θ e Θ: zeros, antissimetria
  - [x] Derivadas analíticas vs numéricas: dθ/dx, dΘ/dx, dθ/dU, dΘ/dU
  - [x] Números quânticos: pares e ímpares
  - [x] Residual: dimensões, valores
  - [x] Jacobiano: dimensões, diagonal, validação numérica (< 1e-10)
  - [x] dF/dU: validação numérica
  - [x] Energia: U=0, dimensões

- [x] **`test/test_nonlinear_solvers.f90`** (302 linhas, 9 testes ✅):
  - [x] Sistema linear: 2×2, identidade, preservação de inputs
  - [x] Jacobiano: validação numérica
  - [x] Newton: Fermi gas (U=0), sistema pequeno (U=4)
  - [x] Convergência: flags, redução de resíduo
  - [x] Line search: eficácia

- [x] **`test/test_continuation.f90`** (198 linhas, 5 testes ✅):
  - [x] `estimate_dxdU`: diferenças finitas simples
  - [x] `sweep_forward`: 3 pontos, convergência total
  - [x] `sweep_backward`: 3 pontos
  - [x] `sweep_bidirectional`: consistência entre métodos

#### 🏆 Conquistas da Fase 1:
- ✅ **31 testes unitários** passando (100% de sucesso)
- ✅ **Jacobiano validado numericamente** (erro < 1e-10)
- ✅ **Continuation method implementado**: predictor-corrector com sweeps bidirecional
- ✅ **Casos especiais tratados**: U=0 (Fermi gas livre)
- ✅ **Newton robusto**: Line search + detecção de estagnação + NaN checking
- ✅ **Código documentado**: Comentários FORD-compliant em todos os módulos
- ✅ **Performance**: Continuation 5-10x mais rápido que soluções independentes

**Duração:** 4 dias
**Linhas de código:** ~1159 (produção) + ~946 (testes)
**Status:** ✅ **FASE 1 COMPLETA!**

---

### Fase 2: Geração e I/O de Tabelas ✅ COMPLETA

#### ✅ Completo (100%):
- [x] **`table_io.f90`** (~400+ linhas, totalmente testado):
  - [x] Tipo `xc_table_t` para armazenar tabelas XC
  - [x] `read_cpp_table()` - Leitura de tabelas ASCII legadas (formato C++)
  - [x] `write_fortran_table()` - Escrita em formato binário nativo Fortran
  - [x] `read_fortran_table()` - Leitura de formato binário (~10x mais rápido que ASCII)
  - [x] `extract_U_from_filename()` - Parser de nome de arquivo `lsda_hub_uX.XX`
  - [x] `deallocate_table()` - Gerenciamento de memória
  - [x] `print_table_info()` - Diagnóstico e debug

- [x] **`convert_tables.f90`** (executável utilitário):
  - [x] Conversão em batch de 25 tabelas C++ → Fortran binário
  - [x] Valores de U: 1.00, 1.10, 2.00, 3.00, 4.00, 4.10, 5.00, 5.90, 6.00, 6.10, 6.90, 7.00, 7.10, 7.90, 8.00, 8.10, 8.90, 9.00, 9.10, 10.00, 12.00, 14.00, 16.00, 18.00, 20.00
  - [x] Argumentos de linha de comando: `fpm run convert_tables -- <input_dir> <output_dir>`
  - [x] Relatório de progresso e estatísticas de conversão

- [x] **`test_table_io.f90`** (274 linhas - 10 testes unitários):
  - [x] Leitura de tabelas C++ ASCII
  - [x] Escrita/leitura de formato binário Fortran
  - [x] Validação de roundtrip (ASCII → binário → memória)
  - [x] Parsing de U a partir do nome do arquivo

- [x] **`bethe_tables.f90`** (325 linhas, 6 testes - totalmente implementado):
  - [x] Tipo `grid_params_t` para configurar grid de densidades
  - [x] `compute_E0()` - Energia não-interagente (Fermi gas livre)
  - [x] `compute_E_xc()` - Energia XC: E_xc = E_BA - E_0
  - [x] `compute_V_xc_numerical()` - Potenciais XC via derivadas de 5 pontos
  - [x] `generate_xc_table()` - Geração completa de tabela para dado U
  - [x] `generate_table_grid()` - Geração flexível de grid com parâmetros customizados
  - [x] Tratamento especial para casos limite:
    - [x] U=0 (Fermi gas livre - retorna E_xc=0)
    - [x] Half-filling (n=1, m=0)
    - [x] Polarizado (m=n)
  - [x] Integração total com módulos `bethe_equations`, `nonlinear_solvers`, `table_io`

- [x] **`test_bethe_tables.f90`** (170 linhas - 6 testes):
  - [x] Teste E0 para half-filling
  - [x] Teste E0 para sistema polarizado
  - [x] Teste E_xc = 0 para U=0
  - [x] Teste simetria V_xc (V_xc_up = V_xc_dn quando n_up = n_dn)
  - [x] Teste parâmetros padrão do grid
  - [x] Teste geração de tabela pequena

#### 🎉 Conquistas da Fase 2:
- ✅ **16 testes unitários** (10 I/O + 6 geração) passando (100%)
- ✅ **Pipeline completo**: Bethe Ansatz → E_xc → V_xc → Tabela → I/O
- ✅ **Derivadas numéricas** de 5 pontos para alta precisão
- ✅ **Grid flexível** com parâmetros configuráveis
- ✅ **Casos especiais** corretamente tratados
- ✅ **Total Fase 2:** 888 linhas produção + 444 linhas testes

**Duração:** ~3 dias
**Status:** ✅ **FASE 2 COMPLETA!**

#### 🔜 Melhorias Futuras (Opcionais):
- [ ] Paralelização OpenMP do grid (n, m) - embaraçosamente paralelo
- [ ] Validação física: comparação quantitativa com tabelas C++ legadas
- [ ] Otimização de performance (profiling, vetorização)
- [ ] `table_manager.f90`: Cache inteligente para múltiplos U

---

### Fase 3: Interpolação de Splines 2D ✅ COMPLETA

**Objetivo:** Implementar interpolação bicúbica 2D para avaliar funcionais XC em pontos arbitrários (n, m).

#### ✅ Completo (100%):
- [x] **`spline2d.f90`** (351 linhas, 5 testes, 100% testado):
  - [x] Tipo `spline2d_t` para armazenar coeficientes da spline em grids irregulares
  - [x] `spline1d_coeff()` - Algoritmo de Thomas para spline 1D (natural/clamped BC)
  - [x] `spline2d_init()` - Construir splines separáveis a partir de tabela 2D
  - [x] `spline2d_eval()` - Avaliar spline em ponto (n, m) via interpolação separável
  - [x] `find_interval()` - Busca binária para localizar intervalo do grid
  - [x] Tratamento de grids irregulares (n_y varia com x)
  - [x] Arrays 0-indexed internos para compatibilidade com algoritmo clássico

- [x] **`xc_lsda.f90`** (335 linhas, 6 testes, 100% testado):
  - [x] Tipo `xc_lsda_t` contendo splines de exc, vxc_up, vxc_dw
  - [x] `xc_lsda_init()` - Carregar tabela e construir 3 splines (exc, vxc_up, vxc_dw)
  - [x] `get_exc(n_up, n_dw)` - Energia XC por partícula via interpolação
  - [x] `get_vxc(n_up, n_dw, v_xc_up, v_xc_dw)` - Potenciais XC para ambos os spins
  - [x] Conversão (n_up, n_dw) ↔ (n, m) via `convert_to_nm()`
  - [x] **4 regiões de simetria física** mapeadas para Region I:
    - [x] Region I (m≥0, n≤1): Identidade
    - [x] Region II (m<0, n≤1): Spin exchange
    - [x] Region III (m<0, n>1): Particle-hole
    - [x] Region IV (m≥0, n>1): Combinada
  - [x] Tratamento especial para U=0, n=0, densidades fora da faixa física

- [x] **`test_spline2d.f90`** (196 linhas, 5 testes ✅):
  - [x] Init/destroy: alocação, inicialização, cleanup
  - [x] Interpolação exata em pontos do grid (erro < 1e-9)
  - [x] Funções lineares: spline exata para f(x,y) = ax + by + c
  - [x] Funções separáveis: f(x,y) = g(x)·h(y) com alta precisão
  - [x] Casos limite: single x point, bounds checking

- [x] **`test_xc_lsda.f90`** (200 linhas, 6 testes ✅):
  - [x] Init/destroy com tabelas reais (U=4.00, U=2.00)
  - [x] Avaliação de exc retorna valores válidos e não-zero para U>0
  - [x] **Simetria de spin**: exc(n_up, n_dw) = exc(n_dw, n_up)
  - [x] **Simetria de potenciais**: V_up(n_up, n_dw) = V_dw(n_dw, n_up)
  - [x] Determinação de regiões (I, II, III, IV)
  - [x] Transformações de simetria corretas

#### 🏆 Conquistas da Fase 3:
- ✅ **11 testes unitários** passando (100% de sucesso)
- ✅ **Spline 2D separável** implementada com grid irregular
- ✅ **Simetrias físicas** mapeando todo o domínio (n, m) para tabela compacta
- ✅ **Integração total** com pipeline Bethe → Tabelas → Splines
- ✅ **Código robusto**: tratamento de U=0, densidades zero, bounds checking
- ✅ **Convenção padronizada**: n_dw (não n_dn) para spin-down
- ✅ **Total Fase 3:** 686 linhas produção + 396 linhas testes

**Estratégia Implementada:**
1. ✅ Spline 1D em cada direção (separável)
2. ✅ Para cada n_i fixo: spline cúbica natural em m
3. ✅ Com valores interpolados: interpolação linear em n
4. ✅ Simetrias físicas reduzem domínio de [0,1]×[-1,1] → [0,1]×[0,n]

**Duração:** ~2 dias
**Status:** ✅ **FASE 3 COMPLETA!**

---

### Fase 4: Potenciais & Sistema de Erros ✅ COMPLETA

**Objetivo:** Implementar sistema de potenciais externos e tratamento centralizado de erros.

#### ✅ Completo (100%):
- [x] **`lsda_errors.f90`** (224 linhas, 13 testes):
  - [x] Códigos de erro organizados por categoria (input 1-99, numerical 100-199, I/O 200-299, memory 300-399)
  - [x] `get_error_message()` - Mensagens legíveis para cada código de erro
  - [x] `error_handler()` - Handler centralizado com opção fatal
  - [x] `check_bounds()`, `check_positive()`, `check_range()` - Utilitários de validação

- [x] **`potential_uniform.f90`** (34 linhas): V(i) = V₀
  - [x] Potencial constante (shift global de energia)

- [x] **`potential_harmonic.f90`** (46 linhas): V(i) = 0.5·k·(i-center)²
  - [x] Armadilha harmônica parabólica
  - [x] Simetria de paridade V(i) = V(L+1-i)
  - [x] Modela optical traps, cria shell structure

- [x] **`potential_impurity.f90`** (191 linhas):
  - [x] `potential_impurity_single()` - Impureza pontual única
  - [x] `potential_impurity_multiple()` - Múltiplas impurezas (com soma se sobrepõem)
  - [x] `potential_impurity_random()` - Impurezas aleatórias com concentração fixa

- [x] **`potential_random.f90`** (152 linhas):
  - [x] `potential_random_uniform()` - Desordem uniforme V(i) ~ U[-W/2, W/2]
  - [x] `potential_random_gaussian()` - Desordem gaussiana V(i) ~ N(0, σ²)
  - [x] Box-Muller transform para geração de normais
  - [x] Modela localização de Anderson

- [x] **`potential_barrier.f90`** (157 linhas):
  - [x] `potential_barrier_single()` - Barreira retangular única
  - [x] `potential_barrier_double()` - Dupla barreira (poço quântico)
  - [x] Tunelamento quântico, ressonâncias Fabry-Pérot

- [x] **`potential_quasiperiodic.f90`** (100 linhas): V(i) = λ·cos(2πβi + φ)
  - [x] Aubry-André-Harper (AAH) model
  - [x] Extended phase (λ < 2): estados deslocalizados
  - [x] Critical phase (λ = 2): funções de onda multifractais
  - [x] Localized phase (λ > 2): estados exponencialmente localizados
  - [x] Modela localização de Anderson sem desordem

- [x] **`potential_factory.f90`** (186 linhas):
  - [x] `create_potential()` - Factory para criar potenciais via string
  - [x] `get_potential_info()` - Informações sobre cada tipo
  - [x] Suporte para 7 tipos: uniform, harmonic, impurity_single, random_uniform, random_gaussian, barrier_single, barrier_double, quasiperiodic

- [x] **`test_potentials.f90`** (585 linhas, 21 testes):
  - [x] Testes com explicações físicas detalhadas nos comentários
  - [x] Uniform: constância, Harmonic: simetria/mínimo central
  - [x] Impurity: posição/bounds/overlap/concentração
  - [x] Random: média zero, distribuições corretas
  - [x] Barrier: largura/bounds/poço quântico/não-sobreposição
  - [x] Quasiperiodic: golden ratio, phase shift, critical point, localization
  - [x] Factory: criação/comparação/tipo inválido

- [x] **`test_lsda_errors.f90`** (284 linhas, 13 testes):
  - [x] Verificação de códigos em intervalos corretos
  - [x] Mensagens para todos os tipos de erro
  - [x] Utilitários de validação (bounds, positive, range)

#### 🏆 Conquistas da Fase 4:
- ✅ **34 testes unitários** passando (100% de sucesso)
- ✅ **7 tipos de potenciais** implementados com física completa
- ✅ **Sistema de erros robusto** para todo o projeto
- ✅ **Factory pattern** para criação dinâmica de potenciais
- ✅ **Documentação física detalhada** em todos os testes
- ✅ **Total Fase 4:** 1090 linhas produção + 869 linhas testes

**Física Implementada:**
- ✅ Armadilha harmônica (optical traps, cold atoms)
- ✅ Localização de Anderson (random disorder & quasiperiodic AAH)
- ✅ Transição metal-isolante (AAH model, λ = 2 critical point)
- ✅ Tunelamento quântico (barriers)
- ✅ Ressonâncias Fabry-Pérot (double barriers)
- ✅ Impurezas magnéticas (random impurities)

**Duração:** ~1-2 dias
**Status:** ✅ **FASE 4 COMPLETA!**

---

### Fase 5: Ciclo Auto-Consistente (3-4 dias) 🔜 TODO

- [ ] `density_calculator.f90`: Ocupação de níveis
- [ ] `convergence_monitor.f90`: Critérios de parada
- [ ] `mixing_schemes.f90`: Linear mixing
- [ ] `ks_cycle.f90`: Loop SCF completo
- [ ] Testes:
- [ ] U=0, BC periódica → Fermi gas
- [ ] Half-filling, U>0 → comparar literatura

**🎉 MILESTONE:** Código funcional end-to-end!

---

### Fase 6: Features Avançadas (1 semana) 🔜 TODO

- [ ] `degeneracy_handler.f90`: Tratamento de níveis degenerados
- [ ] `symmetry.f90`: Exploração de paridade
- [ ] `twisted_bc.f90`: Boundary conditions torcidas
- [ ] Potenciais avançados (impurity, barrier, random, etc)
- [ ] Testes para cada feature

---

### Fase 7: Otimização (ongoing) 🔜 TODO

- [ ] Paralelização OpenMP (Bethe Ansatz + KS loop)
- [ ] Profiling e otimização de hotspots
- [ ] I/O melhorado (HDF5?)
- [ ] Documentação completa (FORD)
- [ ] Benchmarks vs código C++ original

---

## 🧪 Estratégia de Testes

### Pirâmide de Testes

```
         /\
        /  \
       /E2E \       (End-to-End: ciclo completo, casos físicos)
      /------\
     /Integr. \     (Integração: módulos combinados)
    /----------\
   /Unit Tests  \   (Unitários: funções individuais)
  /--------------\
```

### Tipos de Testes

#### 1. Testes Unitários (Fortuno)
```fortran
! test/test_bethe_equations.f90
- test_theta_function()
- test_residual_computation()
- test_jacobian_analytical_vs_numerical()
- test_quantum_numbers_initialization()

! test/test_splines.f90
- test_spline1d_linear_function()
- test_spline1d_quadratic_function()
- test_spline2d_separable_function()
- test_interpolation_accuracy()

! test/test_hamiltonian.f90
- test_tight_binding_matrix_open_bc()
- test_tight_binding_matrix_periodic_bc()
- test_eigenvalues_free_fermions()
```

#### 2. Testes de Integração
```fortran
! test/test_integration.f90
- test_bethe_to_splines_pipeline()
- test_hamiltonian_diagonalization_with_xc()
- test_ks_single_iteration()
```

#### 3. Testes End-to-End
```fortran
! test/test_e2e.f90
- test_fermi_gas_u0()
- test_half_filling_u4()
- test_harmonic_trap_convergence()
```

#### 4. Testes de Regressão
- Comparar saída com resultados do código C++ original
- Armazenar resultados de referência em `data/reference_results/`

### Cobertura de Código

**Meta:** > 80% de cobertura

```bash
# Usando gfortran + gcov
fpm test --flag "-fprofile-arcs -ftest-coverage"
gcov src/**/*.f90
lcov --capture --directory . --output-file coverage.info
genhtml coverage.info --output-directory coverage_html
```

---

## 🔬 Física e Validação

### Casos de Validação Obrigatórios

#### 1. Fermi Gas (U=0)
```
Input:  N=20, L=20, U=0, BC=periodic
Output: E/L = -2·sin(π/2) / π ≈ -0.6366
        n(i) = 1.0 (uniforme)
```

#### 2. Half-Filling (n=1)
```
Input:  N=L, U=4.0, BC=periodic
Output: E/L = função de U (comparar com Essler et al.)
        n(i) = 1.0 (uniforme)
        m = 0 (não-magnético)
```

#### 3. Armadilha Harmônica
```
Input:  N=20, L=40, U=4.0, V(i) = k*(i-L/2)^2
Output: n(i) = perfil gaussiano (shell structure)
        Comparar com Thomas-Fermi para k→0
```

#### 4. Dupla Barreira (Quantum Well)
```
Input:  Dupla barreira, U=4.0
Output: Estados localizados no poço
        Tunelamento ressonante
```

### Benchmarks de Performance

| Caso                  | C++ (original) | Fortran (meta) | Status |
|-----------------------|----------------|----------------|--------|
| Bethe (N=100)         | N/A (tabelas)  | < 1s           | 🔜     |
| Spline interpolation  | ~10μs          | < 5μs          | 🔜     |
| Diagonalização (N=100)| ~50ms (Givens) | < 5ms (LAPACK) | 🔜     |
| Ciclo KS (10 iter)    | ~5s            | < 2s           | 🔜     |

---

## 📋 Parâmetros Técnicos Definidos

### Intervalo de U
- **Range:** 0 ≤ U/t ≤ 100
- **Unidades:** U em unidades de hopping t
- **Casos especiais:**
  - U = 0: Fermi gas livre
  - U → ∞: Limite de forte acoplamento

### Tamanhos de Sistema
- **Sítios (L):** Até 200 sites
- **Elétrons (N):** 0 ≤ N ≤ 2L
- **Densidades:** 0 ≤ n = N/L ≤ 2
  - n = 1: Half-filling (caso especial)
  - n < 1: Dopagem tipo-n
  - n > 1: Dopagem tipo-p

### Precisão Numérica
- **Padrão:** `real64` (double precision, ~16 dígitos)
- **Futuro:** Possível upgrade para `real128` se necessário
- **Tolerância:** TOL = 1.0e-16 (padrão do código original)
- **Threshold para U=0:** SMALL = 1.0e-9

### Build e Compilação
```bash
# Build padrão
fpm build

# Build otimizado
fpm build --profile release --flag "-O3 -march=native"

# Com OpenMP (futuro)
fpm build --flag "-fopenmp"

# Rodar
fpm run

# Testes
fpm test
```

---

## 📚 Referências

### Papers Fundamentais

1. **Lieb & Wu (1968)**  
   "Absence of Mott transition in an exact solution of the short-range, one-band model in one dimension"  
   *Physical Review Letters*, 20(25), 1445.

2. **Essler et al. (2005)**  
   *The One-Dimensional Hubbard Model*  
   Cambridge University Press. (Livro completo sobre Bethe Ansatz)

3. **Capelle & Campo (2013)**  
   "Density functionals and model Hamiltonians: Pillars of many-particle physics"  
   *Physics Reports*, 528(3), 91-159.

4. **Xianlong et al. (2006)**  
   "Lattice density functional theory at finite temperature with strongly density-dependent exchange-correlation potentials"  
   *Physical Review B*, 73(16), 165120.

### Documentação Técnica

- **LAPACK Users' Guide**: https://netlib.org/lapack/lug/
- **Fortran 2018 Standard**: https://j3-fortran.org/
- **fpm Documentation**: https://fpm.fortran-lang.org/
- **Fortuno**: https://github.com/fortuno-repos/fortuno
- **FORD (Fortran Documenter)**: https://github.com/Fortran-FOSS-Programmers/ford

### Código de Referência

- Código C++ original (Vivaldo Campo Jr)
- DMFT solvers (TRIQS, w2dynamics)
- Exact diagonalization codes (ALPS, ITensor)

---

## 📝 Notas de Implementação

### Decisões Tomadas

1. **Números Quânticos**: São `real(dp)` (não `integer`) porque podem ser semi-inteiros quando N é par
2. **Índices**: Arrays 1-indexed para sítios físicos (convenção do problema)
3. **Precisão**: `real64` (double) é suficiente para U ∈ [0, 100]
4. **U → 0**: Tratamento especial em θ e Θ para evitar divisão por zero
5. **Jacobiano**: Implementação analítica (não diferenças finitas) para máxima precisão
6. **Documentação**: Padrão FORD para geração automática de docs

### TODOs e Decisões Pendentes

- [ ] **Grid de tabelas:** Quantos pontos (n,m,U)? Espaçamento uniforme ou adaptativo?
- [x] **Formato de output:** ✅ Binário Fortran nativo (implementado em `table_io.f90`)
- [x] **Formato de input:** ✅ ASCII C++ legado + binário Fortran
- [ ] **Paralelização:** OpenMP apenas ou também MPI para grids grandes?
- [ ] **Precisão:** Float64 suficiente ou Float128 em alguns casos?
- [ ] **Estados excitados:** Implementar? (mudando {I_j}, {J_α})
- [x] **Broyden:** ✅ Implementado apenas Newton (decisão: Newton suficiente para tabelas)
- [x] **Checkpointing:** ✅ Não necessário (continuation rápido o suficiente)

### Perguntas em Aberto

1. Como tratar U < 0 (interação atrativa)? Usar simetria ou resolver separadamente?
2. Implementar TBA (Thermodynamic Bethe Ansatz) para L → ∞?
3. Adicionar temperatura T > 0 (Yang-Yang)?
4. Implementar funcionais GGA além do LDA?

---

## 📊 Status do Projeto

**Versão:** 0.6.0-dev
**Status:** ✅ Fases 1-5 Completas → 🔄 Fase 6 em Progresso (Densidade & SCF - 20% completo)
**Última atualização:** 2025-01-16

### Progresso Geral

```
[████████████████████████████████] 100% Fase 1: Bethe Ansatz Core (COMPLETO ✅)
[████████████████████████████████] 100% Fase 2: Geração de Tabelas XC (COMPLETO ✅)
[████████████████████████████████] 100% Fase 3: Splines 2D (COMPLETO ✅)
[████████████████████████████████] 100% Fase 4: Potenciais & Erros (COMPLETO ✅)
[████████████████████████████████] 100% Fase 5: Hamiltoniano & Diagonalização (COMPLETO ✅)
[██████░░░░░░░░░░░░░░░░░░░░░░░░░░]  20% Fase 6: Densidade & SCF Cycle (EM PROGRESSO 🔄)
[░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░]   0% Fase 7: Otimização
```

### Checklist de Progresso

#### Infraestrutura ✅
- [x] Estrutura fpm
- [x] Tipos básicos (`lsda_types.f90`)
- [x] Constantes (`lsda_constants.f90`)
- [x] Sistema de erros (`lsda_errors.f90`) ✅
- [x] Sistema de testes (Fortuno configurado)
- [ ] CI/CD

#### Core Physics 🔄
- [x] **Fase 1 - Bethe Ansatz Core** (100% ✅):
  - [x] Equações de Lieb-Wu (`bethe_equations.f90`) ✅
  - [x] Solvers Newton-Raphson (`nonlinear_solvers.f90`) ✅
  - [x] Continuation methods (`continuation.f90`) ✅
  - [x] Testes unitários (31 testes, 100% passando) ✅

- [x] **Fase 2 - Geração de Tabelas XC** (100% ✅):
  - [x] I/O de tabelas (`table_io.f90`) ✅
  - [x] Geração de tabelas (`bethe_tables.f90`) ✅
  - [x] Utilitário de conversão (`convert_tables.f90`) ✅
  - [x] Testes unitários (16 testes, 100% passando) ✅

- [x] **Fase 3 - Splines 2D** (100% ✅):
  - [x] Interpolação bicúbica (`spline2d.f90`) ✅
  - [x] Interface XC funcional (`xc_lsda.f90`) ✅
  - [x] Testes unitários (11 testes, 100% passando) ✅

- [x] **Fase 4 - Potenciais & Erros** (100% ✅):
  - [x] Sistema de erros centralizado (`lsda_errors.f90`) ✅
  - [x] 7 tipos de potenciais implementados (incl. quasiperiodic AAH) ✅
  - [x] Factory pattern para potenciais ✅
  - [x] Testes unitários (34 testes, 100% passando) ✅

- [x] **Fase 5 - Hamiltoniano & Diagonalização** (100% ✅):
  - [x] Boundary conditions (`boundary_conditions.f90`) ✅
    - [x] Implementação: BC_OPEN, BC_PERIODIC, BC_TWISTED ✅
    - [x] Validação de parâmetros ✅
    - [x] Eigenvalues analíticos para free particles ✅
    - [x] Testes unitários (17 testes, 100% passando) ✅
  - [x] Construção do Hamiltoniano (`hamiltonian_builder.f90`) ✅
    - [x] `validate_hamiltonian_inputs()`: validação com NaN/Inf ✅
    - [x] `build_hamiltonian()`: H real com BCs ✅
    - [x] `build_hamiltonian_complex()`: H complexo (twisted BC) ✅
    - [x] `build_hamiltonian_free()`: H livre (validação) ✅
    - [x] `compute_effective_potential()`: V_eff = V_ext + V_xc ✅
    - [x] Bug fix: loop de hopping corrigido ✅
    - [x] Testes unitários (18 testes, 100% passando) ✅
  - [x] Wrapper LAPACK para diagonalização (`lapack_wrapper.f90`) ✅
    - [x] `validate_diagonalization_inputs()`: validação de dimensões ✅
    - [x] `diagonalize_symmetric_real()`: DSYEVD para matrizes reais simétricas ✅
    - [x] `diagonalize_symmetric_real_values_only()`: eigenvalues only (mais rápido) ✅
    - [x] `diagonalize_hermitian_complex()`: ZHEEVD para matrizes complexas Hermitianas ✅
    - [x] `diagonalize_hermitian_complex_values_only()`: eigenvalues only ✅
    - [x] Workspace query em duas fases (lwork=-1 → allocate) ✅
    - [x] Interface LAPACK sem bind(C) (convenção Fortran nativa) ✅
    - [x] Testes unitários (18 testes, 100% passando) ✅
  - [x] Tratamento de degenerescências (`degeneracy_handler.f90`) ✅
    - [x] `find_degenerate_subspaces()`: detectar grupos onde |λᵢ - λⱼ| < tol ✅
    - [x] `orthonormalize_degenerate_subspace()`: QR (DGEQRF/DORGQR) para vetores reais ✅
    - [x] `orthonormalize_degenerate_subspace_complex()`: Gram-Schmidt modificado ✅
    - [x] `compute_degeneracy_count()`: contar degenerescências ✅
    - [x] `verify_orthonormality()`: verificar ||V^T V - I|| < tol ✅
    - [x] Bug fix: removido double conjugation (DOT_PRODUCT já conjuga) ✅
    - [x] Bug fix: workspace query separada para DORGQR ✅
    - [x] Testes unitários (13 testes, 100% passando) ✅

- [~] **Fase 6 - Densidade & SCF Cycle** (20% 🔄):
  - [x] Cálculo de densidade (`density_calculator.f90`) ✅
    - [x] `compute_density_spin()`: n_σ(i) = Σⱼ |ψⱼ(i)|² (real/complex overload) ✅
    - [x] `compute_total_density()`: n(i) = n↑(i) + n↓(i) ✅
    - [x] `verify_particle_number()`: Σn(i) = N ✅
    - [x] `check_density_bounds()`: 0 ≤ n_σ(i) ≤ 1, 0 ≤ n(i) ≤ 2 ✅
    - [x] Bug fix: variável i→j em loop de `check_density_bounds` ✅
    - [x] Testes unitários (6 testes, 100% passando) ✅
  - [ ] Mixing SCF (`scf_mixer.f90`) 🔜
    - [ ] Linear mixing: ρ_new = α·ρ_out + (1-α)·ρ_in
    - [ ] Broyden mixing: aceleração de convergência
    - [ ] Anderson mixing: alternativa
  - [ ] Ciclo Kohn-Sham (`ks_cycle.f90`) 🔜
    - [ ] Loop SCF completo: H → diag → density → V_xc → H'
    - [ ] Monitoramento de convergência
    - [ ] Mixing adaptativo

- [ ] **Fase 7**: Otimização (opcional)
  - [ ] Simetria de paridade (`symmetry.f90`)
    - [ ] `check_parity_symmetry()`: detectar V(i) = V(L+1-i)
    - [ ] `block_diagonalize_hamiltonian()`: split H → H_even, H_odd
    - [ ] `reconstruct_eigenstates()`: merge eigenvectors
    - [ ] Speedup 4x para potenciais simétricos
  - [ ] Paralelização OpenMP
  - [ ] Profiling e otimização

#### Features 🔄
- [x] Potenciais (7 tipos completos: uniform, harmonic, impurity, random, barrier, quasiperiodic) ✅
- [x] Boundary Conditions (Open, Periodic, Twisted) ✅
- [x] Diagonalização LAPACK (real simétrico & complexo Hermitiano) ✅
- [x] Tratamento de degenerescências (QR/Gram-Schmidt) ✅
- [ ] Simetria de paridade (próximo 🔜)

#### Qualidade ✅
- [x] Testes unitários Fase 1 (31 testes, 100% passando) ✅
- [x] Testes unitários Fase 2 (16 testes, 100% passando) ✅
- [x] Testes unitários Fase 3 (11 testes, 100% passando) ✅
- [x] Testes unitários Fase 4 (34 testes, 100% passando) ✅
- [x] Testes unitários Fase 5 (66 testes, 100% passando) ✅
- [x] Testes unitários Fase 6 (6 testes, 100% passando) ✅
- [x] **Total: 164 testes, 100% passando** ✅
- [x] Pipeline Bethe → Tabelas → Splines → Potenciais → Hamiltoniano → Diagonalização → Densidade validado ✅
- [ ] Testes E2E (ciclo KS completo)
- [ ] Documentação completa (FORD)
- [ ] Benchmarks de performance

---

## 🎓 Para Novos Desenvolvedores

### Quick Start

```bash
# Clonar repositório
git clone https://github.com/guycanella/lsda-hubbard-fortran.git
cd lsda-hubbard-fortran

# Build
fpm build

# Rodar testes (quando implementados)
fpm test

# Exemplo simples (quando implementado)
fpm run --example harmonic_trap
```

### Onde Começar?

1. **Leia:** `README.md` (uso básico) e este `PROJECT_CONTEXT.md` (contexto técnico)
2. **Entenda:** Leia `src/types/lsda_types.f90` e `lsda_constants.f90` para ver estruturas de dados
3. **Estude:** Revise módulos completos da Fase 1:
   - `src/bethe_ansatz/bethe_equations.f90` - Equações de Lieb-Wu
   - `src/bethe_ansatz/nonlinear_solvers.f90` - Newton-Raphson
   - `src/bethe_ansatz/continuation.f90` - Continuation methods
   - `src/bethe_ansatz/table_io.f90` - I/O de tabelas
4. **Contribua:** Próximo arquivo: `src/bethe_ansatz/bethe_tables.f90` (geração de tabelas XC)

### Recursos de Aprendizado

- **Fortran moderno:** https://fortran-lang.org/learn/
- **Bethe Ansatz:** Essler et al., "The One-Dimensional Hubbard Model"
- **DFT:** Capelle & Campo, "Density functionals and model Hamiltonians"
- **Newton-Raphson:** Numerical Recipes (Press et al.)

---

## 📄 Licença

Este projeto é licenciado sob a [MIT License](LICENSE).

---

## 👨‍💻 Informações de Desenvolvimento

**Mantido por:** Guilherme Canella
**Contato:** guycanella@gmail.com
**Repositório:** https://github.com/guycanella/lsdaks
**Última atualização:** 2025-01-14
**Status:** Fases 1-4 Completas (100%) → Fase 5 em Progresso (67% - Hamiltoniano & BCs completos, falta simetria)

---

## 📅 Histórico de Mudanças

### 2025-01-16 - Fase 6: Cálculo de Densidade Implementado! 🎉
- ✅ **MILESTONE:** Fase 5 completa (100%)! Fase 6 iniciada (densidade de autoestados KS).

  **`density_calculator.f90` implementado** (203 linhas, 6 testes):
  - ✅ `compute_density_spin()`: Cálculo de n_σ(i) = Σⱼ |ψⱼ(i)|² para estados ocupados
    - Interface genérica: overload para eigenvectors reais (OBC/PBC) e complexos (TBC)
    - Para T=0 (ground state): ocupar primeiros N níveis
  - ✅ `compute_total_density()`: n(i) = n↑(i) + n↓(i)
  - ✅ `verify_particle_number()`: Verifica Σn(i) = N dentro de TOL=1e-10
  - ✅ `check_density_bounds()`: Valida bounds físicos
    - 0 ≤ n_σ(i) ≤ 1 para cada spin (máximo 1 férmion por site/spin)
    - 0 ≤ n(i) ≤ 2 para densidade total (Pauli exclusion)
  - ✅ **Bug fix**: Variável `i→j` em loop de `check_density_bounds` (linha 195)

  **Testes implementados** (291 linhas, 6 testes):
  - ✅ `test_single_electron_density`: 1 elétron em caixa 1D (OBC)
    - Verifica n(i) = |ψ₁(i)|², densidade máxima no centro
  - ✅ `test_half_filling_unpolarized`: N=L, PBC, U=0
    - Densidade uniforme n(i)=1, simetria de spin n_up=n_dw=0.5
  - ✅ `test_particle_number_conservation`: Σn(i) = N (3 elétrons, L=8)
  - ✅ `test_density_positivity`: n(i) ≥ 0 sempre (física básica)
  - ✅ `test_physical_bounds`: Verifica 0 ≤ n_σ(i) ≤ 1, 0 ≤ n(i) ≤ 2
  - ✅ `test_density_from_harmonic_trap`: Armadilha harmônica
    - Shell structure: densidade maior no centro, decai nas bordas
    - Perfil qualitativo correto (gaussiano-like)

  **Correções durante implementação:**
  - ✅ `lsda_errors.f90`: Adicionado `ERROR_UNPHYSICAL_DENSITY` aos exports públicos
  - ✅ Interface `apply_potential_harmonic`: Corrigida chamada (center calculado automaticamente)
  - ✅ Interface `diagonalize_symmetric_real`: Ordem de parâmetros corrigida (H, L, ...) não (L, H, ...)

  **Estatísticas Fase 6 (parcial):**
  - ✅ Total: 203 linhas produção + 291 linhas testes (6 testes)
  - ✅ **Pipeline completo:** Bethe → Tables → Splines → Potentials → Hamiltonian → Diagonalization → **Density!**
  - 🔜 Próximo: `scf_mixer.f90` (mixing schemes) + `ks_cycle.f90` (loop SCF)

  **Total do Projeto:** 164 testes, 100% passando! 🎉

---

### 2025-01-15 - Fase 5: Diagonalização LAPACK & Degenerescências! 🎉
- ✅ **MILESTONE:** Diagonalização de matrizes simétricas/Hermitianas completa!

  **`lapack_wrapper.f90` implementado** (347 linhas, 18 testes):
  - ✅ `validate_diagonalization_inputs()`: Validação de dimensões
    - L > 0, size(H) == (L,L), size(eigvals) == L, size(eigvecs) == (L,L)
  - ✅ `diagonalize_symmetric_real()`: Wrapper DSYEVD para matrizes reais simétricas
    - Calcula eigenvalues E eigenvectors
    - Eigenvalues retornados em ordem crescente (ground state = E₁)
    - Eigenvectors normalizados e ortogonais
  - ✅ `diagonalize_symmetric_real_values_only()`: Eigenvalues only (mais rápido ~2x)
  - ✅ `diagonalize_hermitian_complex()`: Wrapper ZHEEVD para Hermitianas complexas
    - Eigenvalues são SEMPRE reais (teorema fundamental QM)
    - Suporta twisted boundary conditions (Hamiltoniano complexo)
  - ✅ `diagonalize_hermitian_complex_values_only()`: Eigenvalues only
  - ✅ Workspace query em duas fases:
    - Fase 1: lwork=-1, query optimal workspace size
    - Fase 2: allocate(work(lwork)), chamar LAPACK novamente
  - ✅ **Bug fix crítico**: Removido `bind(C)` das interfaces LAPACK
    - LAPACK usa convenção Fortran nativa, não C!
    - `bind(C)` causava falha em workspace query no gfortran

  **`degeneracy_handler.f90` implementado** (405 linhas, 13 testes):
  - ✅ `find_degenerate_subspaces()`: Detecta grupos degenerados
    - Varre eigenvalues, identifica grupos onde |λᵢ - λⱼ| < DEG_TOL (1.0e-8)
    - Retorna array 2D: subspaces(n_subspaces, max_deg)
    - Exemplo: eigenvalues [1, 2, 2, 3, 3, 3] → 2 subspaces: [2,3] e [4,5,6]
  - ✅ `orthonormalize_degenerate_subspace()`: QR decomposition para vetores reais
    - Usa DGEQRF (QR factorization) + DORGQR (generate Q)
    - Mais estável numericamente que Gram-Schmidt
    - LAPACK handles workspace automaticamente via query
  - ✅ `orthonormalize_degenerate_subspace_complex()`: Gram-Schmidt modificado
    - Para vetores complexos (e.g., twisted BC)
    - Modified Gram-Schmidt: v_k ⊥ span{v₁,...,v_{k-1}} iterativamente
  - ✅ `compute_degeneracy_count()`: Conta quantos eigenvalues são degenerados com índice dado
  - ✅ `verify_orthonormality()`: Verifica ||V^T V - I||∞ < tol usando DGEMM
  - ✅ **Bug fix crítico 1**: Removido double conjugation
    - `DOT_PRODUCT(a,b)` em Fortran JÁ FAZ `SUM(CONJG(a)*b)` para vetores complexos!
    - Estava fazendo `DOT_PRODUCT(CONJG(a),b)` = `SUM(a*b)` → ERRADO
  - ✅ **Bug fix crítico 2**: Workspace query separada para DORGQR
    - DGEQRF e DORGQR podem precisar workspaces de tamanhos diferentes!
    - Antes: usava lwork do DGEQRF para DORGQR → falha em alguns sistemas
    - Agora: query separada para cada rotina LAPACK

  **Física dos Eigenproblemas:**
  - ✅ **Real simétrico**: H = H^T (open/periodic BC sem campo magnético)
  - ✅ **Complexo Hermitiano**: H = H† (twisted BC, Aharonov-Bohm phase)
  - ✅ **Degenerescências**: Ocorrem por simetrias (translação, paridade, spin)
    - Exemplo: PBC com L=10 → eigenvalues vêm em pares ±k (exceto k=0, L/2)
  - ✅ **Orthonormalização**: LAPACK pode retornar base arbitrária no subespaço degenerado
    - QR/Gram-Schmidt garante base ortonormal canônica

  **Testes implementados** (31 novos testes, 100% passando):
  - ✅ 18 testes `test_lapack_wrapper.f90`:
    - Validação de inputs (6 testes)
    - Diagonalização real (7 testes): identity, diagonal, 2×2 analítico, tridiagonal tight-binding, ordering, normalização, ortogonalidade
    - Eigenvalues only (1 teste)
    - Complexo Hermitiano (4 testes): identity, 2×2 analítico, eigenvalues reais, values only
  - ✅ 13 testes `test_degeneracy_handler.f90`:
    - Detecção de degenerescências (5 testes): nenhuma, par, tripla, múltiplos grupos, todos degenerados
    - Contagem (3 testes): single, par, tripla
    - Orthonormalization (3 testes): real pair/triple, complex pair
    - Verificação (2 testes): identity perfeita, detectar não-ortogonalidade

  **Estatísticas Fase 5:**
  - ✅ Total: 1228 linhas código produção + 1531 linhas testes
  - ✅ 66 testes (100% passando)
  - ✅ 4 módulos completos: boundary_conditions, hamiltonian_builder, lapack_wrapper, degeneracy_handler
  - ✅ Pipeline completo: Bethe Ansatz → Tables → Splines → Potentials → Hamiltonian → **Diagonalization!**
  - 🔜 Próximo: symmetry.f90 (explorar simetria de paridade para 4x speedup)

  **Total do Projeto:** 158 testes, 100% passando! 🎉

---

### 2025-01-14 - Fase 5: Hamiltoniano & Boundary Conditions! 🎉
- ✅ **MILESTONE:** Construção do Hamiltoniano tight-binding completa!

  **`boundary_conditions.f90` implementado** (256 linhas, 17 testes):
  - ✅ Enum com tipos de BC: `BC_OPEN`, `BC_PERIODIC`, `BC_TWISTED`
  - ✅ `validate_bc_parameters()`: Validação completa
    - BC type válido (1, 2, ou 3)
    - Sistema com L > 1 (mínimo 2 sites para tight-binding)
    - Para BC_TWISTED: theta obrigatório e em [0, 2π)
  - ✅ `apply_boundary_conditions()`: BCs para matrizes reais
    - BC_OPEN: sem modificação (H já é tridiagonal)
    - BC_PERIODIC: H(1,L) = H(L,1) = -1 (cria anel)
    - BC_TWISTED: retorna erro (use versão complexa)
  - ✅ `apply_boundary_conditions_complex()`: BCs para matrizes complexas
    - BC_OPEN: sem modificação
    - BC_PERIODIC: H(1,L) = H(L,1) = -1
    - BC_TWISTED: H(1,L) = -exp(iθ), H(L,1) = -exp(-iθ) (efeito Aharonov-Bohm)
  - ✅ `get_free_particle_eigenvalues()`: Eigenvalues analíticos para validação
    - OBC: E_n = -2cos(nπ/(L+1)), n=1,...,L (standing waves)
    - PBC: E_k = -2cos(2πk/L), k=0,...,L-1 (Bloch waves)
    - TBC: E_k(θ) = -2cos((2πk+θ)/L), k=0,...,L-1 (shifted spectrum)

  **Física das Boundary Conditions:**
  - ✅ **OBC**: Hard-wall boundaries, edge states, confinamento quântico
  - ✅ **PBC**: Conservação de momento, propriedades bulk, Bethe Ansatz
  - ✅ **TBC**: Persistent currents, efeito Aharonov-Bohm, flux threading
  - ✅ Antiperiodic BC (θ=π): meio quantum de fluxo, quebra degenerescências

  **`hamiltonian_builder.f90` implementado** (220 linhas, 18 testes):
  - ✅ `validate_hamiltonian_inputs()`: Validação robusta
    - L > 0 (sistema físico)
    - size(V_ext) == size(V_xc) == L
    - NaN/Inf checking usando `ieee_is_finite()` (importado de `ieee_arithmetic`)
  - ✅ `build_hamiltonian()`: Construção H real
    - Diagonal: H(i,i) = V_ext(i) + V_xc(i) (on-site energies)
    - Off-diagonal: H(i,i±1) = -t = -1 (hopping)
    - Aplica BCs via `apply_boundary_conditions()`
  - ✅ `build_hamiltonian_complex()`: Construção H complexo
    - Similar ao real mas com tipo `complex(dp)`
    - Suporta BC_TWISTED com theta
    - Diagonal sempre real (potenciais on-site)
  - ✅ `build_hamiltonian_free()`: H livre (U=0, V=0)
    - Apenas hopping, sem potenciais
    - Útil para validação contra eigenvalues analíticos
  - ✅ `compute_effective_potential()`: V_eff = V_ext + V_xc
    - Helper function simples
    - Validação de size matching

  **Bug crítico corrigido:**
  - ❌ **Problema:** Todas as 3 funções tinham loop incorreto:
    ```fortran
    do i = 1, L - 1
        if (i > 1) then  ! ❌ Pula i=1!
            H(i,i+1) = -1.0_dp
            H(i+1,i) = -1.0_dp
        end if
    end do
    ```
  - ✅ **Solução:** Removido `if (i > 1)` em todas as funções
  - ✅ **Impacto:** Sem o fix, H(1,2) e H(2,1) nunca eram setados → Hamiltoniano incorreto!
  - ✅ Detectado pelos testes `test_build_free_hamiltonian_open_bc` e `test_build_hamiltonian_offdiagonal`

  **Correção no sistema de erros:**
  - ✅ `ERROR_NOT_A_NUMBER` adicionado aos exports públicos de `lsda_errors.f90`
  - ✅ Mensagem de erro adicionada: "Array contains NaN or Inf values"

  **Testes implementados** (866 linhas, 35 testes):
  - ✅ **`test_boundary_conditions.f90`** (433 linhas, 17 testes):
    - Validação de BC: open, periodic, twisted (com/sem theta, ranges)
    - Aplicação de BC: open (no-op), periodic (edges), twisted (complex phases)
    - Antiperiodic (θ=π): H(1,L) = H(L,1) = +1
    - Free particle eigenvalues: OBC, PBC, TBC (validação analítica)
    - Size mismatch detection
  - ✅ **`test_hamiltonian_builder.f90`** (433 linhas, 18 testes):
    - Validação: inputs válidos, L inválido, size mismatches, NaN/Inf detection
    - Effective potential: computação e size mismatch
    - Free Hamiltonian: estrutura tridiagonal, OBC (sem edges), PBC (com edges)
    - Full Hamiltonian: diagonal (V_ext+V_xc), off-diagonal (hopping), BCs, simetria hermitiana
    - Complex Hamiltonian: diagonal real para potenciais on-site
    - Error handling completo

  **Física validada:**
  - ✅ Estrutura tridiagonal do tight-binding
  - ✅ Hermitianidade (H† = H → H simétrica para H real)
  - ✅ Dispersion relation E(k) = -2cos(k) para free particles
  - ✅ Boundary effects: OBC vs PBC vs TBC
  - ✅ Aharonov-Bohm phase em TBC

  **Estatísticas Fase 5 (parcial):**
  - **Código produção:** 476 linhas (2 módulos completos)
  - **Testes:** 866 linhas (35 testes, 100% passando)
  - **Total do projeto:** 110 testes (antes 92 + 18 novos)
  - **Próximo:** `symmetry.f90` para exploração de paridade

  **🎯 Próximo Passo:**
  - `src/hamiltonian/symmetry.f90`: Simetria de paridade
    - Para V(i) = V(L+1-i): H se block-diagonaliza em setores even/odd
    - Cada setor tem dimensão L/2 → speedup 4x na diagonalização
    - Funções: `check_parity_symmetry()`, `block_diagonalize_hamiltonian()`, `reconstruct_eigenstates()`

### 2025-01-13 (Parte 2) - Potencial Quasiperiódico Adicionado! 🎉
- ✅ **Correção e implementação do potencial quasiperiódico (AAH model)**

  **`potential_quasiperiodic.f90` reescrito** (100 linhas):
  - ✅ Corrigido padrão de erro: substituído tipo customizado `ErrorHandler` por `integer :: ierr`
  - ✅ Unificado em uma única subroutine: `apply_potential_quasiperiodic(lambda, beta, phi, L, V, ierr)`
  - ✅ Validação de parâmetros seguindo padrão do projeto
  - ✅ Uso de constantes corretas: `TWOPI` de `lsda_constants.f90`
  - ✅ Fórmula: V(i) = λ·cos(2πβi + φ) com i-1 para indexação física começar em 0

  **Física do Aubry-André-Harper (AAH):**
  - ✅ Extended phase (λ < 2): Estados deslocalizados
  - ✅ Critical phase (λ = 2): Funções de onda multifractais (transição metal-isolante)
  - ✅ Localized phase (λ > 2): Estados exponencialmente localizados
  - ✅ Golden ratio β = (√5-1)/2 para máxima incomensurabilidade
  - ✅ Modela localização de Anderson sem desordem

  **4 Testes quasiperiódicos adicionados** (83 linhas):
  - ✅ `test_quasiperiodic_golden_ratio`: Testa bounds [-λ, λ] e variação do potencial
  - ✅ `test_quasiperiodic_phase_shift`: Verifica que φ = π inverte o potencial
  - ✅ `test_quasiperiodic_critical_point`: Testa λ = 2 (ponto crítico)
  - ✅ `test_quasiperiodic_localization`: Testa λ = 5 (regime localizado)

  **`potential_factory.f90` atualizado** (186 linhas):
  - ✅ Adicionado suporte para "quasiperiodic" com 3 parâmetros: [lambda, beta, phi]
  - ✅ Documentação e info string completas
  - ✅ Factory agora suporta 7 tipos de potenciais

  **Estatísticas da atualização:**
  - **Código produção:** +100 linhas (potential_quasiperiodic) + 13 linhas (factory)
  - **Testes:** +83 linhas (4 novos testes)
  - **Total Fase 4 atualizado:** 1090 linhas produção + 869 linhas testes
  - **Total de testes:** 92 (antes 88 + 4 novos)

### 2025-01-13 (Parte 1) - Fase 4: COMPLETA! 🎉
- ✅ **MILESTONE:** Sistema de potenciais e erros totalmente funcional!

  **Módulo `lsda_errors.f90` implementado** (224 linhas, 13 testes):
  - ✅ Códigos de erro organizados: input (1-99), numerical (100-199), I/O (200-299), memory (300-399)
  - ✅ `get_error_message()` - Mensagens legíveis para cada código
  - ✅ `error_handler()` - Handler centralizado com opção fatal
  - ✅ `check_bounds()`, `check_positive()`, `check_range()` - Utilitários de validação
  - ✅ Integração com todos os módulos de potenciais

  **6 Módulos de potenciais base implementados** (753 linhas produção):
  - ✅ **`potential_uniform.f90`** (34 linhas): V(i) = V₀
  - ✅ **`potential_harmonic.f90`** (46 linhas): V(i) = 0.5·k·(i-center)² (optical traps)
  - ✅ **`potential_impurity.f90`** (191 linhas): single/multiple/random impurities
  - ✅ **`potential_random.f90`** (152 linhas): uniform/Gaussian disorder (Anderson localization)
  - ✅ **`potential_barrier.f90`** (157 linhas): single/double barriers (quantum tunneling)
  - ✅ **`potential_factory.f90`** (173 linhas): Factory pattern para criação dinâmica

  **Testes implementados** (786 linhas, 30 testes):
  - ✅ **`test_potentials.f90`** (502 linhas, 17 testes):
    - Uniform: constância em todos os sites
    - Harmonic: simetria de paridade, mínimo central
    - Impurity: posicionamento, bounds, overlap, concentração aleatória
    - Random: média zero, distribuições corretas (uniform/Gaussian)
    - Barrier: largura, bounds, separação do poço quântico, não-sobreposição
    - Factory: criação via string, comparação com chamadas diretas
  - ✅ **`test_lsda_errors.f90`** (284 linhas, 13 testes):
    - Códigos em intervalos corretos
    - Mensagens para todos os tipos
    - Utilitários de validação

  **Física Implementada:**
  - ✅ Armadilhas harmônicas (cold atoms, optical traps)
  - ✅ Localização de Anderson (random disorder, W/t regime)
  - ✅ Tunelamento quântico (barriers, T ~ exp(-2κw))
  - ✅ Ressonâncias Fabry-Pérot (double barriers, quasi-bound states)
  - ✅ Impurezas magnéticas diluídas

  **Correções Técnicas:**
  - ✅ Renomeados `potential_uniform()` → `apply_potential_uniform()` para evitar conflito de nomes
  - ✅ Renomeados `potential_harmonic()` → `apply_potential_harmonic()` para evitar conflito de nomes
  - ✅ Adicionado parâmetro `ierr` em uniform e harmonic para consistência

  **Estatísticas Fase 4 (inicial):**
  - **Código produção:** 977 linhas (7 módulos base)
  - **Testes:** 786 linhas (30 testes, 100% passando)
  - **Física:** 6 tipos de potenciais com explicações detalhadas nos testes

  **🎉 GRAND TOTAL (Fases 1+2+3+4 - após quasiperiodic):**
  - **15 módulos produção:** 3635 linhas
  - **2 executáveis:** 208 linhas (main.f90 + convert_tables.f90)
  - **9 suítes de testes:** 2665 linhas, 92 testes (100% passando)
  - **Total geral:** ~6508 linhas de código

### 2025-01-12 - Fase 3: COMPLETA! 🎉
- ✅ **MILESTONE:** Pipeline XC totalmente funcional de ponta a ponta!

  **Módulo `spline2d.f90` implementado** (351 linhas, 5 testes):
  - ✅ Tipo `spline2d_t` para grids irregulares 2D (n_y varia com x)
  - ✅ `spline1d_coeff()` - Algoritmo de Thomas para splines cúbicas 1D
  - ✅ `spline2d_init()` - Construção de splines separáveis
  - ✅ `spline2d_eval()` - Avaliação em (x, y) via interpolação separável
  - ✅ `find_interval()` - Busca binária para localização no grid
  - ✅ Tratamento de boundary conditions (natural e clamped)
  - ✅ Arrays 0-indexed internos para compatibilidade com algoritmo clássico

  **Módulo `xc_lsda.f90` implementado** (335 linhas, 6 testes):
  - ✅ Tipo `xc_lsda_t` com 3 splines: exc, vxc_up, vxc_dw
  - ✅ `xc_lsda_init()` - Carregamento de tabela e construção de splines
  - ✅ `get_exc(n_up, n_dw)` - Energia XC por partícula
  - ✅ `get_vxc(n_up, n_dw, v_xc_up, v_xc_dw)` - Potenciais XC
  - ✅ **4 regiões de simetria física:**
    - Region I (m≥0, n≤1): Identidade
    - Region II (m<0, n≤1): Spin exchange
    - Region III (m<0, n>1): Particle-hole
    - Region IV (m≥0, n>1): Combinada
  - ✅ `determine_region()`, `apply_symmetry_transform()`, `convert_to_nm()`
  - ✅ Casos especiais: U=0, n=0, bounds checking
  - ✅ **Padronização de nomenclatura:** n_dw (não n_dn) para spin-down

  **Testes implementados** (396 linhas, 11 testes):
  - ✅ `test_spline2d.f90` (5 testes): Init/destroy, interpolação exata, funções lineares/separáveis
  - ✅ `test_xc_lsda.f90` (6 testes): Init/destroy, simetrias de spin, regiões, transformações

  **Estatísticas Fase 3:**
  - **Código produção:** 686 linhas (spline2d + xc_lsda)
  - **Testes:** 396 linhas (11 testes, 100% passando)
  - **Pipeline completo:** Bethe Ansatz → Tabelas → Splines → XC funcional pronto para KS!

  **🎉 GRAND TOTAL (Fases 1+2+3):**
  - **7 módulos produção:** 2545 linhas
  - **2 executáveis:** 208 linhas (main.f90 + convert_tables.f90)
  - **7 suítes de testes:** 1796 linhas, 58 testes (100% passando)
  - **Total geral:** ~4549 linhas de código

### 2025-01-11 - Fase 2: COMPLETA! 🎉
- ✅ **MILESTONE:** Geração de tabelas XC totalmente funcional!

  **Módulo `bethe_tables.f90` implementado** (325 linhas, 6 testes):
  - ✅ Tipo `grid_params_t` para configurar grid de densidades
  - ✅ `compute_E0()` - Energia cinética não-interagente (Fermi gas livre)
  - ✅ `compute_E_xc()` - Energia XC via Bethe Ansatz: E_xc = E_BA - E_0
  - ✅ `compute_V_xc_numerical()` - Potenciais XC via derivadas de 5 pontos
  - ✅ `generate_xc_table()` - Geração completa de tabela para dado U
  - ✅ `generate_table_grid()` - Grid flexível com parâmetros customizados
  - ✅ Casos especiais: U=0, half-filling, polarizado

  **Testes `test_bethe_tables.f90`** (170 linhas, 6 testes):
  - ✅ E0 para half-filling e sistema polarizado
  - ✅ E_xc = 0 para U=0 (validação física)
  - ✅ Simetria V_xc (V_up = V_dn quando n_up = n_dn)
  - ✅ Parâmetros padrão do grid
  - ✅ Geração de tabela pequena

  **Estatísticas Fase 2:**
  - **Código produção:** 888 linhas (table_io + bethe_tables + convert_tables)
  - **Testes:** 444 linhas (16 testes, 100% passando)
  - **Pipeline completo:** Bethe Ansatz → E_xc → V_xc → Tabela → I/O binário

### 2025-01-09 - Fase 2: I/O de Tabelas Completo ✅
- ✅ Sistema de I/O de tabelas totalmente funcional!
  - **`table_io.f90`** (364 linhas, 10 testes): Leitura/escrita ASCII/binário
  - **`convert_tables.f90`** (199 linhas): Utilitário de conversão batch
  - **`test_table_io.f90`** (274 linhas): Validação roundtrip completa

### 2025-11-07 - Fase 1: COMPLETA ✅
- ✅ **MILESTONE:** Fase 1 totalmente concluída!
  - **3 módulos completos:** `bethe_equations.f90`, `nonlinear_solvers.f90`, `continuation.f90`
  - **1159 linhas de código produção** (487 + 303 + 369)
  - **946 linhas de testes** (446 + 302 + 198)
  - **31 testes unitários passando** (17 + 9 + 5)

- ✅ **`continuation.f90`** (369 linhas, 5 testes):
  - Implementado predictor-corrector method
  - Sweep forward, backward e bidirectional
  - Speedup típico de 5-10x vs soluções independentes
  - Validado com testes de consistência

- ✅ **Refinamentos finais:**
  - Derivadas dθ/dU e dΘ/dU implementadas
  - `compute_dFdU()` para continuation method
  - `compute_energy()` para cálculo de E = -2·Σcos(k)
  - Documentação FORD-compliant completa

### 2025-11-06 - Fase 1: Continuation + Newton ✅
- ✅ **`nonlinear_solvers.f90`** (303 linhas, 9 testes):
  - Newton-Raphson com line search
  - Wrapper LAPACK DGESV
  - Tratamento robusto: U=0, NaN, matrizes singulares
  - 100% dos testes passando

### 2025-11-05 - Fase 1: Bethe Equations ✅
- ✅ **`bethe_equations.f90`** (487 linhas, 17 testes):
  - Funções θ e Θ completas
  - Jacobiano analítico validado (< 1e-10)
  - Números quânticos inicializados
  - 100% dos testes passando

### 2025-11-03 - Fase 0: Infraestrutura ✅
- ✅ Criada estrutura completa do projeto com fpm
- ✅ Módulos base: `lsda_types`, `lsda_constants`
- ✅ Fortuno configurado e funcionando

---

**Este documento é vivo e deve ser atualizado conforme o projeto evolui!** 🚀