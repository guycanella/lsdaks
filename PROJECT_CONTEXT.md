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
│   │   └── lsda_constants.f90 # Constantes físicas e numéricas
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
│   ├── xc_functional/         # 🔜 TODO
│   │   ├── spline2d.f90       # Interpolação bicúbica 2D
│   │   └── xc_lsda.f90        # Interface exc, Vxc_up, Vxc_dn
│   │
│   ├── potentials/            # 🔜 TODO
│   │   ├── potential_base.f90      # Classe abstrata base
│   │   ├── potential_uniform.f90   # Potencial uniforme
│   │   ├── potential_harmonic.f90  # Armadilha harmônica
│   │   ├── potential_impurity.f90  # Impurezas
│   │   ├── potential_random.f90    # Aleatório
│   │   ├── potential_barrier.f90   # Barreiras periódicas/duplas
│   │   └── potential_factory.f90   # Factory pattern
│   │
│   ├── hamiltonian/           # 🔜 TODO
│   │   ├── hamiltonian_builder.f90 # Tight-binding com Veff
│   │   ├── boundary_conditions.f90 # Open, periodic, twisted
│   │   └── symmetry.f90            # Exploração de simetria de paridade
│   │
│   ├── diagonalization/       # 🔜 TODO
│   │   ├── lapack_wrapper.f90      # Interface para DSYEV/DSYEVD
│   │   └── degeneracy_handler.f90  # Tratamento de níveis degenerados
│   │
│   ├── density/               # 🔜 TODO
│   │   ├── density_calculator.f90  # Ocupação de níveis
│   │   └── fermi_distribution.f90  # Distribuição de Fermi
│   │
│   ├── convergence/           # 🔜 TODO
│   │   ├── convergence_monitor.f90 # Critérios de parada
│   │   └── mixing_schemes.f90      # Mixing linear, Broyden, etc
│   │
│   └── kohn_sham/             # 🔜 TODO
│       └── ks_cycle.f90       # Loop SCF completo
│
├── app/                        # 🔄 EM PROGRESSO
│   ├── main.f90               # Ponto de entrada (placeholder)
│   └── convert_tables.f90     # ✅ COMPLETO - Utilitário conversão tabelas
│
├── test/                       # 🔄 EM PROGRESSO (47 testes, 100% passando)
│   ├── test_bethe_equations.f90      # ✅ COMPLETO - 17 testes
│   ├── test_nonlinear_solvers.f90    # ✅ COMPLETO - 9 testes
│   ├── test_continuation.f90         # ✅ COMPLETO - 5 testes
│   ├── test_table_io.f90             # ✅ COMPLETO - 10 testes
│   ├── test_bethe_tables.f90         # ✅ COMPLETO - 6 testes (NEW!)
│   ├── test_splines.f90              # 🔜 TODO
│   ├── test_potentials.f90           # 🔜 TODO
│   ├── test_hamiltonian.f90          # 🔜 TODO
│   └── test_ks_cycle.f90             # 🔜 TODO
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

### Fase 3: Interpolação de Splines 2D (3-4 dias) 🔜 PRÓXIMA PRIORIDADE

**Objetivo:** Implementar interpolação bicúbica 2D para avaliar funcionais XC em pontos arbitrários (n, m).

#### 📋 Tarefas:
- [ ] **`spline2d.f90`**: Interpolação bicúbica separável
  - [ ] Tipo `spline2d_t` para armazenar coeficientes da spline
  - [ ] `spline2d_init()` - Construir spline a partir de tabela
  - [ ] `spline2d_eval()` - Avaliar spline em ponto (n, m)
  - [ ] Derivadas analíticas da spline (para forças atômicas)
  - [ ] Tratamento de pontos fora do grid (extrapolação ou erro)

- [ ] **`xc_lsda.f90`**: Interface de alto nível para funcional XC
  - [ ] Tipo `xc_lsda_t` contendo splines de exc, vxc_up, vxc_dn
  - [ ] `xc_lsda_init()` - Carregar tabela e construir splines
  - [ ] `get_exc(n_up, n_dn)` - Energia XC por partícula
  - [ ] `get_vxc(n_up, n_dn)` - Potenciais XC para ambos os spins
  - [ ] Conversão (n_up, n_dn) ↔ (n, m)
  - [ ] Cache de último ponto avaliado (otimização)

- [ ] **`test_spline2d.f90`**: Testes de interpolação
  - [ ] Interpolação exata em pontos do grid
  - [ ] Continuidade C² (derivadas suaves)
  - [ ] Comparação com Bethe Ansatz em pontos intermediários
  - [ ] Teste de performance (tempo de eval)

- [ ] **`test_xc_lsda.f90`**: Testes de interface
  - [ ] Carregar tabela real e interpolar
  - [ ] Verificar simetrias físicas (troca de spin)
  - [ ] Casos limite: U→0, n→0, m→±n

**Estratégia de Implementação:**
1. Spline 1D em cada direção (separável)
2. Para cada n fixo: spline cúbica em m
3. Com valores interpolados: spline cúbica em n
4. Derivadas analíticas via regra da cadeia

**Validação:**
- Erro de interpolação < 1e-6 comparado com Bethe Ansatz exato
- Benchmark: eval < 1μs por ponto

---

### Fase 4: Hamiltoniano Básico (2-3 dias) 🔜 TODO

- [ ] `potential_uniform.f90`, `potential_harmonic.f90`
- [ ] `hamiltonian_builder.f90`: Tight-binding + Veff
- [ ] `lapack_wrapper.f90`: Interface DSYEVD
- [ ] `boundary_conditions.f90`: Open, periodic
- [ ] Teste end-to-end simples (1 iteração KS)

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

**Versão:** 0.2.0-dev
**Status:** ✅ Fases 1 & 2 Completas → 🔜 Iniciando Fase 3 (Splines 2D)
**Última atualização:** 2025-01-11

### Progresso Geral

```
[████████████████████████████████] 100% Fase 1: Bethe Ansatz Core (COMPLETO ✅)
[████████████████████████████████] 100% Fase 2: Geração de Tabelas XC (COMPLETO ✅)
[░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░]   0% Fase 3: Splines 2D (PRÓXIMA 🔜)
[░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░]   0% Fase 4: Hamiltoniano
[░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░]   0% Fase 5: Ciclo KS
[░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░]   0% Fase 6: Features
[░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░]   0% Fase 7: Otimização
```

### Checklist de Progresso

#### Infraestrutura ✅
- [x] Estrutura fpm
- [x] Tipos básicos (`lsda_types.f90`)
- [x] Constantes (`lsda_constants.f90`)
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

- [ ] **Fase 3 - Splines 2D** (0% 🔜):
  - [ ] Interpolação bicúbica (`spline2d.f90`)
  - [ ] Interface XC funcional (`xc_lsda.f90`)

- [ ] **Fases 4-7**: Hamiltoniano, Ciclo KS, Features, Otimização

#### Features 🔜
- [ ] Potenciais
- [ ] Simetria
- [ ] Twisted BC
- [ ] Degenerescências

#### Qualidade ✅
- [x] Testes unitários Fase 1 (31 testes, 100% passando) ✅
- [x] Testes unitários Fase 2 (16 testes, 100% passando) ✅
- [x] **Total: 47 testes, 100% passando** ✅
- [x] Pipeline Bethe → Tabelas validado ✅
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
**Última atualização:** 2025-01-11
**Status:** Fases 1 & 2 Completas (100%) → Iniciando Fase 3 (Splines 2D)

---

## 📅 Histórico de Mudanças

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