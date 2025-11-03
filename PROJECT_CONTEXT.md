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
│   ├── types/                  # Tipos derivados e estruturas de dados
│   │   ├── lsda_types.f90     # Tipos principais (SystemParams, State, etc)
│   │   └── lsda_constants.f90 # Constantes físicas e numéricas
│   │
│   ├── io/                     # Entrada/Saída
│   │   ├── input_parser.f90   # Parse de argumentos e arquivos
│   │   ├── output_writer.f90  # Escrita de resultados
│   │   └── logger.f90         # Sistema de logging
│   │
│   ├── bethe_ansatz/          # Solução das equações de Lieb-Wu
│   │   ├── bethe_equations.f90      # Definição de F(x)=0 e Jacobiano
│   │   ├── nonlinear_solvers.f90    # Newton, Broyden, híbridos
│   │   ├── continuation.f90         # Sweep em U com preditor-corretor
│   │   └── bethe_tables.f90         # Geração/cache de tabelas
│   │
│   ├── xc_functional/         # Funcional de troca-correlação
│   │   ├── spline2d.f90       # Interpolação bicúbica 2D
│   │   └── xc_lsda.f90        # Interface exc, Vxc_up, Vxc_dn
│   │
│   ├── potentials/            # Potenciais externos
│   │   ├── potential_base.f90      # Classe abstrata base
│   │   ├── potential_uniform.f90   # Potencial uniforme
│   │   ├── potential_harmonic.f90  # Armadilha harmônica
│   │   ├── potential_impurity.f90  # Impurezas
│   │   ├── potential_random.f90    # Aleatório
│   │   ├── potential_barrier.f90   # Barreiras periódicas/duplas
│   │   └── potential_factory.f90   # Factory pattern
│   │
│   ├── hamiltonian/           # Construção de Hamiltoniano
│   │   ├── hamiltonian_builder.f90 # Tight-binding com Veff
│   │   ├── boundary_conditions.f90 # Open, periodic, twisted
│   │   └── symmetry.f90            # Exploração de simetria de paridade
│   │
│   ├── diagonalization/       # Solvers de autovalores
│   │   ├── lapack_wrapper.f90      # Interface para DSYEV/DSYEVD
│   │   └── degeneracy_handler.f90  # Tratamento de níveis degenerados
│   │
│   ├── density/               # Cálculo de densidades eletrônicas
│   │   ├── density_calculator.f90  # Ocupação de níveis
│   │   └── fermi_distribution.f90  # Distribuição de Fermi
│   │
│   ├── convergence/           # Gerenciamento de convergência
│   │   ├── convergence_monitor.f90 # Critérios de parada
│   │   └── mixing_schemes.f90      # Mixing linear, Broyden, etc
│   │
│   └── kohn_sham/             # Ciclo auto-consistente principal
│       └── ks_cycle.f90       # Loop SCF completo
│
├── test/                       # Testes unitários (Fortuno)
│   ├── test_types.f90
│   ├── test_bethe_ansatz.f90
│   ├── test_splines.f90
│   ├── test_potentials.f90
│   ├── test_hamiltonian.f90
│   └── test_ks_cycle.f90
│
├── examples/                   # Exemplos de uso
│   ├── harmonic_trap.f90
│   ├── double_barrier.f90
│   └── half_filling.f90
│
└── data/                       # Dados auxiliares
    ├── potential_params/       # Parâmetros de potenciais
    └── reference_results/      # Resultados de referência (validação)
```

---

## 🧮 Física do Problema

### Modelo de Hubbard 1D

Hamiltoniano:
```
H = -t Σᵢⱼσ (cᵢσ† cⱼσ + h.c.) + U Σᵢ nᵢ↑ nᵢ↓ + Σᵢσ Vᵢˢᵗ nᵢσ
```

- **t = 1**: Hopping (unidade de energia)
- **U**: Interação on-site (Hubbard U)
- **Vᵢˢᵗ**: Potencial externo (armadilha, impurezas, etc)

### Teoria do Funcional da Densidade (DFT)

Kohn-Sham equations:
```
[-∇² + Vₑₓₜ(r) + VH(r) + Vxc(r)] ψᵢ(r) = εᵢ ψᵢ(r)
```

Para o Hubbard 1D:
```
Hₖₛ = H₀ + Vₑₓₜ + U·n₋σ + Vxc
```

### Bethe Ansatz de Lieb-Wu

Solução exata para o estado fundamental via Ansatz de Bethe:

**Equações integrais:**
```
k_j = (2π/L)·Iⱼ + (1/L) Σ_α θ(k_j - Λ_α)
(2π/L)·Jα = Σⱼ θ(Λα - k_j) - Σ_β Θ(Λα - Λ_β)
```

Onde:
- `θ(x) = 2·atan(2x/U)` (espalhamento carga-spin)
- `Θ(x) = 2·atan(x/U)` (espalhamento spin-spin)
- `{Iⱼ}`: números quânticos de carga (N↑ inteiros/semi-inteiros)
- `{Jα}`: números quânticos de spin (N↓-1 inteiros/semi-inteiros)

**Energia:**
```
E = -2 Σⱼ cos(k_j)
```

**Funcional XC:**
```
Exc = E_BA[n↑, n↓] - E₀[n↑, n↓]
Vxc↑ = ∂Exc/∂n↑
Vxc↓ = ∂Exc/∂n↓
```

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
author = "Your Name"
maintainer = "your.email@example.com"

[build]
auto-executables = true
auto-tests = true
auto-examples = true

[dependencies]
fortuno = "*"

[dev-dependencies]

[library]
source-dir = "src"

[executable]
name = "lsda"
source-dir = "src"
main = "lsda_main.f90"
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
1. Para (n_up, n_down) → calcular (n=n↑+n↓, m=n↑-n↓)
2. Para cada `n_i` fixo: interpolação spline 1D em `m`
3. Com valores interpolados: segunda interpolação spline 1D em `n`

**Vantagens:**
- Código limpo e moderno (sem dependências pesadas)
- Facilmente testável
- Totalmente controlado
- C² contínuo

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
- Arrays 0-indexed quando possível (padrão Fortran moderno)
- 1-indexed apenas quando fisicamente motivado (sítios da rede)

**Precisão:**
```fortran
use, intrinsic :: iso_fortran_env, only: real64
real(real64) :: variable  ! Float64 (double precision)
```

**Nomes:**
- Módulos: `snake_case` (ex: `bethe_ansatz`)
- Tipos: `snake_case_t` (ex: `system_params_t`)
- Funções/subrotinas: `snake_case` (ex: `solve_newton`)
- Constantes: `UPPER_SNAKE_CASE` (ex: `MAX_ITER`)

**Documentação:**
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

**Variáveis:** `x = [k₁, k₂, ..., k_N↑, Λ₁, Λ₂, ..., Λ_M]` onde M = N↓

**Sistema não-linear:** `F(x) = 0`

```fortran
F_j^k = k_j - (2π/L)·I_j - (1/L) Σ_α θ(k_j - Λ_α)

F_α^Λ = (2π/L)·J_α - Σ_j θ(Λ_α - k_j) + Σ_β Θ(Λ_α - Λ_β)
```

**Jacobiano analítico:**
```
J = [ ∂F^k/∂k    ∂F^k/∂Λ  ]
    [ ∂F^Λ/∂k    ∂F^Λ/∂Λ  ]
```

Com derivadas:
```fortran
θ'(x) = 4U / (U² + 4x²)
Θ'(x) = 2U / (U² + x²)
```

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

**Ideia:** Aproximar J⁻¹ iterativamente sem recalcular Jacobiano

**Atualização de Broyden:**
```
B_{k+1} = B_k + (Δx - B_k·ΔF) ⊗ ΔF^T / (ΔF^T·ΔF)
```

Onde `B ≈ J⁻¹` (inversa aproximada).

**Vantagens:**
- Não precisa calcular J a cada iteração
- Custo O(n²) por iteração (vs O(n³) do Newton)

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

**Problema:** k ∈ [-π,π], Λ ~ O(U) → mal-condicionado para U grande

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

Para N↑=5: I = [-2, -1, 0, 1, 2]  
Para N↓=3: J = [-1, 0, 1]

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

Grid (n, m, U) é **embaraçosamente paralelo**:

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

### Fase 0: Infraestrutura (1-2 dias) ✅
- [x] Criar estrutura fpm
- [x] Módulo de tipos (`lsda_types.f90`)
- [x] Módulo de constantes (`lsda_constants.f90`)
- [x] Configurar Fortuno

### Fase 1: Bethe Ansatz (1-2 semanas) 🔄
- [ ] `bethe_equations.f90`: F(x)=0 e Jacobiano analítico
- [ ] `nonlinear_solvers.f90`: Newton + line search
- [ ] `nonlinear_solvers.f90`: Broyden (quasi-Newton)
- [ ] `continuation.f90`: Sweep em U com preditor-corretor
- [ ] Testes unitários extensivos
  - [ ] U=0 (Fermi gas)
  - [ ] U→∞ (forte acoplamento)
  - [ ] Half-filling
  - [ ] Comparação com literatura (Essler, Lieb-Wu)
- [ ] Geração de tabelas de teste (n, m, U)

### Fase 2: Splines 2D (3-4 dias)
- [ ] `spline2d.f90`: Interpolação bicúbica
- [ ] `xc_lsda.f90`: Interface exc, Vxc_up, Vxc_dn
- [ ] Testes de interpolação vs valores exatos
- [ ] Benchmark de performance

### Fase 3: Hamiltoniano Básico (2-3 dias)
- [ ] `potential_uniform.f90`, `potential_harmonic.f90`
- [ ] `hamiltonian_builder.f90`: Tight-binding + Veff
- [ ] `lapack_wrapper.f90`: Interface DSYEVD
- [ ] `boundary_conditions.f90`: Open, periodic
- [ ] Teste end-to-end simples (1 iteração KS)

### Fase 4: Ciclo Auto-Consistente (3-4 dias)
- [ ] `density_calculator.f90`: Ocupação de níveis
- [ ] `convergence_monitor.f90`: Critérios de parada
- [ ] `mixing_schemes.f90`: Linear mixing
- [ ] `ks_cycle.f90`: Loop SCF completo
- [ ] Testes:
  - [ ] U=0, BC periódica → Fermi gas
  - [ ] Half-filling, U>0 → comparar literatura

**🎉 MILESTONE:** Código funcional end-to-end!

### Fase 5: Features Avançadas (1 semana)
- [ ] `degeneracy_handler.f90`: Tratamento de níveis degenerados
- [ ] `symmetry.f90`: Exploração de paridade
- [ ] `twisted_bc.f90`: Boundary conditions torcidas
- [ ] Potenciais avançados (impurity, barrier, random, etc)
- [ ] Testes para cada feature

### Fase 6: Otimização (ongoing)
- [ ] Paralelização OpenMP (Bethe Ansatz + KS loop)
- [ ] Profiling e otimização de hotspots
- [ ] I/O melhorado (HDF5?)
- [ ] Documentação completa (Doxygen/FORD)
- [ ] Benchmarks vs código C++ original

---

## 🧪 Estratégia de Testes

### Pirâmide de Testes

```
        /\
       /E2E\        (End-to-End: ciclo completo, casos físicos)
      /------\
     /Integr.\     (Integração: módulos combinados)
    /----------\
   /Unit Tests \   (Unitários: funções individuais)
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
gcov src/*.f90
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

| Caso                  | C++ (original) | Fortran (meta) |
|-----------------------|----------------|----------------|
| Bethe (N=100)         | ~1.0s          | < 0.5s         |
| Spline interpolation  | ~10μs          | < 5μs          |
| Diagonalização (N=100)| ~50ms (Givens) | < 5ms (LAPACK) |
| Ciclo KS (10 iter)    | ~5s            | < 2s           |

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

### Código de Referência

- Código C++ original (este projeto)
- DMFT solvers (TRIQS, w2dynamics)
- Exact diagonalization codes (ALPS, ITensor)

---

## 🤝 Contribuindo

### Workflow

1. **Branch por feature**
   ```bash
   git checkout -b feature/bethe-ansatz-solver
   ```

2. **Commits atômicos**
   ```bash
   git commit -m "feat: add Newton solver with line search"
   git commit -m "test: validate against Lieb-Wu U=0 limit"
   ```

3. **Pull Request**
   - Descrição clara da mudança
   - Testes passando
   - Cobertura mantida/aumentada

### Convenções de Commit

```
feat:     Nova funcionalidade
fix:      Correção de bug
test:     Adiciona/modifica testes
docs:     Documentação
perf:     Otimização de performance
refactor: Refatoração (sem mudança de comportamento)
style:    Formatação, lint
chore:    Tarefas de manutenção
```

---

## 📝 Notas de Implementação

### TODOs e Decisões Pendentes

- [ ] **Grid de tabelas:** Quantos pontos (n,m,U)? Espaçamento uniforme ou adaptativo?
- [ ] **Formato de output:** ASCII, HDF5, NetCDF?
- [ ] **Paralelização:** OpenMP apenas ou também MPI para grids grandes?
- [ ] **Precisão:** Float64 suficiente ou Float128 em alguns casos?
- [ ] **Estados excitados:** Implementar? (mudando {I_j}, {J_α})

### Perguntas em Aberto

1. Como tratar U < 0 (interação atrativa)? Usar simetria ou resolver separadamente?
2. Implementar TBA (Thermodynamic Bethe Ansatz) para L → ∞?
3. Adicionar temperatura T > 0 (Yang-Yang)?

---

## 📊 Status do Projeto

**Versão:** 0.1.0-dev  
**Status:** 🔄 Em desenvolvimento (Fase 0)  
**Última atualização:** 2025-01-XX

### Checklist de Progresso

#### Infraestrutura
- [ ] Estrutura fpm
- [ ] Tipos básicos
- [ ] Sistema de testes
- [ ] CI/CD

#### Core Physics
- [ ] Bethe Ansatz solver
- [ ] Splines 2D
- [ ] XC functional
- [ ] Hamiltoniano
- [ ] Ciclo KS

#### Features
- [ ] Potenciais
- [ ] Simetria
- [ ] Twisted BC
- [ ] Degenerescências

#### Qualidade
- [ ] Testes unitários (>80% coverage)
- [ ] Testes de integração
- [ ] Testes E2E
- [ ] Documentação completa
- [ ] Benchmarks

---

## 🎓 Para Novos Desenvolvedores

### Quick Start

```bash
# Clonar repositório
git clone https://github.com/seu-usuario/lsda-hubbard-fortran.git
cd lsda-hubbard-fortran

# Build
fpm build

# Rodar testes
fpm test

# Exemplo simples
fpm run --example harmonic_trap
```

### Onde Começar?

1. **Leia:** `README.md` (uso básico) e este `PROJECT_CONTEXT.md` (contexto técnico)
2. **Entenda:** Leia `src/types/lsda_types.f90` para ver estruturas de dados
3. **Explore:** Rode exemplos em `examples/`
4. **Contribua:** Escolha uma issue com label `good-first-issue`

### Recursos de Aprendizado

- **Fortran moderno:** https://fortran-lang.org/learn/
- **Bethe Ansatz:** Essler et al., "The One-Dimensional Hubbard Model"
- **DFT:** Capelle & Campo, "Density functionals and model Hamiltonians"

---

## 📄 Licença

Este projeto é licenciado sob a [MIT License](LICENSE).

---

**Última atualização:** 2025-11-03
**Mantenedores:** Guilherme Canella
**Contato:** guycanella@gmail.com