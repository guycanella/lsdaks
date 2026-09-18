!> Closing sweep of the generated XC tables against the C++ reference set.
!!
!! VALIDATION-ONLY tool: this program exists to license the deletion of
!! `original/` and `data/tables/`.  It regenerates a table with the Fortran
!! generator for every reference `U`, feeds it through the path the SCF really
!! uses (`xc_lsda_init` on the generated file, then `get_exc` / `get_vxc`) and
!! reports the worst deviation at the nodes of the reference table.  Delete it
!! together with the reference tables; it is deliberately NOT part of
!! `fpm test`, because it costs tens of minutes.
!!
!! The `n = 1, m = 1` corner is excluded: `xc_lsda` short-circuits it to zero
!! (`is_corner_shortcut`, a literal transcription of the C++ reference), so it
!! says nothing about the generator.
!!
!! Usage: sweep_reference_tables <scratch-dir> <reference table> [...]
!! Build: scripts/sweep_reference_tables.sh
program sweep_reference_tables
    use lsda_constants, only: dp
    use bethe_tables, only: generate_xc_table, grid_params_t
    use table_io, only: xc_table_t, read_fortran_table, write_fortran_table, &
                        deallocate_table
    use xc_lsda, only: xc_lsda_t, xc_lsda_init, get_exc, get_vxc, xc_lsda_destroy
    use lsda_errors, only: ERROR_SUCCESS
    implicit none

    real(dp), parameter :: TOL_EXC = 1.0e-6_dp, TOL_VXC = 1.0e-4_dp
    real(dp), parameter :: CORNER_TOL = 1.0e-14_dp

    character(len=512) :: scratch, path, gen_path
    type(grid_params_t) :: params
    type(xc_table_t) :: ref, gen
    type(xc_lsda_t) :: xc
    integer :: nargs, k, i, j, ierr, n_eval, n_bad_exc, n_bad_vxc, n_skipped
    integer :: i_worst, j_worst, i_worst_v, j_worst_v
    integer(8) :: c0, c1, rate
    real(dp) :: U, n, m, n_up, n_dw, exc, vup, vdn, gen_wall
    real(dp) :: d_exc, d_up, d_dn, w_exc, w_up, w_dn, n_worst, m_worst
    real(dp) :: n_worst_v, m_worst_v
    logical :: failed_any

    nargs = command_argument_count()
    if (nargs < 2) then
        print '(A)', "Usage: sweep_reference_tables <scratch-dir> <reference table> [...]"
        stop 1
    end if

    call get_command_argument(1, scratch)
    failed_any = .false.

    print '(A)', "| U | worst |d exc| | worst |d Vxc_up| | worst |d Vxc_down| " // &
                 "| at (n, m) | nodes > tol (exc / Vxc) | nodes | gen wall [s] |"
    print '(A)', "|---|---------------|-------------------|---------------------" // &
                 "|-----------|-------------------------|-------|--------------|"

    do k = 2, nargs
        call get_command_argument(k, path)

        call read_fortran_table(trim(path), ref, ierr)
        if (ierr /= ERROR_SUCCESS) then
            print '(A,A)', "ERROR: cannot read ", trim(path)
            stop 1
        end if
        U = ref%U

        call system_clock(c0, rate)
        call generate_xc_table(U, params, gen, ierr)
        call system_clock(c1)
        gen_wall = real(c1 - c0, dp) / real(max(rate, 1_8), dp)
        if (ierr /= ERROR_SUCCESS) then
            print '(A,F0.2)', "ERROR: generation failed at U = ", U
            stop 1
        end if

        write(gen_path, '(A,A,F0.2,A)') trim(scratch), '/gen_u', U, '.dat'
        call write_fortran_table(trim(gen_path), gen, ierr)
        if (ierr /= ERROR_SUCCESS) then
            print '(A,A)', "ERROR: cannot write ", trim(gen_path)
            stop 1
        end if

        call xc_lsda_init(xc, trim(gen_path), ierr)
        if (ierr /= ERROR_SUCCESS) then
            print '(A,A)', "ERROR: xc_lsda_init failed on ", trim(gen_path)
            stop 1
        end if

        w_exc = 0.0_dp
        w_up = 0.0_dp
        w_dn = 0.0_dp
        n_eval = 0
        n_bad_exc = 0
        n_bad_vxc = 0
        n_skipped = 0
        i_worst = 1
        j_worst = 1
        i_worst_v = 1
        j_worst_v = 1

        do i = 1, ref%n_points_n
            n = ref%n_grid(i)
            do j = 1, ref%n_points_m
                m = ref%m_grid(j, i)
                if (m > n) cycle
                ! The (1, 1) corner is a shortcut in xc_lsda, not a solve.
                if (n >= 1.0_dp - CORNER_TOL .and. m >= 1.0_dp - CORNER_TOL) cycle

                n_up = 0.5_dp * (n + m)
                n_dw = 0.5_dp * (n - m)

                call get_exc(xc, n_up, n_dw, exc, ierr)
                if (ierr /= ERROR_SUCCESS) then
                    n_skipped = n_skipped + 1
                    cycle
                end if
                call get_vxc(xc, n_up, n_dw, vup, vdn, ierr)
                if (ierr /= ERROR_SUCCESS) then
                    n_skipped = n_skipped + 1
                    cycle
                end if

                d_exc = abs(exc - ref%exc(j, i))
                d_up = abs(vup - ref%vxc_up(j, i))
                d_dn = abs(vdn - ref%vxc_down(j, i))

                n_eval = n_eval + 1
                if (d_exc > TOL_EXC) n_bad_exc = n_bad_exc + 1
                if (max(d_up, d_dn) > TOL_VXC) n_bad_vxc = n_bad_vxc + 1

                if (d_exc > w_exc) then
                    w_exc = d_exc
                    i_worst = i
                    j_worst = j
                end if
                if (max(d_up, d_dn) > max(w_up, w_dn)) then
                    i_worst_v = i
                    j_worst_v = j
                end if
                w_up = max(w_up, d_up)
                w_dn = max(w_dn, d_dn)
            end do
        end do

        n_worst = ref%n_grid(i_worst)
        m_worst = ref%m_grid(j_worst, i_worst)
        n_worst_v = ref%n_grid(i_worst_v)
        m_worst_v = ref%m_grid(j_worst_v, i_worst_v)

        print '(A,F6.2,A,ES10.3,A,ES10.3,A,ES10.3,A,F6.4,A,ES9.2,A,F6.4,A,ES9.2, &
               &A,I0,A,I0,A,I0,A,F8.2,A)', &
            "| ", U, " | ", w_exc, " | ", w_up, " | ", w_dn, &
            " | exc (", n_worst, ", ", m_worst, ") / Vxc (", n_worst_v, ", ", &
            m_worst_v, ") | ", n_bad_exc, " / ", n_bad_vxc, " | ", n_eval, &
            " | ", gen_wall, " |"
        print '(A,I0)', "  skipped nodes: ", n_skipped
        if (n_bad_exc > 0 .or. n_bad_vxc > 0 .or. n_skipped > 0) failed_any = .true.

        call xc_lsda_destroy(xc)
        call deallocate_table(gen)
        call deallocate_table(ref)
    end do

    if (failed_any) then
        print '(A)', "Validation failed: tolerance violations or skipped nodes were found."
        stop 1
    end if

end program sweep_reference_tables
