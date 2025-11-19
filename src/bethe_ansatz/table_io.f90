!> Input/output routines for XC functional tables
module table_io
    use, intrinsic :: iso_fortran_env, only: real64, int32, error_unit
    use lsda_constants, only: dp
    use lsda_errors, only: ERROR_SUCCESS, ERROR_FILE_NOT_FOUND, &
                                    ERROR_FILE_READ, ERROR_FILE_WRITE, &
                                    ERROR_INVALID_INPUT
    implicit none
    private

    public :: xc_table_t
    public :: read_cpp_table, write_fortran_table, read_fortran_table
    public :: deallocate_table, print_table_info
    public :: extract_U_from_filename

    integer, parameter :: MAX_LINE_LEN = 256

    type :: xc_table_t
        real(dp) :: U
        integer :: n_points_n
        integer :: n_points_m
        real(dp), allocatable :: n_grid(:)
        real(dp), allocatable :: m_grid(:,:)
        real(dp), allocatable :: exc(:,:)
        real(dp), allocatable :: vxc_up(:,:)
        real(dp), allocatable :: vxc_down(:,:)
    end type xc_table_t

contains

    subroutine read_cpp_table(filename, table, ierr)
        character(len=*), intent(in) :: filename
        type(xc_table_t), intent(out) :: table
        integer, intent(out) :: ierr
        integer :: unit, io_stat
        character(len=MAX_LINE_LEN) :: line
        integer :: n_blocks, n_mag_points
        integer :: i_block, i_mag
        real(dp) :: n_val, m_val, exc_val, vxc_up_val, vxc_down_val

        call extract_U_from_filename(filename, table%U, ierr)
        if (ierr /= ERROR_SUCCESS) return

        call count_blocks_and_points(filename, n_blocks, n_mag_points, ierr)
        if (ierr /= ERROR_SUCCESS) return

        table%n_points_n = n_blocks
        table%n_points_m = n_mag_points

        allocate(table%n_grid(n_blocks))
        allocate(table%m_grid(n_mag_points, n_blocks))
        allocate(table%exc(n_mag_points, n_blocks))
        allocate(table%vxc_up(n_mag_points, n_blocks))
        allocate(table%vxc_down(n_mag_points, n_blocks))

        open(newunit=unit, file=filename, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            ierr = ERROR_FILE_NOT_FOUND
            return
        end if

        i_block = 0
        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) then
                ierr = ERROR_FILE_READ
                close(unit)
                exit
            end if

            if (line(1:2) == 'n:') then
                i_block = i_block + 1
                read(line(3:), *, iostat=io_stat) n_val
                if (io_stat /= 0) then
                    ierr = ERROR_FILE_READ
                    close(unit)
                    return
                end if
                table%n_grid(i_block) = n_val

                read(unit, '(A)', iostat=io_stat) line  ! Skip header
                if (io_stat /= 0) then
                    ierr = ERROR_FILE_READ
                    close(unit)
                    return
                end if

                do i_mag = 1, n_mag_points
                    read(unit, '(A)', iostat=io_stat) line
                    
                    if (io_stat /= 0 .or. line(1:2) == 'n:') then
                        if (i_mag > 1) then
                            table%m_grid(i_mag:n_mag_points, i_block) = &
                                table%m_grid(i_mag-1, i_block)
                            table%exc(i_mag:n_mag_points, i_block) = &
                                table%exc(i_mag-1, i_block)
                            table%vxc_up(i_mag:n_mag_points, i_block) = &
                                table%vxc_up(i_mag-1, i_block)
                            table%vxc_down(i_mag:n_mag_points, i_block) = &
                                table%vxc_down(i_mag-1, i_block)
                        end if
                        
                        if (line(1:2) == 'n:') backspace(unit)
                        exit
                    end if
                    
                    read(line, *, iostat=io_stat) m_val, exc_val, vxc_up_val, vxc_down_val
                    if (io_stat /= 0) then
                        ierr = ERROR_FILE_READ
                        close(unit)
                        return
                    end if

                    table%m_grid(i_mag, i_block) = m_val
                    table%exc(i_mag, i_block) = exc_val
                    table%vxc_up(i_mag, i_block) = vxc_up_val
                    table%vxc_down(i_mag, i_block) = vxc_down_val
                end do
            end if
        end do

        close(unit)
        ierr = ERROR_SUCCESS
    end subroutine read_cpp_table

    subroutine count_blocks_and_points(filename, n_blocks, n_mag_points, ierr)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: n_blocks, n_mag_points, ierr
        integer :: unit, io_stat
        character(len=MAX_LINE_LEN) :: line
        integer :: points_in_block, max_points
        logical :: in_data_section

        open(newunit=unit, file=filename, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            ierr = ERROR_FILE_NOT_FOUND
            return
        end if

        n_blocks = 0
        max_points = 0
        points_in_block = 0
        in_data_section = .false.

        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) then
                ierr = ERROR_FILE_READ
                close(unit)
                exit
            end if

            if (line(1:2) == 'n:') then
                if (n_blocks > 0) then
                    max_points = max(max_points, points_in_block)
                end if
                
                n_blocks = n_blocks + 1
                points_in_block = 0
                in_data_section = .false.
                
            else if (line(1:4) == '#mag') then
                in_data_section = .true.
                
            else if (in_data_section .and. len_trim(line) > 0) then
                points_in_block = points_in_block + 1
            end if
        end do

        if (points_in_block > 0) then
            max_points = max(max_points, points_in_block)
        end if
        
        n_mag_points = max_points
        close(unit)

        if (n_blocks == 0 .or. n_mag_points == 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if
        
        ierr = ERROR_SUCCESS
    end subroutine count_blocks_and_points

    subroutine extract_U_from_filename(filename, U, ierr)
        character(len=*), intent(in) :: filename
        real(dp), intent(out) :: U
        integer, intent(out) :: ierr
        integer :: pos_u, pos_end, i, io_stat
        character(len=32) :: u_string
        logical :: found_dot

        pos_u = index(filename, '_u')
        if (pos_u == 0) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        u_string = filename(pos_u+2:)
        
        pos_end = 0
        found_dot = .false.
        
        do i = 1, len_trim(u_string)
            if (u_string(i:i) == '.') then
                if (found_dot) then
                    pos_end = i - 1
                    exit
                else
                    found_dot = .true.
                end if
            else if (.not. ((u_string(i:i) >= '0' .and. u_string(i:i) <= '9') .or. &
                            u_string(i:i) == '-' .or. u_string(i:i) == '+' .or. &
                            u_string(i:i) == 'e' .or. u_string(i:i) == 'E')) then
                pos_end = i - 1
                exit
            end if
        end do
        
        if (pos_end == 0) pos_end = len_trim(u_string)
        
        read(u_string(1:pos_end), *, iostat=io_stat) U
        
        if (io_stat /= 0) then
            ierr = ERROR_FILE_READ
            return
        end if
        
        ierr = ERROR_SUCCESS

    end subroutine extract_U_from_filename
    
    subroutine write_fortran_table(filename, table, ierr)
        character(len=*), intent(in) :: filename
        type(xc_table_t), intent(in) :: table
        integer, intent(out) :: ierr
        integer :: unit, io_stat

        if (.not. allocated(table%n_grid)) then
            ierr = ERROR_INVALID_INPUT
            return
        end if

        open(newunit=unit, file=filename, status='replace', action='write', &
             form='unformatted', access='stream', iostat=io_stat)
        if (io_stat /= 0) then
            ierr = ERROR_FILE_WRITE
            return
        end if

        write(unit, iostat=io_stat) table%U
        if (io_stat /= 0) goto 100
        write(unit, iostat=io_stat) table%n_points_n, table%n_points_m
        if (io_stat /= 0) goto 100
        write(unit, iostat=io_stat) table%n_grid
        if (io_stat /= 0) goto 100
        write(unit, iostat=io_stat) table%m_grid
        if (io_stat /= 0) goto 100
        write(unit, iostat=io_stat) table%exc
        if (io_stat /= 0) goto 100
        write(unit, iostat=io_stat) table%vxc_up
        if (io_stat /= 0) goto 100
        write(unit, iostat=io_stat) table%vxc_down
        if (io_stat /= 0) goto 100

        close(unit)
        ierr = ERROR_SUCCESS
        return

100     ierr = ERROR_FILE_WRITE
        close(unit)
    end subroutine write_fortran_table


    subroutine read_fortran_table(filename, table, ierr)
        character(len=*), intent(in) :: filename
        type(xc_table_t), intent(out) :: table
        integer, intent(out) :: ierr
        integer :: unit, io_stat

        call deallocate_table(table)

        open(newunit=unit, file=filename, status='old', action='read', &
             form='unformatted', access='stream', iostat=io_stat)
        if (io_stat /= 0) then
            ierr = ERROR_FILE_NOT_FOUND
            return
        end if

        read(unit, iostat=io_stat) table%U
        if (io_stat /= 0) goto 200
        read(unit, iostat=io_stat) table%n_points_n, table%n_points_m
        if (io_stat /= 0) goto 200

        if (table%n_points_n <= 0 .or. table%n_points_m <= 0) then
            ierr = ERROR_INVALID_INPUT
            close(unit)
            return
        end if

        allocate(table%n_grid(table%n_points_n))
        allocate(table%m_grid(table%n_points_m, table%n_points_n))
        allocate(table%exc(table%n_points_m, table%n_points_n))
        allocate(table%vxc_up(table%n_points_m, table%n_points_n))
        allocate(table%vxc_down(table%n_points_m, table%n_points_n))

        read(unit, iostat=io_stat) table%n_grid
        if (io_stat /= 0) goto 200
        read(unit, iostat=io_stat) table%m_grid
        if (io_stat /= 0) goto 200
        read(unit, iostat=io_stat) table%exc
        if (io_stat /= 0) goto 200
        read(unit, iostat=io_stat) table%vxc_up
        if (io_stat /= 0) goto 200
        read(unit, iostat=io_stat) table%vxc_down
        if (io_stat /= 0) goto 200

        close(unit)
        ierr = ERROR_SUCCESS
        return

200     ierr = ERROR_FILE_READ
        close(unit)
        call deallocate_table(table)
    end subroutine read_fortran_table


    subroutine deallocate_table(table)
        type(xc_table_t), intent(inout) :: table
        if (allocated(table%n_grid)) deallocate(table%n_grid)
        if (allocated(table%m_grid)) deallocate(table%m_grid)
        if (allocated(table%exc)) deallocate(table%exc)
        if (allocated(table%vxc_up)) deallocate(table%vxc_up)
        if (allocated(table%vxc_down)) deallocate(table%vxc_down)
    end subroutine deallocate_table


    subroutine print_table_info(table, unit)
        type(xc_table_t), intent(in) :: table
        integer, intent(in), optional :: unit
        integer :: out_unit
        real(dp) :: n_min, n_max, m_min, m_max, exc_min, exc_max, vxc_min, vxc_max

        out_unit = 6
        if (present(unit)) out_unit = unit

        if (.not. allocated(table%n_grid)) then
            write(out_unit, '(A)') "XC Table: [UNALLOCATED]"
            return
        end if

        n_min = minval(table%n_grid)
        n_max = maxval(table%n_grid)
        m_min = minval(table%m_grid)
        m_max = maxval(table%m_grid)
        exc_min = minval(table%exc)
        exc_max = maxval(table%exc)
        vxc_min = min(minval(table%vxc_up), minval(table%vxc_down))
        vxc_max = max(maxval(table%vxc_up), maxval(table%vxc_down))

        write(out_unit, '(A)') "========================================="
        write(out_unit, '(A)') "XC Table Summary"
        write(out_unit, '(A)') "========================================="
        write(out_unit, '(A,F8.4)') "  Hubbard U:           ", table%U
        write(out_unit, '(A,I6)')   "  Density points:      ", table%n_points_n
        write(out_unit, '(A,I6)')   "  Magnetization points:", table%n_points_m
        write(out_unit, '(A)')      "-----------------------------------------"
        write(out_unit, '(A,F10.6,A,F10.6)') "  Density range:       ", n_min, " to ", n_max
        write(out_unit, '(A,ES12.5,A,ES12.5)') "  Magnetization range: ", m_min, " to ", m_max
        write(out_unit, '(A)')      "-----------------------------------------"
        write(out_unit, '(A,ES12.5,A,ES12.5)') "  E_xc range:          ", exc_min, " to ", exc_max
        write(out_unit, '(A,ES12.5,A,ES12.5)') "  V_xc range:          ", vxc_min, " to ", vxc_max
        write(out_unit, '(A)') "========================================="
    end subroutine print_table_info

end module table_io