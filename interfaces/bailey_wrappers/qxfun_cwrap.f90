module qxfun_cwrap
    use iso_c_binding
    use qxmodule
    implicit none
contains

subroutine qx_add(a,b,c) bind(C,name="qxadd_")
    real(qxknd), intent(in)  :: a, b
    real(qxknd), intent(out) :: c
    c = a + b
end subroutine

subroutine qx_sub(a,b,c) bind(C,name="qxsub_")
    real(qxknd), intent(in)  :: a, b
    real(qxknd), intent(out) :: c
    c = a - b
end subroutine

subroutine qx_mul(a,b,c) bind(C,name="qxmul_")
    real(qxknd), intent(in)  :: a, b
    real(qxknd), intent(out) :: c
    c = a * b
end subroutine

subroutine qx_div(a,b,c) bind(C,name="qxdiv_")
    real(qxknd), intent(in)  :: a, b
    real(qxknd), intent(out) :: c
    c = a / b
end subroutine

subroutine qx_fromdbl(d,a) bind(C,name="qxdqd_")
    real(c_double), intent(in)  :: d
    real(qxknd), intent(out) :: a
    a = real(d, qxknd)
end subroutine

subroutine qx_sqrt(a,b) bind(C,name="qxsqrt_")
    real(qxknd), intent(in)  :: a
    real(qxknd), intent(out) :: b
    b = sqrt(a)
end subroutine

subroutine qx_abs(a,b) bind(C,name="qxabs_")
    real(qxknd), intent(in)  :: a
    real(qxknd), intent(out) :: b
    b = abs(a)
end subroutine

subroutine qx_tostr(a, nd, s, str_len) bind(C,name="qx_to_string")
    real(qxknd), intent(in) :: a
    integer(c_int), intent(in) :: nd
    character(c_char), intent(out) :: s(str_len)
    integer(c_int), value :: str_len
    character(len=:), allocatable :: tmp
    character(len=64) :: fmt
    integer :: i, c_len, nd_i, w_eff, d_eff

    nd_i = int(nd)
    d_eff = max(1, min(nd_i, 33))
    w_eff = max(d_eff + 10, 20)
    allocate(character(len=w_eff) :: tmp)
    write(fmt, '(A,I0,A,I0,A)') '(ES', w_eff, '.', d_eff, 'E3)'
    write(tmp, fmt) a
    c_len = min(len_trim(tmp), str_len - 1)
    do i = 1, c_len
        s(i) = tmp(i:i)
    end do
    s(c_len + 1) = c_null_char
end subroutine qx_tostr

subroutine qx_fromstr(str, strlen, a) bind(C,name="qxfromstr_")
    character(c_char), intent(in) :: str(*)
    integer(c_int), value :: strlen
    real(qxknd), intent(out) :: a
    character(len=256) :: buffer

    call copy_c_string(str, strlen, buffer)
    read(buffer, *, err=100, end=100) a
    return
100 continue
    a = 0.0_qxknd
end subroutine qx_fromstr

subroutine qx_read_line(str, strlen, value) bind(C,name="qx_read_line")
    character(c_char), intent(in) :: str(*)
    integer(c_int), value :: strlen
    real(c_long_double), intent(out) :: value

    real(qxknd) :: tmp
    character(len=256) :: buffer
    integer :: ios

    call copy_c_string(str, strlen, buffer)
    read(buffer, *, iostat=ios) tmp
    if (ios /= 0) then
        value = 0.0_c_long_double
    else
        value = real(tmp, c_long_double)
    end if
end subroutine qx_read_line

subroutine copy_c_string(str, strlen, buffer)
    character(c_char), intent(in) :: str(*)
    integer(c_int), value :: strlen
    character(*), intent(out) :: buffer
    integer :: i, copy_len

    buffer = ' '
    copy_len = min(strlen, len(buffer))
    do i = 1, copy_len
        if (str(i) == c_null_char) exit
        buffer(i:i) = str(i)
    end do
end subroutine copy_c_string

end module
