module ddfun_cwrap
    use iso_c_binding
    use ddmodule
    use ddfuna, only: ddinpc, ddcpr
    implicit none
contains

subroutine dd_add(a,b,c) bind(C,name="ddadd_")
    real(c_double), intent(in)  :: a(2), b(2)
    real(c_double), intent(out) :: c(2)
    type(dd_real) :: da, db, dc
    da%ddr = a
    db%ddr = b
    dc = da + db
    c = dc%ddr
end subroutine

subroutine dd_sub(a,b,c) bind(C,name="ddsub_")
    real(c_double), intent(in)  :: a(2), b(2)
    real(c_double), intent(out) :: c(2)
    type(dd_real) :: da, db, dc
    da%ddr = a
    db%ddr = b
    dc = da - db
    c = dc%ddr
end subroutine

subroutine dd_mul(a,b,c) bind(C,name="ddmul_")
    real(c_double), intent(in)  :: a(2), b(2)
    real(c_double), intent(out) :: c(2)
    type(dd_real) :: da, db, dc
    da%ddr = a
    db%ddr = b
    dc = da * db
    c = dc%ddr
end subroutine

subroutine dd_div(a,b,c) bind(C,name="dddiv_")
    real(c_double), intent(in)  :: a(2), b(2)
    real(c_double), intent(out) :: c(2)
    type(dd_real) :: da, db, dc
    da%ddr = a
    db%ddr = b
    dc = da / db
    c = dc%ddr
end subroutine

subroutine dd_fromdbl(d,a) bind(C,name="dddqd_")
    real(c_double), intent(in)  :: d
    real(c_double), intent(out) :: a(2)
    type(dd_real) :: da
    integer :: nd_i
    da%ddr(1) = d
    da%ddr(2) = 0.0d0
    a = da%ddr
end subroutine

subroutine dd_sqrt(a,b) bind(C,name="ddsqrt_")
    real(c_double), intent(in)  :: a(2)
    real(c_double), intent(out) :: b(2)
    type(dd_real) :: da, db
    da%ddr = a
    db = sqrt(da)
    b = db%ddr
end subroutine

subroutine dd_abs(a,b) bind(C,name="ddabs_")
    real(c_double), intent(in)  :: a(2)
    real(c_double), intent(out) :: b(2)
    type(dd_real) :: da, db
    da%ddr = a
    db = abs(da)
    b = db%ddr
end subroutine

subroutine dd_compare(a,b,ic) bind(C,name="ddcpr_")
    real(c_double), intent(in)  :: a(2), b(2)
    integer(c_int), intent(out) :: ic
    type(dd_real) :: da, db
    integer :: ic_local

    da%ddr = a
    db%ddr = b
    call ddcpr(da%ddr, db%ddr, ic_local)
    ic = ic_local
end subroutine

subroutine dd_tostr(a,nd,s,str_len) bind(C,name="dd_to_string")
    real(c_double), intent(in) :: a(2)
    integer(c_int), intent(in) :: nd
    character(c_char), intent(out) :: s(str_len)
    integer(c_int), value :: str_len
    type(dd_real) :: da
    integer :: i, c_len, nd_i, w_eff, d_eff
    character(1) :: buf(256)

    da%ddr = a
    nd_i = int(nd)
    d_eff = max(1, min(nd_i, 33))
    w_eff = max(d_eff + 10, 20)
    if (w_eff > size(buf)) w_eff = size(buf)

    call ddeform(da, w_eff, d_eff, buf)

    c_len = min(w_eff, str_len-1)
    do i = 1, c_len
        s(i) = buf(i)
    end do
    s(c_len+1) = c_null_char
end subroutine

subroutine dd_fromstr(str, strlen, a) bind(C,name="ddfromstr_")
    character(c_char), intent(in) :: str(*)
    integer(c_int), value :: strlen
    real(c_double), intent(out) :: a(2)

    type(dd_real) :: tmp
    character(len=120) :: buffer

    call copy_c_string(str, strlen, buffer)
    call ddinpc(buffer, tmp%ddr)
    a = tmp%ddr
end subroutine

subroutine dd_read_line(str, strlen, value) bind(C,name="dd_read_line")
    character(c_char), intent(in) :: str(*)
    integer(c_int), value :: strlen
    real(c_double), intent(out) :: value(2)

    type(dd_real) :: tmp
    character(len=120) :: buffer
    integer :: unit, ios

    call copy_c_string(str, strlen, buffer)
    open(newunit=unit, status='scratch', action='readwrite', form='formatted', iostat=ios)
    if (ios /= 0) then
        value = 0.0d0
        return
    end if

    write(unit, '(a)') trim(buffer)
    rewind(unit)
    call ddread(unit, tmp)
    close(unit)

    value = tmp%ddr
end subroutine

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
