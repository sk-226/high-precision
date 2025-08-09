module ddfun_cwrap
    use iso_c_binding
    use ddmodule
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

subroutine dd_tostr(a,nd,s,str_len) bind(C,name="ddtoqd_")
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

end module