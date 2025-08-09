module dqfun_cwrap
    use iso_c_binding
    use dqmodule
    implicit none
contains

subroutine dq_add(a,b,c) bind(C,name="dqadd_")
    real(c_long_double), intent(in)  :: a(2), b(2)
    real(c_long_double), intent(out) :: c(2)
    type(dq_real) :: da, db, dc
    da%dqr(1) = real(a(1), dqknd)
    da%dqr(2) = real(a(2), dqknd)
    db%dqr(1) = real(b(1), dqknd)
    db%dqr(2) = real(b(2), dqknd)
    dc = da + db
    c(1) = real(dc%dqr(1), c_long_double)
    c(2) = real(dc%dqr(2), c_long_double)
end subroutine

subroutine dq_sub(a,b,c) bind(C,name="dqsub_")
    real(c_long_double), intent(in)  :: a(2), b(2)
    real(c_long_double), intent(out) :: c(2)
    type(dq_real) :: da, db, dc
    da%dqr(1) = real(a(1), dqknd)
    da%dqr(2) = real(a(2), dqknd)
    db%dqr(1) = real(b(1), dqknd)
    db%dqr(2) = real(b(2), dqknd)
    dc = da - db
    c(1) = real(dc%dqr(1), c_long_double)
    c(2) = real(dc%dqr(2), c_long_double)
end subroutine

subroutine dq_mul(a,b,c) bind(C,name="dqmul_")
    real(c_long_double), intent(in)  :: a(2), b(2)
    real(c_long_double), intent(out) :: c(2)
    type(dq_real) :: da, db, dc
    da%dqr(1) = real(a(1), dqknd)
    da%dqr(2) = real(a(2), dqknd)
    db%dqr(1) = real(b(1), dqknd)
    db%dqr(2) = real(b(2), dqknd)
    dc = da * db
    c(1) = real(dc%dqr(1), c_long_double)
    c(2) = real(dc%dqr(2), c_long_double)
end subroutine

subroutine dq_div(a,b,c) bind(C,name="dqdiv_")
    real(c_long_double), intent(in)  :: a(2), b(2)
    real(c_long_double), intent(out) :: c(2)
    type(dq_real) :: da, db, dc
    da%dqr(1) = real(a(1), dqknd)
    da%dqr(2) = real(a(2), dqknd)
    db%dqr(1) = real(b(1), dqknd)
    db%dqr(2) = real(b(2), dqknd)
    dc = da / db
    c(1) = real(dc%dqr(1), c_long_double)
    c(2) = real(dc%dqr(2), c_long_double)
end subroutine

subroutine dq_fromdbl(d,a) bind(C,name="dqdqd_")
    real(c_double), intent(in)  :: d
    real(c_long_double), intent(out) :: a(2)
    type(dq_real) :: da
    da%dqr(1) = real(d, dqknd)
    da%dqr(2) = 0.0_dqknd
    a(1) = real(da%dqr(1), c_long_double)
    a(2) = real(da%dqr(2), c_long_double)
end subroutine

subroutine dq_sqrt(a,b) bind(C,name="dqsqrt_")
    real(c_long_double), intent(in)  :: a(2)
    real(c_long_double), intent(out) :: b(2)
    type(dq_real) :: da, db
    da%dqr(1) = real(a(1), dqknd)
    da%dqr(2) = real(a(2), dqknd)
    db = sqrt(da)
    b(1) = real(db%dqr(1), c_long_double)
    b(2) = real(db%dqr(2), c_long_double)
end subroutine

subroutine dq_tostr(a,nd,s,str_len) bind(C,name="dqtoqd_")
    real(c_long_double), intent(in) :: a(2)
    integer(c_int), intent(in) :: nd
    character(c_char), intent(out) :: s(str_len)
    integer(c_int), value :: str_len
    type(dq_real) :: da
    integer :: i, c_len, nd_i, w_eff, d_eff
    character(1) :: buf(512)

    da%dqr(1) = real(a(1), dqknd)
    da%dqr(2) = real(a(2), dqknd)
    nd_i = int(nd)
    d_eff = max(1, min(nd_i, 64))
    w_eff = max(d_eff + 10, 30)
    if (w_eff > size(buf)) w_eff = size(buf)

    call dqeform(da, w_eff, d_eff, buf)

    c_len = min(w_eff, str_len-1)
    do i = 1, c_len
        s(i) = buf(i)
    end do
    s(c_len+1) = c_null_char
end subroutine

end module