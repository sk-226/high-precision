
!  QXFUN: A quad precision package with special functions

!  High-level language interface module (QXMODULE).

!  Revision date:  18 Mar 2023

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired)
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package enhances an IEEE quad precision floating-point facility (approx.
!    33 digit accuracy) to include a library of numerous special functions, all by
!    making only very minor changes to existing Fortran programs.  The package should
!    run correctly on any Unix-based system supporting a Fortran-2008 compiler and
!    IEEE 128-bit floating-point arithmetic, in hardware or software (for example, the
!    GNU gfortran compiler and the Intel ifort compiler). Note however that results
!    are NOT guaranteed to the last bit.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.

!    Three related software packages by the same author are DDFUN (double-double real
!    and complex), DQFUN (double-quad real and complex) and MPFUN2020 (arbitrary
!    precision real and complex). They are available from the same site as this suite:
!    http://www.davidhbailey.com/dhbsoftware/

!  DOCUMENTATION:
!    See the README-qxfun.txt file in the main QXFUN directory.

!  DESCRIPTION OF THIS MODULE (QXFUNE):
!    This module contains all special functions.

module qxfune

integer, public, parameter:: qxknd = selected_real_kind (33, 4931), qxldb = 6, &
  qxnbt = 113, qxnwx = 1
real (qxknd), public, parameter:: qxdpw = 34.016389510029875059152495103867712_qxknd, &
  qxlogb = 38.816242111356937327364998801657888_qxknd, qxrdfz = 2.e0_qxknd**(-110)
real (qxknd), public, parameter:: &
  qxegammacon = 5.7721566490153286060651209008240247066034e-1_qxknd, &
  qxlog2con = 0.6931471805599453094172321214581765680755e0_qxknd, &
  qxpicon = 3.1415926535897932384626433832795027974791e0_qxknd

contains

!  These are basic utility routines.

subroutine qxabrt
implicit none
stop
end subroutine

subroutine qxegamc (ra)
implicit none
real (qxknd), intent(out)::ra
ra = qxegammacon
return
end

subroutine qxlog2c (ra)
implicit none
real (qxknd), intent(out)::ra
ra = qxlog2con
return
end

subroutine qxpic (ra)
implicit none
real (qxknd), intent(out)::ra
ra = qxpicon
return
end

integer function qxsgn (ra)
implicit none
real (qxknd), intent(in):: ra
integer ia
if (ra == 0.e0_qxknd) then
  qxsgn = 0
elseif (ra > 0.e0_qxknd) then
  qxsgn = 1
else
  qxsgn = -1
endif
return
end function qxsgn

!   Special functions start here.

subroutine qxagmr (a, b, c)

!   This performs the arithmetic-geometric mean (AGM) iterations on A and B.
!   The AGM algorithm is as follows: Set a_0 = a and b_0 = b, then iterate

!    a_{k+1} = (a_k + b_k)/2
!    b_{k+1} = sqrt (a_k * b_k)

!   until convergence (i.e., until a_k = b_k to available precision).
!   The result is returned in C.

implicit none
real (qxknd), intent(in):: a, b
real (qxknd), intent(out):: c
integer, parameter:: itrmax = 100, ibz = 2
integer j
real (qxknd) eps, s0, s1, s2, s3, tc1, tc2

if (qxsgn (a) <= 0 .or. qxsgn (b) <= 0) then
  write (qxldb, 1)
1 format ('*** QXAGMR: argument <= 0')
  call qxabrt
endif

eps = 2.e0_qxknd ** (ibz - qxnwx*qxnbt)
s1 = a
s2 = b

do j = 1, itrmax
  s0 = s1 + s2
  s3 = s0 * 0.5e0_qxknd
  s0 = s1 * s2
  s2 = sqrt (s0)
  s1 = s3

!   Check for convergence.

  s0 = s1 - s2
  tc1 = abs (s0) - eps * abs (s1)
  if (qxsgn (tc1) < 0) goto 100
enddo

write (qxldb, 2)
2 format ('*** QXAGMR: Iteration limit exceeded.')
call qxabrt

100 continue

c = s1

return
end subroutine qxagmr

subroutine qxberner (nb2, berne)

!   This returns the array berne, containing Bernoulli numbers indexed 2*k for
!   k = 1 to n, using a polynomial Newton iteration scheme as described in the
!   main documentation paper.

implicit none
integer, intent(in):: nb2
real (qxknd), intent(out):: berne(nb2)
integer, parameter:: ibz = 2, itrmax = 10000
integer i, i1, ic1, j, kn, n, n1, nn1
real (qxknd) d1, dd1, dd2, dd3
real (qxknd) c1(0:nb2), cp2, p1(0:nb2), p2(0:nb2), q(0:nb2), q1, &
  r(0:nb2), s(0:nb2), t1, t2, t3, t4, eps, tc1

!  End of declaration

n = nb2

cp2 = qxpicon**2
c1(0) = 1.e0_qxknd
p1(0) = 1.e0_qxknd
p2(0) = 1.e0_qxknd
q(0) = 1.e0_qxknd

!   Construct numerator and denominator polynomials.

do i = 1, n
  c1(i) = 0.e0_qxknd
  dd1 = 2.e0_qxknd * (i + 1) - 3.e0_qxknd
  dd2 = dd1 + 1.e0_qxknd
  dd3 = dd2 + 1.e0_qxknd
  t1 = cp2 * p1(i-1)
  p1(i) = t1 / (dd1 * dd2)
  t1 = cp2 * p2(i-1)
  p2(i) = t1 / (dd2 * dd3)
  q(i) = 0.e0_qxknd
enddo

kn = min (4, nb2)
eps = 2.e0_qxknd ** (ibz - qxnwx*qxnbt)
q1 = 0.e0_qxknd

!   Perform Newton iterations with dynamic precision levels, using an
!   iteration formula similar to that used to evaluate reciprocals.

do j = 1, itrmax
  call qxpolymul (kn, p2, q, r)
  call qxpolysub (kn, c1, r, s)
  call qxpolymul (kn, s, q, r)
  call qxpolyadd (kn, q, r, q)
  t1 = q(kn) - q1
  t2 = abs (t1)
  tc1 = t2 - eps
  if (qxsgn (tc1) < 0) then
    if (kn == n) goto 100
    if (kn < n) then
      kn = min (2 * kn, n)
      q1 = 0.e0_qxknd
    endif
  else
    q1 = q(kn)
  endif
enddo

write (qxldb, 8)
8 format ('*** QXBERNER: Loop end error')
call qxabrt

100 continue

!   Multiply numerator polynomial by reciprocal of denominator polynomial.

call qxpolymul (n, p1, q, r)

!   Apply formula to produce Bernoulli numbers.

t1 = -2.e0_qxknd
t2 = 1.e0_qxknd

do i = 1, n
  d1 = - real (2*i-1, qxknd) * real (2*i, qxknd)
  t1 = t1 * d1
  t3 = 4.e0_qxknd * cp2
  t2 = t3 * t2
  t4 = 0.5e0_qxknd * t1 / t2
  berne(i) = t4 * abs (r(i))
enddo

return
end subroutine qxberner

subroutine qxpolyadd (n, a, b, c)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: n
real (qxknd), intent(in):: a(0:n), b(0:n)
real (qxknd), intent(out):: c(0:n)
integer k
real (qxknd) t1, t2

!  End of declaration

do k = 0, n
  t1 = a(k)
  t2 = b(k)
  c(k) = t1 + t2
enddo

return
end subroutine qxpolyadd

subroutine qxpolysub (n, a, b, c)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: n
real (qxknd), intent(in):: a(0:n), b(0:n)
real (qxknd), intent(out):: c(0:n)
integer k
real (qxknd) t1, t2

!  End of declaration

do k = 0, n
  t1 = a(k)
  t2 = b(k)
  c(k) = t1 - t2
enddo

return
end subroutine qxpolysub

subroutine qxpolymul (n, a, b, c)

!   This adds two polynomials (ignoring high-order terms), as is required
!   by mpberne. The output array C may NOT be the same as A or B.

implicit none
integer, intent(in):: n
real (qxknd), intent(in):: a(0:n), b(0:n)
real (qxknd), intent(out):: c(0:n)
integer j, k
real (qxknd) t0, t1, t2, t3

!  End of declaration

do k = 0, n
  t0 = 0.e0_qxknd

  do j = 0, k
    t1 = a(j)
    t2 = b(k-j)
    t3 = t1 * t2
    t0 = t0 + t3
  enddo

  c(k) = t0
enddo

return
end subroutine qxpolymul

subroutine qxbesselinr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselI (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.25.2 for modest RR,
!   and DLMF 10.40.1 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (qxknd), intent(in):: rr
real (qxknd), intent(out):: ss
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dfrac = 1.5e0_qxknd
integer ic1, k, nua, n1
real (qxknd) d1
real (qxknd) f1, f2, sum, td, tn, t1, t2, t3, t4, rra, tc1, tc2, tc3, eps

!  End of declaration

!   Check for RR = 0.

if (qxsgn (rr) == 0) then
  write (qxldb, 2)
2 format ('*** QXBESSELINR: Second argument is zero')
  call qxabrt
endif

nua = abs (nu)
eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
rra = abs (rr)
d1 = real (rra, qxknd)

if (d1 < dfrac * qxdpw) then
  tn = 1.e0_qxknd
  f1 = 1.e0_qxknd
  f2 = 1.e0_qxknd
  t1 = 0.25e0_qxknd * rra**2

  do k = 1, nua
    f2 = real (k, qxknd) * f2
  enddo

  td = f1 * f2
  t2 = tn / td
  sum = t2

  do k = 1, itrmax
    f1 = f1 * real (k, qxknd)
    f2 = f2 * real (k + nua, qxknd)
    tn = t1 * tn
    td = f1 * f2
    t2 = tn / td
    sum = sum + t2
    tc1 = abs (t2) - eps * abs (sum)
    if (qxsgn (tc1) < 0) goto 100
  enddo

  write (qxldb, 4)
  4 format ('*** QXBESSELINR: Loop end error 1')
  call qxabrt

100 continue

  t3 = sum * (0.5e0_qxknd * rra)**nua
else
  sum = 1.e0_qxknd
  t1 = real (2 * nua, qxknd)**2
  tn = 1.e0_qxknd
  td = 1.e0_qxknd

  do k = 1, itrmax
    t2 = t1 - real (2*k - 1, qxknd)**2
    tn = -tn * t2
    td = td * real (8*k, qxknd) * rra
    t4 = tn / td
    sum = sum + t4
    tc2 = abs (t4) - eps * abs (sum)
    if (qxsgn (tc2) < 0) goto 110
  enddo

write (qxldb, 5)
5 format ('*** QXBESSELINR: Loop end error 2')
call qxabrt

110 continue

  t1 = exp (rra)
  t2 = 2.e0_qxknd * qxpicon * rra
  t4 = sqrt (t2)
  t2 = t1 / t4
  t3 = t2 * sum
endif

if (qxsgn (rr) < 0 .and. mod (nu, 2) /= 0) t3 = - t3
ss = t3

return
end subroutine qxbesselinr

subroutine qxbesselir (qq, rr, ss)

!   This evaluates the modified Bessel function BesselI (QQ,RR) for QQ and RR
!   both real. The algorithm is DLMF formula 10.25.2 for modest RR, and
!   DLMF 10.40.1 for large RR, relative to precision.

implicit none
real (qxknd), intent(in):: qq, rr
real (qxknd), intent(out):: ss
integer ic1, i0, i1, i2, k, n1
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dfrac = 1.5e0_qxknd
real (qxknd) d1
real (qxknd) f1, f2, sum, td, tn, t1, t2, t3, t4, rra, tc1, tc2, tc3, eps

!  End of declaration

tc1 = 2.e0_qxknd
eps = 2.e0_qxknd ** (-qxnwx*qxnbt)

!   If QQ is integer, call mpbesselinr; if qq < 0 and rr <= 0, then error.

i0 = qxsgn (qq - anint (qq))
i1 = qxsgn (qq)
i2 = qxsgn (rr)
if (i0 == 0) then
  n1 = nint (qq)
  call qxbesselinr (n1, rr, t3)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (qxldb, 3)
3 format ('*** QXBESSELIR: First argument < 0 and second argument <= 0')
  call qxabrt
endif

rra = abs (rr)
d1 = real (rra, qxknd)

if (d1 < dfrac * qxdpw) then
  tn = 1.e0_qxknd
  f1 = 1.e0_qxknd
  t1 = qq + f1
  f2 = gamma (t1)
  t1 = 0.25e0_qxknd * rra**2
  td = f1 * f2
  sum = tn / td

  do k = 1, itrmax
    f1 = f1 * real (k, qxknd)
    t3 = real (k, qxknd)
    t4 = qq + t3
    f2 = f2 * t4
    tn = t1 * tn
    td = f1 * f2
    t2 = tn / td
    sum = sum + t2
    tc1 = abs (t2) - eps * abs (sum)
    if (qxsgn (tc1) < 0) goto 100
  enddo

  write (qxldb, 4)
  4 format ('*** QXBESSELIR: Loop end error 1')
  call qxabrt

100 continue

  t1 = 0.5e0_qxknd * rr
  t2 = t1 ** qq
  t3 = sum * t2
else
  sum = 1.e0_qxknd
  t1 = 4.e0_qxknd * qq**2
  tn = 1.e0_qxknd
  td = 1.e0_qxknd

  do k = 1, itrmax
    t2 = real (2*k - 1, qxknd)
    t2 = t1 - t2**2
    t3 = tn * t2
    t2 = real (8*k, qxknd) * rra
    td = td * t2
    t4 = tn / td
    sum = sum + t4
    tc2 = abs (t4) - eps * abs (sum)
    if (qxsgn (tc2) < 0) goto 110
 enddo

write (qxldb, 5)
5 format ('*** QXBESSELIR: Loop end error 2')
call qxabrt

110 continue

  t1 = exp (rra)
  t2 = 2.e0_qxknd * qxpicon * rra
  t4 = sqrt (t2)
  t3 = t2 * t1 / t4
endif

120 continue

ss = t3

return
end subroutine qxbesselir

subroutine qxbesseljnr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselJ (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (qxknd), intent(in):: rr
real (qxknd), intent(out):: ss
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dfrac = 1.5e0_qxknd
integer ic1, ic2, k, nua, n1
real (qxknd) d1, d2
real (qxknd) f1, f2, sum1, sum2, td1, td2, tn1, tn2, t1, t2, t3, &
  t41, t42, t5, rra, rr2, tc1, tc2, tc3, eps

!  End of declaration

!   Check for RR = 0.

if (qxsgn (rr) == 0) then
  write (qxldb, 2)
2 format ('*** QXBESSELJNR: Second argument is zero')
  call qxabrt
endif

nua = abs (nu)
eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
d1 = real (rr, qxknd)

if (d1 < dfrac * qxdpw) then
  rra = abs (rr)
  tn1 = 1.e0_qxknd
  f1 = 1.e0_qxknd
  f2 = 1.e0_qxknd
  t1 = 0.25e0_qxknd * rra**2

  do k = 1, nua
    f2 = real (k, qxknd) * f2
  enddo

  td1 = f1 * f2
  sum1 = tn1 / td1

  do k = 1, itrmax
    f1 = f1 * real (k, qxknd)
    f2 = f2 * real (k + nua, qxknd)
    tn1 = - t1 * tn1
    td1 = f1 * f2
    t2 = tn1 / td1
    sum1 = sum1 + t2
    tc1 = abs (t2) - eps * abs (sum1)
    if (qxsgn (tc1) < 0) goto 100
  enddo

  write (qxldb, 4)
4 format ('*** QXBESSELJNR: Loop end error 1')
  call qxabrt

100 continue

  t1 = 0.5e0_qxknd * rra
  t3 = sum1 * t1 ** nua
else
  rra = abs (rr)
  rr2 = rra**2
  t1 = real (2 * nua, qxknd)**2
  tn1 = 1.e0_qxknd
  tn2 = (t1 - tn1 ) / 8.e0_qxknd
  td1 = 1.e0_qxknd
  td2 = rra
  sum1 = tn1 / td1
  sum2 = tn2 / td2

  do k = 1, itrmax
    d1 = real (4 * k - 3, qxknd)**2
    d2 = real (4 * k - 1, qxknd)**2
    t2 = d1
    t3 = t1 - t2
    t2 = d2
    t5 = t1 - t2
    t2 = t3 * t5
    tn1 = - tn1 * t2
    d1 = real (2*k - 1, qxknd) * real (2*k, qxknd) * 64.e0_qxknd
    t2 = td1 * d1
    td1 = t2 * rr2
    t41 = tn1 / td1
    sum1 = sum1 + t41

    d1 = (2 * (2*k) - 1)**2
    d2 = (2 * (2*k + 1) - 1)**2
    t2 = d1
    t3 = t1 - t2
    t2 = d2
    t5 = t1 - t2
    t2 = t3 * t5
    tn2 = - tn2 * t2
    d1 = real (2*k, qxknd) * real (2*k + 1, qxknd) * 64.e0_qxknd
    td2 = (td2 * d1) * rr2
    t42 = tn2 / td2
    sum2 = sum2 + t42
    tc1 = abs (t41) - eps * abs (sum1)
    tc2 = abs (t42) - eps * abs (sum2)
    if (qxsgn (tc1) < 0 .and. qxsgn (tc2) < 0) goto 110
  enddo

  write (qxldb, 5)
5 format ('*** QXBESSELJNR: Loop end error 2')
  call qxabrt

110 continue

  t1 = qxpicon * 0.5e0_qxknd * real (nua, qxknd)
  t2 = rra - t1
  t1 = qxpicon * 0.25e0_qxknd
  t3 = t2 - t1
  t41 = cos (t3)
  t42 = sin (t3)
  t1 = t41 * sum1
  t2 = t42 * sum2
  t5 = t1 - t2
  t1 = qxpicon * rra
  t2 = 2.e0_qxknd
  t3 = t2 / t1
  t1 = sqrt (t3)
  t3 = t1 * t5
endif

if (mod (nu, 2) /= 0) then
  if (nu < 0 .and. qxsgn (rr) > 0 .or. nu > 0 .and. qxsgn (rr) < 0) t3 = - t3
endif

ss = t3

return
end subroutine qxbesseljnr

subroutine qxbesseljr (qq, rr, ss)

!   This evaluates the modified Bessel function BesselJ (QQ,RR) for QQ and RR
!   both MPR. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
real (qxknd), intent(in):: qq, rr
real (qxknd), intent(out):: ss
integer ic1, ic2, i1, i2, k, n1
real (qxknd) d1, d2
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dfrac = 1.5e0_qxknd
real (qxknd) f1, f2, sum1, sum2, td1, td2, tn1, tn2, t1, t2, t3, &
  t4, t41, t42, t5, rra, rr2, tc1, tc2, tc3, eps

!  End of declaration

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)

!   If QQ is integer, call mpbesseljnr; if RR <= 0, then error.

t2 = qq - anint (qq)
if (qxsgn (t2) == 0) then
  n1 = nint (qq)
  call qxbesseljnr (n1, rr, t3)
  goto 120
elseif (qxsgn (rr) <= 0) then
  write (qxldb, 3)
3 format ('***QXBESSELJR: Second argument <= 0')
  call qxabrt
endif

d1 = real (rr, qxknd)

if (d1 < dfrac * qxdpw) then
  rra = abs (rr)
  tn1 = 1.e0_qxknd
  f1 = 1.e0_qxknd
  t1 = qq + f1
  f2 = gamma (t1)
  t2 = rra**2
  t1 = t2 * 0.25e0_qxknd
  td1 = f1 * f2
  t2 = tn1 / td1
  sum1 = t2

  do k = 1, itrmax
    f1 = f1 * real (k, qxknd)
    t3 = real (k, qxknd)
    t4 = qq + t3
    f2 = f2 * t4
    tn1 = - t1 * tn1
    td1 = f1 * f2
    t2 = tn1 / td1
    sum1 = sum1 + t2
    tc1 = abs (t2) - eps * abs (sum1)
    if (qxsgn (tc1) < 0) goto 100
  enddo

  write (qxldb, 4)
4 format ('*** QXBESSELJR: Loop end error 1')
  call qxabrt

100 continue

  t1 = rr * 0.5e0_qxknd
  t2 = t1 ** qq
  t3 = sum1 * t2
else
  rra = abs (rr)
  rr2 = rra**2
  t2 = qq**2
  t1 = rr2 * 4.e0_qxknd
  tn1 = 1.e0_qxknd
  t2 = t1 - tn1
  tn2 = t2 / 8.e0_qxknd
  td1 = 1.e0_qxknd
  td2 = rra
  sum1 = tn1 / td1
  sum2 = tn2 / td2

  do k = 1, itrmax
    d1 = real (4 * k - 3, qxknd)**2
    d2 = real (4 * k - 1, qxknd)**2
    t2 = d1
    t3 = t1 - t2
    t2 = d2
    t5 = t1 - t2
    t2 = t3 * t5
    tn1 = - tn1 * t2
    d1 = real (2*k - 1, qxknd) * real (2*k, qxknd) * 64.e0_qxknd
    t2 = td1 * d1
    td1 = t2 * rr2
    t41 = tn1 / td1
    sum1 = sum1 + t41

    d1 = (2 * (2*k) - 1)**2
    d2 = (2 * (2*k + 1) - 1)**2
    t2 = d1
    t3 = t1 - t2
    t2 = d2
    t5 = t1 - t2
    t2 = t3 * t5
    tn2 = - tn2 * t2
    d1 = real (2*k, qxknd) * real (2*k + 1, qxknd) * 64.e0_qxknd
    t2 = td2 * d1
    td2 = t2 * rr2
    t42 = tn2 / td2
    sum2 = sum2 + t42
    tc1 = abs (t41) - eps * abs (sum1)
    tc2 = abs (t42) - eps * abs (sum2)
    if (qxsgn (tc1) < 0 .and. qxsgn (tc2) < 0) goto 110
  enddo

  write (qxldb, 5)
5 format ('*** QXBESSELJR: Loop end error 2')
  call qxabrt

110 continue

  t2 = qxpicon * qq
  t1 = t2 * 0.5e0_qxknd
  t2 = rra - t1
  t1 = qxpicon * 0.25e0_qxknd
  t3 = t2 - t1
  t41 = cos (t2)
  t42 = sin (t2)
  t1 = t41 * sum1
  t2 = t42 * sum2
  t5 = t1 - t2
  t1 = qxpicon * rra
  t2 = 2.e0_qxknd
  t3 = t2 / t1
  t1 = sqrt (t3)
  t3 = t1 * t5
endif

! if (mod (nu, 2) /= 0) then
!   if (nu < 0 .and. rr > 0.d0 .or. nu > 0 .and. rr < 0.d0) t3 = - t3
! endif

120 continue

ss = t3

return
end subroutine qxbesseljr

subroutine qxbesselknr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselK (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.31.1 for modest RR,
!   and DLMF 10.40.2 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (qxknd), intent(in):: rr
real (qxknd), intent(out):: ss
integer, parameter:: itrmax = 1000000
integer ic1, k, nua, n1
real (qxknd) d1
real (qxknd), parameter:: dfrac = 1.5e0_qxknd
real (qxknd) f1, f2, f3, f4, f5, sum1, sum2, sum3, td, &
  tn, t1, t2, t3, t4, rra, tc1, tc2, tc3, eps

!  End of declaration

!   Check for RR = 0.

if (qxsgn (rr) == 0) then
  write (qxldb, 2)
2 format ('*** QXBESSELKNR: Second argument is zero')
  call qxabrt
endif

nua = abs (nu)
eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
rra = abs (rr)
d1 = rra

if (d1 < dfrac * qxdpw) then
  t2 = rra**2
  t1 = t2 * 0.25e0_qxknd
  f1 = 1.e0_qxknd
  f2 = 1.e0_qxknd
  f3 = 1.e0_qxknd
  sum1 = 0.e0_qxknd

  do k = 1, nua - 1
    f1 = f1 * real (k, qxknd)
  enddo

  do k = 0, nua - 1
    if (k > 0) then
      f1 = f1 / real (nua - k, qxknd)
      f2 = -t1 * f2
      f3 = f3 * real (k, qxknd)
    endif
    t3 = f1 * f2
    t2 = t3 / f3
    sum1 = sum1 + t2
  enddo

  t2 = sum1 * 0.5e0_qxknd
  t3 = rra * 0.5e0_qxknd
  t4 = t3 ** nua
  sum1 = t2 / t4
  t3 = log (rra * 0.5e0_qxknd)
  d1 = (-1.e0_qxknd) ** (nua + 1)
  t2 = t3 * d1
  call qxbesselinr (nua, rra, t3)
  sum2 = t2 * t3
  f1 = - qxegammacon
  f2 = f1
  f3 = 1.e0_qxknd
  f4 = 1.e0_qxknd
  f5 = 1.e0_qxknd

  do k = 1, nua
    t2 = 1.e0_qxknd
    t3 = t2 / real (k, qxknd)
    f2 = f2 + t3
    f5 = f5 * real (k, qxknd)
  enddo

  t2 = f1 + f2
  t3 = t2 * f3
  t4 = f4 * f5
  sum3 = t3 / t4

  do k = 1, itrmax
    t2 = 1.e0_qxknd
    t3 = t2 / real (k, qxknd)
    f1 = f1 + t3
    t3 = t2 / real (nua + k, qxknd)
    f2 = f2 + t3
    f3 = t1 * f3
    f4 = f4 * real (k, qxknd)
    f5 = f5 * real (nua + k, qxknd)
    t2 = f1 + f2
    t3 = t2 * f3
    t4 = f4 * f5
    t2 = t3 / t4
    sum3 = sum3 + t2
    tc1 = abs (t2) - eps * abs (sum3)
    if (qxsgn (tc1) < 0) goto 100
  enddo

  write (qxldb, 5)
5 format ('*** QXBESSELKNR: Loop end error 1')
  call qxabrt

100 continue

  t2 = rra * 0.5e0_qxknd
  t3 = t2 ** nua
  d1 = (-1.e0_qxknd)**nua * 0.5e0_qxknd
  t4 = t3 * d1
  sum3 = t4 * sum3
  t2 = sum1 + sum2
  t3 = t2 + sum3
else
  sum1 = 1.e0_qxknd
  d1 = 4.e0_qxknd * real (nua, qxknd)**2
  t1 = d1
  tn = 1.e0_qxknd
  td = 1.e0_qxknd

  do k = 1, itrmax
    d1 = real (2*k - 1, qxknd)
    t2 = d1
    t3 = t2**2
    t2 = t1 - t3
    tn = tn * t2
    t2 = rra * real (8*k, qxknd)
    td = td * t2
    t4 = tn / td
    sum1 = sum1 + t4
    tc2 = abs (t4) - eps * abs (sum1)
    if (qxsgn (tc2) < 0) goto 110
  enddo

write (qxldb, 6)
6 format ('*** QXBESSELKNR: Loop end error 2')
call qxabrt

110 continue

  t1 = exp (rra)
  t2 = rra * 2.e0_qxknd
  t3 = qxpicon / t2
  t4 = sqrt (t3)
  t2 = t4 / t1
  t3 = t2 * sum1
endif

if (qxsgn (rr) < 0 .and. mod (nu, 2) /= 0) t3 = - t3
ss = t3

return
end subroutine qxbesselknr

subroutine qxbesselkr (qq, rr, ss)

!   This evaluates the Bessel function BesselK (QQ,RR) for QQ and RR
!   both MPR. This uses DLMF formula 10.27.4.

implicit none
real (qxknd), intent(in):: qq, rr
real (qxknd), intent(out):: ss
integer i0, i1, i2, n1
real (qxknd) d1
real (qxknd) t1, t2, t3, t4

!  End of declaration

!   If QQ is integer, call mpbesselknr; if qq < 0 and rr <= 0, then error.

i0 = qxsgn (qq - anint (qq))
i1 = qxsgn (qq)
i2 = qxsgn (rr)
if (i0 == 0) then
  n1 = nint (qq)
  call qxbesselknr (n1, rr, t1)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (qxldb, 2)
2 format ('*** QXBESSELKR: First argument < 0 and second argument <= 0')
  call qxabrt
endif

t1 = - qq
call qxbesselir (t1, rr, t2)
call qxbesselir (qq, rr, t3)
t4 = t2 - t3
t1 = qq * qxpicon
t3 = sin (t1)
t2 = t4 / t3
t3 = qxpicon * t2
t1 = t3 * 0.5e0_qxknd

120 continue

ss = t1
return
end subroutine qxbesselkr

subroutine qxbesselynr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselY (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.8.1 for modest RR,
!   and DLMF 10.17.4 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (qxknd), intent(in):: rr
real (qxknd), intent(out):: ss
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dfrac = 1.5e0_qxknd
integer ic1, ic2, k, nua, n1
real (qxknd) d1, d2
real (qxknd) f1, f2, f3, f4, f5, rra, rr2, sum1, sum2, sum3, td1, td2, &
  tn1, tn2, t1, t2, t3, t4, t41, t42, t5, tc1, tc2, tc3, eps

!  End of declaration
!   Check for RR = 0.

if (qxsgn (rr) == 0) then
  write (qxldb, 2)
2 format ('*** QXBESSELYNR: argument is negative or too large')
  call qxabrt
endif

nua = abs (nu)
eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
d1 = rr

if (d1 < dfrac * qxdpw) then
  rra = abs (rr)
  t2 = rra**2
  t1 = t2 * 0.25e0_qxknd
  f1 = 1.e0_qxknd
  f2 = 1.e0_qxknd
  f3 = 1.e0_qxknd
  sum1 = 0.e0_qxknd

  do k = 1, nua - 1
    f1 = f1 * real (k, qxknd)
  enddo

  do k = 0, nua - 1
    if (k > 0) then
      f1 = f1 / real (nua - k, qxknd)
      f2 = t1 * f2
      f3 = f3 * real (k, qxknd)
    endif
    t3 = f1 * f2
    t2 = t3 / f3
    sum1 = sum1 + t2
  enddo

  t3 = rra * 0.5e0_qxknd
  t4 = t3 ** nua
  sum1 = - sum1 / t4
  t2 = rra * 0.5e0_qxknd
  t3 = log (t2)
  t2 = t3 * 2.e0_qxknd
  call qxbesseljnr (nua, rra, t3)
  sum2 = t2 * t3

  f1 = -qxegammacon
  f2 = f1
  f3 = 1.e0_qxknd
  f4 = 1.e0_qxknd
  f5 = 1.e0_qxknd

  do k = 1, nua
    t2 = 1.e0_qxknd
    t3 = t2 / real (k, qxknd)
    f2 = f2 + t3
    f5 = f5 * real (k, qxknd)
  enddo

  t2 = f1 + f2
  t3 = t2 * f3
  t4 = f4 * f5
  sum3 = t3 / t4

  do k = 1, itrmax
    t2 = 1.e0_qxknd
    t3 = t2 / real (k, qxknd)
    f1 = f1 + t3
    t3 = t2 / real (nua + k, qxknd)
    f2 = f2 + t3
    f3 = - t1 * f3
    f4 = f4 * real (k, qxknd)
    f5 = f5 * real (nua + k, qxknd)
    t2 = f1 + f2
    t3 = t2 * f3
    t4 = f4 * f5
    t2 = t3 / t4
    sum3 = sum3 + t2
    tc1 = abs (t2) - eps * abs (sum3)
    if (qxsgn (tc1) < 0) goto 100
  enddo

  write (qxldb, 6)
6 format ('*** QXBESSELYNR: Loop end error 1')
  call qxabrt

100 continue

  t2 = rra * 0.5e0_qxknd
  t3 = t2 ** nua
  sum3 = -t3 * sum3
  t2 = sum1 + sum2
  t4 = t2 + sum3
  t2 = qxpicon
  t3 = t4 / t2
else
  rra = abs (rr)
  rr2 = rra**2
  t1 = 4.e0_qxknd * real (nua, qxknd)**2
  tn1 = 1.e0_qxknd
  t2 = t1 - tn1
  tn2 = t2 / 8.e0_qxknd
  td1 = 1.e0_qxknd
  td2 = rra
  sum1 = tn1 / td1
  sum2 = tn2 / td2

  do k = 1, itrmax
    d1 = real (4 * k - 3, qxknd)**2
    d2 = real (4 * k - 1, qxknd)**2
    t2 = d1
    t3 = t1 - t2
    t2 = d2
    t5 = t1 - t2
    t2 = t3 * t5
    tn1 = -tn1 * t2
    d1 = real (2*k - 1, qxknd) * real (2*k, qxknd) * 64.e0_qxknd
    t2 = td1 * d1
    td1 = t2 * rr2
    t41 = tn1 / td1
    sum1 = sum1 + t41

    d1 = real (2 * (2*k) - 1, qxknd)**2
    d2 = real (2 * (2*k + 1) - 1, qxknd)**2
    t2 = d1
    t3 = t1 - t2
    t2 = d2
    t5 = t1 - t2
    t2 = t3 * t5
    tn2 = - tn2 * t2
    d1 = real (2*k, qxknd) * real (2*k + 1, qxknd) * 64.e0_qxknd
    t2 = td2 * d1
    td2 = t2 * rr2
    t42 = tn2 / td2
    sum2 = sum2 + t42
    tc1 = abs (t41) - eps * abs (sum1)
    tc2 = abs (t42) - eps * abs (sum2)
    if (qxsgn (tc1) < 0 .and. qxsgn (tc2) < 0) goto 110
  enddo

  write (qxldb, 5)
5 format ('*** QXBESSELYNR: Loop end error 2')
  call qxabrt

110 continue

  t1 = qxpicon * 0.5e0_qxknd * real (nua, qxknd)
  t2 = rra - t1
  t1 = qxpicon * 0.25e0_qxknd
  t3 = t2 - t1
  t41 = cos (t3)
  t42 = sin (t3)
  t1 = t42 * sum1
  t2 = t41 * sum2
  t5 = t1 + t2
  t1 = qxpicon * rra
  t2 = 2.e0_qxknd
  t3 = t2 / t1
  t1 = sqrt (t3)
  t3 = t1 * t5
endif

if (mod (nu, 2) /= 0) then
  ic1 = qxsgn (rr)
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) t3 = - t3
endif

ss = t3

return
end subroutine qxbesselynr

subroutine qxbesselyr (qq, rr, ss)

!   This evaluates the modified Bessel function BesselY (QQ,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2.

implicit none
real (qxknd), intent(in):: qq, rr
real (qxknd), intent(out):: ss
integer i0, i1, i2, n1
real (qxknd) d1
real (qxknd) t1, t2, t3, t4

!  End of declaration

!   If QQ is integer, call mpbesselynr; if qq < 0 and rr <= 0, then error.

i0 = qxsgn (qq - anint (qq))
i1 = qxsgn (qq)
i2 = qxsgn (rr)
if (i0 == 0) then
  n1 = nint (qq)
  call qxbesselynr (n1, rr, t1)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (qxldb, 2)
2 format ('*** QXBESSELYR: First argument < 0 and second argument <= 0')
  call qxabrt
endif

t1 = qq * qxpicon
t2 = cos (t1)
t3 = sin (t1)
call qxbesseljr (qq, rr, t4)
t1 = t4 * t2
t2 = - qq
call qxbesseljr (t2, rr, t4)
t2 = t1 - t4
t1 = t2 / t3

120 continue

ss = t1

return
end subroutine qxbesselyr

subroutine qxdigammabe (nb2, berne, x, y)

!  This evaluates the digamma function, using asymptotic formula DLMF 5.11.2:
!  dig(x) ~ log(x) - 1/(2*x) - Sum_{k=1}^inf B[2k] / (2*k*x^(2*k)).
!  Before using this formula, the recursion dig(x+1) = dig(x) + 1/x is used
!  to shift the argument up by IQQ, where IQQ is set based on MPNW below.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent (in):: nb2
real (qxknd), intent(in):: berne(nb2), x
real (qxknd), intent(out):: y
real (qxknd), parameter:: dber = 1.4e0_qxknd, dfrac = 0.4e0_qxknd
integer k, i1, i2, ic1, iqq, n1
real (qxknd) d1
real (qxknd) f1, sum1, sum2, t1, t2, t3, t4, t5, xq, tc1, tc2, tc3, eps

!  End of declaration

iqq = dfrac * qxdpw
eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
f1 = 1.e0_qxknd

!   Check if argument is less than or equal to 0 -- undefined.

if (qxsgn (x) <= 0) then
  write (qxldb, 2)
2 format ('*** QXDIGAMMABE: Argument <= 0')
  call qxabrt
endif

!   Check if berne array has been initialized.

d1 = berne(1)
if (abs (d1 - 1.e0_qxknd / 6.e0_qxknd) > qxrdfz .or. nb2 < int (dber * qxdpw)) then
  write (qxldb, 3) int (dber * qxdpw)
3 format ('*** QXDIGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries using QXBERNE or QXBERNER.')
  call qxabrt
endif

sum1 = 0.e0_qxknd
sum2 = 0.e0_qxknd
xq = x + real (iqq, qxknd)

do k = 0, iqq - 1
  t1 = real (k, qxknd)
  t2 = x + t1
  t3 = f1 / t2
  sum1 = sum1 + t3
enddo

t1 = 1.e0_qxknd
t2 = xq**2

do k = 1, nb2
  t1 = t1 * t2
  t4 = t1 * 2.e0_qxknd * real (k, qxknd)
  t3 = berne(k) / t4
  sum2 = sum2 + t3
  tc1 = abs (t3) - eps * abs (sum2)
  if (qxsgn (tc1) < 0) goto 110
enddo

write (qxldb, 4)
4 format ('*** QXDIGAMMABE: Loop end error: Increase NB2')
call qxabrt

110 continue

t1 = - sum1
t2 = log (xq)
t3 = t1 + t2
t4 = xq * 2.e0_qxknd
t5 = f1 / t4
t2 = t3 - t5
y = t2 - sum2

return
end subroutine qxdigammabe

subroutine qxerfr (z, terf)

!   This evaluates the erf function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (z == 0) then
!     erf = 0
!   elseif (z > sqrt(B*log(2))) then
!     erf = 1
!   elseif (z < -sqrt(B*log(2))) then
!     erf = -1
!   elseif (abs(z) < B/dcon + 8) then
!     erf = 2 / (sqrt(pi)*exp(z^2)) * Sum_{k>=0} 2^k * z^(2*k+1)
!             / (1.3....(2*k+1))
!   else
!     erf = 1 - 1 / (sqrt(pi)*exp(z^2))
!             * Sum_{k>=0} (-1)^k * (1.3...(2*k-1)) / (2^k * z^(2*k+1))
!   endif

implicit none
real (qxknd), intent(in):: z
real (qxknd), intent(out):: terf
integer, parameter:: itrmx = 100000
real (qxknd), parameter:: dcon = 100.e0_qxknd
integer ic1, ic2, ic3, k, nbt, n1
real (qxknd) d1, d2

real (qxknd) t1, t2, t3, t4, t5, t6, t7, z2, tc1, tc2, tc3, eps

!  End of declaration

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
nbt = qxnbt
d1 = aint (1.e0_qxknd + sqrt (nbt * log (2.e0_qxknd)))
d2 = aint (nbt / dcon + 8.e0_qxknd)
t1 = d1
t2 = d2

if (qxsgn (z) == 0) then
  terf = 0.e0_qxknd
elseif (qxsgn (z - t1) > 0) then
  terf = 1.e0_qxknd
elseif (qxsgn (z + t1) < 0) then
  terf = -1.e0_qxknd
elseif (qxsgn (z - t2) < 0) then
  z2 = z**2
  t1 = 0.e0_qxknd
  t2 = z
  t3 = 1.e0_qxknd
  t5 = 10.e0_qxknd

  do k = 0, itrmx
    if (k > 0) then
      t6 = z2 * 2.e0_qxknd
      t2 = t6 * t2
      d1 = real (2 * k + 1, qxknd)
      t3 = t3 * d1
    endif

    t4 = t2 / t3
    t1 = t1 + t4
    t6 = t4 / t1
    tc1 = t6 - eps
    tc2 = t5 - t6
    if (qxsgn (tc1) < 0 .or. qxsgn (tc2) < 0) goto 120
    t5 = t6
  enddo

write (qxldb, 3) 1, itrmx
3 format ('*** QXERFR: iteration limit exceeded',2i10)
call qxabrt

120 continue

  t3 = t1 * 2.e0_qxknd
  t4 = sqrt (qxpicon)
  t5 = exp (z2)
  t6 = t4 * t5
  terf = t3 / t6
else
  z2 = z**2
  t1 = 0.e0_qxknd
  t2 = 1.e0_qxknd
  t3 = z
  t5 = 10.e0_qxknd

  do k = 0, itrmx
    if (k > 0) then
      d1 = -(2.e0_qxknd * k - 1.e0_qxknd)
      t2 = t2 * d1
      t3 = t2 * t3
    endif

    t4 = t2 / t3
    t1 = t1 + t4
    t6 = t4 / t1
    tc1 = t6 - eps
    tc2 = t5 - t6
    if (qxsgn (tc1) < 0 .or. qxsgn (tc2) < 0) goto 130
    t5 = t6
  enddo

write (qxldb, 3) 2, itrmx
call qxabrt

130 continue

  t2 = 1.e0_qxknd
  t3 = sqrt (qxpicon)
  t4 = exp (z2)
  t5 = t3 * t4
  t6 = t1 / t5
  terf = t2 - t6
  if (qxsgn (z) < 0) terf = - terf
endif

return
end subroutine qxerfr

subroutine qxerfcr (z, terfc)

!   This evaluates the erfc function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (z == 0) then
!     erfc = 1
!   elseif (z > sqrt(B*log(2))) then
!     erfc = 0
!   elseif (z < -sqrt(B*log(2))) then
!     erfc = 2
!   elseif (abs(z) < B/dcon + 8) then
!     erfc = 1 - 2 / (sqrt(pi)*exp(z^2)) * Sum_{k>=0} 2^k * z^(2*k+1)
!               / (1.3....(2*k+1))
!   else
!     erfc = 1 / (sqrt(pi)*exp(z^2))
!             * Sum_{k>=0} (-1)^k * (1.3...(2*k-1)) / (2^k * z^(2*k+1))
!   endif

implicit none
real (qxknd), intent(in):: z
real (qxknd), intent(out):: terfc
integer, parameter:: itrmx = 100000
real (qxknd), parameter:: dcon = 100.e0_qxknd
integer ic1, ic2, ic3, k, nbt, n1
real (qxknd) d1, d2

real (qxknd) t1, t2, t3, t4, t5, t6, t7, z2, tc1, tc2, tc3, eps

!  End of declaration

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
nbt = qxnbt
d1 = aint (1.e0_qxknd + sqrt (nbt * log (2.e0_qxknd)))
d2 = aint (nbt / dcon + 8.e0_qxknd)
t1 = d1
t2 = d2

if (qxsgn (z) == 0) then
  terfc = 1.e0_qxknd
elseif (qxsgn (z - t1) > 0) then
  terfc = 0.e0_qxknd
elseif (qxsgn (z + t1) < 0) then
  terfc = 2.e0_qxknd
elseif (qxsgn (z - t2) < 0) then
  z2 = z**2
  t1 = 0.e0_qxknd
  t2 = z
  t3 = 1.e0_qxknd
  t5 = 10.e0_qxknd

  do k = 0, itrmx
    if (k > 0) then
      t6 = z2 * 2.e0_qxknd
      t2 = t6 * t2
      d1 = real (2 * k + 1, qxknd)
      t3 = t3 * d1
    endif

    t4 = t2 / t3
    t1 = t1 + t4
    t6 = t4 / t1
    tc1 = t6 - eps
    tc2 = t5 - t6
    if (qxsgn (tc1) < 0 .or. qxsgn (tc2) < 0) goto 120
    t5 = t6
  enddo

write (qxldb, 3) 1, itrmx
3 format ('*** QXERFRC: iteration limit exceeded',2i10)
call qxabrt

120 continue

  t2 = 1.e0_qxknd
  t3 = t1 * 2.e0_qxknd
  t4 = sqrt (qxpicon)
  t5 = exp (z2)
  t6 = t4 * t5
  terfc = t2 - t3 / t6
else
  z2 = z**2
  t1 = 0.e0_qxknd
  t2 = 1.e0_qxknd
  t3 = z
  t5 = 10.e0_qxknd

  do k = 0, itrmx
    if (k > 0) then
      d1 = -(2.e0_qxknd * k - 1.e0_qxknd)
      t2 = t2 * d1
      t3 = t2 * t3
    endif

    t4 = t2 / t3
    t1 = t1 + t4
    t6 = t4 / t1
    tc1 = t6 - eps
    tc2 = t5 - t6
    if (qxsgn (tc1) < 0 .or. qxsgn (tc2) < 0) goto 130
    t5 = t6
  enddo

write (qxldb, 3) 2, itrmx
call qxabrt

130 continue

  t3 = sqrt (qxpicon)
  t4 = exp (z2)
  t5 = t3 * t4
  t6 = t1 / t5
  if (qxsgn (z) < 0) t6 = 2.e0_qxknd - t6
  terfc = t6
endif

return
end subroutine qxerfcr

subroutine qxexpint (x, y)

!   This evaluates the exponential integral function Ei(x):
!   Ei(x) = - incgamma (0, -x)

implicit none
real (qxknd), intent(in):: x
real (qxknd), intent(out):: y
real (qxknd) t1, t2, t3

!  End of declaration

if (qxsgn (x) == 0) then
  write (qxldb, 2)
2 format ('*** QXEXPINT: argument is zero')
  call qxabrt
endif

t1 = 0.e0_qxknd
t2 = -x
call qxincgammar (t1, t2, t3)
y = - t3
return
end subroutine qxexpint

subroutine qxgammar (t, z)

!   This evaluates the gamma function, using an algorithm of R. W. Potter.
!   The argument t must not exceed 10^8 in size (this limit is set below),
!   must not be zero, and if negative must not be integer.

!   In the parameter statement below:
!     itrmx = limit of number of iterations in series; default = 100000.
!     con1 = 1/2 * log (10) to DP accuracy.
!     dmax = maximum size of input argument.

implicit none
real (qxknd), intent(in):: t
real (qxknd), intent(out):: z
integer, parameter:: itrmx = 100000
real (qxknd), parameter:: dmax = 1e8_qxknd
integer i, i1, ic1, j, nt, n1, n2, n3, qxnw
real (qxknd) alpha, al2, d1, d2, d3
real (qxknd) f1, sum1, sum2, tn, t1, t2, t3, t4, t5, t6, tc1, tc2, tc3, eps

!  End of declaration

qxnw = qxnwx
eps = 2.e0_qxknd ** (-qxnw*qxnbt)

ic1 = t - anint (t)
i1 = qxsgn (t)
if (i1 == 0 .or. d1 > dmax .or. (i1 < 0 .and. ic1 == 0)) then
  write (qxldb, 2) dmax
2 format ('*** QXGAMMAR: input argument must have absolute value <=',f10.0,','/ &
  'must not be zero, and if negative must not be an integer.')
  call qxabrt
endif

al2 = qxlog2con
f1 = 1.e0_qxknd

!   Find the integer and fractional parts of t.

t2 = aint (t)
t3 = t - t2

if (qxsgn (t3) == 0) then

!   If t is a positive integer, then apply the usual factorial recursion.

  nt = int (t2)
  t1 = f1

  do i = 2, nt - 1
    t1 = t1 * real (i, qxknd)
  enddo

  z = t1
  goto 120
elseif (qxsgn (t) > 0) then

!   Apply the identity Gamma[t+1] = t * Gamma[t] to reduce the input argument
!   to the unit interval.

  nt = t2
  t1 = f1
  tn = t3

  do i = 1, nt
    t4 = real (i, qxknd)
    t5 = t - t4
    t1 = t1 * t5
  enddo
else

!   Apply the gamma identity to reduce a negative argument to the unit interval.

  t4 = f1 - t
  t3 = aint (t4)
  t5 = t4 - t3
  nt = t3
  t1 = f1
  t2 = f1 - t5
  tn = t2

  do i = 0, nt - 1
    t4 = real (i, qxknd)
    t5 = t + t4
    t1 = t1 / t5
  enddo
endif

!   Calculate alpha = bits of precision * log(2) / 2, then take the next even
!   integer value, so that alpha/2 and alpha^2/4 can be calculated exactly in DP.

alpha = 2.e0_qxknd * aint (0.25e0_qxknd * qxnw * qxnbt * al2 + 1.e0_qxknd)
d2 = 0.25e0_qxknd * alpha**2
t2 = tn
t3 = f1 / t2
sum1 = t3

!   Evaluate the series with t.

do j = 1, itrmx
  t6 = real (j, qxknd)
  t4 = t2 + t6
  t5 = t4 * real (j, qxknd)
  t6 = t3 / t5
  t3 = t6 * d2
  sum1 = sum1 + t3
  tc1 = abs (t3) - eps * abs (sum1)
  if (qxsgn (tc1) < 0) goto 100
enddo

write (qxldb, 3) 1, itrmx
3 format ('*** QXGAMMAR: iteration limit exceeded',2i10)
call qxabrt

100 continue

t2 = - tn
t3 = f1 / t2
sum2 = t3

!   Evaluate the same series with -t.

do j = 1, itrmx
  t6 = real (j, qxknd)
  t4 = t2 + t6
  t5 = t4 * real (j, qxknd)
  t6 = t3 / t5
  t3 = t6 * d2
  sum2 = sum2 + t3
  tc1 = abs (t3) - eps * abs (sum2)
  if (qxsgn (tc1) < 0) goto 110
enddo

write (qxldb, 3) 2, itrmx
call qxabrt

110 continue

!   Compute sqrt (pi * sum1 / (tn * sin (pi * tn) * sum2))
!   and (alpha/2)^tn terms. Also, multiply by the factor t1, from the
!   If block above.

t2 = qxpicon
t3 = t2 * tn
t4 = cos (t3)
t5 = sin (t3)
t6 = t5 * sum2
t5 = tn * t6
t3 = t2 * sum1
t6 = - t3 / t5
t2 = sqrt (t6)
t3 = 0.5e0_qxknd * alpha
t4 = log (t3)
t5 = tn * t4
t6 = exp (t5)
t3 = t2 * t6
t4 = t1 * t3
z = t4

120 continue

return
end subroutine qxgammar

subroutine qxhurwitzzetan (is, aa, zz)

!   This returns the Hurwitz zeta function of IS and AA, using an algorithm from:
!   David H. Bailey and Jonathan M. Borwein, "Crandall's computation of the
!   incomplete gamma function and the Hurwitz zeta function with applications to
!   Dirichlet L-series," Applied Mathematics and Computation, vol. 268C (Oct 2015),
!   pg. 462-477, preprint at:
!   https://www.davidhbailey.com/dhbpapers/lerch.pdf
!   This is limited to IS >= 2 and 0 < AA < 1.

implicit none
integer, intent(in):: is
real (qxknd), intent(in):: aa
real (qxknd), intent(out):: zz
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: pi = 3.1415926535897932385e0_qxknd
integer i1, ic1, ic2, ic3, k, n1
real (qxknd) d1, dk
real (qxknd) gs1, gs2, ss, sum1, sum2, sum3, ss1, ss2, ss3, ss4, s1, s2, s3, &
  t1, t2, t3, t4, t5, t6, t7, t8, tc1, tc2, tc3, eps

!  End of declaration

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)

if (is <= 0) then
  write (qxldb, 3)
3 format ('*** QXHURWITZZETAN: IS less than or equal to 0:')
  call qxabrt
endif

t1 = 0.e0_qxknd
t2 = 1.e0_qxknd
ic1 = qxsgn (aa - t1)
ic2 = qxsgn (aa - t2)
if (ic1 <= 0 .or. ic2 >= 0) then
  write (qxldb, 4)
4 format ('*** QXHURWITZZETAN: AA must be in the range (0, 1)')
  call qxabrt
endif

ss = real (is, qxknd)
ss1 = ss * 0.5e0_qxknd
t1 = 1.e0_qxknd
t2 = t1 + ss
ss2 = t2 * 0.5e0_qxknd
t2 = t1 - ss
ss3 = t2 * 0.5e0_qxknd
ss4 = t1 - ss1
gs1 = gamma (ss1)
gs2 = gamma (ss2)
t2 = aa**2
t1 = qxpicon * t2

call qxincgammar (ss1, t1, t2)
t3 = t2 / gs1
call qxincgammar (ss2, t1, t2)
t4 = t2 / gs2
t2 = t3 + t4
t3 = abs (aa)
t4 = t3 ** is
sum1 = t2 / t4
sum2 = 0.e0_qxknd
sum3 = 0.e0_qxknd

do k = 1, itrmax
  dk = real (k, qxknd)
  t5 = dk
  t6 = t5 + aa
  t5 = t6**2
  t1 = qxpicon * t5
  t5 = - dk
  t6 = t5 + aa
  t7 = t6**2
  t2 = qxpicon * t7
  t6 = t5**2
  t3 = qxpicon * t6
  t5 = qxpicon * (2.e0_qxknd * dk)
  t4 = t5 * aa

  call qxincgammar (ss1, t1, t5)
  t6 = t5 / gs1
  call qxincgammar (ss2, t1, t5)
  t7 = t5 / gs2
  t5 = t6 + t7
  t6 = dk
  t7 = t6 + aa
  t6 = abs (t7)
  t7 = t6 ** is
  s1 = t5 / t7

  call qxincgammar (ss1, t2, t5)
  t6 = t5 / gs1
  call qxincgammar (ss2, t2, t5)
  t7 = t5 / gs2
  t5 = t6 - t7
  t6 = - dk
  t7 = t6 + aa
  t6 = abs (t7)
  t7 = t6 ** is
  s2 = t5 / t7
  sum1 = sum1 + s1
  sum2 = sum2 + s2

  call qxincgammar (ss3, t3, t5)
  t6 = cos (t4)
  t7 = sin (t4)
  t8 = t5 * t6
  t5 = t8 / gs1
  call qxincgammar (ss4, t3, t6)
  t8 = t6 * t7
  t6 = t8 / gs2
  t7 = t5 + t6
  t5 = dk
  i1 = 1 - is
  t6 = t5 ** i1
  s3 = t7 / t6
  sum3 = sum3 + s3
  tc1 = abs (s1) - eps * abs (sum1)
  tc2 = abs (s2) - eps * abs (sum2)
  tc3 = abs (s3) - eps * abs (sum3)
  if (qxsgn (tc1) < 0 .and. qxsgn (tc2) < 0 .and. qxsgn (tc3) < 0) goto 100
enddo

write (qxldb, 5)
5 format ('*** QXHURWITZZETAN: Loop end error')
call qxabrt

100 continue

if (mod (is, 2) == 0) then
  i1 = is / 2
  t2 = qxpicon ** i1
  t3 = 1.e0_qxknd
  t4 = ss - t3
  t5 = gamma (ss1)
  t6 = t4 * t5
  t1 = t2 / t6
else
  i1 = (is - 1) / 2
  t2 = qxpicon ** i1
  t3 = sqrt (qxpicon)
  t4 = t2 * t3
  t2 = 1.e0_qxknd
  t3 = ss - t2
  t5 = gamma (ss1)
  t6 = t3 * t5
  t1 = t4 / t6
endif

t3 = qxpicon ** is
t4 = sqrt (qxpicon)
t2 = t3 / t4
t3 = sum1 * 0.5e0_qxknd
t4 = sum2 * 0.5e0_qxknd
t5 = sum3 * t2
t6 = t1 + t3
t7 = t6 + t4
zz = t7 + t5

return
end subroutine qxhurwitzzetan

subroutine qxhurwitzzetanbe (nb2, berne, iss, aa, zz)

!  This evaluates the Hurwitz zeta function, using the combination of
!  the definition formula (for large iss), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb2, iss
real (qxknd), intent(in):: berne(nb2), aa
real (qxknd), intent(out):: zz
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dber = 1.4e0_qxknd, dcon = 0.6e0_qxknd
integer i1, i2, ic1, iqq, k, n1
real (qxknd) d1, dp
real (qxknd) aq, aq2, s1, s2, s3, s4, t1, t2, t3, t4, t5, t6, eps, f1, &
  tc1, tc2, tc3

!  End of declaration

!   Check if berne array has been initialized.

d1 = berne(1)
if (abs (d1 - 1.e0_qxknd / 6.e0_qxknd) > qxrdfz .or. nb2 < int (dber * qxdpw)) then
  write (qxldb, 3) int (dber * qxdpw)
3 format ('*** QXHURWITZZETANBE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries using QXBERNE or QXBERNER.')
  call qxabrt
endif

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
f1 = 1.e0_qxknd
s1 = 0.e0_qxknd
s2 = 0.e0_qxknd
s3 = 0.e0_qxknd
s4 = 0.e0_qxknd

if (iss <= 0) then
  write (qxldb, 4)
4 format ('*** QXHURWITZZETANBE: ISS <= 0')
  call qxabrt
endif

if (qxsgn (aa) < 0) then
  write (qxldb, 5)
5 format ('*** QXHURWITZZETANBE: AA < 0')
  call qxabrt
endif

dp = anint (qxdpw)

!   If iss > a certain value, then use definition formula.

if (iss > 2.303e0_qxknd * dp / log (2.515e0_qxknd * dp)) then
  do k = 0, itrmax
!    t1 = 1.d0 / (aa + dble (k))**iss
!    s1 = s1 + t1

    t1 = real (k, qxknd)
    t2 = aa + t1
    t3 = t2 ** iss
    t1 = f1 / t3
    s1 = s1 + t1
    tc1 = abs (t1) - eps * abs (s1)
    if (qxsgn (tc1) < 0) goto 110
  enddo

  write (6, 6)
6 format ('*** QXHURWITZZETANBE: Loop end error 1')
  call qxabrt
endif

d1 = aa
iqq = max (dcon * qxdpw - d1, 0.e0_qxknd)

do k = 0, iqq - 1
  t1 = real (k, qxknd)
  t2 = aa + t1
  t3 = t2 ** iss
  t1 = f1 / t3
  s1 = s1 + t1
enddo

t1 = real (iqq, qxknd)
aq = aa + t1
t1 = real (iss - 1, qxknd)
t2 = aq ** (iss - 1)
t3 = t1 * t2
s2 = f1 / t3
t1 = aq ** iss
t2 = t1 * 2.e0_qxknd
s3 = f1 / t2

t1 = real (iss, qxknd)
t2 = 1.e0_qxknd
t3 = aq ** (iss - 1)
aq2 = aq**2

do k = 1, nb2
  if (k > 1) then
    t5 = t1 * real (iss + 2*k - 3, qxknd)
    t1 = t5 * real (iss + 2*k - 2, qxknd)
  endif
  t5 = t2 * real (2 * k - 1, qxknd)
  t2 = t5 * real (2 * k, qxknd)
  t3 = t3 * aq2
  t5 = berne(k) * t1
  t6 = t2 * t3
  t4 = t5 / t6
  s4 = s4 + t4
  tc2 = abs (t4) - eps * abs (s4)
  if (qxsgn (tc2) < 0) goto 110
enddo

write (6, 7)
7 format ('*** QXHURWITZZETANBE: End loop error 2; call qxBERNE with larger NB.')
call qxabrt

110 continue

t5 = s1 + s2
t6 = t5 + s3
zz = t6 + s4

return
end subroutine qxhurwitzzetanbe

subroutine qxhypergeompfq (np, nq, aa, bb, xx, yy)

!  This returns the HypergeometricPFQ function, namely the sum of the infinite series

!  Sum_0^infinity poch(aa(1),n)*poch(aa(2),n)*...*poch(aa(np),n) /
!      poch(bb(1),n)*poch(bb(2),n)*...*poch(bb(nq),n) * xx^n / n!

!  This subroutine evaluates the HypergeometricPFQ function directly according to
!  this definitional formula. The arrays aa and bb must be dimensioned as shown below.
!  NP and NQ are limited to [1,10].

implicit none
integer, intent(in):: np, nq
real (qxknd), intent(in):: aa(np), bb(nq), xx
real (qxknd), intent(out):: yy
integer, parameter:: itrmax = 1000000, npq = 10
integer i1, i2, ic1, j, k
real (qxknd) sum, td, tn, t1, t2, t3, t4, tc1, tc2, tc3, eps

!  End of declaration

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)

if (np < 1 .or. np > npq .or. nq < 1 .or. nq > npq) then
  write (qxldb, 2) npq
2 format ('*** QXHYPERGEOMPFQ: NP and NQ must be between 1 and',i4)
  call qxabrt
endif

sum = 1.e0_qxknd
td = 1.e0_qxknd
tn = 1.e0_qxknd

do k = 1, itrmax
  t1 = real (k - 1, qxknd)

  do j = 1, np
    t2 = t1 + aa(j)
    tn = tn * t2
  enddo

  do j = 1, nq
    t2 = t1 + bb(j)
    td = td * t2
  enddo

  tn = tn * xx
  td = td * real (k, qxknd)
  t1 = tn / td
  sum = sum + t1
  tc1 = abs (t1) - eps * sum
  if (qxsgn (tc1) < 0) goto 100
enddo

write (qxldb, 3) itrmax
3 format ('*** QXHYPERGEOMPFQ: Loop end error',i10)
call qxabrt

100  continue

yy = sum

return
end subroutine qxhypergeompfq

subroutine qxincgammar (s, z, g)

!  This returns the incomplete gamma function, using a combination of formula
!  8.7.3 of the DLMF (for modest-sized z), formula 8.11.2 (for large z),
!  a formula from the Wikipedia page for the case S = 0, and another formula
!  from the Wikipedia page for the case S = negative integer. The formula
!  for the case S = 0 requires increased working precision, up to 2.5X normal,
!  depending on the size of Z.

implicit none
real (qxknd), intent(in):: s, z
real (qxknd), intent(out):: g
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dmax = 0.833e0_qxknd
integer ic1, k, nn, n1, n2
real (qxknd) d1, d2, bits
real (qxknd) t0, t1, t2, t3, t4, t5, f1, tc1, tc2, tc3, eps

!  End of declaration

n1 = qxsgn (s)
n2 = qxsgn (z)
if (n2 == 0 .or. n1 /= 0 .and. n2 < 0) then
  write (qxldb, 2)
2 format ('*** QXINCGAMMAR: The second argument must not be zero,'/ &
    'and must not be negative unless the first is zero.')
  call qxabrt
endif

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
f1 = 1.e0_qxknd
d1 = z
bits = qxnbt

if (abs (d1) < dmax * bits) then

!   This is for modest-sized z.

  t1 = aint (s)
  ic1 = qxsgn (s - t1)
  nn = s

  if (ic1 == 0 .and. nn == 1) then
    t1 = exp (-z)
    goto 200
  elseif (ic1 == 0 .and. nn <= 0) then

!    S is zero or a negative integer -- use a different algorithm. In
!    either event, first compute incgamma for S = 0. For large Z, the
!    working precision must be increased, up to 2.5X times normal.

!    mpnw2 = min (mpnw1 + 1.5d0 * d1 / (dmax * bits) * mpnw, 5*mpnw/2+1.d0)

    t0 = z
    t1 = z
    t2 = 1.e0_qxknd

    do k = 2, itrmax
      if (mod (k, 2) == 1) then
        d1 = real (k, qxknd)
        t3 = f1 / d1
        t2 = t2 + t3
      endif
      t3 = z * t1
      d1 = real (2 * k, qxknd)
      t1 = t3 / d1
      t3 = t1 * t2
      t0 = t0 + t3
      tc1 = abs (t3) - eps * abs (t0)
      if (qxsgn (tc1) < 0) goto 100
    enddo

    write (qxldb, 4)
4   format ('*** QXINCGAMMAR: Loop end error 1')
    call qxabrt

100  continue

    t1 = -qxegammacon
    t3 = abs (z)
    t2 = log (t3)
    t3 = t1 - t2
    t4 = z * (-0.5e0_qxknd)
    t5 = exp (t4)
    t4 = t5 * t0
    t1 = t3 + t4
    if (nn == 0) goto 200

!   S is negative integer (not zero).

    nn = abs (nn)
    t0 = 1.e0_qxknd
    t2 = t0

    do k = 1, nn - 1
      t0 = t0 * real (k, qxknd)
      t2 = t0
    enddo

    t5 = t0 * real (nn, qxknd)

    do k = 1, nn - 1
      t3 = t2 * z
      t4 = t3 / real (nn - k, qxknd)
      t2 = - t4
      t0 = t0 + t2
    enddo

    t2 = exp (z)
    t3 = t0 / t2
    t4 = z**nn
    t2 = t3 / t4

    if (mod (nn, 2) == 0) then
      t3 = t2 + t1
    else
      t3 = t2 - t1
    endif
    t1 = t3 / t5
    goto 200
  endif

  t1 = gamma (s)
  t3 = s * t1
  t2 = f1 / t3
  t0 = t2

  do k = 1, itrmax
    t5 = t2  * z
    t3 = real (k, qxknd)
    t4 = s + t3
    t2 = t5 / t4
    t0 = t0 + t2
    tc2 = abs (t2) - eps * abs (t0)
    if (qxsgn (tc2) < 0) goto 110
  enddo

  write (qxldb, 5) itrmax
5   format ('*** QXINCGAMMAR: Loop end error 1')
  call qxabrt

110 continue

  t2 = z ** s
  t3 = exp (z)
  t4 = t2 / t3
  t5 = t4 * t0
  t2 = f1 - t5
  t1 = t1 * t2
  goto 200
else

!   This is for large z. Note that if S is a positive integer, this loop
!   is finite.

  t0 = 1.e0_qxknd
  t1 = 1.e0_qxknd

  do k = 1, itrmax
    t2 = real (k, qxknd)
    t3 = s - t2
    t4 = t1 * t3
    t1 = t4 / z
    t0 = t0 + t1
    tc3 = abs (t1) - eps * abs (t0)
    if (qxsgn (tc3) < 0) goto 120
  enddo

  write (qxldb, 6)
6 format ('*** QXINCGAMMAR: Loop end error 3')
  call qxabrt

120 continue

  t2 = s - f1
  t3 = z ** t2
  t4 = exp (z)
  t2 = t3 / t4
  t1 = t2 * t0
  goto 200
endif

200 continue

g = t1

return
end subroutine qxincgammar

subroutine qxpolygamma (nn, x, y)

!   This returns polygamma (nn, x) for nn >= 0 and 0 < x < 1, by calling
!   mphurwitzzetan.

implicit none
integer, intent(in):: nn
real (qxknd), intent(in):: x
real (qxknd), intent(out):: y
integer ic1, ic2, k
real (qxknd) t1, t2, t3

!  End of declaration

if (nn <= 0) then
  write (qxldb, 2)
2 format ('*** QXPOLYGAMMA: NN <= 0')
  call qxabrt
endif

t1 = 0.e0_qxknd
t2 = 1.e0_qxknd
ic1 = qxsgn (x - t1)
ic2 = qxsgn (x - t2)

if (ic1 <= 0 .or. ic2 >= 0) then
  write (qxldb, 3)
3 format ('*** QXPOLYGAMMA: X must be in the range (0, 1)')
  call qxabrt
endif

t1 = 1.e0_qxknd

do k = 1, nn
  t1 = t1 * real (k, qxknd)
enddo

if (mod (nn + 1, 2) == 1) t1 = - t1
call qxhurwitzzetan (nn + 1, x, t2)
y = t1 * t2

return
end subroutine qxpolygamma

subroutine qxpolygammabe (nb2, berne, nn, x, y)

!  This returns polygamma (nn, x) for nn >= 0, by calling mphurwitzzetanbe.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb2, nn
real (qxknd), intent(in):: berne(nb2), x
real (qxknd), intent(out):: y
real (qxknd), parameter:: dber = 1.4e0_qxknd
integer i1, i2, k, n1
real (qxknd) d1
real (qxknd) t1, t2, t3

!  End of declaration

!   Check if berne array has been initialized.

d1 = berne(1)
if (abs (d1 - 1.e0_qxknd / 6.e0_qxknd) > qxrdfz .or. nb2 < int (dber * qxdpw)) then
  write (qxldb, 3) int (dber * qxdpw)
3 format ('*** QXPOLYGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries using QXBERNE or QXBERNER.')
  call qxabrt
endif

if (nn <= 0) then
  write (qxldb, 4)
4 format ('*** QXPOLYGAMMABE: NN <= 0')
  call qxabrt
endif

if (qxsgn (x) < 0) then
  write (qxldb, 5)
5 format ('*** QXPOLYGAMMABE: X < 0')
  call qxabrt
endif

t1 = 1.e0_qxknd

do k = 1, nn
  t1 = t1 * real (k, qxknd)
enddo

if (mod (nn + 1, 2) == 1) t1 = - t1
call qxhurwitzzetanbe (nb2, berne, nn + 1, x, t2)
y = t1 * t2

return
end subroutine qxpolygammabe

subroutine qxpolylogini (nn, arr)

!   Initializes the MP array arr with data for mppolylogneg.
!   NN must be in the range (-nmax, -1).

implicit none
integer, intent(in):: nn
real (qxknd), intent(out):: arr(1:abs(nn))
integer, parameter:: nmax = 1000
integer i1, i2, k, n, nna
real (qxknd) aa(2,abs(nn)), t1, t2

!  End of declaration

nna = abs (nn)
i1 = 2
i2 = 1
aa(1,1) = 1.e0_qxknd
aa(2,1) = 1.e0_qxknd

do k = 2, nna
  aa(1,k) = 0.e0_qxknd
  aa(2,k) = 0.e0_qxknd
enddo

do n = 2, nna
  i1 = 3 - i1
  i2 = 3 - i1

  do k = 2, n
    t1 = aa(i1,k-1) * real (n + 1 - k, qxknd)
    t2 = aa(i1,k) * real (k, qxknd)
    aa(i2,k) = t1 + t2
  enddo
enddo

do k = 1, nna
  arr(k) = aa(i2,k)
enddo

return
end subroutine qxpolylogini

subroutine qxpolylogneg (nn, arr, x, y)

!   This returns polylog (nn, x) for the case nn < 0. Before calling this,
!   one must call mppolylognini to initialize the array arr for this NN.
!   The dimensions of arr must be as shown below.
!   NN must be in the range (-nmax, -1).
!   The parameter nmxa is the maximum number of additional words of
!   precision needed to overcome cancelation errors when x is negative,
!   for nmax = 1000.

implicit none
integer, intent(in):: nn
real (qxknd), intent(in):: arr(1:abs(nn)), x
real (qxknd), intent(out):: y
integer, parameter:: nmax = 1000
integer i1, i2, k, n1, n2, nna
real (qxknd) d1, d2
real (qxknd) t1, t2, t3, t4

!  End of declaration

nna = abs (nn)
d1 = arr(1)
d2 = arr(nna)

if (d1 /= 1.e0_qxknd .or. d2 /= 1.e0_qxknd) then
  write (qxldb, 2)
2 format ('*** QXPOLYLOGNEG: Uninitialized or inadequately sized arrays'/ &
  'Call qxpolylogini or polylog_ini to initialize array. See documentation.')
  call qxabrt
endif

t1 = x
t2 = t1

do k = 2, nna
  t1 = x * t1
  t3 = arr(k) * t1
  t2 = t2 + t3
enddo

t3 = 1.e0_qxknd
t4 = t3 - x
t3 = t4 ** (nna + 1)
y = t2 / t3

return
end subroutine qxpolylogneg

subroutine qxpolylogpos (nn, x, y)

!   This returns polylog (nn, x) for the case nn >= 0.

implicit none
integer, intent(in):: nn
real (qxknd), intent(in):: x
real (qxknd), intent(out):: y
integer, parameter:: itrmax = 1000000
integer ic1, k
real (qxknd) t1, t2, t3, t4, t5, tc1, tc2, tc3, eps

!  End of declaration

if (nn < 0) then
  write (qxldb, 1)
1 format ('*** QXPOLYLOGPOS: N is less than zero.'/ &
  'For negative n, call qxpolylogneg or polylog_neg. See documentation.')
  call qxabrt
endif

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
t1 = abs (x)
t2 = 1.e0_qxknd
ic1 = qxsgn (t1 - t2)

if (ic1 >= 0) then
  write (qxldb, 3)
3 format ('*** QXPOLYLOGPOS: |X| must be less than one.')
  call qxabrt
endif

if (nn == 0) then
  t1 = 1.e0_qxknd
  t2 = t1 - x
  y = x / t2
else
  t1 = x
  t2 = x

  do k = 2, itrmax
    t2 = x * t2
    t3 = real (k, qxknd)
    t4 = t3 ** nn
    t3 = t2 / t4
    t1 = t1 + t3
    tc1 = abs (t3) - eps * abs (t1)
    if (qxsgn (tc1) < 0) goto 100
  enddo

  write (qxldb, 4)
4 format ('*** QXPOLYLOGPOS: Loop end error')
  call qxabrt

100 continue

  y = t1
endif

return
end subroutine qxpolylogpos

subroutine qxstruvehn (nu, ss, zz)

!   This returns the StruveH function with integer arg NU and MPFR argument SS.

implicit none
integer, intent(in):: nu
real (qxknd), intent(in):: ss
real (qxknd), intent(out):: zz
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dmax = 1000.e0_qxknd
integer ic1, k, n1
real (qxknd) d1
real (qxknd) sum, td1, td2, tn1, tnm1, t1, t2, tc1, tc2, tc3, eps

!  End of declaration

if (nu < 0) then
  write (qxldb, 2)
2 format ('*** QXSTRUVEHN: NU < 0')
  call qxabrt
endif

d1 = abs (ss)
if (d1 > dmax) then
  write (qxldb, 3) dmax
3 format ('*** QXSTRUVEHN: ABS(SS) >',f8.2)
  call qxabrt
endif

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
tn1 = 1.e0_qxknd
t1 = ss**2
tnm1 = t1 * (-0.25e0_qxknd)
td1 = sqrt (qxpicon) * 0.5e0_qxknd
td2 = td1

do k = 1, nu
  td2 = td2 * (real (k, qxknd) + 0.5e0_qxknd)
enddo

sum = tn1 / (td1 * td2)

do k = 1, itrmax
  tn1 = tn1 * tnm1
  td1 = td1 * (real (k, qxknd) + 0.5e0_qxknd)

  td2 = td2 * (real (nu + k, qxknd) + 0.5e0_qxknd)
  t2 = td1 * td2
  t1 = tn1 / t2
  sum = sum + t1
  tc1 = abs (t1) - eps * abs (sum)
  if (qxsgn (tc1) < 0) goto 100
enddo

write (qxldb, 5)
5 format ('*** QXSTRUVEHN: Loop end error')
call qxabrt

100 continue

t1 = ss * 0.5e0_qxknd
t2 = t1 ** (nu + 1)
zz = t2 * sum
 
return
end subroutine qxstruvehn

subroutine qxzetar (ss, zz)

!   This returns the zeta function of an MPR argument SS using an algorithm
!   due to Peter Borwein.

implicit none
real (qxknd), intent(in):: ss
real (qxknd), intent(out):: zz
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dfrac = 1.e0_qxknd+qxdpw
integer i, ic1, iss, j, n, n1, n2
real (qxknd) d1, d2
real (qxknd) f1, s, t1, t2, t3, t4, t5, tn, tt, tc1, tc2, tc3, eps
real (qxknd) sgn

!  End of declaration

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
f1 = 1.e0_qxknd
ic1 = qxsgn (ss - f1)
t2 = ss - anint (ss)

if (ic1 == 0) then
  write (qxldb, 2)
2 format ('*** QXZETAR: argument is 1')
  call qxabrt
elseif (qxsgn (t2) == 0) then

!   The argument is an integer value. Call mpzetaintr instead.

  iss = ss
  call qxzetaintr (iss, t1)
  goto 200
elseif (qxsgn (ss) < 0) then

!   If arg < 0, compute zeta(1-ss), and later apply Riemann's formula.

  tt = f1 - ss
else
  tt = ss
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = qxnbt * log (2.e0_qxknd) / log (2.e0_qxknd * qxnbt / 3.e0_qxknd)
d2 = tt

if (d2 > d1) then

!   Evaluate the infinite series.

  t1 = 1.e0_qxknd

  do i = 2, itrmax
    t4 = real (i, qxknd)
    t2 = t4 ** tt
    t3 = f1 / t2
    t1 = t1 + t3
    tc1 = abs (t3) - eps * abs (t1)
    if (qxsgn (tc1) < 0) goto 200
  enddo

  write (qxldb, 3) 1, itrmax
3 format ('*** QXZETAR: iteration limit exceeded',2i10)
  call qxabrt
endif

n = dfrac
t1 = 2.e0_qxknd
tn = t1 ** n
t1 = - tn
t2 = 0.e0_qxknd
s = 0.e0_qxknd
sgn = 1.e0_qxknd

do j = 0, 2 * n - 1
  t4 = real (j + 1, qxknd)
  t3 = t4 ** tt
  t4 = t1 / t3
  t5 = t4 * sgn
  t4 = s + t5
  s = t4
  sgn = - sgn

  if (j < n - 1) then
    t2 = 0.e0_qxknd
  elseif (j == n - 1) then
    t2 = 1.e0_qxknd
  else
    t3 = t2 * real (2 * n - j, qxknd)
    t2 = t3 / real (j + 1 - n, qxknd)
  endif
  t1 = t1 + t2
enddo

t3 = 1 - tt
t2 = 2.e0_qxknd
t4 = t2 ** t3
t2 = f1 - t4
t3 = tn * t2
t1 = - s / t3

!   If original argument was negative, apply Riemann's formula.

if (qxsgn (ss) < 0) then
  t3 = gamma (tt)
  t2 = t1 * t3
  t1 = qxpicon * tt
  t3 = t1 * 0.5e0_qxknd
  t4 = cos (t3)
  t5 = sin (t3)
  t1 = t2 * t4
  t2 = qxpicon * 2.e0_qxknd
  t3 = t2 ** tt
  t1 = t1 / t3 * 2.e0_qxknd
endif

200 continue

zz = t1

return
end subroutine qxzetar

subroutine qxzetaintr (iss, zz)

!   This returns the zeta function of the integer argument ISS using an algorithm
!   due to Peter Borwein.

implicit none
integer, intent(in):: iss
real (qxknd), intent(out):: zz
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dfrac = 1.e0_qxknd+qxdpw, pi = 3.1415926535897932385e0_qxknd
integer i, ic1, j, n, n1, itt
real (qxknd) d1, sgn
real (qxknd) f1, s, t1, t2, t3, t4, t5, tn, tc1, tc2, tc3, eps

!  End of declaration

eps = 2.e0_qxknd ** (-qxnwx*qxnbt)
f1 = 1.e0_qxknd

if (iss == 1) then
  write (qxldb, 2)
2 format ('*** QXZETAINTR: argument is 1')
  call qxabrt
elseif (iss == 0) then

!   Argument is zero -- result is -1/2.

  t1 = -0.5e0_qxknd
  goto 200
elseif (iss < 0) then

!   If argument is a negative even integer, the result is zero.

  if (mod (iss, 2) == 0) then
    t1 = 0.e0_qxknd
    goto 200
  endif

!   Otherwise if arg < 0, compute zeta(1-is), and later apply Riemann's formula.

  itt = 1 - iss
else
  itt = iss
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = qxnbt * log (2.e0_qxknd) / log (2.e0_qxknd * qxnbt / 3.e0_qxknd)

if (itt > d1) then

!   Evaluate the infinite series.

  t1 = 1.e0_qxknd

  do i = 2, itrmax
    t4 = real (i, qxknd)
    t2 = t4 ** itt
    t3 = f1 / t2
    t1 = t1 + t3
    tc1 = abs (t3) - eps * abs (t1)
    if (qxsgn (tc1) < 0) goto 200
  enddo

  write (qxldb, 3) 1, itrmax
3 format ('*** QXZETAINTR: iteration limit exceeded',2i10)
  call qxabrt
endif

n = dfrac
t1 = 2.e0_qxknd
tn = t1 ** n
t1 = - tn
t2 = 0.e0_qxknd
s = 0.e0_qxknd
sgn = 1.e0_qxknd

do j = 0, 2 * n - 1
  t4 = real (j + 1, qxknd)
  t3 = t4 ** itt
  t4 = t1 / t3
  t5 = t4 * sgn
  s = s + t5
  sgn = - sgn

  if (j < n - 1) then
    t2 = 0.e0_qxknd
  elseif (j == n - 1) then
    t2 = 1.e0_qxknd
  else
    t3 = t2 * real (2 * n - j, qxknd)
    t2 = t3 / real (j + 1 - n, qxknd)
  endif

  t1 = t1 + t2
enddo

t2 = 2.e0_qxknd
t4 = t2 ** (1 - itt)
t2 = f1 - t4
t3 = tn * t2
t1 = - s / t3

!   If original argument was negative, apply Riemann's formula.

if (iss < 0) then
  t3 = 1.e0_qxknd

  do i = 1, itt - 1
    t3 = t3 * real (i, qxknd)
  enddo

  t2 = t1 * t3
  t1 = qxpicon * real (itt, qxknd)
  t3 = t1 * 0.5e0_qxknd
  t4 = cos (t3)
  t5 = sin (t3)
  t1 = t2 * t4
  t2 = qxpicon * 2.e0_qxknd
  t3 = t2 ** itt
  t1 = t1 / t3 * 2.e0_qxknd
endif

200 continue

zz = t1

return
end subroutine qxzetaintr

subroutine qxzetabe (nb2, berne, s, z)

!  This evaluates the Riemann zeta function, using the combination of
!  the definition formula (for large s), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb2
real (qxknd), intent(in):: berne(nb2), s
real (qxknd), intent(out):: z
integer, parameter:: itrmax = 1000000
real (qxknd), parameter:: dber = 1.5e0_qxknd, dfrac = 0.6e0_qxknd
integer i, i1, i2, ic1, k, n1, n2, nn
real (qxknd) d1, d2
real (qxknd) t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, tt, f1, tc1, tc2, tc3, eps

!  End of declaration

!   Check if berne array has been initialized.

d1 = berne(1)
if (abs (d1 - 1.e0_qxknd / 6.e0_qxknd) > qxrdfz .or. nb2 < int (dber * qxdpw)) then
  write (qxldb, 3) int (dber * qxdpw)
3 format ('*** QXZETABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries.')
  call qxabrt
endif

i = 0
k = 0
eps = 2.e0_qxknd ** (-qxnwx*qxnbt)

!   Check if argument is 1 -- undefined.

t0 = 1.e0_qxknd
ic1 = qxsgn (s - t0)

if (ic1 == 0) then
  write (qxldb, 2)
2 format ('*** QXZETABE: argument is 1')
  call qxabrt
endif

f1 = 1.e0_qxknd

!   Check if argument is zero. If so, result is - 1/2.

if (qxsgn (s) == 0) then
  t1 = -0.5e0_qxknd
  goto 200
endif

!   Check if argument is negative.

if (qxsgn (s) < 0) then

!   Check if argument is a negative even integer. If so, the result is zero.

  t1 = s * 0.5e0_qxknd
  t3 = t1 - anint (t1)
  if (qxsgn (t3) == 0) then
    t1 = 0.e0_qxknd
    goto 200
  endif

!   Otherwise compute zeta(1-s), and later apply the reflection formula.

  tt = f1 - s
else
  tt = s
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = qxlogb / log (32.e0_qxknd)
d2 = tt

if (d2 > d1) then
  t1 = 1.e0_qxknd

  do i = 2, itrmax
    t4 = real (i, qxknd)
    t2 = t4 ** tt
    t3 = f1 / t2
    t1 = t1 + t3
    tc1 = abs (t3) - eps * abs (t1)
    if (qxsgn (tc1) < 0) goto 200
  enddo

  write (qxldb, 4) 1, itrmax
4 format ('*** QXZETABE: iteration limit exceeded',2i10)
  call qxabrt
endif

t0 = 1.e0_qxknd
nn = dfrac * qxdpw

do k = 2, nn
 t2 = real (k, qxknd)
 t1 = t2 ** tt
 t2 = f1 / t1
 t0 = t0 + t2
enddo

t2 = real (nn, qxknd)
t3 = tt - f1
t4 = t1 * t3
t3 = t2 / t4
t2 = t0 + t3
t3 = 0.5e0_qxknd
t4 = t3 / t1
t0 = t2 - t4

t3 = tt
d1 = 12.e0_qxknd * real (nn, qxknd)
t4 = t1 * d1
t2 = t3 / t4
t5 = t1 * real (nn, qxknd)
t6 = real (nn, qxknd)
t9 = t6**2

do k = 2, min (nb2, itrmax)
  t4 = real (2*k - 2, qxknd)
  t6 = tt + t4
  t7 = real (2*k - 3, qxknd)
  t8 = tt + t7
  t7 = t6 * t8
  t4 = t3 * t7
  t6 = real (2*k - 1, qxknd)
  t7 = real (2*k - 2, qxknd)
  t8 = t6 * t7
  t3 = t4 / t8
  t5 = t5 * t9
  t4 = t3 * berne(k)
  t6 = t5 * real (2*k, qxknd)
  t7 = t4 / t6
  t2 = t2 + t7
  tc2 = abs (t7) - eps * abs (t2)
  if (qxsgn (tc2) < 0) goto 110
enddo

write (qxldb, 4) 2, min (nb2, itrmax)
call qxabrt

110 continue

t1 = t0 + t2

!   If original argument was negative, apply the reflection formula.

if (qxsgn (s) < 0) then
  t3 = gamma (tt)
  t2 = t1 * t3
  t1 = qxpicon * tt
  t3 = t1 * 0.5e0_qxknd
  t4 = cos (t3)
  t5 = sin (t3)
  t1 = t2 * t4
  t2 = qxpicon * 2.e0_qxknd
  t3 = t2 ** t3
  t2 = t1 / t3
  t1 = t2 * 2.e0_qxknd
endif

200 continue

z = t1

return
end subroutine qxzetabe

end module qxfune
