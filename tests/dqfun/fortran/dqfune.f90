!*****************************************************************************

!  DQFUN: A double-quad precision package with special functions
!  Special functions module (module DQFUNE)

!  Revision date:  27 Feb 2023

!  AUTHOR:
!    David H. Bailey
!    Lawrence Berkeley National Lab (retired)
!    Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs. All basic arithmetic
!    operations and transcendental functions are supported, together with numerous
!    special functions.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:

!    David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package,"
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2020.pdf.

!  DESCRIPTION OF THIS MODULE (DQFUNE):
!    This module contains subroutines to perform special functions. Additional
!    functions will be added as they are completed.

!  NOTE ON PROGRAMMING CONVENTION FOR THIS MODULE:
!    This module is designed to facilitate easy translation (using the program
!    convmpfr.f90), for use in the MPFUN-MPFR package, and also for translation
!    (using the program convmpdq.f90), for use in the DQFUN double-quad package.

module dqfune
use dqfuna

contains

!  These routines perform simple operations on the MP data structure. Those
!  listed between !> and !>> are for MPFUN20-Fort; those between !>> and !>>>
!  are for MPFUN20-MPFR. Those between !>>> and !>>>> are for DQFUN. The
!  translator selects the proper set.
!>
! 
! subroutine mpinitwds (ra, mpnw)
! implicit none
! integer (mpiknd), intent(out):: ra(0:)
! integer, intent(in):: mpnw
! 
! ra(0) = mpnw + 6
! ra(1) = mpnw
! ra(2) = 0; ra(3) = 0; ra(4) = 0
! return
! end subroutine mpinitwds
! 
! integer function mpwprecr (ra)
! implicit none
! integer (mpiknd), intent(in):: ra(0:)
! 
! mpwprecr = ra(1)
! return
! end function mpwprecr
! 
! integer function mpspacer (ra)
! implicit none
! integer (mpiknd), intent(in):: ra(0:)
! 
! mpspacer = ra(0)
! return
! end function mpspacer
! 
!>>

! subroutine mpabrt (ier)
! implicit none
! integer, intent(in):: ier
! write (mpldb, 1) ier
! 1 format ('*** MPABRT: Execution terminated, error code =',i4)
! stop
! end subroutine mpabrt

! subroutine mpinitwds (ra, mpnw)
! implicit none
! integer (mpiknd), intent(out):: ra(0:)
! integer, intent(in):: mpnw
!
! ra(0) = mpnw + 6
! ra(1) = mpnw * mpnbt
! ra(2) = 1
! ra(3) = mpnan
! ra(4) = loc (ra(4)) + 8
! ra(mpnw+5) = 0
! return
! end subroutine mpinitwds

! subroutine mpfixlocr (ia)
! implicit none
! integer (mpiknd), intent(out):: ia(0:)
!
! ia(4) = loc (ia(4)) + 8
! return
! end subroutine

! integer function mpwprecr (ra)
! implicit none
! integer (mpiknd), intent(in):: ra(0:)
!
! mpwprecr = ra(1) / mpnbt
! return
! end function mpwprecr

! integer function mpspacer (ra)
! implicit none
! integer (mpiknd), intent(in):: ra(0:)
!
! mpspacer = ra(0)
! return
! end function mpspacer

!>>>

integer function dqwprecr (ra)
implicit none
real (dqknd), intent(in):: ra(2)

dqwprecr = 2
return
end function dqwprecr

integer function dqspacer (ra)
implicit none
real (dqknd), intent(in):: ra(2)

dqspacer = 2
return
end function dqspacer

!>>>>
! integer function ddwprecr (ra)
! implicit none
! real (ddknd), intent(in):: ra(2)
!
! ddwprecr = 2
! return
! end function ddwprecr

! integer function ddspacer (ra)
! implicit none
! real (ddknd), intent(in):: ra(2)
!
! ddspacer = 2
! return
! end function ddspacer
!>>>>>

subroutine dqberner (nb2, berne)

!   This returns the array berne, containing Bernoulli numbers indexed 2*k for
!   k = 1 to n, to mpnw words precision. This is done by first computing
!   zeta(2*k), based on the following known formulas:

!   coth (pi*x) = cosh (pi*x) / sinh (pi*x)

!            1      1 + (pi*x)^2/2! + (pi*x)^4/4! + ...
!        =  ---- * -------------------------------------
!           pi*x    1 + (pi*x)^2/3! + (pi*x)^4/5! + ...

!        = 1/(pi*x) * (1 + (pi*x)^2/3 - (pi*x)^4/45 + 2*(pi*x)^6/945 - ...)

!        = 2/(pi*x) * Sum_{k >= 1} (-1)^(k+1) * zeta(2*k) * x^{2*k}

!   The strategy is to calculate the coefficients of the series by polynomial
!   operations. Polynomial division is performed by computing the reciprocal
!   of the denominator polynomial, by a polynomial Newton iteration, as follows.
!   Let N(x) be the polynomial approximation to the numerator series; let D(x) be
!   a polynomial approximation to the numerator numerator series; and let Q_k(x)
!   be polynomial approximations to R(x) = 1/D(x). Then iterate:

!   Q_{k+1} = Q_k(x) + [1 - D(x)*Q_k(x)]*Q_k(x)

!   In these iterations, both the degree of the polynomial Q_k(x) and the
!   precision level in words are initially set to 4. When convergence is
!   achieved at this precision level, the degree is doubled, and iterations are
!   continued, etc., until the final desired degree is achieved. Then the
!   precision level is doubled and iterations are performed in a similar way,
!   until the final desired precision level is achieved. The reciprocal polynomial
!   R(x) produced by this process is then multiplied by the numerator polynomial
!   N(x) to yield an approximation to the quotient series. The even zeta values
!   are then the coefficients of this series, scaled according to the formula above.

!   Once the even integer zeta values have been computed in this way, the even
!   Bernoulli numbers are computed via the formula (for n > 0):

!   B(2*n) = (-1)^(n-1) * 2 * (2*n)! * zeta(2*n) / (2*pi)^(2*n)

!   Note: The notation in the description above is not the same as in the code below.

implicit none
integer, intent(in):: nb2
real (dqknd), intent(out):: berne(1:2,nb2)
integer, parameter:: ibz = 6, idb = 0
real (dqknd), parameter:: alog102 = 0.30102999566398119q0, pi = 3.1415926535897932385q0
integer i, i1, ic1, j, kn, n, n1, nn1
real (dqknd) d1, dd1, dd2, dd3
real (dqknd) c1(1:2,0:nb2), cp2(1:2), p1(1:2,0:nb2), &
  p2(1:2,0:nb2), q(1:2,0:nb2), q1(1:2), &
  r(1:2,0:nb2), s(1:2,0:nb2), t1(1:2), t2(1:2), &
  t3(1:2), t4(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

n = nb2
dqnw1 = min (dqnw + 1, dqnwx)
dqnw2 = dqnw1

if (idb > 0) write (dqldb, 3) n
3 format ('Even Bernoulli number calculation; n =',2i6)


do i = 0, nb2
enddo

call dqmul (dqpicon, dqpicon, cp2)
call dqdmc (1.q0, 0, c1(1:2,0))
call dqdmc (1.q0, 0, p1(1:2,0))
call dqdmc (1.q0, 0, p2(1:2,0))
call dqdmc (1.q0, 0, q(1:2,0))

!   Construct numerator and denominator polynomials.

do i = 1, n
  call dqdmc (0.q0, 0, c1(1:2,i))
  dd1 = 2.q0 * (i + 1) - 3.q0
  dd2 = dd1 + 1.q0
  dd3 = dd2 + 1.q0
  call dqmul (cp2, p1(1:2,i-1), t1)
  call dqdivd (t1, dd1 * dd2, p1(1:2,i))
  call dqmul (cp2, p2(1:2,i-1), t1)
  call dqdivd (t1, dd2 * dd3, p2(1:2,i))
  call dqdmc (0.q0, 0, q(1:2,i))
enddo

kn = 4
dqnw2 = min (4, dqnwx)

call dqdmc (2.q0, 0, t1)
call dqnpwr (t1, ibz - dqnw2 * dqnbt, eps)
if (idb > 0) then
  call dqmdc (eps, dd1, nn1)
  write (dqldb, 4) dqnw2, nint (alog102*nn1)
4 format ('mpnw2, log10eps =',2i6)
endif
call dqdmc (0.q0, 0, q1)

!   Perform Newton iterations with dynamic precision levels, using an
!   iteration formula similar to that used to evaluate reciprocals.

do j = 1, 10000
  if (idb > 0) write (dqldb, 5) j, kn
5 format ('j, kn =',3i6)

  call dqpolymul (kn, p2, q, r)
  call dqpolysub (kn, c1, r, s)
  call dqpolymul (kn, s, q, r)
  call dqpolyadd (kn, q, r, q)
  call dqsub (q(1:2,kn), q1, t1)

  if (idb > 0) then
    call dqmdc (t1, dd1, nn1)
    if (dd1 == 0.q0) then
      write (dqldb, 6)
6     format ('Newton error = 0')
    else
      write (dqldb, 7) nint (alog102*nn1)
7     format ('Newton error = 10^',i6)
    endif
  endif

  call dqabs (t1, t2)
  call dqcpr (t2, eps, ic1)
  if (ic1 < 0) then
    if (kn == n .and. dqnw2 == dqnw1) goto 100
    if (kn < n) then
      kn = min (2 * kn, n)
      call dqdmc (0.q0, 0, q1)
    elseif (dqnw2 < dqnw1) then
      dqnw2 = min (2 * dqnw2, dqnw1)
      call dqdmc (2.q0, 0, t1)
      call dqnpwr (t1, ibz - dqnw2 * dqnbt, eps)
      call dqdmc (0.q0, 0, q1)
      if (idb > 0) then
        call dqmdc (eps, dd1, nn1)
        write (dqldb, 4) dqnw2, nint (alog102*nn1)
      endif
    endif
  else
    call dqeq (q(1:2,kn), q1)
  endif
enddo

write (dqldb, 8)
8 format ('*** DQBERNER: Loop end error')
call dqabrt

100 continue

if (idb > 0) write (dqldb, 9)
9 format ('Even zeta computation complete')

!   Multiply numerator polynomial by reciprocal of denominator polynomial.

call dqpolymul (n, p1, q, r)

!   Apply formula to produce Bernoulli numbers.

call dqdmc (-2.q0, 0, t1)
call dqdmc (1.q0, 0, t2)

do i = 1, n
  d1 = - real (2*i-1, dqknd) * real (2*i, dqknd)
  call dqmuld (t1, d1, t3)
  call dqeq (t3, t1)
  call dqmuld (cp2, 4.q0, t3)
  call dqmul (t3, t2, t4)
  call dqeq (t4, t2)
  call dqmuld (t1, 0.5q0, t3)
  call dqdiv (t3, t2, t4)
  call dqabs (r(1:2,i), t3)
  call dqmul (t4, t3, berne(1:2,i))
enddo

if (idb > 0) write (dqldb, 10)
10 format ('Bernoulli number computation complete')
return
end subroutine dqberner

subroutine dqpolyadd (n, a, b, c)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: n
real (dqknd), intent(in):: a(1:2,0:n), b(1:2,0:n)
real (dqknd), intent(out):: c(1:2,0:n)
integer k
real (dqknd) t1(1:2), t2(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx


do k = 0, n
  call dqeq (a(1:2,k), t1)
  call dqeq (b(1:2,k), t2)
  call dqadd (t1, t2, c(1:2,k))
enddo

return
end subroutine dqpolyadd

subroutine dqpolysub (n, a, b, c)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: n
real (dqknd), intent(in):: a(1:2,0:n), b(1:2,0:n)
real (dqknd), intent(out):: c(1:2,0:n)
integer k
real (dqknd) t1(1:2), t2(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx


do k = 0, n
  call dqeq (a(1:2,k), t1)
  call dqeq (b(1:2,k), t2)
  call dqsub (t1, t2, c(1:2,k))
enddo

return
end subroutine dqpolysub

subroutine dqpolymul (n, a, b, c)

!   This adds two polynomials (ignoring high-order terms), as is required
!   by mpberne. The output array C may not be the same as A or B.

implicit none
integer, intent(in):: n
real (dqknd), intent(in):: a(1:2,0:n), b(1:2,0:n)
real (dqknd), intent(out):: c(1:2,0:n)
integer j, k
real (dqknd) t0(1:2), t1(1:2), t2(1:2), t3(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx


do k = 0, n
  call dqdmc (0.q0, 0, t0)

  do j = 0, k
    call dqeq (a(1:2,j), t1)
    call dqeq (b(1:2,k-j), t2)
    call dqmul (t1, t2, t3)
    call dqadd (t0, t3, t2)
    call dqeq (t2, t0)
  enddo

  call dqeq (t0, c(1:2,k))
enddo

return
end subroutine dqpolymul

subroutine dqbesselinr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselI (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.25.2 for modest RR,
!   and DLMF 10.40.1 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (dqknd), intent(in):: rr(1:2)
real (dqknd), intent(out):: ss(1:2)
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dfrac = 1.5q0, pi = 3.1415926535897932385q0
integer ic1, k, nua, n1
real (dqknd) d1
real (dqknd) f1(1:2), f2(1:2), sum(1:2), td(1:2), &
  tn(1:2), t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  rra(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

!   Check for RR = 0.

if (dqsgn (rr) == 0) then
  write (dqldb, 2)
2 format ('*** DQBESSELINR: Second argument is zero')
  call dqabrt
endif

dqnw1 = min (dqnw + 1, dqnwx)
nua = abs (nu)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)
call dqabs (rr, rra)
call dqmdc (rra, d1, n1)
d1 = d1 * 2.q0 ** n1

if (d1 < dfrac * dqnw1 * dqdpw) then
  call dqdmc (1.q0, 0, tn)
  call dqdmc (1.q0, 0, f1)
  call dqdmc (1.q0, 0, f2)
  call dqmul (rra, rra, t2)
  call dqmuld (t2, 0.25q0, t1)

  do k = 1, nua
    call dqmuld (f2, real (k, dqknd), t2)
    call dqeq (t2, f2)
  enddo

  call dqmul (f1, f2, td)
  call dqdiv (tn, td, t2)
  call dqeq (t2, sum)

  do k = 1, itrmax
    call dqmuld (f1, real (k, dqknd), t2)
    call dqeq (t2, f1)
    call dqmuld (f2, real (k + nua, dqknd), t2)
    call dqeq (t2, f2)
    call dqmul (t1, tn, t2)
    call dqeq (t2, tn)
    call dqmul (f1, f2, td)
    call dqdiv (tn, td, t2)
    call dqadd (sum, t2, t3)
    call dqeq (t3, sum)

    call dqabs (t2, tc1)
    call dqmul (eps, sum, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (dqldb, 4)
  4 format ('*** DQBESSELINR: Loop end error 1')
  call dqabrt

100 continue

  call dqmuld (rra, 0.5q0, t1)
  call dqnpwr (t1, nua, t2)
  call dqmul (sum, t2, t3)
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

  call dqdmc (1.q0, 0, sum)
  d1 = 4.q0 * real (nua, dqknd)**2
  call dqdmc (d1, 0, t1)
  call dqdmc (1.q0, 0, tn)
  call dqdmc (1.q0, 0, td)

  do k = 1, itrmax
!  t2 = -t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.q0 * k - 1.q0
    call dqdmc (d1, 0, t2)
    call dqmul (t2, t2, t3)
    call dqsub (t1, t3, t2)
    call dqmul (tn, t2, t3)
    call dqneg (t3, tn)
    call dqmuld (rra, 8.q0 * k, t2)
    call dqmul (td, t2, t3)
    call dqeq (t3, td)
    call dqdiv (tn, td, t4)
    call dqadd (sum, t4, t3)
    call dqeq (t3, sum)

!   if (abs (t4) / abs (sum1) < eps) goto 110

    call dqabs (t4, tc1)
    call dqmul (eps, sum, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 110
  enddo

write (dqldb, 5)
5 format ('*** DQBESSELINR: Loop end error 2')
call dqabrt

110 continue

! t1 = exp (xa) / sqrt (2.d0 * mppi (mpnw) * xa)
! besseli = t1 * sum1

  call dqexp (rra, t1)
  call dqmuld (dqpicon, 2.q0, t2)
  call dqmul (t2, rra, t3)
  call dqsqrt (t3, t4)
  call dqdiv (t1, t4, t2)
  call dqmul (t2, sum, t3)
endif

! if (x < 0.d0 .and. mod (nu, 2) /= 0) besseli = - besseli

if (dqsgn (rr) < 0 .and. mod (nu, 2) /= 0) then
  call dqneg (t3, t4)
  call dqeq (t4, t3)
endif

call dqeq (t3, ss)

return
end subroutine dqbesselinr

subroutine dqbesselir (qq, rr, ss)

!   This evaluates the modified Bessel function BesselI (QQ,RR) for QQ and RR
!   both MPR. The algorithm is DLMF formula 10.25.2 for modest RR, and
!   DLMF 10.40.1 for large RR, relative to precision.

implicit none
real (dqknd), intent(in):: qq(1:2), rr(1:2)
real (dqknd), intent(out):: ss(1:2)
integer ic1, i1, i2, k, n1
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dfrac = 1.5q0, pi = 3.1415926535897932385q0
real (dqknd) d1
real (dqknd) f1(1:2), f2(1:2), sum(1:2), td(1:2), &
  tn(1:2), t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  rra(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

!   If QQ is integer, call mpbesselinr; if qq < 0 and rr <= 0, then error.

call dqinfr (qq, t1, t2)
i1 = dqsgn (qq)
i2 = dqsgn (rr)
if (dqsgn (t2) == 0) then
  call dqmdc (qq, d1, n1)
  d1 = d1 * 2.q0**n1
  n1 = nint (d1)
  call dqbesselinr (n1, rr, t3)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (dqldb, 3)
3 format ('*** DQBESSELIR: First argument < 0 and second argument <= 0')
  call dqabrt
endif

call dqabs (rr, rra)
call dqmdc (rra, d1, n1)
d1 = d1 * 2.q0 ** n1

if (d1 < dfrac * dqnw1 * dqdpw) then
  call dqdmc (1.q0, 0, tn)
  call dqdmc (1.q0, 0, f1)
  call dqadd (qq, f1, t1)
  call dqgammar (t1, f2)
  call dqmul (rra, rra, t2)
  call dqmuld (t2, 0.25q0, t1)

  call dqmul (f1, f2, td)
  call dqdiv (tn, td, t2)
  call dqeq (t2, sum)

  do k = 1, itrmax
    call dqmuld (f1, real (k, dqknd), t2)
    call dqeq (t2, f1)
    call dqdmc (real (k, dqknd), 0, t3)
    call dqadd (qq, t3, t4)
    call dqmul (f2, t4, t3)
    call dqeq (t3, f2)
    call dqmul (t1, tn, t2)
    call dqeq (t2, tn)
    call dqmul (f1, f2, td)
    call dqdiv (tn, td, t2)
    call dqadd (sum, t2, t3)
    call dqeq (t3, sum)

    call dqabs (t2, tc1)
    call dqmul (eps, sum, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (dqldb, 4)
  4 format ('*** DQBESSELIR: Loop end error 1')
  call dqabrt

100 continue

  call dqmuld (rra, 0.5q0, t1)
  call dqpower (t1, qq, t2)
  call dqmul (sum, t2, t3)
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

  call dqdmc (1.q0, 0, sum)
  call dqmul (qq, qq, t2)
  call dqmuld (t2, 4.q0, t1)
  call dqdmc (1.q0, 0, tn)
  call dqdmc (1.q0, 0, td)

  do k = 1, itrmax
!  t2 = -t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.q0 * k - 1.q0
    call dqdmc (d1, 0, t2)
    call dqmul (t2, t2, t3)
    call dqsub (t1, t3, t2)
    call dqmul (tn, t2, t3)
    call dqneg (t3, tn)
    call dqmuld (rra, 8.q0 * k, t2)
    call dqmul (td, t2, t3)
    call dqeq (t3, td)
    call dqdiv (tn, td, t4)
    call dqadd (sum, t4, t3)
    call dqeq (t3, sum)

!   if (abs (t4) / abs (sum1) < eps) goto 110

    call dqabs (t4, tc1)
    call dqmul (eps, sum, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 110
  enddo

write (dqldb, 5)
5 format ('*** DQBESSELIR: Loop end error 2')
call dqabrt

110 continue

! t1 = exp (xa) / sqrt (2.d0 * mppi (mpnw) * xa)
! besseli = t1 * sum1

  call dqexp (rra, t1)
  call dqmuld (dqpicon, 2.q0, t2)
  call dqmul (t2, rra, t3)
  call dqsqrt (t3, t4)
  call dqdiv (t1, t4, t2)
  call dqmul (t2, sum, t3)
endif

120 continue

call dqeq (t3, ss)

return
end subroutine dqbesselir

subroutine dqbesseljnr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselJ (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (dqknd), intent(in):: rr(1:2)
real (dqknd), intent(out):: ss(1:2)
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dfrac = 1.5q0, pi = 3.1415926535897932385q0
integer ic1, ic2, k, nua, n1
real (dqknd) d1, d2
real (dqknd) f1(1:2), f2(1:2), sum1(1:2), &
  sum2(1:2), td1(1:2), td2(1:2), tn1(1:2), &
  tn2(1:2), t1(1:2), t2(1:2), t3(1:2), &
  t41(1:2), t42(1:2), t5(1:2), rra(1:2), &
  rr2(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

!   Check for RR = 0.

if (dqsgn (rr) == 0) then
  write (dqldb, 2)
2 format ('*** DQBESSELJNR: Second argument is zero')
  call dqabrt
endif

dqnw1 = min (2 * dqnw + 1, dqnwx)
nua = abs (nu)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw*dqnbt, eps)
call dqmdc (rr, d1, n1)
d1 = abs (d1) * 2.q0 ** n1

if (d1 < dfrac * dqnw1 * dqdpw) then
  dqnw1 = min (dqnw + nint (d1 / (dfrac * dqdpw)), 2 * dqnw + 1, dqnwx)
  call dqabs (rr, rra)
  call dqdmc (1.q0, 0, tn1)
  call dqdmc (1.q0, 0, f1)
  call dqdmc (1.q0, 0, f2)
  call dqmul (rra, rra, t2)
  call dqmuld (t2, 0.25q0, t1)

  do k = 1, nua
    call dqmuld (f2, real (k, dqknd), t2)
    call dqeq (t2, f2)
  enddo

  call dqmul (f1, f2, td1)
  call dqdiv (tn1, td1, t2)
  call dqeq (t2, sum1)

  do k = 1, itrmax
    call dqmuld (f1, real (k, dqknd), t2)
    call dqeq (t2, f1)
    call dqmuld (f2, real (k + nua, dqknd), t2)
    call dqeq (t2, f2)
    call dqmul (t1, tn1, t2)
    call dqneg (t2, tn1)
    call dqmul (f1, f2, td1)
    call dqdiv (tn1, td1, t2)
    call dqadd (sum1, t2, t3)
    call dqeq (t3, sum1)

    call dqabs (t2, tc1)
    call dqmul (eps, sum1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (dqldb, 4)
4 format ('*** DQBESSELJNR: Loop end error 1')
  call dqabrt

100 continue

  call dqmuld (rra, 0.5q0, t1)
  call dqnpwr (t1, nua, t2)
  call dqmul (sum1, t2, t3)
else
! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  dqnw1 = min (dqnw + 1, dqnwx)
  call dqabs (rr, rra)
  call dqmul (rra, rra, rr2)
  d1 = 4.q0 * real (nua, dqknd)**2
  call dqdmc (d1, 0, t1)
  call dqdmc (1.q0, 0, tn1)
  call dqsub (t1, tn1, t2)
  call dqdivd (t2, 8.q0, tn2)
  call dqdmc (1.q0, 0, td1)
  call dqeq (rra, td2)
  call dqdiv (tn1, td1, sum1)
  call dqdiv (tn2, td2, sum2)

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.q0*(2.q0*k-1.q0) - 1.q0)**2
    d2 = (2.q0*(2.q0*k) - 1.q0)**2
    call dqdmc (d1, 0, t2)
    call dqsub (t1, t2, t3)
    call dqdmc (d2, 0, t2)
    call dqsub (t1, t2, t5)
    call dqmul (t3, t5, t2)
    call dqmul (tn1, t2, t3)
    call dqneg (t3, tn1)

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = real (2*k-1, dqknd) * real (2*k, dqknd) * 64.q0
    call dqmuld (td1, d1, t2)
    call dqmul (t2, rr2, td1)

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

    call dqdiv (tn1, td1, t41)
    call dqadd (sum1, t41, t2)
    call dqeq (t2, sum1)

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.q0*(2.q0*k) - 1.q0)**2
    d2 = (2.q0*(2.q0*k+1.q0) - 1.q0)**2
    call dqdmc (d1, 0, t2)
    call dqsub (t1, t2, t3)
    call dqdmc (d2, 0, t2)
    call dqsub (t1, t2, t5)
    call dqmul (t3, t5, t2)
    call dqmul (tn2, t2, t3)
    call dqneg (t3, tn2)

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = real (2*k, dqknd) * real (2*k+1, dqknd) * 64.q0
    call dqmuld (td2, d1, t2)
    call dqmul (t2, rr2, td2)

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

    call dqdiv (tn2, td2, t42)
    call dqadd (sum2, t42, t2)
    call dqeq (t2, sum2)
    
!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

    call dqabs (t41, tc1)
    call dqmul (eps, sum1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    call dqabs (t42, tc1)
    call dqmul (eps, sum2, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic2)
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (dqldb, 5)
5 format ('*** DQBESSELJNR: Loop end error 2')
  call dqabrt

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

  call dqmuld (dqpicon, 0.5q0 * nua, t1)
  call dqsub (rra, t1, t2)
  call dqmuld (dqpicon, 0.25q0, t1)
  call dqsub (t2, t1, t3)
  call dqcssnr (t3, t41, t42)
  call dqmul (t41, sum1, t1)
  call dqmul (t42, sum2, t2)
  call dqsub (t1, t2, t5)
  call dqmul (dqpicon, rra, t1)
  call dqdmc (2.q0, 0, t2)
  call dqdiv (t2, t1, t3)
  call dqsqrt (t3, t1)
  call dqmul (t1, t5, t3)
endif

if (mod (nu, 2) /= 0) then
!  if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) besselj = - besselj

  ic1 = dqsgn (rr)
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
    call dqneg (t3, t2)
    call dqeq (t2, t3)
  endif
endif

call dqeq (t3, ss)

return
end subroutine dqbesseljnr

subroutine dqbesseljr (qq, rr, ss)

!   This evaluates the modified Bessel function BesselJ (QQ,RR) for QQ and RR
!   both MPR. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
real (dqknd), intent(in):: qq(1:2), rr(1:2)
real (dqknd), intent(out):: ss(1:2)
integer ic1, ic2, i1, i2, k, n1
real (dqknd) d1, d2
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dfrac = 1.5q0, pi = 3.1415926535897932385q0
real (dqknd) f1(1:2), f2(1:2), sum1(1:2), &
  sum2(1:2), td1(1:2), td2(1:2), tn1(1:2), &
  tn2(1:2), t1(1:2), t2(1:2), t3(1:2), &
  t4(1:2), t41(1:2), t42(1:2), t5(1:2), &
  rra(1:2), rr2(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (2 * dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw*dqnbt, eps)

!   If QQ is integer, call mpbesseljnr; if RR <= 0, then error.

call dqinfr (qq, t1, t2)
i1 = dqsgn (qq)
i2 = dqsgn (rr)
if (dqsgn (t2) == 0) then
  call dqmdc (qq, d1, n1)
  d1 = d1 * 2.q0**n1
  n1 = nint (d1)
  call dqbesseljnr (n1, rr, t3)
  goto 120
elseif (i2 <= 0) then
  write (dqldb, 3)
3 format ('*** DQBESSELJR: Second argument <= 0')
  call dqabrt
endif

call dqmdc (rr, d1, n1)
d1 = abs (d1) * 2.q0 ** n1

if (d1 < dfrac * dqnw1 * dqdpw) then
  dqnw1 = min (dqnw + nint (d1 / (dfrac * dqdpw)), 2 * dqnw + 1, dqnwx)
  call dqabs (rr, rra)
  call dqdmc (1.q0, 0, tn1)
  call dqdmc (1.q0, 0, f1)
  call dqadd (qq, f1, t1)
  call dqgammar (t1, f2)
  call dqmul (rra, rra, t2)
  call dqmuld (t2, 0.25q0, t1)

  call dqmul (f1, f2, td1)
  call dqdiv (tn1, td1, t2)
  call dqeq (t2, sum1)

  do k = 1, itrmax
    call dqmuld (f1, real (k, dqknd), t2)
    call dqeq (t2, f1)
    call dqdmc (real (k, dqknd), 0, t3)
    call dqadd (qq, t3, t4)
    call dqmul (f2, t4, t3)
    call dqeq (t3, f2)
    call dqmul (t1, tn1, t2)
    call dqneg (t2, tn1)
    call dqmul (f1, f2, td1)
    call dqdiv (tn1, td1, t2)
    call dqadd (sum1, t2, t3)
    call dqeq (t3, sum1)

    call dqabs (t2, tc1)
    call dqmul (eps, sum1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (dqldb, 4)
4 format ('*** DQBESSELJR: Loop end error 1')
  call dqabrt

100 continue

  call dqmuld (rr, 0.5q0, t1)
  call dqpower (t1, qq, t2)
  call dqmul (sum1, t2, t3)
else
! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  dqnw1 = min (dqnw + 1, dqnwx)
  call dqabs (rr, rra)
  call dqmul (rra, rra, rr2)
  call dqmul (qq, qq, t2)
  call dqmuld (t2, 4.q0, t1)
  call dqdmc (1.q0, 0, tn1)
  call dqsub (t1, tn1, t2)
  call dqdivd (t2, 8.q0, tn2)
  call dqdmc (1.q0, 0, td1)
  call dqeq (rra, td2)
  call dqdiv (tn1, td1, sum1)
  call dqdiv (tn2, td2, sum2)

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.q0*(2.q0*k-1.q0) - 1.q0)**2
    d2 = (2.q0*(2.q0*k) - 1.q0)**2
    call dqdmc (d1, 0, t2)
    call dqsub (t1, t2, t3)
    call dqdmc (d2, 0, t2)
    call dqsub (t1, t2, t5)
    call dqmul (t3, t5, t2)
    call dqmul (tn1, t2, t3)
    call dqneg (t3, tn1)

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = real (2*k-1, dqknd) * real (2*k, dqknd) * 64.q0
    call dqmuld (td1, d1, t2)
    call dqmul (t2, rr2, td1)

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

    call dqdiv (tn1, td1, t41)
    call dqadd (sum1, t41, t2)
    call dqeq (t2, sum1)

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.q0*(2.q0*k) - 1.q0)**2
    d2 = (2.q0*(2.q0*k+1.q0) - 1.q0)**2
    call dqdmc (d1, 0, t2)
    call dqsub (t1, t2, t3)
    call dqdmc (d2, 0, t2)
    call dqsub (t1, t2, t5)
    call dqmul (t3, t5, t2)
    call dqmul (tn2, t2, t3)
    call dqneg (t3, tn2)

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = real (2*k, dqknd) * real (2*k+1, dqknd) * 64.q0
    call dqmuld (td2, d1, t2)
    call dqmul (t2, rr2, td2)

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

    call dqdiv (tn2, td2, t42)
    call dqadd (sum2, t42, t2)
    call dqeq (t2, sum2)

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

    call dqabs (t41, tc1)
    call dqmul (eps, sum1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    call dqabs (t42, tc1)
    call dqmul (eps, sum2, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic2)
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (dqldb, 5)
5 format ('*** DQBESSELJR: Loop end error 2')
  call dqabrt

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

!   call mpmuld (mppicon, 0.5d0 * nua, t1, mpnw1)

  call dqmul (dqpicon, qq, t2)
  call dqmuld (t2, 0.5q0, t1)
  call dqsub (rra, t1, t2)
  call dqmuld (dqpicon, 0.25q0, t1)
  call dqsub (t2, t1, t3)
  call dqcssnr (t3, t41, t42)
  call dqmul (t41, sum1, t1)
  call dqmul (t42, sum2, t2)
  call dqsub (t1, t2, t5)
  call dqmul (dqpicon, rra, t1)
  call dqdmc (2.q0, 0, t2)
  call dqdiv (t2, t1, t3)
  call dqsqrt (t3, t1)
  call dqmul (t1, t5, t3)
endif

! if (mod (nu, 2) /= 0) then
!  if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) besselj = - besselj

!   ic1 = mpsgn (rr)
!   if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
!     call mpneg (t3, t2, mpnw1)
!     call mpeq (t2, t3, mpnw1)
!   endif
! endif

120 continue

call dqeq (t3, ss)

return
end subroutine dqbesseljr

subroutine dqbesselknr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselK (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.31.1 for modest RR,
!   and DLMF 10.40.2 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (dqknd), intent(in):: rr(1:2)
real (dqknd), intent(out):: ss(1:2)
integer, parameter:: itrmax = 1000000
integer ic1, k, nua, n1
real (dqknd) d1
real (dqknd), parameter:: dfrac = 1.5q0, egam = 0.5772156649015328606q0, &
  pi = 3.1415926535897932385q0
real (dqknd) f1(1:2), f2(1:2), f3(1:2), f4(1:2), &
  f5(1:2), sum1(1:2), sum2(1:2), sum3(1:2), td(1:2), &
  tn(1:2), t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  rra(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

!   Check for RR = 0.

if (dqsgn (rr) == 0) then
  write (dqldb, 2)
2 format ('*** DQBESSELKNR: Second argument is zero')
  call dqabrt
endif

dqnw1 = min (dqnw + 1, dqnwx)
nua = abs (nu)
dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)
call dqabs (rr, rra)
call dqmdc (rra, d1, n1)
d1 = d1 * 2.q0 ** n1

if (d1 < dfrac * dqnw1 * dqdpw) then
  call dqmul (rra, rra, t2)
  call dqmuld (t2, 0.25q0, t1)
  call dqdmc (1.q0, 0, f1)
  call dqdmc (1.q0, 0, f2)
  call dqdmc (1.q0, 0, f3)
  call dqdmc (0.q0, 0,  sum1)

  do k = 1, nua - 1
    call dqmuld (f1, real (k, dqknd), t2)
    call dqeq (t2, f1)
  enddo

  do k = 0, nua - 1
    if (k > 0) then
      call dqdivd (f1, real (nua - k, dqknd), t2)
      call dqeq (t2, f1)
      call dqmul (t1, f2, t2)
      call dqneg (t2, f2)
      call dqmuld (f3, real (k, dqknd), t2)
      call dqeq (t2, f3)
    endif
    call dqmul (f1, f2, t3)
    call dqdiv (t3, f3, t2)
    call dqadd (sum1, t2, t3)
    call dqeq (t3, sum1)
  enddo

  call dqmuld (sum1, 0.5q0, t2)
  call dqmuld (rra, 0.5q0, t3)
  call dqnpwr (t3, nua, t4)
  call dqdiv (t2, t4, sum1)

  call dqmuld (rra, 0.5q0, t2)
  call dqlog (t2, t3)
  d1 = (-1.q0) ** (nua + 1)
  call dqmuld (t3, d1, t2)
  call dqbesselinr (nua, rra, t3)
  call dqmul (t2, t3, sum2)

  call dqneg (dqegammacon, f1)
  call dqeq (f1, f2)
  call dqdmc (1.q0, 0, f3)
  call dqdmc (1.q0, 0, f4)
  call dqdmc (1.q0, 0, f5)

  do k = 1, nua
    call dqdmc (1.q0, 0, t2)
    call dqdivd (t2, real (k, dqknd), t3)
    call dqadd (f2, t3, t4)
    call dqeq (t4, f2)
    call dqmuld (f5, real (k, dqknd), t2)
    call dqeq (t2, f5)
  enddo

  call dqadd (f1, f2, t2)
  call dqmul (t2, f3, t3)
  call dqmul (f4, f5, t4)
  call dqdiv (t3, t4, sum3)

  do k = 1, itrmax
    call dqdmc (1.q0, 0, t2)
    call dqdivd (t2, real (k, dqknd), t3)
    call dqadd (f1, t3, t4)
    call dqeq (t4, f1)
    call dqdivd (t2, real (nua + k, dqknd), t3)
    call dqadd (f2, t3, t4)
    call dqeq (t4, f2)
    call dqmul (t1, f3, t2)
    call dqeq (t2, f3)
    call dqmuld (f4, real (k, dqknd), t2)
    call dqeq (t2, f4)
    call dqmuld (f5, real (nua + k, dqknd), t2)
    call dqeq (t2, f5)
    call dqadd (f1, f2, t2)
    call dqmul (t2, f3, t3)
    call dqmul (f4, f5, t4)
    call dqdiv (t3, t4, t2)
    call dqadd (sum3, t2, t3)
    call dqeq (t3, sum3)

    call dqabs (t2, tc1)
    call dqmul (eps, sum3, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (dqldb, 5)
5 format ('*** DQBESSELKNR: Loop end error 1')
  call dqabrt

100 continue

  call dqmuld (rra, 0.5q0, t2)
  call dqnpwr (t2, nua, t3)
  d1 = (-1.q0)**nua * 0.5q0
  call dqmuld (t3, d1, t4)
  call dqmul (t4, sum3, t2)
  call dqeq (t2, sum3)
  call dqadd (sum1, sum2, t2)
  call dqadd (t2, sum3, t3)
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

  call dqdmc (1.q0, 0, sum1)
  d1 = 4.q0 * real (nua, dqknd)**2
  call dqdmc (d1, 0, t1)
  call dqdmc (1.q0, 0, tn)
  call dqdmc (1.q0, 0, td)

  do k = 1, itrmax
!  t2 = t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.q0 * k - 1.q0
    call dqdmc (d1, 0, t2)
    call dqmul (t2, t2, t3)
    call dqsub (t1, t3, t2)
    call dqmul (tn, t2, t3)
    call dqeq (t3, tn)
    call dqmuld (rra, 8.q0 * k, t2)
    call dqmul (td, t2, t3)
    call dqeq (t3, td)
    call dqdiv (tn, td, t4)
    call dqadd (sum1, t4, t3)
    call dqeq (t3, sum1)

!   if (abs (t4) / abs (sum1) < eps) goto 110

    call dqabs (t4, tc1)
    call dqmul (eps, sum1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 110
  enddo

write (dqldb, 6)
6 format ('*** DQBESSELKNR: Loop end error 2')
call dqabrt

110 continue

! t1 = sqrt (mppi (mpnw) / (2.d0 * xa)) / exp (xa)
! besseli = t1 * sum1

  call dqexp (rra, t1)
  call dqmuld (rra, 2.q0, t2)
  call dqdiv (dqpicon, t2, t3)
  call dqsqrt (t3, t4)
  call dqdiv (t4, t1, t2)
  call dqmul (t2, sum1, t3)
endif

! if (x < 0.d0 .and. mod (nu, 2) /= 0) besselk = - besselk

if (dqsgn (rr) < 0 .and. mod (nu, 2) /= 0) then
  call dqneg (t3, t4)
  call dqeq (t4, t3)
endif
call dqeq (t3, ss)
return
end subroutine dqbesselknr

subroutine dqbesselkr (qq, rr, ss)

!   This evaluates the Bessel function BesselK (QQ,RR) for QQ and RR
!   both MPR. This uses DLMF formula 10.27.4.

implicit none
real (dqknd), intent(in):: qq(1:2), rr(1:2)
real (dqknd), intent(out):: ss(1:2)
integer i1, i2, n1
real (dqknd) d1
real (dqknd) t1(1:2), t2(1:2), t3(1:2), t4(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)

!   If QQ is integer, call mpbesselknr; if qq < 0 and rr <= 0, then error.

call dqinfr (qq, t1, t2)
i1 = dqsgn (qq)
i2 = dqsgn (rr)
if (dqsgn (t2) == 0) then
  call dqmdc (qq, d1, n1)
  d1 = d1 * 2.q0**n1
  n1 = nint (d1)
  call dqbesselknr (n1, rr, t1)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (dqldb, 2)
2 format ('*** DQBESSELKR: First argument < 0 and second argument <= 0')
  call dqabrt
endif

call dqneg (qq, t1)
call dqbesselir (t1, rr, t2)
call dqbesselir (qq, rr, t3)
call dqsub (t2, t3, t4)
call dqmul (qq, dqpicon, t1)
call dqcssnr (t1, t2, t3)
call dqdiv (t4, t3, t2)
call dqmul (dqpicon, t2, t3)
call dqmuld (t3, 0.5q0, t1)

120 continue

call dqeq (t1, ss)
return
end subroutine dqbesselkr

subroutine dqbesselynr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselY (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.8.1 for modest RR,
!   and DLMF 10.17.4 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (dqknd), intent(in):: rr(1:2)
real (dqknd), intent(out):: ss(1:2)
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dfrac = 1.5q0, egam = 0.5772156649015328606q0, &
  pi = 3.1415926535897932385q0
integer ic1, ic2, k, nua, n1
real (dqknd) d1, d2
real (dqknd) f1(1:2), f2(1:2), f3(1:2), f4(1:2), &
  f5(1:2), rra(1:2), rr2(1:2), sum1(1:2), &
  sum2(1:2), sum3(1:2), td1(1:2), td2(1:2), &
  tn1(1:2), tn2(1:2), t1(1:2), t2(1:2), &
  t3(1:2), t4(1:2), t41(1:2), t42(1:2), &
  t5(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

!   Check for RR = 0.

if (dqsgn (rr) == 0) then
  write (dqldb, 2)
2 format ('*** DQBESSELYNR: argument is negative or too large')
  call dqabrt
endif

dqnw1 = min (2 * dqnw + 1, dqnwx)
nua = abs (nu)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw*dqnbt, eps)
call dqmdc (rr, d1, n1)
d1 = abs (d1) * 2.q0 ** n1

if (d1 < dfrac * dqnw1 * dqdpw) then
  dqnw1 = min (dqnw + nint (d1 / (dfrac * dqdpw)), 2 * dqnw + 1, dqnwx)
  call dqabs (rr, rra)
  call dqmul (rra, rra, t2)
  call dqmuld (t2, 0.25q0, t1)
  call dqdmc (1.q0, 0, f1)
  call dqdmc (1.q0, 0, f2)
  call dqdmc (1.q0, 0, f3)
  call dqdmc (0.q0, 0,  sum1)

  do k = 1, nua - 1
    call dqmuld (f1, real (k, dqknd), t2)
    call dqeq (t2, f1)
  enddo

  do k = 0, nua - 1
    if (k > 0) then
      call dqdivd (f1, real (nua - k, dqknd), t2)
      call dqeq (t2, f1)
      call dqmul (t1, f2, t2)
      call dqeq (t2, f2)
      call dqmuld (f3, real (k, dqknd), t2)
      call dqeq (t2, f3)
    endif
    call dqmul (f1, f2, t3)
    call dqdiv (t3, f3, t2)
    call dqadd (sum1, t2, t3)
    call dqeq (t3, sum1)
  enddo

  call dqmuld (rra, 0.5q0, t3)
  call dqnpwr (t3, nua, t4)
  call dqdiv (sum1, t4, t3)
  call dqneg (t3, sum1)

  call dqmuld (rra, 0.5q0, t2)
  call dqlog (t2, t3)
  call dqmuld (t3, 2.q0, t2)
  call dqbesseljnr (nua, rra, t3)
  call dqmul (t2, t3, sum2)

  call dqneg (dqegammacon, f1)
  call dqeq (f1, f2)
  call dqdmc (1.q0, 0, f3)
  call dqdmc (1.q0, 0, f4)
  call dqdmc (1.q0, 0, f5)

  do k = 1, nua
    call dqdmc (1.q0, 0, t2)
    call dqdivd (t2, real (k, dqknd), t3)
    call dqadd (f2, t3, t4)
    call dqeq (t4, f2)
    call dqmuld (f5, real (k, dqknd), t2)
    call dqeq (t2, f5)
  enddo

  call dqadd (f1, f2, t2)
  call dqmul (t2, f3, t3)
  call dqmul (f4, f5, t4)
  call dqdiv (t3, t4, sum3)

  do k = 1, itrmax
    call dqdmc (1.q0, 0, t2)
    call dqdivd (t2, real (k, dqknd), t3)
    call dqadd (f1, t3, t4)
    call dqeq (t4, f1)
    call dqdivd (t2, real (nua + k, dqknd), t3)
    call dqadd (f2, t3, t4)
    call dqeq (t4, f2)
    call dqmul (t1, f3, t2)
    call dqneg (t2, f3)
    call dqmuld (f4, real (k, dqknd), t2)
    call dqeq (t2, f4)
    call dqmuld (f5, real (nua + k, dqknd), t2)
    call dqeq (t2, f5)
    call dqadd (f1, f2, t2)
    call dqmul (t2, f3, t3)
    call dqmul (f4, f5, t4)
    call dqdiv (t3, t4, t2)
    call dqadd (sum3, t2, t3)
    call dqeq (t3, sum3)

    call dqabs (t2, tc1)
    call dqmul (eps, sum3, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (dqldb, 6)
6 format ('*** DQBESSELYNR: Loop end error 1')
  call dqabrt

100 continue

  call dqmuld (rra, 0.5q0, t2)
  call dqnpwr (t2, nua, t3)
  call dqmul (t3, sum3, t2)
  call dqneg (t2, sum3)

  call dqadd (sum1, sum2, t2)
  call dqadd (t2, sum3, t4)
  call dqeq (dqpicon, t2)
  call dqdiv (t4, t2, t3)
else

! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  dqnw1 = min (dqnw + 1, dqnwx)
  call dqabs (rr, rra)
  call dqmul (rra, rra, rr2)
  d1 = 4.q0 * real (nua, dqknd)**2
  call dqdmc (d1, 0, t1)
  call dqdmc (1.q0, 0, tn1)
  call dqsub (t1, tn1, t2)
  call dqdivd (t2, 8.q0, tn2)
  call dqdmc (1.q0, 0, td1)
  call dqeq (rra, td2)
  call dqdiv (tn1, td1, sum1)
  call dqdiv (tn2, td2, sum2)

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.q0*(2.q0*k-1.q0) - 1.q0)**2
    d2 = (2.q0*(2.q0*k) - 1.q0)**2
    call dqdmc (d1, 0, t2)
    call dqsub (t1, t2, t3)
    call dqdmc (d2, 0, t2)
    call dqsub (t1, t2, t5)
    call dqmul (t3, t5, t2)
    call dqmul (tn1, t2, t3)
    call dqneg (t3, tn1)

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = real (2*k-1, dqknd) * real (2*k, dqknd) * 64.q0
    call dqmuld (td1, d1, t2)
    call dqmul (t2, rr2, td1)

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

    call dqdiv (tn1, td1, t41)
    call dqadd (sum1, t41, t2)
    call dqeq (t2, sum1)

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.q0*(2.q0*k) - 1.q0)**2
    d2 = (2.q0*(2.q0*k+1.q0) - 1.q0)**2
    call dqdmc (d1, 0, t2)
    call dqsub (t1, t2, t3)
    call dqdmc (d2, 0, t2)
    call dqsub (t1, t2, t5)
    call dqmul (t3, t5, t2)
    call dqmul (tn2, t2, t3)
    call dqneg (t3, tn2)

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = real (2*k, dqknd) * real (2*k+1, dqknd) * 64.q0
    call dqmuld (td2, d1, t2)
    call dqmul (t2, rr2, td2)

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

    call dqdiv (tn2, td2, t42)
    call dqadd (sum2, t42, t2)
    call dqeq (t2, sum2)

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

    call dqabs (t41, tc1)
    call dqmul (eps, sum1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    call dqabs (t42, tc1)
    call dqmul (eps, sum2, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic2)
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (dqldb, 5)
5 format ('*** DQBESSELYNR: Loop end error 2')
  call dqabrt

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

  call dqmuld (dqpicon, 0.5q0 * nua, t1)
  call dqsub (rra, t1, t2)
  call dqmuld (dqpicon, 0.25q0, t1)
  call dqsub (t2, t1, t3)
  call dqcssnr (t3, t41, t42)
  call dqmul (t42, sum1, t1)
  call dqmul (t41, sum2, t2)
  call dqadd (t1, t2, t5)
  call dqmul (dqpicon, rra, t1)
  call dqdmc (2.q0, 0, t2)
  call dqdiv (t2, t1, t3)
  call dqsqrt (t3, t1)
  call dqmul (t1, t5, t3)
endif

if (mod (nu, 2) /= 0) then
!   if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) bessely = - bessely

  ic1 = dqsgn (rr)
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
    call dqneg (t3, t4)
    call dqeq (t4, t3)
  endif
endif

call dqeq (t3, ss)
return
end subroutine dqbesselynr

subroutine dqbesselyr (qq, rr, ss)

!   This evaluates the modified Bessel function BesselY (QQ,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2.

implicit none
real (dqknd), intent(in):: qq(1:2), rr(1:2)
real (dqknd), intent(out):: ss(1:2)
integer i1, i2, n1
real (dqknd) d1
real (dqknd) t1(1:2), t2(1:2), t3(1:2), t4(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)

!   If QQ is integer, call mpbesselynr; if qq < 0 and rr <= 0, then error.

call dqinfr (qq, t1, t2)
i1 = dqsgn (qq)
i2 = dqsgn (rr)
if (dqsgn (t2) == 0) then
  call dqmdc (qq, d1, n1)
  d1 = d1 * 2.q0**n1
  n1 = nint (d1)
  call dqbesselynr (n1, rr, t1)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (dqldb, 2)
2 format ('*** DQBESSELYR: First argument < 0 and second argument <= 0')
  call dqabrt
endif

call dqmul (qq, dqpicon, t1)
call dqcssnr (t1, t2, t3)
call dqbesseljr (qq, rr, t4)
call dqmul (t4, t2, t1)
call dqneg (qq, t2)
call dqbesseljr (t2, rr, t4)
call dqsub (t1, t4, t2)
call dqdiv (t2, t3, t1)

120 continue

call dqeq (t1, ss)
return
end subroutine dqbesselyr

subroutine dqdigammabe (nb2, berne, x, y)

!  This evaluates the digamma function, using asymptotic formula DLMF 5.11.2:
!  dig(x) ~ log(x) - 1/(2*x) - Sum_{k=1}^inf B[2k] / (2*k*x^(2*k)).
!  Before using this formula, the recursion dig(x+1) = dig(x) + 1/x is used
!  to shift the argument up by IQQ, where IQQ is set based on MPNW below.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent (in):: nb2
real (dqknd), intent(in):: berne(1:2,nb2), x(1:2)
real (dqknd), intent(out):: y(1:2)
real (dqknd), parameter:: dber = 1.4q0, dfrac = 0.4q0
integer k, i1, i2, ic1, iqq, n1
real (dqknd) d1
real (dqknd) f1(1:2), sum1(1:2), sum2(1:2), &
  t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  t5(1:2), xq(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)
iqq = dfrac * dqnw1 * dqdpw
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw*dqnbt, eps)
call dqdmc (1.q0, 0, f1)

!   Check if argument is less than or equal to 0 -- undefined.

if (dqsgn (x) <= 0) then
  write (dqldb, 2)
2 format ('*** DQDIGAMMABE: Argument is less than or equal to 0')
  call dqabrt
endif

!   Check if berne array has been initialized.

call dqmdc (berne(1:2,1), d1, n1)
d1 = d1 * 2.q0 ** n1
if (dqwprecr (berne(1:2,1)) < dqnw .or. &
  abs (d1 - 1.q0 / 6.q0) > dqrdfz .or. nb2 < int (dber * dqdpw * dqnw)) then
  write (dqldb, 3) int (dber * dqdpw * dqnw)
3 format ('*** DQDIGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries using DQBERNE or DQBERNER.')
  call dqabrt
endif

! sum1 = mpreal (0.d0, nwds)
! sum2 = mpreal (0.d0, nwds)
! xq = x + dble (iqq)

call dqdmc (0.q0, 0, sum1)
call dqdmc (0.q0, 0, sum2)
call dqdmc (real (iqq, dqknd), 0, t1)
call dqadd (x, t1, t2)
call dqeq (t2, xq)

do k = 0, iqq - 1
!   sum1 = sum1 + f1 / (x + dble (k))

  call dqdmc (real (k, dqknd), 0, t1)
  call dqadd (x, t1, t2)
  call dqdiv (f1, t2, t3)
  call dqadd (sum1, t3, t1)
  call dqeq (t1, sum1)
enddo

! t1 = mpreal (1.d0, nwds)
! t2 = xq ** 2

call dqdmc (1.q0, 0, t1)
call dqmul (xq, xq, t2)

do k = 1, nb2
!  t1 = t1 * t2
!  t3 = bb(k) / (2.d0 * dble (k) * t1)
!  sum2 = sum2 + t3

  call dqmul (t1, t2, t3)
  call dqeq (t3, t1)
  call dqmuld (t1, 2.q0 * real (k, dqknd), t4)
  call dqdiv (berne(1:2,k), t4, t3)
  call dqadd (sum2, t3, t4)
  call dqeq (t4, sum2)

!  if (abs (t3 / sum2) < eps) goto 100

  call dqabs (t3, tc1)
  call dqmul (eps, sum2, tc3)
  call dqabs (tc3, tc2)
  call dqcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 110
enddo

write (dqldb, 4)
4 format ('*** DQDIGAMMABE: Loop end error: Increase NB2')
call dqabrt

110 continue

! digammax = -sum1 + log (xq) - 1.d0 / (2.d0 * xq) - sum2

call dqneg (sum1, t1)
call dqlog (xq, t2)
call dqadd (t1, t2, t3)
call dqmuld (xq, 2.q0, t4)
call dqdiv (f1, t4, t5)
call dqsub (t3, t5, t2)
call dqsub (t2, sum2, t1)
call dqeq (t1, y)
return
end subroutine dqdigammabe

subroutine dqerfr (z, terf)

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
real (dqknd), intent(in):: z(1:2)
real (dqknd), intent(out):: terf(1:2)
integer, parameter:: itrmx = 100000
real (dqknd), parameter:: dcon = 100.q0, pi = 3.1415926535897932385q0
integer ic1, ic2, ic3, k, nbt, n1
real (dqknd) d1, d2

real (dqknd) t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  t5(1:2), t6(1:2), t7(1:2), z2(1:2), tc1(1:2), &
  tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

nbt = dqnw * dqnbt
d1 = aint (1.q0 + sqrt (nbt * log (2.q0)))
d2 = aint (nbt / dcon + 8.q0)
call dqdmc (d1, 0, t1)
call dqdmc (d2, 0, t2)
call dqcpr (z, t1, ic1)
! t1(2) = - t1(2)
call dqneg (t1, t3)
call dqeq (t3, t1)
call dqcpr (z, t1, ic2)
call dqcpr (z, t2, ic3)

if (dqsgn (z) == 0) then
  call dqdmc (0.q0, 0, terf)
elseif (ic1 > 0) then
  call dqdmc (1.q0, 0, terf)
elseif (ic2 < 0) then
  call dqdmc (-1.q0, 0, terf)
elseif (ic3 < 0) then
  call dqmul (z, z, z2)
  call dqdmc (0.q0, 0, t1)
  call dqeq (z, t2)
  call dqdmc (1.q0, 0, t3)
  call dqdmc (1.q10, 0, t5)

  do k = 0, itrmx
    if (k > 0) then
      call dqmuld (z2, 2.q0, t6)
      call dqmul (t6, t2, t7)
      call dqeq (t7, t2)
      d1 = 2.q0 * k + 1.q0
      call dqmuld (t3, d1, t6)
      call dqeq (t6, t3)
    endif

    call dqdiv (t2, t3, t4)
    call dqadd (t1, t4, t6)
    call dqeq (t6, t1)
    call dqdiv (t4, t1, t6)
    call dqcpr (t6, eps, ic1)
    call dqcpr (t6, t5, ic2)
    if (ic1 <= 0 .or. ic2 >= 0) goto 120
    call dqeq (t6, t5)
  enddo

write (dqldb, 3) 1, itrmx
3 format ('*** DQERFR: iteration limit exceeded',2i10)
call dqabrt

120 continue

  call dqmuld (t1, 2.q0, t3)
  call dqsqrt (dqpicon, t4)
  call dqexp (z2, t5)
  call dqmul (t4, t5, t6)
  call dqdiv (t3, t6, t7)
  call dqeq (t7, terf)
else
  call dqmul (z, z, z2)
  call dqdmc (0.q0, 0, t1)
  call dqdmc (1.q0, 0, t2)
  call dqabs (z, t3)
  call dqdmc (1.q10, 0, t5)

  do k = 0, itrmx
    if (k > 0) then
      d1 = -(2.q0 * k - 1.q0)
      call dqmuld (t2, d1, t6)
      call dqeq (t6, t2)
      call dqmul (t2, t3, t6)
      call dqeq (t6, t3)
    endif

    call dqdiv (t2, t3, t4)
    call dqadd (t1, t4, t6)
    call dqeq (t6, t1)
    call dqdiv (t4, t1, t6)
    call dqcpr (t6, eps, ic1)
    call dqcpr (t6, t5, ic2)
    if (ic1 <= 0 .or. ic2 >= 0) goto 130
    call dqeq (t6, t5)
  enddo

write (dqldb, 3) 2, itrmx
call dqabrt

130 continue

  call dqdmc (1.q0, 0, t2)
  call dqsqrt (dqpicon, t3)
  call dqexp (z2, t4)
  call dqmul (t3, t4, t5)
  call dqdiv (t1, t5, t6)
  call dqsub (t2, t6, t7)
  call dqeq (t7, terf)
  if (dqsgn (z) < 0) then
    call dqneg (terf, t6)
    call dqeq (t6, terf)
  endif
endif

return
end subroutine dqerfr

subroutine dqerfcr (z, terfc)

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
real (dqknd), intent(in):: z(1:2)
real (dqknd), intent(out):: terfc(1:2)
integer, parameter:: itrmx = 100000
real (dqknd), parameter:: dcon = 100.q0, pi = 3.1415926535897932385q0
integer ic1, ic2, ic3, k, nbt, n1
real (dqknd) d1, d2
real (dqknd) t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  t5(1:2), t6(1:2), t7(1:2), z2(1:2), tc1(1:2), &
  tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

nbt = dqnw * dqnbt
d1 = aint (1.q0 + sqrt (nbt * log (2.q0)))
d2 = aint (nbt / dcon + 8.q0)
call dqdmc (d1, 0, t1)
call dqdmc (d2, 0, t2)
call dqcpr (z, t1, ic1)
call dqneg (t1, t3)
call dqeq (t3, t1)
call dqcpr (z, t1, ic2)
call dqcpr (z, t2, ic3)

if (dqsgn (z) == 0) then
  call dqdmc (1.q0, 0, terfc)
elseif (ic1 > 0) then
  call dqdmc (0.q0, 0, terfc)
elseif (ic2 < 0) then
  call dqdmc (2.q0, 0, terfc)
elseif (ic3 < 0) then
  call dqmul (z, z, z2)
  call dqdmc (0.q0, 0, t1)
  call dqeq (z, t2)
  call dqdmc (1.q0, 0, t3)
  call dqdmc (1.q10, 0, t5)

  do k = 0, itrmx
    if (k > 0) then
      call dqmuld (z2, 2.q0, t6)
      call dqmul (t6, t2, t7)
      call dqeq (t7, t2)
      d1 = 2.q0 * k + 1.q0
      call dqmuld (t3, d1, t6)
      call dqeq (t6, t3)
    endif

    call dqdiv (t2, t3, t4)
    call dqadd (t1, t4, t6)
    call dqeq (t6, t1)
    call dqdiv (t4, t1, t6)
    call dqcpr (t6, eps, ic1)
    call dqcpr (t6, t5, ic2)
    if (ic1 <= 0 .or. ic2 >= 0) goto 120
    call dqeq (t6, t5)
  enddo

write (dqldb, 3) 1, itrmx
3 format ('*** DQERFCR: iteration limit exceeded',2i10)
call dqabrt

120 continue

  call dqdmc (1.q0, 0, t2)
  call dqmuld (t1, 2.q0, t3)
  call dqsqrt (dqpicon, t4)
  call dqexp (z2, t5)
  call dqmul (t4, t5, t6)
  call dqdiv (t3, t6, t7)
  call dqsub (t2, t7, t6)
  call dqeq (t6, terfc)
else
  call dqmul (z, z, z2)
  call dqdmc (0.q0, 0, t1)
  call dqdmc (1.q0, 0, t2)
  call dqabs (z, t3)
  call dqdmc (1.q10, 0, t5)

  do k = 0, itrmx
    if (k > 0) then
      d1 = -(2.q0 * k - 1.q0)
      call dqmuld (t2, d1, t6)
      call dqeq (t6, t2)
      call dqmul (t2, t3, t6)
      call dqeq (t6, t3)
    endif

    call dqdiv (t2, t3, t4)
    call dqadd (t1, t4, t6)
    call dqeq (t6, t1)
    call dqdiv (t4, t1, t6)
    call dqcpr (t6, eps, ic1)
    call dqcpr (t6, t5, ic2)
    if (ic1 <= 0 .or. ic2 >= 0) goto 130
    call dqeq (t6, t5)
  enddo

write (dqldb, 3) 2, itrmx
call dqabrt

130 continue

  call dqsqrt (dqpicon, t3)
  call dqexp (z2, t4)
  call dqmul (t3, t4, t5)
  call dqdiv (t1, t5, t6)
  if (dqsgn (z) < 0) then
    call dqdmc (2.q0, 0, t2)
    call dqsub (t2, t6, t7)
    call dqeq (t7, t6)
  endif

  call dqeq (t6, terfc)
endif

return
end subroutine dqerfcr

subroutine dqexpint (x, y)

!   This evaluates the exponential integral function Ei(x):
!   Ei(x) = - incgamma (0, -x)

implicit none
real (dqknd), intent(in):: x(1:2)
real (dqknd), intent(out):: y(1:2)
real (dqknd) t1(1:2), t2(1:2), t3(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

if (dqsgn (x) == 0) then
  write (dqldb, 2)
2 format ('*** DQEXPINT: argument is zero')
  call dqabrt
endif

call dqdmc (0.q0, 0, t1)
call dqneg (x, t2)
call dqincgammar (t1, t2, t3)
call dqneg (t3, y)
return
end subroutine dqexpint

subroutine dqgammar (t, z)

!   This evaluates the gamma function, using an algorithm of R. W. Potter.
!   The argument t must not exceed 10^8 in size (this limit is set below),
!   must not be zero, and if negative must not be integer.

!   In the parameter statement below:
!     itrmx = limit of number of iterations in series; default = 100000.
!     con1 = 1/2 * log (10) to DP accuracy.
!     dmax = maximum size of input argument.

implicit none
real (dqknd), intent(in):: t(1:2)
real (dqknd), intent(out):: z(1:2)
integer, parameter:: itrmx = 100000
real (dqknd), parameter:: al2 = 0.69314718055994530942q0, dmax = 1d8, &
  pi = 3.1415926535897932385q0
integer i, i1, ic1, j, nt, n1, n2, n3
real (dqknd) alpha, d1, d2, d3
real (dqknd) f1(1:2), sum1(1:2), sum2(1:2), tn(1:2), &
  t1(1:2), t2(1:2), t3(1:2), t4(1:2), t5(1:2), &
  t6(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

call dqmdc (t, d1, n1)
d1 = d1 * 2.q0**n1
call dqnint (t, t1)
call dqcpr (t, t1, ic1)
i1 = dqsgn (t)
if (i1 == 0 .or. d1 > dmax .or. (i1 < 0 .and. ic1 == 0)) then
  write (dqldb, 2) dmax
2 format ('*** DQGAMMAR: input argument must have absolute value <=',f10.0,','/ &
  'must not be zero, and if negative must not be an integer.')
  call dqabrt
endif

call dqdmc (1.q0, 0, f1)

!   Find the integer and fractional parts of t.

call dqinfr (t, t2, t3)

if (dqsgn (t3) == 0) then

!   If t is a positive integer, then apply the usual factorial recursion.

  call dqmdc (t2, d2, n2)
  nt = d2 * 2.q0 ** n2
  call dqeq (f1, t1)

  do i = 2, nt - 1
    call dqmuld (t1, real (i, dqknd), t2)
    call dqeq (t2, t1)
  enddo

  call dqeq (t1, z)
  goto 120
elseif (dqsgn (t) > 0) then

!   Apply the identity Gamma[t+1] = t * Gamma[t] to reduce the input argument
!   to the unit interval.

  call dqmdc (t2, d2, n2)
  nt = d2 * 2.q0 ** n2
  call dqeq (f1, t1)
  call dqeq (t3, tn)

  do i = 1, nt
    call dqdmc (real (i, dqknd), 0, t4)
    call dqsub (t, t4, t5)
    call dqmul (t1, t5, t6)
    call dqeq (t6, t1)
  enddo
else

!   Apply the gamma identity to reduce a negative argument to the unit interval.

  call dqsub (f1, t, t4)
  call dqinfr (t4, t3, t5)
  call dqmdc (t3, d3, n3)
  nt = d3 * 2.q0 ** n3

  call dqeq (f1, t1)
  call dqsub (f1, t5, t2)
  call dqeq (t2, tn)

  do i = 0, nt - 1
    call dqdmc (real (i, dqknd), 0, t4)
    call dqadd (t, t4, t5)
    call dqdiv (t1, t5, t6)
    call dqeq (t6, t1)
  enddo
endif

!   Calculate alpha = bits of precision * log(2) / 2, then take the next even
!   integer value, so that alpha/2 and alpha^2/4 can be calculated exactly in DP.

alpha = 2.q0 * aint (0.25q0 * (dqnw1 + 1) * dqnbt * al2 + 1.q0)
d2 = 0.25q0 * alpha**2
call dqeq (tn, t2)
call dqdiv (f1, t2, t3)
call dqeq (t3, sum1)

!   Evaluate the series with t.

do j = 1, itrmx
  call dqdmc (real (j, dqknd), 0, t6)
  call dqadd (t2, t6, t4)
  call dqmuld (t4, real (j, dqknd), t5)
  call dqdiv (t3, t5, t6)
  call dqmuld (t6, d2, t3)
  call dqadd (sum1, t3, t4)
  call dqeq (t4, sum1)

  call dqabs (t3, tc1)
  call dqmul (eps, sum1, tc3)
  call dqabs (tc3, tc2)
  call dqcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 100
enddo

write (dqldb, 3) 1, itrmx
3 format ('*** DQGAMMAR: iteration limit exceeded',2i10)
call dqabrt

100 continue

call dqneg (tn, t2)
call dqdiv (f1, t2, t3)
call dqeq (t3, sum2)

!   Evaluate the same series with -t.

do j = 1, itrmx
  call dqdmc (real (j, dqknd), 0, t6)
  call dqadd (t2, t6, t4)
  call dqmuld (t4, real (j, dqknd), t5)
  call dqdiv (t3, t5, t6)
  call dqmuld (t6, d2, t3)
  call dqadd (sum2, t3, t4)
  call dqeq (t4, sum2)

  call dqabs (t3, tc1)
  call dqmul (eps, sum2, tc3)
  call dqabs (tc3, tc2)
  call dqcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 110
enddo

write (dqldb, 3) 2, itrmx
call dqabrt

110 continue

!   Compute sqrt (pi * sum1 / (tn * sin (pi * tn) * sum2))
!   and (alpha/2)^tn terms. Also, multiply by the factor t1, from the
!   If block above.

call dqeq (dqpicon, t2)
call dqmul (t2, tn, t3)
call dqcssnr (t3, t4, t5)
call dqmul (t5, sum2, t6)
call dqmul (tn, t6, t5)
call dqmul (t2, sum1, t3)
call dqdiv (t3, t5, t6)
call dqneg (t6, t4)
call dqeq (t4, t6)
call dqsqrt (t6, t2)
call dqdmc (0.5q0 * alpha, 0, t3)
call dqlog (t3, t4)
call dqmul (tn, t4, t5)
call dqexp (t5, t6)
call dqmul (t2, t6, t3)
call dqmul (t1, t3, t4)

!   Round to mpnw words precision.

call dqeq (t4, z)

120 continue

return
end subroutine dqgammar

subroutine dqhurwitzzetan (is, aa, zz)

!   This returns the Hurwitz zeta function of IS and AA, using an algorithm from:
!   David H. Bailey and Jonathan M. Borwein, "Crandall's computation of the
!   incomplete gamma function and the Hurwitz zeta function with applications to
!   Dirichlet L-series," Applied Mathematics and Computation, vol. 268C (Oct 2015),
!   pg. 462-477, preprint at:
!   https://www.davidhbailey.com/dhbpapers/lerch.pdf
!   This is limited to IS >= 2 and 0 < AA < 1.

implicit none
integer, intent(in):: is
real (dqknd), intent(in):: aa(1:2)
real (dqknd), intent(out):: zz(1:2)
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: pi = 3.1415926535897932385q0
integer i1, ic1, ic2, ic3, k, n1
real (dqknd) d1, dk
real (dqknd) gs1(1:2), gs2(1:2), ss(1:2), &
  sum1(1:2), sum2(1:2), sum3(1:2), ss1(1:2), ss2(1:2), &
  ss3(1:2), ss4(1:2), s1(1:2), s2(1:2), s3(1:2), &
  t1(1:2), t2(1:2), t3(1:2), t4(1:2), t5(1:2), &
  t6(1:2), t7(1:2), t8(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

if (is <= 0) then
  write (dqldb, 3)
3 format ('*** DQHURWITZZETAN: IS less than or equal to 0:')
  call dqabrt
endif

call dqdmc (0.q0, 0, t1)
call dqdmc (1.q0, 0, t2)
call dqcpr (aa, t1, ic1)
call dqcpr (aa, t2, ic2)
if (ic1 <= 0 .or. ic2 >= 0) then
  write (dqldb, 4)
4 format ('*** DQHURWITZZETAN: AA must be in the range (0, 1)')
  call dqabrt
endif

! ss = mpreal (dble (is), mpnw)
! ss1 = 0.5d0 * ss
! ss2 = 0.5d0 * (ss + 1.d0)
! ss3 = 0.5d0 * (1.d0 - ss)
! ss4 = 1.d0 - 0.5d0 * ss

call dqdmc (real (is, dqknd), 0, ss)
call dqmuld (ss, 0.5q0, ss1)
call dqdmc (1.q0, 0, t1)
call dqadd (t1, ss, t2)
call dqmuld (t2, 0.5q0, ss2)
call dqsub (t1, ss, t2)
call dqmuld (t2, 0.5q0, ss3)
call dqsub (t1, ss1, ss4)

! gs1 = gamma (ss1)
! gs2 = gamma (ss2)
! t1 = pi * aa ** 2

call dqgammar (ss1, gs1)
call dqgammar (ss2, gs2)
call dqmul (aa, aa, t2)
call dqmul (dqpicon, t2, t1)

! sum1 = (incgamma (ss1, t1) / gs1 + incgamma (ss2, t1) / gs2) / abs (aa)**is
! sum2 = mpreal (0.d0, mpnw)
! sum3 = mpreal (0.d0, mpnw)

call dqincgammar (ss1, t1, t2)
call dqdiv (t2, gs1, t3)
call dqincgammar (ss2, t1, t2)
call dqdiv (t2, gs2, t4)
call dqadd (t3, t4, t2)
call dqabs (aa, t3)
call dqnpwr (t3, is, t4)
call dqdiv (t2, t4, sum1)
call dqdmc (0.q0, 0, sum2)
call dqdmc (0.q0, 0, sum3)

do k = 1, itrmax
  dk = real (k, dqknd)

!  t1 = pi * (dk + aa)**2
!  t2 = pi * (-dk + aa)**2
!  t3 = dk**2 * pi
!  t4 = 2.d0 * pi * dk * aa

  call dqdmc (dk, 0, t5)
  call dqadd (t5, aa, t6)
  call dqmul (t6, t6, t5)
  call dqmul (dqpicon, t5, t1)
  call dqdmc (-dk, 0, t5)
  call dqadd (t5, aa, t6)
  call dqmul (t6, t6, t7)
  call dqmul (dqpicon, t7, t2)
  call dqmul (t5, t5, t6)
  call dqmul (dqpicon, t6, t3)
  call dqmuld (dqpicon, 2.q0 * dk, t5)
  call dqmul (t5, aa, t4)

!  s1 = (incgamma (ss1, t1) / gs1 + incgamma (ss2, t1) / gs2) / abs (dk + aa)**is

  call dqincgammar (ss1, t1, t5)
  call dqdiv (t5, gs1, t6)
  call dqincgammar (ss2, t1, t5)
  call dqdiv (t5, gs2, t7)
  call dqadd (t6, t7, t5)
  call dqdmc (dk, 0, t6)
  call dqadd (t6, aa, t7)
  call dqabs (t7, t6)
  call dqnpwr (t6, is, t7)
  call dqdiv (t5, t7, s1)

!  s2 = (incgamma (ss1, t2) / gs1 - incgamma (ss2, t2) / gs2) / abs (-dk + aa)**is

  call dqincgammar (ss1, t2, t5)
  call dqdiv (t5, gs1, t6)
  call dqincgammar (ss2, t2, t5)
  call dqdiv (t5, gs2, t7)
  call dqsub (t6, t7, t5)
  call dqdmc (-dk, 0, t6)
  call dqadd (t6, aa, t7)
  call dqabs (t7, t6)
  call dqnpwr (t6, is, t7)
  call dqdiv (t5, t7, s2)

!  sum1 = sum1 + s1
!  sum2 = sum2 + s2

  call dqadd (sum1, s1, t5)
  call dqeq (t5, sum1)
  call dqadd (sum2, s2, t5)
  call dqeq (t5, sum2)

!  s3 = (incgamma (ss3, t3) * cos (t4)/ gs1 + incgamma (ss4, t3) * sin (t4) / gs2) &
!    / mpreal (dk, mpnw)**(1-is)

  call dqincgammar (ss3, t3, t5)
  call dqcssnr (t4, t6, t7)
  call dqmul (t5, t6, t8)
  call dqdiv (t8, gs1, t5)
  call dqincgammar (ss4, t3, t6)
  call dqmul (t6, t7, t8)
  call dqdiv (t8, gs2, t6)
  call dqadd (t5, t6, t7)
  call dqdmc (dk, 0, t5)
  i1 = 1 - is
  call dqnpwr (t5, i1, t6)
  call dqdiv (t7, t6, s3)

!  sum3 = sum3 + s3

  call dqadd (sum3, s3, t5)
  call dqeq (t5, sum3)

!  if (abs (s1) < eps * abs (sum1) .and. abs (s2) < eps * abs (sum2) .and. &
!    abs (s3) < eps * abs (sum3)) goto 100

  call dqabs (s1, tc1)
  call dqmul (eps, sum1, tc3)
  call dqabs (tc3, tc2)
  call dqcpr (tc1, tc2, ic1)
  call dqabs (s2, tc1)
  call dqmul (eps, sum2, tc3)
  call dqabs (tc3, tc2)
  call dqcpr (tc1, tc2, ic2)
  call dqabs (s3, tc1)
  call dqmul (eps, sum3, tc3)
  call dqabs (tc3, tc2)
  call dqcpr (tc1, tc2, ic3)
  if (ic1 <= 0 .and. ic2 <= 0 .and. ic3 <= 0) goto 100
enddo

write (dqldb, 5)
5 format ('*** DQHURWITZZETAN: Loop end error')
call dqabrt

100 continue

if (mod (is, 2) == 0) then
!  t1 = pi ** (is / 2) / ((ss - 1.d0) * gamma (ss1))

  i1 = is / 2
  call dqnpwr (dqpicon, i1, t2)
  call dqdmc (1.q0, 0, t3)
  call dqsub (ss, t3, t4)
  call dqgammar (ss1, t5)
  call dqmul (t4, t5, t6)
  call dqdiv (t2, t6, t1)
else
!  t1 = sqrt (pi) * pi ** ((is - 1) / 2) / ((ss - 1.d0) * gamma (ss1))

  i1 = (is - 1) / 2
  call dqnpwr (dqpicon, i1, t2)
  call dqsqrt (dqpicon, t3)
  call dqmul (t2, t3, t4)
  call dqdmc (1.q0, 0, t2)
  call dqsub (ss, t2, t3)
  call dqgammar (ss1, t5)
  call dqmul (t3, t5, t6)
  call dqdiv (t4, t6, t1)
endif

!t2 = pi ** is / sqrt (pi)

call dqnpwr (dqpicon, is, t3)
call dqsqrt (dqpicon, t4)
call dqdiv (t3, t4, t2)

! hurwitzzetan = t1 + 0.5d0 * sum1 + 0.5d0 * sum2 + t2 * sum3

call dqmuld (sum1, 0.5q0, t3)
call dqmuld (sum2, 0.5q0, t4)
call dqmul (sum3, t2, t5)
call dqadd (t1, t3, t6)
call dqadd (t6, t4, t7)
call dqadd (t7, t5, t1)
call dqeq (t1, zz)

return
end subroutine dqhurwitzzetan

subroutine dqhurwitzzetanbe (nb2, berne, iss, aa, zz)

!  This evaluates the Hurwitz zeta function, using the combination of
!  the definition formula (for large iss), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb2, iss     
real (dqknd), intent(in):: berne(1:2,nb2), aa(1:2)
real (dqknd), intent(out):: zz(1:2)
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dber = 1.4q0, dcon = 0.6q0
integer i1, i2, ic1, iqq, k, n1
real (dqknd) d1, dp
real (dqknd) aq(1:2), aq2(1:2), s1(1:2), s2(1:2), &
  s3(1:2), s4(1:2), t1(1:2), t2(1:2), t3(1:2), &
  t4(1:2), t5(1:2), t6(1:2), eps(1:2), f1(1:2), tc1(1:2), &
  tc2(1:2), tc3(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

!   Check if berne array has been initialized.

call dqmdc (berne(1:2,1), d1, n1)
d1 = d1 * 2.q0 ** n1
if (dqwprecr (berne(1:2,1)) < dqnw .or. &
  abs (d1 - 1.q0 / 6.q0) > dqrdfz .or. nb2 < int (dber * dqdpw * dqnw)) then
  write (dqldb, 3) int (dber * dqdpw * dqnw)
3 format ('*** DQHURWITZZETANBE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries by calling DQBERNE or DQBERNER.')
  call dqabrt
endif

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)
call dqdmc (1.q0, 0, f1)
call dqdmc (0.q0, 0, s1)
call dqdmc (0.q0, 0, s2)
call dqdmc (0.q0, 0, s3)
call dqdmc (0.q0, 0, s4)

if (iss <= 0) then
  write (dqldb, 4)
4 format ('*** DQHURWITZZETANBE: ISS <= 0')
  call dqabrt
endif

if (dqsgn (aa) < 0) then
  write (dqldb, 5)
5 format ('*** DQHURWITZZETANBE: AA < 0')
  call dqabrt
endif

dp = anint (dqnw1 * dqdpw)

!   If iss > a certain value, then use definition formula.

if (iss > 2.303q0 * dp / log (2.515q0 * dp)) then
  do k = 0, itrmax
!    t1 = 1.d0 / (aa + dble (k))**iss
!    s1 = s1 + t1

    call dqdmc (real (k, dqknd), 0, t1)
    call dqadd (aa, t1, t2)
    call dqnpwr (t2, iss, t3)
    call dqdiv (f1, t3, t1)
    call dqadd (s1, t1, t2)
    call dqeq (t2, s1)

!    if (abs (t1 / s1) < eps) goto 110

    call dqabs (t1, tc1)
    call dqmul (eps, s1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 110
  enddo

  write (6, 6)
6 format ('*** DQHURWITZZETANBE: Loop end error 1')
  call dqabrt
endif

call dqmdc (aa, d1, n1)
d1 = d1 * 2.q0**n1
iqq = max (dcon * dqnw1 * dqdpw - d1, 0.q0)

do k = 0, iqq - 1
!  s1 = s1 + 1.d0 / (aa + dble (k))**iss

  call dqdmc (real (k, dqknd), 0, t1)
  call dqadd (aa, t1, t2)
  call dqnpwr (t2, iss, t3)
  call dqdiv (f1, t3, t1)
  call dqadd (s1, t1, t2)
  call dqeq (t2, s1)
enddo

! aq = aa + dble (iqq)

call dqdmc (real (iqq, dqknd), 0, t1)
call dqadd (aa, t1, aq)

! s2 = 1.d0 / (dble (iss - 1) * aq**(iss -  1))

call dqdmc (real (iss - 1, dqknd), 0, t1)
call dqnpwr (aq, iss - 1, t2)
call dqmul (t1, t2, t3)
call dqdiv (f1, t3, s2)

! s3 = 1.d0 / (2.d0 * aq**iss)

call dqnpwr (aq, iss, t1)
call dqmuld (t1, 2.q0, t2)
call dqdiv (f1, t2, s3)

! t1 = mpreal (dble (iss), nwds)
! t2 = mpreal (1.d0, nwds)
! t3 = aq**(iss - 1)
! aq2 = aq**2

call dqdmc (real (iss, dqknd), 0, t1)
call dqdmc (1.q0, 0, t2)
call dqnpwr (aq, iss - 1, t3)
call dqmul (aq, aq, aq2)

do k = 1, nb2
!  if (k > 1) t1 = t1 * dble (iss + 2*k - 3) * dble (iss + 2*k - 2)

  if (k > 1) then
    call dqmuld (t1, real (iss + 2*k - 3, dqknd), t5)
    call dqmuld (t5, real (iss + 2*k - 2, dqknd), t1)
  endif

!  t2 = t2 * dble (2 * k - 1) * dble (2 * k)

  call dqmuld (t2, real (2 * k - 1, dqknd), t5)
  call dqmuld (t5, real (2 * k, dqknd), t2)

!  t3 = t3 * aq2
!  t4 = rb(k) * t1 / (t2 * t3)
!  s4 = s4 + t4

  call dqmul (t3, aq2, t5)
  call dqeq (t5, t3)
  call dqmul (berne(1:2,k), t1, t5)
  call dqmul (t2, t3, t6)
  call dqdiv (t5, t6, t4)
  call dqadd (s4, t4, t5)
  call dqeq (t5, s4)

!  if (abs (t4) < eps) goto 110

  call dqabs (t4, tc1)
  call dqmul (eps, s4, tc3)
  call dqabs (tc3, tc2)
  call dqcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 110
enddo

write (6, 7)
7 format ('*** DQHURWITZZETANBE: End loop error 2; call DQBERNE with larger NB.')
call dqabrt

110 continue

! hurwitz_be = s1 + s2 + s3 + s4

call dqadd (s1, s2, t5)
call dqadd (t5, s3, t6)
call dqadd (t6, s4, s1)
call dqeq (s1, zz)

return
end subroutine dqhurwitzzetanbe

subroutine dqhypergeompfq (np, nq, aa, bb, xx, yy)

!  This returns the HypergeometricPFQ function, namely the sum of the infinite series

!  Sum_0^infinity poch(aa(1),n)*poch(aa(2),n)*...*poch(aa(np),n) /
!      poch(bb(1),n)*poch(bb(2),n)*...*poch(bb(nq),n) * xx^n / n!

!  This subroutine evaluates the HypergeometricPFQ function directly according to
!  this definitional formula. The arrays aa and bb must be dimensioned as shown below.
!  NP and NQ are limited to [1,10].

implicit none
integer, intent(in):: np, nq
real (dqknd), intent(in):: aa(1:2,np), bb(1:2,nq), xx(1:2)
real (dqknd), intent(out):: yy(1:2)
integer, parameter:: itrmax = 1000000, npq = 10
integer i1, i2, ic1, j, k
real (dqknd) sum(1:2), td(1:2), tn(1:2), t1(1:2), &
  t2(1:2), t3(1:2), t4(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

if (np < 1 .or. np > npq .or. nq < 1 .or. nq > npq) then
  write (dqldb, 2) npq
2 format ('*** DQHYPERGEOMPFQ: NP and NQ must be between 1 and',i4)
  call dqabrt
endif

call dqdmc (1.q0, 0, sum)
call dqdmc (1.q0, 0, td)
call dqdmc (1.q0, 0, tn)

do k = 1, itrmax
  call dqdmc (real (k - 1, dqknd), 0, t1)

  do j = 1, np
    call dqadd (t1, aa(1:2,j), t2)
    call dqmul (tn, t2, t3)
    call dqeq (t3, tn)
  enddo

  do j = 1, nq
    call dqadd (t1, bb(1:2,j), t2)
    call dqmul (td, t2, t3)
    call dqeq (t3, td)
  enddo

  call dqmul (tn, xx, t2)
  call dqeq (t2, tn)
  call dqmuld (td, real (k, dqknd), t3)
  call dqeq (t3, td)
  call dqdiv (tn, td, t1)
  call dqadd (sum, t1, t2)
  call dqeq (t2, sum)

  call dqabs (t1, tc1)
  call dqmul (eps, sum, tc3)
  call dqabs (tc3, tc2)
  call dqcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 100
enddo

    write (dqldb, 3) itrmax
3   format ('*** DQHYPERGEOMPFQ: Loop end error',i10)
    call dqabrt

100  continue

call dqeq (sum, yy)
return
end subroutine dqhypergeompfq

subroutine dqincgammar (s, z, g)

!  This returns the incomplete gamma function, using a combination of formula
!  8.7.3 of the DLMF (for modest-sized z), formula 8.11.2 (for large z),
!  a formula from the Wikipedia page for the case S = 0, and another formula
!  from the Wikipedia page for the case S = negative integer. The formula
!  for the case S = 0 requires increased working precision, up to 2.5X normal,
!  depending on the size of Z.

implicit none
real (dqknd), intent(in):: s(1:2), z(1:2)
real (dqknd), intent(out):: g(1:2)
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dmax = 0.833q0, egam = 0.5772156649015328606q0
integer ic1, k, nn, n1, n2
real (dqknd) d1, d2, bits
real (dqknd) t0(1:2), t1(1:2), t2(1:2), &
  t3(1:2), t4(1:2), t5(1:2), f1(1:2), &
  tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

n1 = dqsgn (s)
n2 = dqsgn (z)
if (n2 == 0 .or. n1 /= 0 .and. n2 < 0) then
  write (dqldb, 2)
2 format ('*** DQINCGAMMAR: The second argument must not be zero,'/ &
    'and must not be negative unless the first is zero.')
  call dqabrt
endif

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

call dqdmc (1.q0, 0, f1)
call dqmdc (z, d1, n1)
d1 = d1 * 2.q0 ** n1
bits = dqnw1 * dqnbt

if (abs (d1) < dmax * bits) then

!   This is for modest-sized z.

  call dqinfr (s, t1, t2)
  call dqcpr (s, t1, ic1)
  call dqmdc (s, d2, n2)
  nn = d2 * 2.q0**n2

  if (ic1 == 0 .and. nn == 1) then

!   S = 1; result is exp (-z).

    call dqneg (z, t0)
    call dqexp (t0, t1)
    goto 200
  elseif (ic1 == 0 .and. nn <= 0) then

!    S is zero or a negative integer -- use a different algorithm. In
!    either event, first compute incgamma for S = 0. For large Z, the
!    working precision must be increased, up to 2.5X times normal.

!    mpnw2 = min (mpnw1 + 1.5d0 * d1 / (dmax * bits) * mpnw, 5*mpnw/2+1.d0)

    dqnw2 = dqnw1
    call dqeq (z, t0)
    call dqeq (z, t1)
    call dqdmc (1.q0, 0, t2)

    do k = 2, itrmax
      if (mod (k, 2) == 1) then
        d1 = real (k, dqknd)
        call dqdivd (f1, d1, t3)
        call dqadd (t2, t3, t4)
        call dqeq (t4, t2)
      endif
      call dqmul (z, t1, t3)
      d1 = 2.q0 * real (k, dqknd)
      call dqdivd (t3, d1, t1)
      call dqmul (t1, t2, t3)
      call dqadd (t0, t3, t4)
      call dqeq (t4, t0)

      call dqabs (t3, tc1)
      call dqmul (eps, t0, tc3)
      call dqabs (tc3, tc2)
      call dqcpr (tc1, tc2, ic1)
      if (ic1 <= 0) goto 100
    enddo

    write (dqldb, 4)
4   format ('*** DQINCGAMMAR: Loop end error 1')
    call dqabrt

100  continue

    call dqneg (dqegammacon, t1)
    call dqabs (z, t3)
    call dqlog (t3, t2)
    call dqsub (t1, t2, t3)
    call dqmuld (z, -0.5q0, t4)
    call dqexp (t4, t5)
    call dqmul (t5, t0, t4)
    call dqadd (t3, t4, t1)
    if (nn == 0) goto 200

!   S is negative integer (not zero).

    nn = abs (nn)
    call dqdmc (1.q0, 0, t0)
    call dqeq (t0, t2)

    do k = 1, nn - 1
      call dqmuld (t0, real (k, dqknd), t2)
      call dqeq (t2, t0)
    enddo

    call dqmuld (t0, real (nn, dqknd), t5)

    do k = 1, nn - 1
      call dqmul (t2, z, t3)
      call dqdivd (t3, real (nn - k, dqknd), t4)
      call dqneg (t4, t2)
      call dqadd (t0, t2, t3)
      call dqeq (t3, t0)
    enddo
    
    call dqexp (z, t2)
    call dqdiv (t0, t2, t3)
    call dqnpwr (z, nn, t4)
    call dqdiv (t3, t4, t2)

    if (mod (nn, 2) == 0) then
      call dqadd (t2, t1, t3)
    else
      call dqsub (t2, t1, t3)
    endif
    call dqdiv (t3, t5, t1)
    goto 200
  endif

  call dqgammar (s, t1)
  call dqmul (s, t1, t3)
  call dqdiv (f1, t3, t2)
  call dqeq (t2, t0)

  do k = 1, itrmax
    call dqmul (t2, z, t5)
    call dqdmc (real (k, dqknd), 0, t3)
    call dqadd (s, t3, t4)
    call dqdiv (t5, t4, t2)
    call dqadd (t0, t2, t3)
    call dqeq (t3, t0)

    call dqabs (t2, tc1)
    call dqmul (eps, t0, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 110
  enddo

  write (dqldb, 5) itrmax
5   format ('*** DQINCGAMMAR: Loop end error 1')
  call dqabrt

110 continue

  call dqpower (z, s, t2)
  call dqexp (z, t3)
  call dqdiv (t2, t3, t4)
  call dqmul (t4, t0, t5)
  call dqsub (f1, t5, t2)
  call dqmul (t1, t2, t3)
  call dqeq (t3, t1)
  goto 200
else

!   This is for large z. Note that if S is a positive integer, this loop
!   is finite.

  call dqdmc (1.q0, 0, t0)
  call dqdmc (1.q0, 0, t1)

  do k = 1, itrmax
    call dqdmc (real (k, dqknd), 0, t2)
    call dqsub (s, t2, t3)
    call dqmul (t1, t3, t4)
    call dqdiv (t4, z, t1)
    call dqadd (t0, t1, t2)
    call dqeq (t2, t0)

    call dqabs (t1, tc1)
    call dqmul (eps, t0, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 120
  enddo

  write (dqldb, 6)
6 format ('*** DQINCGAMMAR: Loop end error 3')
  call dqabrt

120 continue

   call dqsub (s, f1, t2)
   call dqpower (z, t2, t3)
   call dqexp (z, t4)
   call dqdiv (t3, t4, t2)
   call dqmul (t2, t0, t1)
   goto 200
endif

200 continue

call dqeq (t1, g)

return
end subroutine dqincgammar

subroutine dqpolygamma (nn, x, y)

!   This returns polygamma (nn, x) for nn >= 0 and 0 < x < 1, by calling
!   mphurwitzzetan.

implicit none
integer, intent(in):: nn
real (dqknd), intent(in):: x(1:2)
real (dqknd), intent(out):: y(1:2)
integer ic1, ic2, k
real (dqknd) t1(1:2), t2(1:2), t3(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw, dqnwx)

if (nn <= 0) then
  write (dqldb, 2)
2 format ('*** DQPOLYGAMMA: NN <= 0')
  call dqabrt
endif

call dqdmc (0.q0, 0, t1)
call dqdmc (1.q0, 0, t2)
call dqcpr (x, t1, ic1)
call dqcpr (x, t2, ic2)
if (ic1 <= 0 .or. ic2 >= 0) then
  write (dqldb, 3)
3 format ('*** DQPOLYGAMMA: X must be in the range (0, 1)')
  call dqabrt
endif

call dqdmc (1.q0, 0, t1)

do k = 1, nn
  call dqmuld (t1, real (k, dqknd), t2)
  call dqeq (t2, t1)
enddo

if (mod (nn + 1, 2) == 1) then
  call dqneg (t1, t2)
  call dqeq (t2, t1)
endif
call dqhurwitzzetan (nn + 1, x, t2)
call dqmul (t1, t2, t3)
call dqeq (t3, y)

return
end subroutine dqpolygamma

subroutine dqpolygammabe (nb2, berne, nn, x, y)

!  This returns polygamma (nn, x) for nn >= 0, by calling mphurwitzzetanbe.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb2, nn
real (dqknd), intent(in):: berne(1:2,nb2), x(1:2)
real (dqknd), intent(out):: y(1:2)
real (dqknd), parameter:: dber = 1.4q0
integer i1, i2, k, n1
real (dqknd) d1
real (dqknd) t1(1:2), t2(1:2), t3(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

!   Check if berne array has sufficient entries.

call dqmdc (berne(1:2,1), d1, n1)
d1 = d1 * 2.q0 ** n1
if (dqwprecr (berne(1:2,1)) < dqnw .or. &
  abs (d1 - 1.q0 / 6.q0) > dqrdfz .or. nb2 < int (dber * dqdpw * dqnw)) then
  write (dqldb, 3) int (dber * dqdpw * dqnw)
3 format ('*** DQPOLYGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries by calling DQBERNE or DQBERER.')
  call dqabrt
endif

dqnw1 = min (dqnw, dqnwx)

if (nn <= 0) then
  write (dqldb, 4)
4 format ('*** DQPOLYGAMMABE: NN <= 0')
  call dqabrt
endif

if (dqsgn (x) < 0) then
  write (dqldb, 5)
5 format ('*** DQPOLYGAMMABE: X < 0')
  call dqabrt
endif

call dqdmc (1.q0, 0, t1)

do k = 1, nn
  call dqmuld (t1, real (k, dqknd), t2)
  call dqeq (t2, t1)
enddo

if (mod (nn + 1, 2) == 1) then
  call dqneg (t1, t2)
  call dqeq (t2, t1)
endif
call dqhurwitzzetanbe (nb2, berne, nn + 1, x, t2)
call dqmul (t1, t2, t3)
call dqeq (t3, y)

return
end subroutine dqpolygammabe

subroutine dqpolylogini (nn, arr)

!   Initializes the MP array arr with data for mppolylogneg.
!   NN must be in the range (-nmax, -1).

implicit none
integer, intent(in):: nn
real (dqknd), intent(out):: arr(1:2,1:abs(nn))
integer, parameter:: nmax = 1000
integer i1, i2, k, n, nna
real (dqknd) aa(1:2,2,abs(nn)), t1(1:2), t2(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

nna = abs (nn)
dqnw1 = min (dqnw + 1, dqnwx)
i1 = 2
i2 = 1
call dqdmc (1.q0, 0, aa(1:2,1,1))
call dqdmc (1.q0, 0, aa(1:2,2,1))

do k = 2, nna
  call dqdmc (0.q0, 0, aa(1:2,1,k))
  call dqdmc (0.q0, 0, aa(1:2,2,k))
enddo

do n = 2, nna
  i1 = 3 - i1
  i2 = 3 - i1

  do k = 2, n
    call dqmuld (aa(1:2,i1,k-1), real (n + 1 - k, dqknd), t1)
    call dqmuld (aa(1:2,i1,k), real (k, dqknd), t2)
    call dqadd (t1, t2, aa(1:2,i2,k))
  enddo
enddo

do k = 1, nna
  call dqeq (aa(1:2,i2,k), arr(1:2,k))
enddo

return
end subroutine dqpolylogini

subroutine dqpolylogneg (nn, arr, x, y)

!   This returns polylog (nn, x) for the case nn < 0. Before calling this,
!   one must call mppolylognini to initialize the array arr for this NN.
!   The dimensions of arr must be as shown below.
!   NN must be in the range (-nmax, -1).
!   The parameter nmxa is the maximum number of additional words of
!   precision needed to overcome cancelation errors when x is negative,
!   for nmax = 1000.

implicit none
integer, intent(in):: nn
real (dqknd), intent(in):: arr(1:2,1:abs(nn)), x(1:2)
real (dqknd), intent(out):: y(1:2)
integer, parameter:: nmax = 1000, nmxa = 8525 / dqnbt + 1
integer i1, i2, k, n1, n2, nna
real (dqknd) d1, d2
real (dqknd) t1(1:2+nmxa), t2(1:2+nmxa), t3(1:2+nmxa), &
  t4(1:2+nmxa)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

nna = abs (nn)
call dqmdc (arr(1:2,1), d1, n1)
d1 = d1 * 2.q0 ** n1
call dqmdc (arr(1:2,nna), d2, n2)
d2 = d2 * 2.q0 ** n2

if (d1 /= 1.q0 .or. d2 /= 1.q0) then
  write (dqldb, 2)
2 format ('*** DQPOLYLOGNEG: Uninitialized or inadequately sized arrays'/ &
  'Call dqpolylogini or polylog_ini to initialize array. See documentation.')
  call dqabrt
endif

dqnw1 = min (dqnw + 1, dqnwx)

if (dqsgn (x) < 0) then
  i1 = (nna + 1) / 2
  call dqmdc (arr(1:2,i1), d1, n1)
  dqnw1 = min (dqnw1 + (n1 + 1) / dqnbt + 1, dqnwx)
endif

call dqeq (x, t1)
call dqeq (t1, t2)

do k = 2, nna
  call dqmul (x, t1, t3)
  call dqeq (t3, t1)
  call dqmul (arr(1:2,k), t1, t3)
  call dqadd (t2, t3, t4)
  call dqeq (t4, t2)
enddo

call dqdmc (1.q0, 0, t3)
call dqsub (t3, x, t4)
call dqnpwr (t4, nna + 1, t3)
call dqdiv (t2, t3, t4)
call dqeq (t4, y)

return
end subroutine dqpolylogneg

subroutine dqpolylogpos (nn, x, y)

!   This returns polylog (nn, x) for the case nn >= 0.

implicit none
integer, intent(in):: nn
real (dqknd), intent(in):: x(1:2)
real (dqknd), intent(out):: y(1:2)
integer, parameter:: itrmax = 1000000
integer ic1, k
real (dqknd) t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  t5(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps (1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

if (nn < 0) then
  write (dqldb, 1)
1 format ('*** DQPOLYLOGPOS: N is less than zero.'/ &
  'For negative n, call dqpolylogneg or polylog_neg. See documentation.')
  call dqabrt
endif

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

call dqabs (x, t1)
call dqdmc (1.q0, 0, t2)
call dqcpr (t1, t2, ic1)
if (ic1 >= 0) then
  write (dqldb, 3)
3 format ('*** DQPOLYLOGPOS: |X| must be less than one.')
  call dqabrt
endif

if (nn == 0) then
  call dqdmc (1.q0, 0, t1)
  call dqsub (t1, x, t2)
  call dqdiv (x, t2, t3)
  call dqeq (t3, y)
else
  call dqeq (x, t1)
  call dqeq (x, t2)

  do k = 2, itrmax
    call dqmul (x, t2, t3)
    call dqeq (t3, t2)
    call dqdmc (real (k, dqknd), 0, t3)
    call dqnpwr (t3, nn, t4)
    call dqdiv (t2, t4, t3)
    call dqadd (t1, t3, t4)
    call dqeq (t4, t1)

    call dqabs (t3, tc1)
    call dqmul (eps, t1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (dqldb, 4)
4 format ('*** DQPOLYLOGPOS: Loop end error')
  call dqabrt

100 continue

  call dqeq (t1, y)
endif

return
end subroutine dqpolylogpos

subroutine dqstruvehn (nu, ss, zz)

!   This returns the StruveH function with integer arg NU and MPFR argument SS.

implicit none
integer, intent(in):: nu
real (dqknd), intent(in):: ss(1:2)
real (dqknd), intent(out):: zz(1:2)
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dmax = 1000.q0, pi = 3.1415926535897932385q0
integer ic1, k, n1
real (dqknd) d1
real (dqknd) sum(1:2), td1(1:2), td2(1:2), tn1(1:2), &
  tnm1(1:2), t1(1:2), t2(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

if (nu < 0) then
  write (dqldb, 2)
2 format ('*** DQSTRUVEHN: NU < 0')
  call dqabrt
endif

call dqmdc (ss, d1, n1)
d1 = abs (d1) * 2.q0**n1
if (d1 > dmax) then
  write (dqldb, 3)
3 format ('*** DQSTRUVEHN: ABS(SS) >',f8.2)
  call dqabrt
endif

dqnw1 = min (dqnw * (1.q0 + d1 / dmax), real (dqnwx, dqknd))

call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw*dqnbt, eps)

! tn1 = mpreal (1.d0, nwds1)
! tnm1 = -0.25d0 * mpreal (ss, nwds1)**2
! td1 = 0.5d0 * sqrt (mppi (nwds1))
! td2 = td1

call dqdmc (1.q0, 0, tn1)
call dqmul (ss, ss, t1)
call dqmuld (t1, -0.25q0, tnm1)
call dqsqrt (dqpicon, t1)
call dqmuld (t1, 0.5q0, td1)
call dqeq (td1, td2)

! do k = 1, nu
!  td2 = (k + 0.5d0) * td2
! enddo

do k = 1, nu
  d1 = k + 0.5q0
  call dqmuld (td2, d1, t1)
  call dqeq (t1, td2)
enddo

! sum = tn1 / (td1 * td2)

call dqmul (td1, td2, t2)
call dqdiv (tn1, t2, sum)

do k = 1, itrmax

!  tn1 = tnm1 * tn1
!  td1 = (k + 0.5d0) * td1
!  td2 = (nu + k + 0.5d0) * td2
!  t1 = tn1 / (td1 * td2)
!  sum = sum + t1

  call dqmul (tnm1, tn1, t1)
  call dqeq (t1, tn1)
  d1 = k + 0.5q0
  call dqmuld (td1, d1, t1)
  call dqeq (t1, td1)
  d1 = nu + k + 0.5q0
  call dqmuld (td2, d1, t1)
  call dqeq (t1, td2)
  call dqmul (td1, td2, t2)
  call dqdiv (tn1, t2, t1)
  call dqadd (sum, t1, t2)
  call dqeq (t2, sum)

!  if (abs (t1) < eps) goto 100

  call dqabs (t1, tc1)
  call dqmul (eps, sum, tc3)
  call dqabs (tc3, tc2)
  call dqcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 100
enddo

write (dqldb, 5)
5 format ('*** DQSTRUVEHN: Loop end error')
call dqabrt

100 continue

! struvehn = (0.5d0 * ss)**(nu + 1) * sum

call dqmuld (ss, 0.5q0, t1)
n1 = nu + 1
call dqnpwr (t1, n1, t2)
call dqmul (t2, sum, t1)
call dqeq (t1, zz)

return
end subroutine dqstruvehn

subroutine dqzetar (ss, zz)

!   This returns the zeta function of an MPR argument SS using an algorithm
!   due to Peter Borwein.

implicit none
real (dqknd), intent(in):: ss(1:2)
real (dqknd), intent(out):: zz(1:2)
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dfrac = 1.q0+dqdpw, pi = 3.1415926535897932385q0
integer i, ic1, iss, j, n, n1, n2
real (dqknd) d1, d2
real (dqknd) f1(1:2), s(1:2), t1(1:2), t2(1:2), &
  t3(1:2), t4(1:2), t5(1:2), tn(1:2), tt(1:2), &
  tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)
real (dqknd) sgn

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

call dqdmc (1.q0, 0, f1)
call dqcpr (ss, f1, ic1)
call dqinfr (ss, t1, t2)

if (ic1 == 0) then
  write (dqldb, 2)
2 format ('*** DQZETAR: argument is 1')
  call dqabrt
elseif (dqsgn (t2) == 0) then

!   The argument is an integer value. Call mpzetaintr instead.

  call dqmdc (ss, d1, n1)
  iss = d1 * 2.q0**n1
  call dqzetaintr (iss, t1)
  goto 200
elseif (dqsgn (ss) < 0) then

!   If arg < 0, compute zeta(1-ss), and later apply Riemann's formula.

  call dqsub (f1, ss, tt)
else
  call dqeq (ss, tt)
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = dqnbt * dqnw * log (2.q0) / log (2.q0 * dqnbt * dqnw / 3.q0)
call dqmdc (tt, d2, n2)
d2 = d2 * 2.q0 ** n2

if (d2 > d1) then

!   Evaluate the infinite series.

  call dqdmc (1.q0, 0, t1)

  do i = 2, itrmax
    call dqdmc (real (i, dqknd), 0, t4)
    call dqpower (t4, tt, t2)
    call dqdiv (f1, t2, t3)
    call dqadd (t1, t3, t2)
    call dqeq (t2, t1)

    call dqabs (t3, tc1)
    call dqmul (eps, t1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 200
  enddo

  write (dqldb, 3) 1, itrmax
3 format ('*** DQZETAR: iteration limit exceeded',2i10)
  call dqabrt
endif

n = dfrac * dqnw1
call dqdmc (2.q0, 0, t1)
call dqnpwr (t1, n, tn)
call dqneg (tn, t1)
call dqdmc (0.q0, 0, t2)
call dqdmc (0.q0, 0, s)
sgn = 1.q0

do j = 0, 2 * n - 1
  call dqdmc (real (j + 1, dqknd), 0, t4)
  call dqpower (t4, tt, t3)
  call dqdiv (t1, t3, t4)
  call dqmuld (t4, sgn, t5)
  call dqadd (s, t5, t4)
  call dqeq (t4, s)
  sgn = - sgn

  if (j < n - 1) then
    call dqdmc (0.q0, 0, t2)
  elseif (j == n - 1) then
    call dqdmc (1.q0, 0, t2)
  else
    call dqmuld (t2, real (2 * n - j, dqknd), t3)
    call dqdivd (t3, real (j + 1 - n, dqknd), t2)
  endif
  call dqadd (t1, t2, t3)
  call dqeq (t3, t1)
enddo

call dqsub (f1, tt, t3)
call dqdmc (2.q0, 0, t2)
call dqpower (t2, t3, t4)
call dqsub (f1, t4, t2)
call dqmul (tn, t2, t3)
call dqdiv (s, t3, t1)
call dqneg (t1, t2)
call dqeq (t2, t1)

!   If original argument was negative, apply Riemann's formula.

if (dqsgn (ss) < 0) then
  call dqgammar (tt, t3)
  call dqmul (t1, t3, t2)
  call dqmul (dqpicon, tt, t1)
  call dqmuld (t1, 0.5q0, t3)
  call dqcssnr (t3, t4, t5)
  call dqmul (t2, t4, t1)
  call dqmuld (dqpicon, 2.q0, t2)
  call dqpower (t2, tt, t3)
  call dqdiv (t1, t3, t2)
  call dqmuld (t2, 2.q0, t1)
endif

200 continue

call dqeq (t1, zz)
return
end subroutine dqzetar

subroutine dqzetaintr (iss, zz)

!   This returns the zeta function of the integer argument ISS using an algorithm
!   due to Peter Borwein.

implicit none
integer, intent(in):: iss
real (dqknd), intent(out):: zz(1:2)
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dfrac = 1.q0+dqdpw, pi = 3.1415926535897932385q0
integer i, ic1, j, n, n1, itt
real (dqknd) d1, sgn
real (dqknd) f1(1:2), s(1:2), t1(1:2), t2(1:2), &
  t3(1:2), t4(1:2), t5(1:2), tn(1:2), &
  tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

dqnw1 = min (dqnw + 1, dqnwx)
call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

call dqdmc (1.q0, 0, f1)

if (iss == 1) then
  write (dqldb, 2)
2 format ('*** DQZETAINTR: argument is 1')
  call dqabrt
elseif (iss == 0) then

!   Argument is zero -- result is -1/2.

  call dqdmc (-0.5q0, 0, t1)
  goto 200
elseif (iss < 0) then

!   If argument is a negative even integer, the result is zero.

  if (mod (iss, 2) == 0) then
    call dqdmc (0.q0, 0, t1)
    goto 200
  endif

!   Otherwise if arg < 0, compute zeta(1-is), and later apply Riemann's formula.

  itt = 1 - iss
else
  itt = iss
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = dqnbt * dqnw * log (2.q0) / log (2.q0 * dqnbt * dqnw / 3.q0)

if (itt > d1) then

!   Evaluate the infinite series.

  call dqdmc (1.q0, 0, t1)

  do i = 2, itrmax
    call dqdmc (real (i, dqknd), 0, t4)
    call dqnpwr (t4, itt, t2)
    call dqdiv (f1, t2, t3)
    call dqadd (t1, t3, t2)
    call dqeq (t2, t1)

    call dqabs (t3, tc1)
    call dqmul (eps, t1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 200
  enddo

  write (dqldb, 3) 1, itrmax
3 format ('*** DQZETAINTR: iteration limit exceeded',2i10)
  call dqabrt
endif

n = dfrac * dqnw1
call dqdmc (2.q0, 0, t1)
call dqnpwr (t1, n, tn)
call dqneg (tn, t1)
call dqdmc (0.q0, 0, t2)
call dqdmc (0.q0, 0, s)

sgn = 1.q0

do j = 0, 2 * n - 1
  call dqdmc (real (j + 1, dqknd), 0, t4)
  call dqnpwr (t4, itt, t3)
  call dqdiv (t1, t3, t4)
  call dqmuld (t4, sgn, t5)
  call dqadd (s, t5, t4)
  call dqeq (t4, s)
  sgn = - sgn

  if (j < n - 1) then
    call dqdmc (0.q0, 0, t2)
  elseif (j == n - 1) then
    call dqdmc (1.q0, 0, t2)
  else
    call dqmuld (t2, real (2 * n - j, dqknd), t3)
    call dqdivd (t3, real (j + 1 - n, dqknd), t2)
  endif

  call dqadd (t1, t2, t3)
  call dqeq (t3, t1)
enddo

call dqdmc (2.q0, 0, t2)
call dqnpwr (t2, 1 - itt, t4)
call dqsub (f1, t4, t2)
call dqmul (tn, t2, t3)
call dqdiv (s, t3, t1)
call dqneg (t1, t2)
call dqeq (t2, t1)

!   If original argument was negative, apply Riemann's formula.

if (iss < 0) then
  call dqdmc (1.q0, 0, t3)
  do i = 1, itt - 1
    call dqmuld (t3, real (i, dqknd), t4)
    call dqeq (t4, t3)
  enddo

  call dqmul (t1, t3, t2)
  call dqmuld (dqpicon, real (itt, dqknd), t1)
  call dqmuld (t1, 0.5q0, t3)
  call dqcssnr (t3, t4, t5)
  call dqmul (t2, t4, t1)
  call dqmuld (dqpicon, 2.q0, t2)
  call dqnpwr (t2, itt, t3)
  call dqdiv (t1, t3, t2)
  call dqmuld (t2, 2.q0, t1)
endif

200 continue

call dqeq (t1, zz)
return
end subroutine dqzetaintr

subroutine dqzetabe (nb2, berne, s, z)

!  This evaluates the Riemann zeta function, using the combination of
!  the definition formula (for large s), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb2
real (dqknd), intent(in):: berne(1:2,nb2), s(1:2)
real (dqknd), intent(out):: z(1:2)
integer, parameter:: itrmax = 1000000
real (dqknd), parameter:: dber = 1.5q0, dfrac = 0.6q0, &
  pi = 3.1415926535897932385q0
integer i, i1, i2, ic1, k, n1, n2, nn
real (dqknd) d1, d2
real (dqknd) t0(1:2), t1(1:2), t2(1:2), t3(1:2), &
  t4(1:2), t5(1:2), t6(1:2), t7(1:2), t8(1:2), &
  t9(1:2), tt(1:2), f1(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer dqnw, dqnw1, dqnw2



dqnw = dqnwx

!   Check if berne array has been initialized.

call dqmdc (berne(1:2,1), d1, n1)
d1 = d1 * 2.q0 ** n1
if (dqwprecr (berne(1:2,1)) < dqnw .or. &
  abs (d1 - 1.q0 / 6.q0) > dqrdfz .or. nb2 < int (dber * dqdpw * dqnw)) then
  write (dqldb, 3) int (dber * dqdpw * dqnw)
3 format ('*** DQZETABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries.')
  call dqabrt
endif

i = 0
k = 0
dqnw1 = min (dqnw + 1, dqnwx)

call dqdmc (2.q0, 0, tc1)
call dqnpwr (tc1, -dqnw1*dqnbt, eps)

!   Check if argument is 1 -- undefined.

call dqdmc (1.q0, 0, t0)
call dqcpr (s, t0, ic1)
if (ic1 == 0) then
  write (dqldb, 2)
2 format ('*** DQZETABE: argument is 1')
  call dqabrt
endif

call dqdmc (1.q0, 0, f1)

!   Check if argument is zero. If so, result is - 1/2.

if (dqsgn (s) == 0) then
  call dqdmc (-0.5q0, 0, t1)
  goto 200
endif

!   Check if argument is negative.

if (dqsgn (s) < 0) then

!   Check if argument is a negative even integer. If so, the result is zero.

  call dqmuld (s, 0.5q0, t1)
  call dqinfr (t1, t2, t3)
  if (dqsgn (t3) == 0) then
    call dqdmc (0.q0, 0, t1)
    goto 200
  endif

!   Otherwise compute zeta(1-s), and later apply the reflection formula.

  call dqsub (f1, s, tt)
else
  call dqeq (s, tt)
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = dqlogb * dqnw1 / log (32.q0 * dqnw1)
call dqmdc (tt, d2, n2)
d2 = d2 * 2.q0**n2
if (d2 > d1) then
  call dqdmc (1.q0, 0, t1)

  do i = 2, itrmax
    call dqdmc (real (i, dqknd), 0, t4)
    call dqpower (t4, tt, t2)
    call dqdiv (f1, t2, t3)
    call dqadd (t1, t3, t2)
    call dqeq (t2, t1)

    call dqabs (t3, tc1)
    call dqmul (eps, t1, tc3)
    call dqabs (tc3, tc2)
    call dqcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 200
  enddo

  write (dqldb, 4) 1, itrmax
4 format ('*** DQZETABE: iteration limit exceeded',2i10)
  call dqabrt
endif

call dqdmc (1.q0, 0, t0)
nn = dfrac * dqdpw * dqnw1

do k = 2, nn
  call dqdmc (real (k, dqknd), 0, t2)
  call dqpower (t2, tt, t1)
  call dqdiv (f1, t1, t2)
  call dqadd (t0, t2, t3)
  call dqeq (t3, t0)
enddo

call dqdmc (real (nn, dqknd), 0, t2)
call dqsub (tt, f1, t3)
call dqmul (t1, t3, t4)
call dqdiv (t2, t4, t3)
call dqadd (t0, t3, t2)
call dqdmc (0.5q0, 0, t3)
call dqdiv (t3, t1, t4)
call dqsub (t2, t4, t0)

call dqeq (tt, t3)
d1 = 12.q0 * real (nn, dqknd)
call dqmuld (t1, d1, t4)
call dqdiv (t3, t4, t2)
call dqmuld (t1, real (nn, dqknd), t5)
call dqdmc (real (nn, dqknd), 0, t6)
call dqmul (t6, t6, t9)

do k = 2, min (nb2, itrmax)
  call dqdmc (real (2 * k - 2, dqknd), 0, t4)
  call dqadd (tt, t4, t6)
  call dqdmc (real (2 * k - 3, dqknd), 0, t7)
  call dqadd (tt, t7, t8)
  call dqmul (t6, t8, t7)
  call dqmul (t3, t7, t4)
  call dqdmc (real (2 * k - 1, dqknd), 0, t6)
  call dqdmc (real (2 * k - 2, dqknd), 0, t7)
  call dqmul (t6, t7, t8)
  call dqdiv (t4, t8, t3)
  call dqmul (t5, t9, t6)
  call dqeq (t6, t5)
  call dqmul (t3, berne(1:2,k), t4)
  call dqmuld (t5, real (2 * k, dqknd), t6)
  call dqdiv (t4, t6, t7)
  call dqadd (t2, t7, t4)
  call dqeq (t4, t2)

  call dqabs (t7, tc1)
  call dqmul (eps, t2, tc3)
  call dqabs (tc3, tc2)
  call dqcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 110
enddo

write (dqldb, 4) 2, min (nb2, itrmax)
call dqabrt

110 continue

call dqadd (t0, t2, t1)

!   If original argument was negative, apply the reflection formula.

if (dqsgn (s) < 0) then
  call dqgammar (tt, t3)
  call dqmul (t1, t3, t2)
  call dqmul (dqpicon, tt, t1)
  call dqmuld (t1, 0.5q0, t3)
  call dqcssnr (t3, t4, t5)
  call dqmul (t2, t4, t1)
  call dqmuld (dqpicon, 2.q0, t2)
  call dqpower (t2, tt, t3)
  call dqdiv (t1, t3, t2)
  call dqmuld (t2, 2.q0, t1)
endif

200 continue

call dqeq (t1, z)

return
end subroutine dqzetabe

end module dqfune
