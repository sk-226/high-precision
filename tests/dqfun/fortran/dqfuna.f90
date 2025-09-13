
!  DQFUN: A double-quad precision package with special functions

!  Computational routine module (DQFUNA).

!  Revision date:  16 Mar 2023

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired)
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to double-quad precision (approximately 65 digits), by making only
!    relatively minor changes to existing Fortran programs. All basic arithmetic
!    operations and transcendental functions are supported, together with numerous
!    special functions. The package should run correctly on any Unix-based system
!    supporting a Fortran-2008 compiler and IEEE 128-bit floating-point arithmetic
!    (in hardware or software). Note however that results are NOT guaranteed to the
!    last bit.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.

!    Two related software packages by the same author are DDFUN (double-double real
!    and complex) and MPFUN-2020 (arbitrary precision real and complex). They are
!    available from the same site as this package:
!    http://www.davidhbailey.com/dhbsoftware/

!  DOCUMENTATION:
!    See the README-dqfun.txt file in the main DQFUN directory.

!  DESCRIPTION OF THIS MODULE (DQFUNA):
!    This module contains most lower-level computational routines.

!  The following notational scheme is used to designate datatypes below:

!  QP  Quad precision [i.e. REAL (KIND (0.Q0))] -- approx. 34 digit accuracy
!  QC  Quad complex [i.e. COMPLEX (KIND (0.Q0))] -- approx. 34 digit accuracy
!  DQR Double-quad real -- approx. 66 decimal digit accuracy
!  DQC Double-quad complex -- approx. 66 decimal digit accuracy 

module dqfuna
integer, public, parameter :: dqknd = selected_real_kind (33, 4931), dqldb = 6, &
  dqnbt = 113, dqnwx = 2
real (dqknd), public, parameter :: dqdpw = 34.016389510029875059152495103867712q0, &
  dqlogb = 77.632484222713874654729997603315776q0, dqrdfz = 1.q0 / 2.q0**110
real (dqknd), public, parameter :: dqpicon(1:2) = [ &
 3.1415926535897932384626433832795027974791q+00, &
 8.6718101301237810247970440260433494722008q-35]
real (dqknd), public, parameter :: dqegammacon(1:2) = [ &
    5.7721566490153286060651209008240247066034q-01, &
   -3.9618179631972089534583347132391012525411q-35]

contains

subroutine dqabrt
implicit none

!   This permits one to insert a call to a vendor-specific traceback routine.

stop
end subroutine dqabrt

subroutine dqabs (a, b)

!   This sets b = abs (a).

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2)

if (a(1) >= 0.q0) then
  b(1) = a(1); b(2) = a(2)
else
  b(1) = - a(1); b(2) = - a(2)
endif

return
end

subroutine dqadd (dqa, dqb, dqc)

!   This subroutine computes dqc = dqa + dqb, where dqa, dqb and dqc are type DQR.

implicit none
real (dqknd), intent(in):: dqa(2), dqb(2)
real (dqknd), intent(out):: dqc(2)
real (dqknd) e, t1, t2

!   Compute dqa + dqb using Knuth's trick.

t1 = dqa(1) + dqb(1)
e = t1 - dqa(1)
t2 = ((dqb(1) - e) + (dqa(1) - (t1 - e))) + dqa(2) + dqb(2)

!   The result is t1 + t2, after normalization.

dqc(1) = t1 + t2
dqc(2) = t2 - (dqc(1) - t1)
return
end subroutine dqadd

subroutine dqacosh (a, b)

!   This computes the inverse hyperbolic cosine of a, using the standard formula.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2)
real (dqknd) f1(2), t1(2), t2(2)

!   Check if a < 1; error.

if (a(1) < 1.d0) then
  write (6, 1)
1 format ('DQACOSH: Argument is < 1.')
  call dqabrt
endif

f1(1) = 1.d0; f1(2) = 0.d0
call dqmul (a, a, t1)
call dqsub (t1, f1, t2)
call dqsqrt (t2, t1)
call dqadd (a, t1, t2)
call dqlog (t2, b)
end subroutine dqacosh

subroutine dqagmr (a, b, c)

!   This performs the arithmetic-geometric mean (AGM) iterations on A and B.
!   The AGM algorithm is as follows: Set a_0 = a and b_0 = b, then iterate

!    a_{k+1} = (a_k + b_k)/2
!    b_{k+1} = sqrt (a_k * b_k)

!   until convergence (i.e., until a_k = b_k to available precision).
!   The result is returned in C.

implicit none
real (dqknd), intent(in):: a(1:2), b(1:2)
real (dqknd), intent(out):: c(1:2)
integer, parameter:: itrmx = 100
integer j
real (dqknd) eps, s0(1:2), s1(1:2), s2(1:2), s3(1:2)

eps = 2.q0 ** (2 - dqnwx*dqnbt)
call dqeq (a, s1)
call dqeq (b, s2)

do j = 1, itrmx
  call dqadd (s1, s2, s0)
  call dqmuld (s0, 0.5q0, s3)
  call dqmul (s1, s2, s0)
  call dqsqrt (s0, s2)
  call dqeq (s3, s1)

!   Check for convergence.

  call dqsub (s1, s2, s0)
  if (s0(1) == 0.q0 .or. s0(1) / s1(1) < eps) goto 100
enddo

write (dqldb, 2)
2 format ('*** DQAGMR: Iteration limit exceeded.')
call dqabrt

100 continue

call dqeq (s1, c)

return
end subroutine dqagmr

subroutine dqang (x, y, a)

!   This computes the DQR angle A subtended by the DQR pair (X, Y) considered as
!   a point in the x-y plane. This is more useful than an arctan or arcsin
!   routine, since it places the result correctly in the full circle, i.e.
!   -Pi < A <= Pi.

!   The Taylor series for Sin converges much more slowly than that of Arcsin.
!   Thus this routine does not employ Taylor series, but instead computes
!   Arccos or Arcsin by solving Cos (a) = x or Sin (a) = y using one of the
!   following Newton iterations, both of which converge to a:

!           z_{k+1} = z_k - [x - Cos (z_k)] / Sin (z_k)
!           z_{k+1} = z_k + [y - Sin (z_k)] / Cos (z_k)

!   The first is selected if Abs (x) <= Abs (y); otherwise the second is used.

implicit none
real (dqknd), intent(in):: x(2), y(2)
real (dqknd), intent(out):: a(2)
real (dqknd) t1, t2, t3
integer i, ix, iy, k, kk, nx, ny
real (dqknd) s0(2), s1(2), s2(2), s3(2), s4(2)
real (dqknd), parameter:: pi(1:2) = &
  [3.1415926535897932384626433832795027974791q+00, &
  8.6718101301237810247970440260433494722008q-35]

!   Check if both X and Y are zero.

if (x(1) == 0.q0 .and. y(1) == 0.q0) then
  write (dqldb, 1)
1 format ('*** DQANG: Both arguments are zero.')
  call dqabrt
  return
endif

!   Check if one of X or Y is zero.

if (x(1) == 0.q0) then
  if (y(1) > 0.q0) then
    call dqmuld (pi, 0.5q0, a)
  else
    call dqmuld (pi, -0.5q0, a)
  endif
  goto 120
elseif (y(1) == 0.q0) then
  if (x(1) > 0.q0) then
      a(1) = 0.q0
      a(2) = 0.q0
  else
    a(1) = pi(1)
    a(2) = pi(2)
  endif
  goto 120
endif

!   Normalize x and y so that x^2 + y^2 = 1.

call dqmul (x, x, s0)
call dqmul (y, y, s1)
call dqadd (s0, s1, s2)
call dqsqrt (s2, s3)
call dqdiv (x, s3, s1)
call dqdiv (y, s3, s2)

!   Compute initial approximation of the angle.

call dqdqdpc (s1, t1)
call dqdqdpc (s2, t2)
t3 = atan2 (t2, t1)
a(1) = t3
a(2) = 0.q0

!   The smaller of x or y will be used from now on to measure convergence.
!   This selects the Newton iteration (of the two listed above) that has the
!   largest denominator.

if (abs (t1) <= abs (t2)) then
  kk = 1
  s0(1) = s1(1)
  s0(2) = s1(2)
else
  kk = 2
  s0(1) = s2(1)
  s0(2) = s2(2)
endif

!   Perform the Newton-Raphson iteration described.

do k = 1, 3
  call dqcssnr (a, s1, s2)
  if (kk == 1) then
    call dqsub (s0, s1, s3)
    call dqdiv (s3, s2, s4)
    call dqsub (a, s4, s1)
  else
    call dqsub (s0, s2, s3)
    call dqdiv (s3, s1, s4)
    call dqadd (a, s4, s1)
  endif
  a(1) = s1(1)
  a(2) = s1(2)
enddo

 120  continue

return
end subroutine dqang

subroutine dqasinh (a, b)

!   This computes the inverse hyperbolic sine of a, using the standard formula.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2)
real (dqknd) f1(2), t1(2), t2(2)

f1(1) = 1.d0; f1(2) = 0.d0
call dqmul (a, a, t1)
call dqadd (t1, f1, t2)
call dqsqrt (t2, t1)
call dqadd (a, t1, t2)
call dqlog (t2, b)
end subroutine dqasinh

subroutine dqatanh (a, b)

!   This computes the inverse hyperbolic tangent of a, using the standard formula.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2)
real (dqknd) f1(2), t1(2), t2(2), t3(2)

!   Check if a <= -1 or a >= 1; error.

if (abs (a(1)) >= 1.d0) then
  write (6, 1)
1 format ('DQATANH: Argument is <= -1 or >= 1.')
  call dqabrt
endif

f1(1) = 1.d0; f1(2) = 0.d0
call dqadd (f1, a, t1)
call dqsub (f1, a, t2)
call dqdiv (t1, t2, t3)
call dqlog (t3, t1)
call dqmuld (t1, 0.5q0, b)
end subroutine dqatanh

subroutine dqcadd (a, b, c)

!   This computes the sum of the DQC numbers A and B and returns the DQC
!   result in C.

implicit none
real (dqknd) a(4), b(4), c(4)

call dqadd (a, b, c)
call dqadd (a(3), b(3), c(3))

return
end subroutine dqcadd

subroutine dqcdiv (a, b, c)

!   This routine divides the DQC numbers A and B to yield the DQC quotient C.
!   This routine employs the formula described in DQCMUL to save multiprecision
!   multiplications.

implicit none
real (dqknd), intent(in):: a(4), b(4)
real (dqknd), intent(out):: c(4)
real (dqknd) f(2), s0(2), s1(2), s2(2), s3(2), s4(2)

if (b(1) == 0.q0 .and. b(3) == 0.q0) then
  write (dqldb, 1)
1 format ('*** DQCDIV: Divisor is zero.')
  call dqabrt
  return
endif

f(1) = 1.q0
f(2) = 0.q0
call dqmul (a, b, s0)
call dqmul (a(3), b(3), s1)
call dqadd (s0, s1, s2)
call dqsub (s0, s1, s3)
call dqadd (a, a(3), s0)
call dqsub (b, b(3), s1)
call dqmul (s0, s1, s4)
call dqsub (s4, s3, s1)
call dqmul (b, b, s0)
call dqmul (b(3), b(3), s3)
call dqadd (s0, s3, s4)
call dqdiv (f, s4, s0)
call dqmul (s2, s0, c)
call dqmul (s1, s0, c(3))

return
end subroutine dqcdiv

subroutine dqceq (a, b)

!   This sets the DQC number B equal to the DQC number A.

implicit none
real (dqknd), intent(in):: a(4)
real (dqknd), intent(out):: b(4)

b(1) = a(1)
b(2) = a(2)
b(3) = a(3)
b(4) = a(4)

return
end subroutine dqceq

subroutine dqcmul (a, b, c)

!   This routine multiplies the DQC numbers A and B to yield the DQC result.

implicit none
real (dqknd), intent(in):: a(4), b(4)
real (dqknd), intent(out):: c(4)
real (dqknd) s0(2), s1(2), s2(2), s3(2)

call dqmul (a, b, s0)
call dqmul (a(3), b(3), s1)
call dqmul (a, b(3), s2)
call dqmul (a(3), b, s3)
call dqsub (s0, s1, c)
call dqadd (s2, s3, c(3))

return
end subroutine dqcmul

subroutine dqcpr (a, b, ic)

!   This routine compares the DQR numbers A and B and returns in IC the value
!   -1, 0, or 1 depending on whether A < B, A = B, or A > B. It is faster
!   than merely subtracting A and B and looking at the sign of the result.

implicit none
real (dqknd), intent(in):: a(2), b(2)
integer, intent(out):: ic

if (a(1) < b(1)) then
  ic = -1
elseif (a(1) == b(1)) then
  if (a(2) < b(2)) then
    ic = -1
  elseif (a(2) == b(2)) then
    ic = 0
  else
    ic = 1
  endif
else
  ic = 1
endif

return
end subroutine dqcpr

subroutine dqcpwr (a, n, b)

!   This computes the N-th power of the DQC number A and returns the DQC
!   result C in B. When N is zero, 1 is returned. When N is negative, the
!   reciprocal of A ^ |N| is returned.

!   This routine employs the binary method for exponentiation.

implicit none
real (dqknd), intent(in):: a(4)
integer, intent(in):: n
real (dqknd), intent(out):: b(4)
real (dqknd), parameter:: cl2 = 1.4426950408889634073599246810018921374q0
integer j, kk, kn, l1, mn, na1, na2, nn
real (dqknd) t1
real (dqknd) s0(4), s1(4), s2(4), s3(4)

if (a(1) == 0.q0 .and. a(3) == 0.q0) then
  if (n >= 0) then
    b(1) = 0.q0
    b(2) = 0.q0
    b(3) = 0.q0
    b(4) = 0.q0
    goto 120
  else
    write (dqldb, 1)
1   format ('*** DQCPWR: Argument is zero and N is negative or zero.')
    call dqabrt
    return
  endif
endif

nn = abs (n)
if (nn == 0) then
  s2(1) = 1.q0
  s2(2) = 0.q0
  s2(3) = 0.q0
  s2(4) = 0.q0
  goto 120
elseif (nn == 1) then
  s2(1) = a(1)
  s2(2) = a(2)
  s2(3) = a(3)
  s2(4) = a(4)
  goto 110
elseif (nn == 2) then
  call dqcmul (a, a, s2)
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN > NN.

t1 = nn
mn = cl2 * log (t1) + 1.q0 + 1.q-14

s0(1) = a(1)
s0(2) = a(2)
s0(3) = a(3)
s0(4) = a(4)
s2(1) = 1.q0
s2(2) = 0.q0
s2(3) = 0.q0
s2(4) = 0.q0
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn /= 2 * kk) then
    call dqcmul (s2, s0, s1)
    s2(1) = s1(1)
    s2(2) = s1(2)
    s2(3) = s1(3)
    s2(4) = s1(4)
  endif
    kn = kk
  if (j < mn) then
    call dqcmul (s0, s0, s1)
    s0(1) = s1(1)
    s0(2) = s1(2)
    s0(3) = s1(3)
    s0(4) = s1(4)
  endif
enddo

!   Compute reciprocal if N is negative.

110  continue

if (n < 0) then
  s1(1) = 1.q0
  s1(2) = 0.q0
  s1(3) = 0.q0
  s1(4) = 0.q0
  call dqcdiv (s1, s2, s0)
  s2(1) = s0(1)
  s2(2) = s0(2)
  s2(3) = s0(3)
  s2(4) = s0(4)
endif

b(1) = s2(1)
b(2) = s2(2)
b(3) = s2(3)
b(4) = s2(4)

120  continue
return
end subroutine dqcpwr

subroutine dqcsqrt (a, b)

!   This routine computes the complex square root of the DQC number C.
!   This routine uses the following formula, where A1 and A2 are the real and
!   imaginary parts of A, and where R = Sqrt [A1 ^ 2 + A2 ^2]:

!      B = Sqrt [(R + A1) / 2] + I Sqrt [(R - A1) / 2]

!   If the imaginary part of A is < 0, then the imaginary part of B is also
!   set to be < 0.

implicit none
real (dqknd), intent(in):: a(4)
real (dqknd), intent(out):: b(4)
real (dqknd) s0(2), s1(2), s2(2)

if (a(1) == 0.q0 .and. a(3) == 0.q0) then
  b(1) = 0.q0
  b(2) = 0.q0
  b(3) = 0.q0
  b(4) = 0.q0
  goto 100
endif

call dqmul (a, a, s0)
call dqmul (a(3), a(3), s1)
call dqadd (s0, s1, s2)
call dqsqrt (s2, s0)

s1(1) = a(1)
s1(2) = a(2)
if (s1(1) < 0.q0) then
  s1(1) = - s1(1)
  s1(2) = - s1(2)
endif
call dqadd (s0, s1, s2)
call dqmuld (s2, 0.5q0, s1)
call dqsqrt (s1, s0)
call dqmuld (s0, 2.q0, s1)
if (a(1) >= 0.q0) then
  b(1) = s0(1)
  b(2) = s0(2)
  call dqdiv (a(3), s1, b(3))
else
  call dqdiv (a(3), s1, b)
  if (b(1) < 0.q0) then
    b(1) = - b(1)
    b(2) = - b(2)
  endif
  b(3) = s0(1)
  b(4) = s0(2)
  if (a(3) < 0.q0) then
    b(3) = - b(3)
    b(4) = - b(4)
  endif
endif

 100  continue
return
end subroutine dqcsqrt

subroutine dqcsshr (a, x, y)

!   This computes the hyperbolic cosine and sine of the DQR number A and
!   returns the two DQR results in X and Y, respectively. 

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: x(2), y(2)
real (dqknd) f(2), s0(2), s1(2), s2(2)

f(1) = 1.q0
f(2) = 0.q0
call dqexp (a, s0)
call dqdiv (f, s0, s1)
call dqadd (s0, s1, s2)
call dqmuld (s2, 0.5q0, x)
call dqsub (s0, s1, s2)
call dqmuld (s2, 0.5q0, y)

return
end subroutine dqcsshr

subroutine dqcssnr (a, x, y)

!   This computes the cosine and sine of the DQR number A and returns the
!   two DQR results in X and Y, respectively.

!   This routine uses the conventional Taylor series for Sin (s):

!   Sin (s) =  s - s^3 / 3! + s^5 / 5! - s^7 / 7! ...

!   where the argument S has been reduced to (-pi, pi). To further accelerate
!   convergence of the series, the reduced argument is divided by 2^NQ. After
!   convergence, the double-angle formulas for cos are applied NQ times.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: x(2), y(2)
integer, parameter:: itrmx = 1000, nq = 5
real (dqknd), parameter:: eps = 1.q-70
real (dqknd) f1(2), f2(2), s0(2), s1(2), s2(2), s3(2), s4(2), s5(2)
integer is, i1, j, na, n1
real (dqknd) d1, t1, t2
real (dqknd), parameter:: pi(1:2) = &
  [3.1415926535897932384626433832795027974791q+00, &
  8.6718101301237810247970440260433494722008q-35]

! End of declaration

na = 2
if (a(1) == 0.d0) na = 0
if (na == 0) then
  x(1) = 1.d0
  x(2) = 0.d0
  y(1) = 0.d0
  y(2) = 0.d0
  goto 120
endif

!   Set f1 = 1 and f2 = 1/2.

f1(1) = 1.d0
f1(2) = 0.d0
f2(1) = 0.5d0
f2(2) = 0.d0

!   Check if argument is too large to compute meaningful cos/sin values.

if (a(1) >= 1.q60) then
  write (dqldb, 3)
3 format ('*** DQCSSNR: argument is too large to compute cos or sin.')
  call dqabrt
endif

!   Reduce to between - Pi and Pi.

call dqmuld (pi, 2.q0, s0)
call dqdiv (a, s0, s1)
call dqnint (s1, s2)
call dqmul (s0, s2, s4)
call dqsub (a, s4, s3)

!   Check if reduced argument is zero. If so then cos = 1 and sin = 0.

if (s3(1) == 0.d0) then
  x(1) = 1.d0
  x(2) = 0.d0
  y(1) = 0.d0
  y(2) = 0.d0
  goto 120
endif

call dqdivd (s3, 2.q0**nq, s0)
s1(1) = s0(1); s1(2) = s0(2)

!   Compute the sin of the reduced argument of s1 using a Taylor series.

call dqmul (s0, s0, s2)
if (s0(1) < 0.d0) then
  is = -1
else
  is = 1
endif

do i1 = 1, itrmx
  t2 = - (2.d0 * i1) * (2.d0 * i1 + 1.d0)
  call dqmul (s2, s1, s3)
  call dqdivd (s3, t2, s1)
  call dqadd (s1, s0, s3)
!  call mpeq (s3, s0, mpnw1)
  s0(1) = s3(1); s0(2) = s3(2)

!   Check for convergence of the series, and adjust working precision
!   for the next term.

  if (abs (s1(1)) < eps) goto 110
enddo

write (dqldb, 4)
4 format ('*** DQCSSNR: Iteration limit exceeded.')
call dqabrt

110 continue

!   Apply the formula cos(2*x) = 2*cos^2(x) - 1 NQ times to produce
!   the cosine of the reduced argument, except that the first iteration is
!   cos(2*x) = 1 - 2*sin^2(x), since we have computed sin(x) above.
!   Note that these calculations are performed as 2 * (cos^2(x) - 1/2) and
!   2 * (1/2 - sin^2(x)), respectively, to avoid loss of precision.

call dqmul (s0, s0, s4)
call dqsub (f2, s4, s5)
call dqmuld (s5, 2.q0, s0)

do j = 2, nq
  call dqmul (s0, s0, s4)
  call dqsub (s4, f2, s5)
  call dqmuld (s5, 2.q0, s0)
enddo

!   Compute sin of result and correct sign.

call dqmul (s0, s0, s4)
call dqsub (f1, s4, s5)
call dqsqrt (s5, s1)

if (is < 1) then
  s1(1) = - s1(1); s1(2) = - s1(2)
endif

115 continue

x(1) = s0(1); x(2) = s0(2)
y(1) = s1(1); y(2) = s1(2)

120 continue

return
end subroutine dqcssnr

subroutine dqcsub (a, b, c)

!   This subracts the DQC numbers A and B and returns the DQC difference in C.

implicit none
real (dqknd) a(4), b(4), c(4)

call dqsub (a, b, c)
call dqsub (a(3), b(3), c(3))

return
end subroutine dqcsub

character(32) function dqdigout (a, n)

!   This converts the double precision input A to a character(32) string of
!   nonblank length N. A must be a whole number, and N must be sufficient
!   to hold it. This is intended for internal use only.

implicit none
integer, intent(in):: n
real (dqknd), intent(in):: a
character(10), parameter:: digits = '0123456789'
real (dqknd) d1, d2
character(32) ca
integer i, k

! End of declaration

ca = ' '
d1 = abs (a)

do i = n, 1, -1
  d2 = aint (d1 / 10.d0)
  k = 1.d0 + (d1 - 10.d0 * d2)
  d1 = d2
  ca(i:i) = digits(k:k)
enddo

dqdigout = ca
return
end function dqdigout

subroutine dqdiv (dqa, dqb, dqc)

!   This divides the DQR number DQA by the DQR number DQB to yield the DQ
!   quotient DQC.

implicit none
real (dqknd), intent(in):: dqa(2), dqb(2)
real (dqknd), intent(out):: dqc(2)
real (dqknd), parameter:: split = 144115188075855873.q0
real (dqknd) a1, a2, b1, b2, cona, conb, c11, c2, c21, e, s1, s2, &
  t1, t2, t11, t12, t21, t22

!   Compute a DQR approximation to the quotient.

s1 = dqa(1) / dqb(1)

!   This splits s1 and dqb(1) into high-order and low-order words.

cona = s1 * split
conb = dqb(1) * split
a1 = cona - (cona - s1)
b1 = conb - (conb - dqb(1))
a2 = s1 - a1
b2 = dqb(1) - b1

!   Multiply s1 * dqb(1) using Dekker's method.

c11 = s1 * dqb(1)
c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2
!>
!   Compute s1 * dqb(2) (only high-order word is needed).

c2 = s1 * dqb(2)

!   Compute (c11, c21) + c2 using Knuth's trick.

t1 = c11 + c2
e = t1 - c11
t2 = ((c2 - e) + (c11 - (t1 - e))) + c21

!   The result is t1 + t2, after normalization.

t12 = t1 + t2
t22 = t2 - (t12 - t1)

!   Compute dqa - (t12, t22) using Knuth's trick.

t11 = dqa(1) - t12
e = t11 - dqa(1)
t21 = ((-t12 - e) + (dqa(1) - (t11 - e))) + dqa(2) - t22

!   Compute high-order word of (t11, t21) and divide by dqb(1).

s2 = (t11 + t21) / dqb(1)

!   The result is s1 + s2, after normalization.

dqc(1) = s1 + s2
dqc(2) = s2 - (dqc(1) - s1)

return
end subroutine dqdiv

subroutine dqdivd (dqa, db, dqc)

!   This routine divides the DQR number A by the QP number B to yield
!   the DQR quotient C. DB must be an integer or exact binary fraction,
!   such as 3., 0.25 or -0.375.

implicit none
real (dqknd), intent(in):: dqa(2), db
real (dqknd), intent(out):: dqc(2)
real (dqknd), parameter:: split = 144115188075855873.q0
real (dqknd) a1, a2, b1, b2, cona, conb, e, t1, t2, t11, t12, t21, t22

!   Compute a DP approximation to the quotient.

t1 = dqa(1) / db
!>
!   On systems with a fused multiply add, such as IBM systems, it is faster to
!   uncomment the next two lines and comment out the following lines until !>.
!   On other systems, do the opposite.

! t12 = t1 * db
! t22 = t1 * db - t12

!   This splits t1 and db into high-order and low-order words.

cona = t1 * split
conb = db * split
a1 = cona - (cona - t1)
b1 = conb - (conb - db)
a2 = t1 - a1
b2 = db - b1

!   Multiply t1 * db using Dekker's method.

t12 = t1 * db
t22 = (((a1 * b1 - t12) + a1 * b2) + a2 * b1) + a2 * b2
!>
!   Compute dqa - (t12, t22) using Knuth's trick.

t11 = dqa(1) - t12
e = t11 - dqa(1)
t21 = ((-t12 - e) + (dqa(1) - (t11 - e))) + dqa(2) - t22

!   Compute high-order word of (t11, t21) and divide by db.

t2 = (t11 + t21) / db

!   The result is t1 + t2, after normalization.

dqc(1) = t1 + t2
dqc(2) = t2 - (dqc(1) - t1)
return
end subroutine dqdivd

subroutine dqdpdqc (a, b)

!   This routine converts the QP number A to DQR form in B. A must
!   be an integer or exact binary fraction, such as 3., 0.25 or -0.375.

implicit none
real(dqknd), intent(in):: a
real (dqknd), intent(out):: b(2)

b(1) = a
b(2) = 0.q0
return
end subroutine dqdpdqc

subroutine dqdqdpc (a, b)

!   This converts the DQR number A to DP.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b

b = a(1)
return
end subroutine dqdqdpc

subroutine dqeq (a, b)

!   This routine sets the DQR number B equal to the DQR number A. 

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2)

b(1) = a(1)
b(2) = a(2)
end subroutine dqeq

real (dqknd) function dqdigin (ca, n)
implicit none
real (dqknd) d1
character(*), intent(in):: ca
integer, intent(in):: n
character(10), parameter:: digits = '0123456789'
integer i, k

d1 = 0.q0

do i = 1, n
  k = index (digits, ca(i:i)) - 1
  if (k < 0) then
    write (dqldb, *) 'DQDIGIN: non-digit in character string'
  elseif (k <= 9) then
    d1 = 10.q0 * d1 + k
  endif
enddo

dqdigin = d1
end function dqdigin

subroutine dqdmc (a, n, b)

!   This converts the QP number A * 2^N to DQR form in B.

!   NOTE however that the product is not fully accurate unless A is an exact
!   binary value.
!   Examples of exact binary values (good): 123456789.d0, 0.25d0, -5.3125d0.
!   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.

implicit none
integer, intent(in):: n
real (dqknd), intent(in):: a
real (dqknd), intent(out):: b(2)
integer i, k, n1, n2
real (dqknd) aa

b(1) = a * 2.q0 ** n
b(2) = 0.q0
return
end subroutine dqdmc

subroutine dqeformat (a, nb, nd, b)

!   Converts the DQR number A into character form in the character(1) array B.
!   NB (input) is the length of the output string, and ND (input) is the
!   number of digits after the decimal point. The format is analogous to
!   Fortran E format. The result is left-justified among the NB cells of B.
!   The condition NB >= ND + 8 must hold or an error message will result.
!   NB cells must be available in array B.

implicit none
real (dqknd), intent(in):: a(2)
integer, intent(in):: nb, nd
character(1), intent(out):: b(nb)
character(10), parameter:: digits = '0123456789'
integer, parameter:: ndpw = 32
real (dqknd), parameter:: d10w = 10.q0**ndpw
integer i, ia, ix, ixp, i1, i2, j, k, na, nexp, nl
character(1) b2(nb+50)
character(32) ca
real (dqknd) aa, an, f(2), s0(2), s1(2), s2(2), s3(2), t1, t2

! End of declaration

if (nb < nd + 8) then
  write (dqldb, 1)
1 format ('*** DQEFORMAT: uninitialized or inadequately sized arrays')
  call dqabrt
endif

ia = sign (1.q0, a(1))
if (a(1) == 0.d0) ia = 0
na = 2
if (ia == 0) na = 0

!   Set f = 10.

f(1) = 10.d0; f(2) = 0.d0

!   Determine power of ten for exponent, and scale input to within 1 and 10.

if (ia < 0) then
  s1(1) = - a(1); s1(2) = - a(2)
else
  s1(1) = a(1); s1(2) = a(2)
endif

if (na > 0) then
  aa = s1(1)
  t1 = log10 (aa)

  if (t1 >= 0.d0) then
    nexp = t1
  else
    nexp = t1 - 1.d0
  endif

  if (nexp == 0) then
  elseif (nexp > 0) then
    call dqnpwr (f, nexp, s0)
    call dqdiv (s1, s0, s2)
    s1(1) = s2(1); s1(2) = s2(2)
  elseif (nexp < 0) then
    call dqnpwr (f, -nexp, s0)
    call dqmul (s1, s0, s2)
    s1(1) = s2(1); s1(2) = s2(2)
  endif

!   If we didn't quite get it exactly right, multiply or divide by 10 to fix.

100 continue

  if (s1(1) < 1.d0) then
    nexp = nexp - 1
    call dqmuld (s1, 10.q0, s0)
    s1(1) = s0(1); s1(2) = s0(2)
    goto 100
  elseif (s1(1) >= 10.d0) then
    nexp = nexp + 1
    call dqdivd (s1, 10.q0, s0)
    s1(1) = s0(1); s1(2) = s0(2)
    goto 100
  endif
else
  nexp = 0
endif

!   Insert sign and first digit.

ix = 0
if (ia == -1) then
  ix = ix + 1
  b2(ix) = '-'
endif
if (na > 0) then
  call dqinfr (s1, s2, s3)
  an = s2(1)
else
  an = 0.d0
endif
ca = dqdigout (an, 1)
ix = ix + 1
b2(ix) = ca(1:1)
ix = ix + 1
b2(ix) = '.'
ixp = ix

!   Set f = an.

f(1) = an
f(2) = 0.d0
call dqsub (s1, f, s0)
call dqmuld (s0, d10w, s1)

!   Calculate the number of remaining chunks.

nl = nd / ndpw + 1

!   Insert the digits of the remaining words.

do j = 1, nl
  if (s1(1) /= 0.d0) then
    call dqinfr (s1, s2, s3)
    an = s2(1)
    f(1) = an
    f(2) = 0.d0
  else
    an = 0.d0
    f(1) = 0.d0
    f(2) = 0.d0
  endif
  
  ca = dqdigout (an, ndpw)

  do i = 1, ndpw
    ix = ix + 1
    if (ix > nb + 50) then
      write (dqldb, 2)
2     format ('DQEFORMAT: Insufficient space in B2 array.')
      call dqabrt
    endif
    b2(ix) = ca(i:i)
  enddo

  call dqsub (s1, f, s0)
  call dqmuld (s0, d10w, s1)
enddo

!   Round the result.

if (ix >= nd + 1) then
  i1 = index (digits, b2(nd+1)) - 1
  if (i1 >= 5) then

!   Perform rounding, beginning at the last digit (position IX). If the rounded
!   digit is 9, set to 0, then repeat at position one digit to left. Continue
!   rounding if necessary until the decimal point is reached.

    do i = ix, ixp + 1, -1
      i2 = index (digits, b2(i)) - 1
      if (i2 <= 8) then
        b2(i) = digits(i2+2:i2+2)
        goto 180
      else
        b2(i) = '0'
      endif
    enddo

!   We have rounded up all digits to the right of the decimal point. If the
!   digit to the left of the decimal point is a 9, then set that digit to 1
!   and increase the exponent by one; otherwise increase that digit by one.

    if (b2(ixp-1) == '9') then
      b2(ixp-1) = '1'
      nexp = nexp + 1
    else
      i1 = index (digits, b2(ixp-1)) - 1
      b2(ixp-1) = digits(i1+2:i1+2)
    endif
  endif
endif

180 continue

!   Done with mantissa. Insert exponent.

ix = nd + 2
if (ia < 0) ix = ix + 1
b2(ix) = 'e'
if (nexp < 0) then
  ix = ix + 1
  b2(ix) = '-'
endif
ca = dqdigout (real (abs (nexp), dqknd), 10)

do k = 1, 10
  if (ca(k:k) /= '0') goto 190
enddo

k = 10

190 continue

do i = k, 10
  ix = ix + 1
  b2(ix) = ca(i:i)
enddo

do i = ix + 1, nb
  b2(i) = ' '
enddo

!   Copy entire b2 array to B.

do i = 1, nb
  b(i) = b2(i)
enddo

return
end subroutine dqeformat

subroutine dqegamc (egam)

!   This returns the value of the Euler gamma constant.

implicit none
real (dqknd) egam(2)

egam(1) = dqegammacon(1); egam(2) = dqegammacon(2)
return
end subroutine dqegamc

subroutine dqexp (a, b)

!   This computes the exponential function of the DQR number A and returns the
!   DQR result in B.

!   This routine uses a modification of the Taylor's series for Exp (t):

!   Exp (t) =  (1 + r + r^2 / 2! + r^3 / 3! + r^4 / 4! ...) ^ q * 2 ^ n

!   where q = 64, r = t' / q, t' = t - n Log(2) and where n is chosen so
!   that -0.5 Log(2) < t' <= 0.5 Log(2). Reducing t mod Log(2) and
!   dividing by 64 insures that -0.004 < r <= 0.004, which accelerates
!   convergence in the above series.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2)
integer, parameter:: nq = 6
integer i, ia, l1, na, nz, n1
real (dqknd) t1, t2
real (dqknd) eps, f(2), s0(2), s1(2), s2(2), s3(2), tl
real (dqknd), parameter:: al2(1:2) = &
  [6.9314718055994530941723212145817657508364q-01, &
  -7.0081394745495851634126620087716238869662q-36]

!   Check for overflows and underflows.

eps = 10.q0 ** (-70)
if (abs (a(1)) >= 11000.q0) then
  if (a(1) > 0.q0) then
    write (dqldb, 1) a(1)
1   format ('*** DQEXP: Argument is too large',f12.6)
    call dqabrt
    return
  else
    call dqdpdqc (0.q0, b)
    goto 130
  endif
endif

f(1) = 1.q0
f(2) = 0.q0

!   Compute the reduced argument A' = A - Log(2) * Nint [A / Log(2)]. Save
!   NZ = Nint [A / Log(2)] for correcting the exponent of the final result.

call dqdiv (a, al2, s0)
call dqnint (s0, s1)
t1 = s1(1)
nz = t1 + sign (1.q-14, t1)
call dqmul (al2, s1, s2)
call dqsub (a, s2, s0)

!   Check if the reduced argument is zero.

if (s0(1) == 0.q0) then
  s0(1) = 1.q0
  s0(2) = 0.q0
  l1 = 0
  goto 120
endif

!   Divide the reduced argument by 2 ^ NQ.

call dqdivd (s0, 2.q0 ** nq, s1)

!   Compute Exp using the usual Taylor series.

s2(1) = 1.q0
s2(2) = 0.q0
s3(1) = 1.q0
s3(2) = 0.q0
l1 = 0

100  l1 = l1 + 1
if (l1 == 100) then
  write (dqldb, 2)
2 format ('*** DQEXP: Iteration limit exceeded.')
  call dqabrt
  return
endif

t2 = l1
call dqmul (s2, s1, s0)
call dqdivd (s0, t2, s2)
call dqadd (s3, s2, s0)
call dqeq (s0, s3)

!   Check for convergence of the series.

if (abs (s2(1)) > eps * abs (s3(1))) goto 100

!   Raise to the (2 ^ NQ)-th power.

do i = 1, nq
  call dqmul (s0, s0, s1)
  s0(1) = s1(1)
  s0(2) = s1(2)
enddo

!  Multiply by 2 ^ NZ.

120  call dqmuld (s0, 2.q0 ** nz, b)

!   Restore original precision level.

 130  continue
return
end subroutine dqexp

subroutine dqfformat (a, nb, nd, b)

!   Converts the DQR number A into character form in the character(1) array B.
!   NB (input) is the length of the output string, and ND (input) is the
!   number of digits after the decimal point. The format is analogous to
!   Fortran F format; the result is right-justified among the NB cells of B.
!   The condition NB >= ND + 8 must hold or an error message will result.
!   However, if it is found during execution that there is not sufficient space,
!   to hold all digits, the entire output field will be filled with asterisks.
!   NB cells of type character(1) must be available in B.

implicit none
integer, intent(in):: nb, nd
real (dqknd), intent(in):: a(2)
character(1), intent(out):: b(nb)
integer i, ixp, i1, i2, i3, j, k, nb2, nb3, nexp
character(1) b2(nb+20)
character(32) ca
real (dqknd) t1

! End of declaration

if (nb < nd + 8) then
  write (dqldb, 1)
1 format ('*** DQFFORMAT: uninitialized or inadequately sized arrays')
  call dqabrt
endif

!   Call dqeformat with sufficient precision.

if (a(1) == 0.q0) then
  nb2 = nd + 11
else
  nb2 = max (log10 (abs (a(1))), 0.q0) + nd + 11
endif
nb3 = nb2 - 8
call dqeformat (a, nb2, nb3, b2)

!   Trim off trailing blanks.

do i = nb2, 1, -1
  if (b2(i) /= ' ') goto 90
enddo

90 continue

nb2 = i

!   Look for the 'e' in B2.

do k = 1, nb2
  if (b2(k) == 'e') goto 100
enddo

write (dqldb, 2)
2 format ('*** DQFFORMAT: Syntax error in output of dqeformat')
call dqabrt

100 continue

!   Check the sign of the exponent.

k = k + 1
if (b2(k) == '-') then
  ixp = -1
  k = k + 1
else
  ixp = 1
endif
j = 0
ca = ' '

!   Copy the exponent into CA.

do i = k, nb2
  j = j + 1
  if (j <= 16) ca(j:j) = b2(i)
enddo

t1 = dqdigin (ca, j)

!   Check if there is enough space in the output array for all digits.

if (ixp == 1 .and. t1 + nd + 3 > nb) then
  do i = 1, nb
    b(i) = '*'
  enddo

  goto 210
endif
nexp = ixp * t1

!   Insert the sign of the number, if any.

i1 = 0
i2 = 0
if (b2(1) == '-') then
  i1 = i1 + 1
  b(i1) = '-'
  i2 = i2 + 1
endif

if (nexp == 0) then

!   Exponent is zero. Copy first digit, period and ND more digits.

  do i = 1, nd + 2
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo

  goto 200
elseif (nexp > 0) then

!   Exponent is positive. Copy first digit, skip the period, then copy
!   nexp digits.

  i1 = i1 + 1
  i2 = i2 + 1
  b(i1) = b2(i2)
  i2 = i2 + 1

  do i = 1, nexp
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo

!   Insert the period.

  i1 = i1 + 1
  b(i1) = '.'

!   Copy nd more digits.

  do i = 1, nd
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo

  goto 200
else

!   Exponent is negative. Insert a zero, then a period, then -1 - nexp
!   zeroes, then the first digit, then the remaining digits up to ND total
!   fractional digits.

  i1 = i1 + 1
  b(i1) = '0'
  i1 = i1 + 1
  b(i1) = '.'
  i3 = min (- nexp - 1, nd - 1)

  do i = 1, i3
    i1 = i1 + 1
    b(i1) = '0'
  enddo

  if (- nexp - 1 < nd) then
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
    i2 = i2 + 1
  endif

  do i = i3 + 2, nd
    i1 = i1 + 1
    i2 = i2 + 1
    b(i1) = b2(i2)
  enddo
endif

200 continue

!   Right-justify in field.

k = nb - i1

do i = 1, i1
  b(nb-i+1) = b(nb-i-k+1)
enddo

do i = 1, k
  b(i) = ' '
enddo

210 continue

return
end subroutine dqfformat

subroutine dqinfr (a, b, c)

!   Sets B to the integer part of the DQR number A and sets C equal to the
!   fractional part of A. Note that if A = -3.3, then B = -3 and C = -0.3.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2), c(2)
real (dqknd), parameter:: t225 = 2.q0 ** 225, t112 = 2.q0 ** 112
integer ic
real (dqknd) con(2), f(2), s0(2), s1(2)
save con
data con / t225, t112/

!   Check if  A  is zero.

if (a(1) == 0.q0)  then
  b(1) = 0.q0
  b(2) = 0.q0
  c(1) = 0.q0
  c(2) = 0.q0
  goto 120
endif

if (a(1) >= t225) then
  write (dqldb, 1)
1 format ('*** DQINFR: Argument is too large.')
  call dqabrt
  return
endif

f(1) = 1.q0
f(2) = 0.q0
if (a(1) > 0.q0) then
  call dqadd (a, con, s0)
  call dqsub (s0, con, b)
  call dqcpr (a, b, ic)
  if (ic >= 0) then
    call dqsub (a, b, c)
  else
    call dqsub (b, f, s1)
    b(1) = s1(1)
    b(2) = s1(2)
    call dqsub (a, b, c)
  endif
else
  call dqsub (a, con, s0)
  call dqadd (s0, con, b)
  call dqcpr (a, b, ic)
  if (ic <= 0) then
    call dqsub (a, b, c)
  else
    call dqadd (b, f, s1)
    b(1) = s1(1)
    b(2) = s1(2)
    call dqsub (a, b, c)
  endif
endif

120  continue
return
end subroutine dqinfr

subroutine dqinp (iu, a)

!   This routine reads the DQR number A from logical unit IU. The input
!   value must be placed on a single line of not more than 120 characters.

implicit none
integer, intent(in):: iu
real (dqknd), intent(out):: a(2)
integer, parameter:: ln = 120
character(120) cs

read (iu, '(a)', end = 100) cs
call dqinpc (cs, a)
goto 110

100 write (dqldb, 1)
1  format ('*** DQINP: End-of-file encountered.')
call dqabrt

110 return
end subroutine dqinp

subroutine dqinpc (a, b)

!   Converts the CHARACTER(*) array A into the DQR number B.

implicit none
character(*), intent(in):: a
real (dqknd), intent(out):: b(2)
character(10), parameter:: dig = '0123456789'
integer i, ib, id, ie, inz, ip, is, ix, k, ln, lnn
real (dqknd) bi
character(1) ai
character(10) ca
real (dqknd) f(2), s0(2), s1(2), s2(2)

id = 0
ip = -1
is = 0
inz = 0
s1(1) = 0.q0
s1(2) = 0.q0
ln = len (a)

do i = ln, 1, -1
  if (a(i:i) /= ' ') goto 90
enddo

90 continue

lnn = i

!   Scan for digits, looking for the period also.

do i = 1, lnn
  ai = a(i:i)
  if (ai == ' ' .and. id == 0) then
  elseif (ai == '.') then
    if (ip >= 0) goto 210
    ip = id
    inz = 1
  elseif (ai == '+') then
    if (id /= 0 .or. ip >= 0 .or. is /= 0) goto 210
    is = 1
  elseif (ai == '-') then
    if (id /= 0 .or. ip >= 0 .or. is /= 0) goto 210
    is = -1
  elseif (ai == 'e' .or. ai == 'E' .or. ai == 'd' .or. ai == 'D') then
    goto 100
  elseif (index (dig, ai) == 0) then
    goto 210
  else
!    read (ai, '(f1.0)') bi
    bi = index (dig, ai) - 1
    if (inz > 0 .or. bi > 0.q0) then
      inz = 1
      id = id + 1
      call dqmuld (s1, 10.q0, s0)
      f(1) = bi
      f(2) = 0.q0
      call dqdpdqc (bi, f)
      call dqadd (s0, f, s1)
    endif
  endif
enddo

100   continue
if (is == -1) then
  s1(1) = - s1(1)
  s1(2) = - s1(2)
endif
k = i
if (ip == -1) ip = id
ie = 0
is = 0
ca = ' '

do i = k + 1, lnn
  ai = a(i:i)
  if (ai == ' ') then
  elseif (ai == '+') then
    if (ie /= 0 .or. is /= 0) goto 210
    is = 1
  elseif (ai == '-') then
    if (ie /= 0 .or. is /= 0) goto 210
    is = -1
  elseif (index (dig, ai) == 0) then
    goto 210
  else
    ie = ie + 1
    if (ie > 3) goto 210
    ca(ie:ie) = ai
  endif
enddo

! read (ca, '(i4)') ie

if (ca == ' ') then
  ie = 0
else
  ie = dqdigin (ca, ie)
endif

if (is == -1) ie = - ie
ie = ie + ip - id
s0(1) = 10.q0
s0(2) = 0.q0
call dqnpwr (s0, ie, s2)
call dqmul (s1, s2, b)
goto 220

210  write (dqldb, 1)
1 format ('*** DQINPC: Syntax error in literal string.')
call dqabrt

220  continue

return
end subroutine dqinpc

subroutine dqlog (a, b)

!   This computes the natural logarithm of the DQR number A and returns the DQ
!   result in B.

!   The Taylor series for Log converges much more slowly than that of Exp.
!   Thus this routine does not employ Taylor series, but instead computes
!   logarithms by solving Exp (b) = a using the following Newton iteration,
!   which converges to b:

!           x_{k+1} = x_k + [a - Exp (x_k)] / Exp (x_k)

implicit none
integer k
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2)
real (dqknd) t1, t2
real (dqknd) s0(2), s1(2), s2(2)

if (a(1) <= 0.q0) then
  write (dqldb, 1)
1 format ('*** DQLOG: Argument is less than or equal to zero.')
  call dqabrt
  return
endif

!   Compute initial approximation of Log (A).

t1 = a(1)
t2 = log (t1)
b(1) = t2
b(2) = 0.q0

!   Perform the Newton-Raphson iteration described above.

do k = 1, 3
  call dqexp (b, s0)
  call dqsub (a, s0, s1)
  call dqdiv (s1, s0, s2)
  call dqadd (b, s2, s1)
  b(1) = s1(1)
  b(2) = s1(2)
enddo

120  continue

return
end subroutine dqlog

subroutine dqlog2c (alog2d)

!   This returns log(2) to DQR precision.

implicit none
real (dqknd), intent(out):: alog2d(2)
real (dqknd) alog2c(2)
save alog2c
data alog2c/ 6.9314718055994530941723212145817657508364q-01, &
  -7.0081394745495851634126620087716238869662q-36/

alog2d(1) = alog2c(1)
alog2d(2) = alog2c(2)
return
end subroutine dqlog2c

subroutine dqqqc (a, b, c)

!   This converts DQR numbers A and B to DQC form in C, i.e. C = A + B i.

implicit none
real (dqknd), intent(in):: a(2), b(2)
real (dqknd), intent(out):: c(4)

c(1) = a(1)
c(2) = a(2)
c(3) = b(1)
c(4) = b(2)
return
end subroutine dqqqc

subroutine dqmdc (a, b, n)

!   This returns a DP approximation the MPR number A in the form B * 2^n.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b
integer, intent(out):: n

b = a(1)
n = 0
return
end subroutine dqmdc

subroutine dqmul (dqa, dqb, dqc)

!   This routine multiplies DQR numbers DQA and DQB to yield the DQR product DQC.

implicit none
real (dqknd), intent(in):: dqa(2), dqb(2)
real (dqknd), intent(out):: dqc(2)
real (dqknd), parameter:: split = 144115188075855873.q0
real (dqknd) a1, a2, b1, b2, cona, conb, c11, c21, c2, e, t1, t2

!   This splits dqa(1) and dqb(1) into high-order and low-order words.

cona = dqa(1) * split
conb = dqb(1) * split
a1 = cona - (cona - dqa(1))
b1 = conb - (conb - dqb(1))
a2 = dqa(1) - a1
b2 = dqb(1) - b1

!   Multilply dqa(1) * dqb(1) using Dekker's method.

c11 = dqa(1) * dqb(1)
c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2
!>
!   Compute dqa(1) * dqb(2) + dqa(2) * dqb(1) (only high-order word is needed).

c2 = dqa(1) * dqb(2) + dqa(2) * dqb(1)

!   Compute (c11, c21) + c2 using Knuth's trick, also adding low-order product.

t1 = c11 + c2
e = t1 - c11
t2 = ((c2 - e) + (c11 - (t1 - e))) + c21 + dqa(2) * dqb(2)

!   The result is t1 + t2, after normalization.

dqc(1) = t1 + t2
dqc(2) = t2 - (dqc(1) - t1)

return
end subroutine dqmul

subroutine dqmuld (dqa, db, dqc)

!   This routine multiplies the DQR number DQA by the QP number DB to yield
!   the DQR product DQC. DB must be an integer or exact binary fraction,
!   such as 3., 0.25 or -0.375.

implicit none
real (dqknd), intent(in):: dqa(2), db
real (dqknd), intent(out):: dqc(2)
real (dqknd), parameter:: split = 144115188075855873.q0
real (dqknd) a1, a2, b1, b2, cona, conb, c11, c21, c2, e, t1, t2

!   This splits dqa(1) and db into high-order and low-order words.

cona = dqa(1) * split
conb = db * split
a1 = cona - (cona - dqa(1))
b1 = conb - (conb - db)
a2 = dqa(1) - a1
b2 = db - b1

!   Multilply dqa(1) * db using Dekker's method.

c11 = dqa(1) * db
c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2
!>
!   Compute dqa(2) * db (only high-order word is needed).

c2 = dqa(2) * db

!   Compute (c11, c21) + c2 using Knuth's trick.

t1 = c11 + c2
e = t1 - c11
t2 = ((c2 - e) + (c11 - (t1 - e))) + c21

!   The result is t1 + t2, after normalization.

dqc(1) = t1 + t2
dqc(2) = t2 - (dqc(1) - t1)
return
end subroutine dqmuld

subroutine dqmuldd (da, db, dqc)

!   This subroutine computes DQC = DA x DB, where DA and DB are quad and DDC is
!   quad-double.

implicit none
real (dqknd), intent(in):: da, db
real (dqknd), intent(out):: dqc(2)
real (dqknd), parameter:: split = 144115188075855873.q0
real (dqknd) a1, a2, b1, b2, cona, conb, s1, s2

!>
!   On systems with a fused multiply add, such as IBM systems, it is faster to
!   uncomment the next two lines and comment out the following lines until !>.
!   On other systems, do the opposite.

! s1 = da * db
! s2 = da * db - s1

!   This splits da and db into high-order and low-order words.

cona = da * split
conb = db * split
a1 = cona - (cona - da)
b1 = conb - (conb - db)
a2 = da - a1
b2 = db - b1

!   Multiply da * db using Dekker's method.

s1 = da * db
s2 = (((a1 * b1 - s1) + a1 * b2) + a2 * b1) + a2 * b2
!>
dqc(1) = s1
dqc(2) = s2

return
end subroutine dqmuldd

subroutine dqmzc (a, b)

!  This converts the DQR real variable A to the DQC variable B.
!  This routine is not intended to be called directly by the user.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(4)

b(1) = a(1)
b(2) = a(2)
b(3) = 0.q0
b(4) = 0.q0
return
end subroutine dqmzc

subroutine dqneg (a, b)

!   This sets b = -a.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2)

b(1) = -a(1); b(2) = -a(2)
return
end

subroutine dqnint (a, b)

!   This sets B equal to the integer (type DQR) nearest to the DQR number A.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2)
real (dqknd), parameter:: t225 = 2.q0 ** 225, t112 = 2.q0 ** 112
real (dqknd) con(2), s0(2)
save con
data con / t225, t112/


!   Check if  A  is zero.

if (a(1) == 0.q0)  then
  b(1) = 0.q0
  b(2) = 0.q0
  goto 120
endif

if (a(1) >= t225) then
  write (dqldb, 1)
1 format ('*** DQNINT: Argument is too large.')
  call dqabrt
  return
endif

if (a(1) > 0.q0) then
  call dqadd (a, con, s0)
  call dqsub (s0, con, b)
else
  call dqsub (a, con, s0)
  call dqadd (s0, con, b)
endif

120  continue
return
end subroutine dqnint

subroutine dqnpwr (a, n, b)

!   This computes the N-th power of the DQR number A and returns the DQR result
!   in B. When N is zero, 1 is returned. When N is negative, the reciprocal
!   of A ^ |N| is returned. 

!   This routine employs the binary method for exponentiation.

implicit none
real (dqknd), intent(in):: a(2)
integer, intent(in):: n
real (dqknd), intent(out):: b(2)
real (dqknd), parameter:: cl2 = 1.4426950408889634073599246810018921374q0
integer j, kk, kn, l1, mn, na1, na2, nn
real (dqknd) t1
real (dqknd) s0(2), s1(2), s2(2), s3(2)

if (a(1) == 0.q0) then
  if (n >= 0) then
    s2(1) = 0.q0
    s2(2) = 0.q0
    goto 120
  else
    write (dqldb, 1)
1   format ('*** DQNPWR: Argument is zero and N is negative or zero.')
    call dqabrt
    return
  endif
endif

nn = abs (n)
if (nn == 0) then
  s2(1) = 1.q0
  s2(2) = 0.q0
  goto 120
elseif (nn == 1) then
  s2(1) = a(1)
  s2(2) = a(2)
  goto 110
elseif (nn == 2) then
  call dqmul (a, a, s2)
  goto 110
endif

!   Determine the least integer MN such that 2 ^ MN > NN.

t1 = nn
mn = cl2 * log (t1) + 1.q0 + 1.d-14
s0(1) = a(1)
s0(2) = a(2)
s2(1) = 1.q0
s2(2) = 0.q0
kn = nn

!   Compute B ^ N using the binary rule for exponentiation.

do j = 1, mn
  kk = kn / 2
  if (kn /= 2 * kk) then
    call dqmul (s2, s0, s1)
    s2(1) = s1(1)
    s2(2) = s1(2)
  endif
  kn = kk
  if (j < mn) then
    call dqmul (s0, s0, s1)
    s0(1) = s1(1)
    s0(2) = s1(2)
  endif
enddo

!   Compute reciprocal if N is negative.

110  continue

if (n < 0) then
  s1(1) = 1.q0
  s1(2) = 0.q0
  call dqdiv (s1, s2, s0)
  s2(1) = s0(1)
  s2(2) = s0(2)
endif

120  continue

b(1) = s2(1)
b(2) = s2(2)
  
return
end subroutine dqnpwr

subroutine dqnrtf (a, n, b)

!   This computes the N-th root of the DQR number A and returns the DQR result
!   in B. N must be at least one.

!   This subroutine employs the following Newton-Raphson iteration, which
!   converges to A ^ (-1/N):

!    X_{k+1} = X_k + (X_k / N) * (1 - A * X_k^N)

!   The reciprocal of the final approximation to A ^ (-1/N) is the N-th root.

implicit none
real (dqknd), intent(in):: a(2)
integer, intent(in):: n
real (dqknd), intent(out):: b(2)
integer i, k
real (dqknd) t1, t2, tn
real (dqknd) f1(2), s0(2), s1(2)

if (a(1) == 0.q0) then
  b(1) = 0.q0
  b(2) = 0.q0
  goto 140
elseif (a(1) < 0.q0) then
  write (dqldb, 1)
1 format ('*** DQNRT: Argument is negative.')
  call dqabrt
  return
endif
if (n <= 0) then
  write (dqldb, 2) n
2 format ('*** DQNRT: Improper value of N',i10)
  call dqabrt
  return
endif

!   Handle cases N = 1 and 2.

if (n == 1) then
  b(1) = a(1)
  b(2) = a(1)
  goto 140
elseif (n == 2) then
  call dqsqrt (a, b)
  goto 140
endif

f1(1) = 1.q0
f1(2) = 0.q0

!   Compute the initial approximation of A ^ (-1/N).

tn = n
t1 = a(1)
t2 = exp (- log (t1) / tn)
b(1) = t2
b(2) = 0.q0

!   Perform the Newton-Raphson iteration described above.

do k = 1, 3
  call dqnpwr (b, n, s0)
  call dqmul (a, s0, s1)
  call dqsub (f1, s1, s0)
  call dqmul (b, s0, s1)
  call dqdivd (s1, tn, s0)
  call dqadd (b, s0, s1)
  b(1) = s1(1)
  b(2) = s1(2)
enddo

!   Take the reciprocal to give final result.

call dqdiv (f1, b, s1)
b(1) = s1(1)
b(2) = s1(2)

140  continue
return
end subroutine dqnrtf

subroutine dqout (iu, nb, nd, a)

!   This routine writes the DQR number A on logical unit IU using a standard
!   Fortran Enb.nd format.

implicit none
integer, intent(in):: iu, nb, nd
real (dqknd), intent(in):: a(2)
integer i
character(1) cs(nb)

call dqeformat (a, nb, nd, cs)
write (iu, '(120a1)') (cs(i), i = 1, nb)

return
end subroutine dqout

subroutine dqpic (pi)

!   This returns Pi to DQR precision.

implicit none
real (dqknd), intent(out):: pi(2)

pi(1) = dqpicon(1); pi(2) = dqpicon(2)
return
end subroutine dqpic

subroutine dqpolyr (n, a, x0, x)

!   This finds the root X, near X0 (input) for the nth-degree DQR polynomial
!   whose coefficients are given in the (n+1)-long vector A. It may be
!   necessary to adjust eps -- default value is 1.q-65.

implicit none
integer, intent(in):: n
real (dqknd), intent(in):: a(2,0:n), x0(2)
real (dqknd), intent(out):: x(2)
real (dqknd), parameter:: eps = 1.q-65
integer i, it
real (dqknd)  ad(2,0:n), t1(2), t2(2), t3(2), t4(2), t5(2), dt1

do i = 0, n - 1
  dt1 = i + 1
  call dqmuld (a(1,i+1), dt1, ad(1,i))
enddo

ad(1,n) = 0.q0
ad(2,n) = 0.q0
x(1) = x0(1)
x(2) = x0(2)

do it = 1, 20
  t1(1) = 0.q0
  t1(2) = 0.q0
  t2(1) = 0.q0
  t2(2) = 0.q0
  t3(1) = 1.q0
  t3(2) = 0.q0

  do i = 0, n
    call dqmul (a(1,i), t3, t4)
    call dqadd (t1, t4, t5)
    t1(1) = t5(1)
    t1(2) = t5(2)
    call dqmul (ad(1,i), t3, t4)
    call dqadd (t2, t4, t5)
    t2(1) = t5(1)
    t2(2) = t5(2)
    call dqmul (t3, x, t4)
    t3(1) = t4(1)
    t3(2) = t4(2)
  enddo

  call dqdiv (t1, t2, t3)
  call dqsub (x, t3, t4)
  x(1) = t4(1)
  x(2) = t4(2)
  if (abs (t3(1)) <= eps) goto 110
enddo

write (dqldb, 1)
1 format ('DQPOLYR: failed to converge.')
call dqabrt

110 continue

return
end subroutine dqpolyr

integer function dqsgn (ra)

!   This function returns 1, 0 or -1, depending on whether ra > 0, ra = 0 or ra < 0.

implicit none
real (dqknd), intent(in):: ra(2)
integer ia

if (ra(1) == 0.q0) then
  dqsgn = 0
elseif (ra(1) > 0.q0) then
  dqsgn = 1
else
  dqsgn = -1
endif
return
end function dqsgn

subroutine dqpower (a, b, c)

!   This computes C = A^B.

implicit none
real (dqknd), intent(in):: a(2), b(2)
real (dqknd), intent(out):: c(2)
real (dqknd) t1(2), t2(2)

if (a(1) <= 0.q0) then
  write (6, 1)
1 format ('DQPOWER: A <= 0')
  call dqabrt
endif

call dqlog (a, t1)
call dqmul (t1, b, t2)
call dqexp (t2, c)
return
end subroutine dqpower

subroutine dqsqrt (a, b)

!   This computes the square root of the DQR number A and returns the DQR result
!   in B.

!   This subroutine employs the following formula (due to Alan Karp):

!          Sqrt(A) = (A * X) + 0.5 * [A - (A * X)^2] * X  (approx.)

!   where X is a double precision approximation to the reciprocal square root,
!   and where the multiplications A * X and [] * X are performed with only
!   double precision.

implicit none
real (dqknd), intent(in):: a(2)
real (dqknd), intent(out):: b(2)
real (dqknd) t1, t2, t3
real (dqknd) f(2), s0(2), s1(2)

if (a(1) == 0.q0) then
  b(1) = 0.q0
  b(2) = 0.q0
  goto 100
endif
t1 = 1.q0 / sqrt (a(1))
t2 = a(1) * t1
call dqmuldd (t2, t2, s0)
call dqsub (a, s0, s1)
t3 = 0.5q0 * s1(1) * t1
s0(1) = t2
s0(2) = 0.q0
s1(1) = t3
s1(2) = 0.q0
call dqadd (s0, s1, b)

100 continue

return
end subroutine dqsqrt

subroutine dqsub (dqa, dqb, dqc)

!   This subroutine computes DDC = DDA - DDB, where all args are DQR.

implicit none
real (dqknd), intent(in):: dqa(2), dqb(2)
real (dqknd), intent(out):: dqc(2)
real (dqknd) e, t1, t2

!   Compute dqa + dqb using Knuth's trick.

t1 = dqa(1) - dqb(1)
e = t1 - dqa(1)
t2 = ((-dqb(1) - e) + (dqa(1) - (t1 - e))) + dqa(2) - dqb(2)

!   The result is t1 + t2, after normalization.

dqc(1) = t1 + t2
dqc(2) = t2 - (dqc(1) - t1)
return
end subroutine dqsub

subroutine dqxzc (a, b)

!  This converts the DC variable A to the DQC variable B.
!  This routine is not intended to be called directly by the user.

implicit none
complex (dqknd), intent(in):: a
real (dqknd), intent(out):: b(4)

b(1) = a
b(2) = 0.q0
b(3) = aimag (a)
b(4) = 0.q0
return
end subroutine dqxzc

end module dqfuna
