program tpphixdq

!   This computes the Poisson phi function, as described in paper by Bailey,
!   Borwein, Crandall and Zucker, and finds the associated minimal polynomial.
!   This version uses the DQFUN double-quad precision package.

!   David H Bailey   23 Feb 2023

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!   Parameters set in parameter statements below:
!     ipal  0: The palindromic technique (for even kd) is NOT implemented.
!           1: The palindromic technique (for even kd) IS implemented.
!           NOTE: If ipal = 1, the parameter n below = 1 + degree of half
!           polynomial. In other words, n = size of relation in call to PSLQM1.
!           The full output polynomial has degree 2*n-2.
!     mpi   Multiplier of pi on RHS of conjectured identity; default = 8.
!     kd    Denominator of arguments of phi_2; default = 12.
!     kp    Numerator of first argument; default = 1.
!     kq    Numerator of second argument; default = 3.
!     lcx   Size of character(1) array chr1; default = 64.
!     lmx   Size of line1; default = 1024.

!   PSLQM1 parameters:
!     idb   Debug level (0 - 3); default = 2.
!     n     Integer relation vector length = 1 + polynomial degree; default = 9.
!     ndp   Full precision level in digits; default = 70.
!     nep   Log10 of full precision epsilon for detections; default = 3 - ndp.
!           ***Must not be smaller than the accuracy of input data. 
!     rb    Log10 of max size (Euclidean norm) of acceptable relation;
!           default = 100. Run will abort if this is exceeded.

use dqmodule
implicit none
integer, parameter:: dpknd = kind (0.d0), idb = 2, ipal = 1, mpi = 8, kd = 12, &
  kp = 1, kq = 3, n = 9, ndp = 70, m1 = 40, m2 = 2, nep = 3 - ndp, lcx = 64, &
  lmx = 1024
real (dpknd), parameter:: rb = 100.d0
integer i, iq, i1, i2, j, j1, k, l1, m, n1, nn, lnm(n), lnm2(2*n)
character(1) chr1(lcx)
character(64) form4, nam(n), nam2(2*n), namx
character(lmx) line1
real (dpknd) d1, d2, second, tm, tm0, tm1, tm2, tm3
external second
type (dq_real) alpha, beta, eps, pi, qq, xx, yy
type (dq_real) epsm
type (dq_complex) zz, c1, c2, c3, c4, theta1, theta2, theta3, theta4
external theta1, theta2, theta3, theta4
type (dq_real) r(n), r2(2*n), x(n)
save

write (6, 3) ipal, kd, kp, kq, mpi, n, nn, ndp, nep
3 format ('Poisson Phi_2 computation and analysis:'/ &
  'ipal =', i4,'; kd =',i4,'; kp =',i4,'; kq =',i4,'; mpi =',i4,'; n =',i4/ &
  'ndp =',i6,'; nep =',i6)

nn = n
pi = dqpi ()
eps = dqreal (10.q0) ** nep
xx = dqreal (real (kp, dqknd)) / dqreal (real (kd, dqknd))
yy = dqreal (real (kq, dqknd)) / dqreal (real (kd, dqknd))

!   Evaluate formula (43) of Section 5.

tm0 = second ()
qq = exp (- pi)
zz = 0.5q0 * pi * dqcmplx (yy, xx)
c1 = theta1 (zz, qq, eps)
c2 = theta2 (zz, qq, eps)
c3 = theta3 (zz, qq, eps)
c4 = theta4 (zz, qq, eps)
alpha = abs (c2 * c4 / (c1 * c3)) ** (mpi / 2)
tm1 = second ()
write (6, 4) tm1 - tm0
4 format ('Alpha CPU time =',f12.2/'Alpha =')
call dqwrite (6, ndp + 8, ndp, alpha)

!   Construct input x vector for PSLQM4.

if (ipal == 0) then

!   Case: palindromic property is NOT used.

  x(1) = 1.q0

  do i = 2, nn
    x(i) = alpha * x(i-1)
  enddo

  nam(1) = '1'
  lnm(1) = 1

  do i = 2, nn
    if (i <= 10) then
      write (namx, '("al^",i1)') i - 1
      nam(i) = namx(1:4)
      lnm(i) = 4
    elseif (i <= 100) then
      write (namx, '("al^",i2)') i - 1
      nam(i) = namx(1:5)
      lnm(i) = 5
    else
      stop
    endif
  enddo

  do i = 1, nn
    lnm2(i) = lnm(i)
    nam2(i) = nam(i)
  enddo
elseif (ipal == 1) then

!   Case: palindromic property IS used.

  beta = alpha + 1.q0 / alpha
  x(1) = 1.q0

  do i = 2, nn
    x(i) = beta * x(i-1)
  enddo

  nam(1) = '1'
  lnm(1) = 1

  do i = 2, nn
    if (i <= 10) then
      write (namx, '("beta^",i1)') i - 1
      nam(i) = namx(1:6)
      lnm(i) = 6
    elseif (i <= 100) then
      write (namx, '("beta^",i2)') i - 1
      nam(i) = namx(1:7)
      lnm(i) = 7
    else
      stop
    endif
  enddo

  nam2(1) = '1'
  lnm2(1) = 1

  do i = 2, 2 * nn
    if (i <= 10) then
      write (namx, '("al^",i1)') i - 1
      nam2(i) = namx(1:4)
      lnm2(i) = 4
    elseif (i <= 100) then
      write (namx, '("al^",i2)') i - 1
      nam2(i) = namx(1:5)
      lnm2(i) = 5
    else
      stop
    endif
  enddo
endif

!   Perform relation search.

tm2 = second ()
call pslqm1 (idb, n, rb, eps, x, iq, r)
tm3 = second ()
write (6, 5) tm3 - tm2, tm3 - tm0
5 format ('PSLQM1 CPU time =',f12.2/'Total CPU time =',f12.2)

!   Output relation in two formats.

if (iq == 1) then

!   Produce format used below for output.

  form4 = '(64a1'
  i1 = 5

  do i = 65, lcx, 64
    form4(i1+1:i1+9) =  '/64a1'
    i1 = i1 + 5
  enddo

  form4(i1+1:i1+9) = '," * ",a)'
  i1 = i1 + 9

!  If palindromic property is used, output recovered half relation.

  if (ipal == 1) then
    write (6, 6)
6   format (/'Recovered half relation: 0 =')

    do i = 1, nn
      if (r(i) /= 0.q0) then
        call dqfform (r(i), lcx, 0, chr1)

        do j = 1, lcx
          if (chr1(j) /= ' ') goto 110
        enddo

110     continue

        j1 = j
        if (chr1(j1) /= '-') then
          if (j1 > 1) then
            j1 = j1 - 1
            chr1(j1) = '+'
          else
            do j = 1, lcx
              chr1(j) = '*'
            enddo
          endif
        endif
        write (6, form4) (chr1(j), j = 1, lcx), nam(i)(1:lnm(i))
      endif
    enddo
  endif

!    If palindromic property is not used, we are done; if it is used, then
!    call doublep to expand to full polynomial.

  if (ipal == 0) then
    n1 = nn

    do i = 1, nn
      r2(i) = r(i)
    enddo
  elseif (ipal == 1) then
    n1 = 2 * nn - 1
    call doublep (2*nn - 1, r, r2)
  endif

  write (6, 7)
7 format (/'Recovered full relation: 0 =')
  l1 = 0

!   Output full polynomial in Mathematica form.

  do i = 1, n1
    if (r2(i) /= 0.q0) then
      call dqfform (r2(i), lcx, 0, chr1)

      do j = 1, lcx
        if (chr1(j) /= ' ') goto 120
      enddo

120   continue

      j1 = j
      if (chr1(j1) /= '-') then
        if (j1 > 1) then
          j1 = j1 - 1
          chr1(j1) = '+'
        else
          do j = 1, lcx
            chr1(j) = '*'
          enddo
        endif
      endif
      write (6, form4) (chr1(j), j = 1, lcx), nam2(i)(1:lnm2(i))

      if (l1 + 100 > lmx) then
        write (6, '("Mathematica polynomial is too long; increase lmx =",i8)') lmx
        stop
      endif
      line1(l1+1:l1+1) = ' '
      l1 = l1 + 1
      k = lcx - j1

      do j = 1, k
        line1(l1+j:l1+j) = chr1(j+j1-1)
      enddo

      l1 = l1 + k
      line1(l1+1:l1+lnm2(i)+1) = '*' // nam2(i)(1:lnm2(i))
      l1 = l1 + lnm2(i) + 1
    endif
  enddo

  i1 = 1
  write (6, 8)
8 format ('Output polynomial in Mathematica notation:')

130 continue

  if (i1 + 64 > l1) goto 140
  i2 = index (line1(i1+65:l1), ' ')
  if (i2 == 0) goto 140
  i2 = i1 + 64 + i2
!  write (6, '(a)') line1(i1:i2) // '\ '
  write (6, '(a)') line1(i1:i2)
  i1 = i2 + 1
  goto 130

140 continue

  write (6, '(a)') line1(i1:l1)
endif

if (iq == 1) then
  write (6, '(a)') 'TEST PASSED'
else
  write (6, '(a)') 'TEST FAILED'
endif
stop
end program tpphixdq

subroutine doublep (n, polyh, poly)
use dqmodule
implicit none
integer, intent(in):: n
type (dq_real), intent(in):: polyh(0:n/2)
type (dq_real), intent(out):: poly(-n/2:n/2)
type (dq_real) x(-n/2-1:n/2+1), y(-n/2-1:n/2+1)
integer i, j, k

do i = -n/2-1, n/2+1
  x(i) = dqreal (0.q0)
enddo

x(0) = polyh(n/2)

do j = 0, n / 2 - 1
  y(-j-1) = x(-j)
  y(j+1) = x(j)
  
  do k = -j, j
    y(k) = x(k-1) + x(k+1)
  enddo
  
  y(0) = y(0) + polyh(n/2-j-1)

  do k = -j-1, j+1
    x(k) = y(k)
  enddo
enddo

do k = -n/2, n/2
  poly(k) = x(k)
enddo

return
end subroutine doublep

type (dq_complex) function theta1 (zz, qq, eps)
use dqmodule
implicit none
type (dq_complex), intent(in):: zz
type (dq_real), intent(in):: qq, eps
integer k1, k2, n
real (dqknd) ds
type (dq_real) r1
type (dq_complex) t0, t1, t2, t3, t4

t0 = dqcmplx (cmplx (0.q0, 0.q0, dqknd))
r1 = sqrt (sqrt (qq))
ds = -1.q0

do k1 = 1, 10001, 2
  ds = - ds
  k2 = k1 ** 2
  t2 = ds * r1 ** k2 * sin (real (k1, dqknd) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, *) 'theta1: loop end error'
stop

100 continue

theta1 = 2.q0 * t0
return
end function theta1

type (dq_complex) function theta2 (zz, qq, eps)
use dqmodule
implicit none
type (dq_complex), intent(in):: zz
type (dq_real), intent(in):: qq, eps
integer k1, k2, n
real (dqknd) ds
type (dq_real) r1
type (dq_complex) t0, t1, t2, t3, t4

t0 = dqcmplx (cmplx (0.q0, 0.q0, dqknd))
r1 = sqrt (sqrt (qq))

do k1 = 1, 10001, 2
  k2 = k1 ** 2
  t2 = r1 ** k2 * cos (real (k1, dqknd) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, *) 'theta2: loop end error'
stop

100 continue

theta2 = 2.q0 * t0
return
end function theta2

type (dq_complex) function theta3 (zz, qq, eps)
use dqmodule
implicit none
type (dq_complex), intent(in):: zz
type (dq_real), intent(in):: qq, eps
integer k1, k2, n
real (dqknd) ds
type (dq_real) r1
type (dq_complex) t0, t1, t2, t3, t4

t0 = dqcmplx (cmplx (0.q0, 0.q0, dqknd))
r1 = qq

do k1 = 1, 10000
  k2 = k1 ** 2
  t2 = r1 ** k2 * cos (2.q0 * real (k1, dqknd) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, *) 'theta3: loop end error'
stop

100 continue

theta3 = cmplx (1.q0, 0.q0, dqknd) + 2.q0 * t0
return
end function theta3

type (dq_complex) function theta4 (zz, qq, eps)
use dqmodule
implicit none
type (dq_complex), intent(in):: zz
type (dq_real), intent(in):: qq, eps
integer k1, k2, n
real (dqknd) ds
type (dq_real) r1
type (dq_complex) t0, t1, t2, t3, t4

t0 = dqcmplx (cmplx (0.q0, 0.q0, dqknd))
r1 = qq
ds = 1.q0

do k1 = 1, 10000
  ds = - ds
  k2 = k1 ** 2
  t2 = ds * r1 ** k2 * cos (2.q0 * real (k1, dqknd) * zz)
  t0 = t0 + t2
  if (abs (t2) < eps) goto 100
enddo

write (6, *) 'theta4: loop end error'
stop

100 continue

theta4 = cmplx (1.q0, 0.q0, dqknd) + 2.q0 * t0
return
end function theta4

!------------------------------

!   The following code performs the one-level, multi-pair PSLQ algorithm.
!   David H. Bailey   23 Feb 2023

subroutine pslqm1 (idb, n, rb, eps, x, iq, r)

!   Arguments are as follows:
!     Name  Type    Description
!     idb   int     Debug flag (0-3); increasing idb produces more output.
!     n     int     Length of input vector x and output relation r.
!     rb    DP      Log10 of max size (Euclidean norm) of acceptable relation.
!                   NOTE: Run is aborted when this is exceeded.
!     eps   dq_real Tolerance for full precision relation detection. 
!     x     dq_real Input dq_real vector.
!     iq    int     Output flag: 0 (unsuccessful) or 1 (successful).
!     r     dq_real Output integer relation vector, if successful.

!   The following parameters are set in this routine:
!     ipi   int     Iteration print interval when idb >= 2; default = 25.
!     ipm   int     Iteration  check interval for iterations; default = 100.
!     itm   int     Maximum iteration count; default = 10^5.
!                   NOTE: Run is aborted when this is exceeded. If itm >= 10^7,
!                     change all "i7" in formats to "i8" or as high as needed.
!     nsq   int     Size of tables used in iterdp and itermp; default = 8.

use dqmodule
implicit none
integer, intent(in):: idb, n
integer, parameter:: dpknd = kind (0.d0)
real (dpknd), intent(in):: rb
type (dq_real), intent(in):: eps, x(n)
integer, intent(out):: iq
type (dq_real), intent(out):: r(n)
integer, parameter:: ipi = 25, ipm = 100, itm = 100000, nsq = 8
integer i, imq, it, izm, j, j1, n1, n2, n3, n4
real (dpknd) d1, d2, d3, d4, tm0, tm1, times(2)
type (dq_real) b(n,n), h(n,n), s(n), syq(n,nsq), y(n), rn, t1, t2, t3, t4
type (dq_real), external:: bound
real (dpknd), external:: dplog10, second

!   Initialize.

if (idb >= 2) write (6, 1) n
1 format ('PSLQM1 integer relation detection: n =',i5)
iq = 0
it = 0
imq = 0
rn = dqreal (0.q0)

do i = 1, 2
  times(i) = 0.q0
enddo

if (idb >= 2) write (6, 2) it
2 format ('Iteration',i7,3x,'Initialization')
tm0 = second ()
call initmp (idb, n, nsq, b, h, syq, x, y)
tm1 = second ()
times(1) = tm1 - tm0

!   Main iterations.

if (idb >= 2) write (6, 3) it
3 format ('Iteration',i7,3x,'Start iterations')

100 continue

!   Perform one iteration.

it = it + 1
if (idb == 3 .or. idb >= 2 .and. mod (it, ipi) == 0) write (6, 4) it
4 format ('Iteration',i7)
tm0 = second ()
call itermp (idb, it, n, nsq, eps, b, h, syq, y, imq, izm)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

!   Test conditions on itermp output flag izm:
!   0: Update was uneventful; periodically output min, max, norm; continue.
!   1: A small value was found in update; go output relation.
!   2: Precision is exhausted; quit.

if (izm == 0) then

!   Periodically output min, max and norm.

  if (mod (it, ipm) == 0) then

!   Find min and max absolute value of y vector.

    call minmax (n, y, t1, t2)
    if (idb >= 2) then
      call decmd (t1, d1, n1)
      call decmd (t2, d2, n2)
      write (6, 5) it, d1, n1, d2, n2
5     format ('Iteration',i7,3x,'Min, max of y =',f11.6,'e',i6,f11.6,'e',i6)
    endif

!   Compute norm bound.

    t1 = bound (n, h)
    rn = max (rn, t1)
    if (idb >= 2) then
      call decmd (t1, d1, n1)
      call decmd (rn, d2, n2)
      write (6, 6) it, d1, n1, d2, n2
6     format ('Iteration',i7,3x,'Norm bound =',f11.6,'e',i5,4x,'Max. bound =', &
        f11.6,'e',i5)
    endif

!   Check if iteration limit or norm bound limit is exceeded; if so, quit.

    if (it > itm) then
      if (idb >= 1) write (6, 7) itm
7     format ('Iteration limit exceeded',i7)
      goto 120
    endif
    if (dplog10 (rn) > rb) then
      if (idb >= 1) write (6, 8) rb
8     format ('Norm bound limit exceeded.',1p,d15.6)
      goto 120
    endif
  endif
  goto 100
elseif (izm == 1) then
  goto 110
elseif (izm == 2) then
  goto 120
endif

110 continue

!   A relation has been detected.

tm0 = second ()
t1 = dqreal (1.q300)
t2 = dqreal (0.q0)
t3 = dqreal (0.q0)

!   Select the relation corresponding to the smallest y entry and compute norm.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

do i = 1, n
  r(i) = b(j1,i)
  t3 = t3 + r(i) ** 2
enddo

t3 = sqrt (t3)
call decmd (t3, d3, n3)
d3 = d3 * 10.q0 ** n3

!   Output the final norm bound and other info.

if (idb >= 1) then
  t4 = bound (n, h)
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  call decmd (t3, d3, n3)
  call decmd (t4, d4, n4)
  write (6, 9) it, d1, n1, d2, n2, d4, n4
9 format ('Iteration',i7,3x,'Relation detected'/ &
  'Min, max of y =',0p,f11.6,'e',i5,f11.6,'e',i5/'Max. bound =',f11.6,'e',i5)
  write (6, 10) j1, d3, n3, d1, n1
10 format ('Index of relation =',i4,3x,'Norm =',f11.6,'e',i5,3x, &
  'Residual =',f11.6,'e',i5)
endif

!   If run was successful, set iq = 1.

if (dplog10 (t3) <= rb) then
  iq = 1
else
  if (idb >= 2) write (6, 11)
11 format ('Relation is too large.')
endif

120 continue

!   Output CPU run times and return.

if (idb >= 2) write (6, 12) times
12 format ('CPU times:'/(5f12.2))

return
end subroutine pslqm1

!------------------------------

!   First-level subroutines.

subroutine minmax (n, y, y1, y2)

!   This returns min|y_k| and max|y_k| using full precision.
!   Input: n, y.
!   Output: y1, y2.

use dqmodule
implicit none
integer, intent(in):: n
type (dq_real), intent(in):: y(n)
type (dq_real), intent(out):: y1, y2
integer i
type (dq_real) t1, t2, t3

t1 = dqreal (1.q300)
t2 = dqreal (0.q0)

!   Find the min and max absolute value in the y vector.

do i = 1, n
  t3 = abs (y(i))
  t1 = min (t1, t3)
  t2 = max (t2, t3)
enddo

y1 = t1
y2 = t2
return
end subroutine minmax

subroutine initmp (idb, n, nsq, b, h, syq, x, y)

!   This initializes DQ arrays at the beginning.
!   This is performed in full precision.
!   Input: idb, n, nsq, x.
!   Output: b, h, syq, y.

use dqmodule
implicit none
integer, intent(in):: idb, n, nsq
type (dq_real), intent(in):: x(n)
type (dq_real), intent(out):: b(n,n), h(n,n), syq(n,nsq), y(n)
integer, parameter:: dpknd = kind (0.d0)
integer i, i1, j
real (dpknd) d1
type (dq_real) s(n), t1, t2

if (idb >= 3) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, x)
endif

!   Set b to the identity matrix.

do j = 1, n
  do i = 1, n
    b(i,j) = dqreal (0.q0)
  enddo

  b(j,j) = dqreal (1.q0)
enddo

t1 = dqreal (0.q0)

!   Compute the s vector, the square root of the partial sum of squares of x,
!   and the y vector, which is the normalized x vector.

do i = n, 1, -1
  t1 = t1 + x(i) ** 2
  s(i) = sqrt (t1)
enddo

t1 = 1.q0 / s(1)

do i = 1, n
  y(i) = t1 * x(i)
  s(i) = t1 * s(i)
enddo

!   Compute the initial h matrix.

do j = 1, n - 1
  do i = 1, j - 1
    h(i,j) = dqreal (0.q0)
  enddo

  h(j,j) = s(j+1) / s(j)
  t1 = y(j) / (s(j) * s(j+1))

  do i = j + 1, n
    h(i,j) = - y(i) * t1
  enddo
enddo

!   Zero the syq array.

do j = 1, nsq
  do i = 1, n
    syq(i,j) = dqreal (0.q0)
  enddo
enddo

if (idb >= 3) then
  write (6, 2)
2 format ('initmp: Initial y vector:')
  call matoutmd (1, n, y)
  write (6, 3)
3 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, h)
endif

return
end subroutine initmp

subroutine itermp (idb, it, n, nsq, eps, b, h, syq, y, imq, izm)

!   This performs one iteration of the PSLQM algorithm using DQ arithmetic.
!   Input: idb, it, n, nsq, eps, b, h, syq, imq.
!   Output: b, h, syq, y, imq, izm.

use dqmodule
implicit none
integer, intent(in):: idb, it, n, nsq
integer, intent(inout):: imq
type (dq_real), intent(in):: eps
type (dq_real), intent(inout):: b(n,n), h(n,n), syq(n,nsq)
type (dq_real), intent(out):: y(n)
integer, intent(out):: izm
integer, parameter:: dpknd = kind (0.d0), ntl = 5
integer i, ii, ij, im, im1, j, j1, j2, k, mpr, mq, n1
real (dpknd) d1, d2
integer ip(n), ir(n), is(n)
type (dq_real) q(n), t(n,n)
type (dq_real) gam, t1, t2, t3, t4, teps

teps = 2.q0 ** ntl * eps
izm = 0
mpr = nint (0.4q0 * n)
gam = sqrt (dqreal (4.q0) / dqreal (3.q0))

!   Compute q vector = {gam^i * |h(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  q(i) = gam ** i * abs (h(i,i))
enddo

call qsortmp (n - 1, q, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest q(i).

do i = 1, n
  is(i) = 0
enddo

if (imq == 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii == 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) /= 0 .or. is(j2) /= 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of y, and rows of b and h.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = y(im)
  y(im) = y(im1)
  y(im1) = t1

  do i = 1, n
    t1 = b(im,i)
    b(im,i) = b(im1,i)
    b(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = h(im,i)
    h(im,i) = h(im1,i)
    h(im1,i) = t1
  enddo
enddo

!   Eliminate the "corners" produced by the above permutation in h.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im <= n - 2) then
    t1 = h(im,im)
    t2 = h(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = h(i,im)
      t4 = h(i,im1)
      h(i,im) = t1 * t3 + t2 * t4
      h(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo

!   Perform reduction on h, using the diagonal scheme.  Multipliers are
!   saved in the t array.

do i = 2, n
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      h(ij,j) = h(ij,j) - t(ij,k) * h(k,j)
    enddo

    t(ij,j) = anint (h(ij,j) / h(j,j))
    h(ij,j) = h(ij,j) - t(ij,j) * h(j,j)
  enddo
enddo

!   Update y, using the t array.  Find min absolute value of y.

t1 = abs (y(n))
j1 = n

do j = 1, n - 1
  do i = j + 1, n
    y(j) = y(j) + t(i,j) * y(i)
  enddo

  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
enddo

!   Update b, using the t array.

do k = 1, n
  do j = 1, n - 1
    do i = j + 1, n
      b(j,k) = b(j,k) + t(i,j) * b(i,k)
    enddo
  enddo
enddo

!  Find the largest entry of b in the same row as the smallest y.

t2 = dqreal (0.q0)

do i = 1, n
  t2 = max (t2, abs (b(j1,i)))
enddo

if (t1 <= t2 * teps) then
  if (idb >= 2) then
    call decmd (t1, d1, n1) 
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'itermp: Small value in y =',f11.6,'e',i5)
  endif
  if (t1 <= t2 * eps) then
    izm = 1
  else
    if (idb >= 1) write (6, 2) it
2   format ('Iteration',i7,3x,'itermp: Precision exhausted.')
    izm = 2
  endif
endif

!   Compare the y vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = dqreal (0.q0)

  do i = 1, n
    t1 = max (t1, abs (y(i) - syq(i,j)))
  enddo

  if (t1 <= t2 * teps) then
    if (idb >= 2) write (6, 3) it, j
 3  format ('Iteration',i7,3x,'itermp: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector y in the table syq.

120 continue

k = 1 + mod (it, nsq)

do i = 1, n
  syq(i,k) = y(i)
enddo

if (idb >= 3) then
  write (6, 4)
4 format ('itermp: Updated y:')
!  call matoutmd (1, n, y)
  call matoutmp (1, n, y)
  write (6, 5)
5 format ('itermp: Updated b matrix:')
  call matoutmd (n, n, b)
  write (6, 6)
6 format ('itermp: Updated h matrix:')
  call matoutmd (n, n - 1, h)
endif

return
end subroutine itermp

!------------------------------

!   Second- and third-level subroutines.

type (dq_real) function bound (n, h)

!   This computes the norm bound using DQ arithmetic.

use dqmodule
implicit none
integer, intent(in):: n
type (dq_real), intent(in):: h(n,n)
integer i
type (dq_real) t1, t2

t1 = dqreal (0.q0)

do i = 1, n - 1
  t1 = max (t1, abs (h(i,i)))
enddo

bound = 1.q0 / t1

return
end function bound

real (kind (0.d0)) function dplog10 (a)

!   For input DQ value a, this routine returns a DP approximation to log10 (a).

use dqmodule
implicit none
type (dq_real), intent(in):: a
integer, parameter:: dpknd = kind (0.d0)
real (dpknd) da, t1

da = qreal (a)
if (da == 0.d0) then
  dplog10 = -999999.d0
else
  dplog10 = log10 (abs (da))
endif

100 continue
return
end function dplog10

subroutine decmd (a, b, ib)

!   For input DQ value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use dqmodule
implicit none
integer ia, ib
integer, parameter:: dpknd = kind (0.d0)
real (dpknd) da, b, t1
type (dq_real) a

da = qreal (a)
if (da /= 0.d0) then
  t1 = log10 (abs (da))
  ib = t1
  if (t1 < 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end subroutine decmd

subroutine matoutmd (n1, n2, a)

!   This outputs the DQ matrix a as a DP matrix.

use dqmodule
implicit none
integer, intent(in):: n1, n2
type (dq_real), intent(in):: a(n1,n2)
integer i, j, ix(n2)
integer, parameter:: dpknd = kind (0.d0)
real (dpknd) dx(n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call decmd (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'e',i5))
enddo

return
end subroutine matoutmd

subroutine matoutmp (n1, n2, a)

!   This outputs the DQ matrix a.  It may be used in place of calls to matoutmd
!   in the code above if greater accuracy is desired in debug output.

use dqmodule
implicit none
integer, intent(in):: n1, n2
type (dq_real), intent(in):: a(n1,n2)
integer i, j
integer, parameter:: m1 = 80, m2 = 70

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call dqwrite (6, m1, m2, a(i,j))
  enddo
enddo

return
end subroutine matoutmp

subroutine qsortmp (n, a, ip)

!   This routine sorts the entries of the N-long DQ vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.
!   Input: n, a.
!   Output: ip.

use dqmodule
implicit none
integer, intent(in):: n
type (dq_real), intent(in):: a(n)
integer, intent(out):: ip(n)

integer i, iq, it, j, jq, jz, k, l
integer ik(50), jk(50)
type (dq_real) s0, s1, s2

do i = 1, n
  ip(i) = i
enddo

if (n == 1) return

k = 1
ik(1) = 1
jk(1) = n

130 continue

i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 < a(ip(l))) goto 160
enddo

i = j
goto 190

160 continue

i = l

do l = j, i, -1
  if (s0 > a(ip(l))) goto 180
enddo

j = i
goto 190

180 continue

j = l
if (i >= j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue

if (s0 >= a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 continue

k = k - 1
jz = 0
if (j == iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 continue

i = i + 1
if (i == jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz == 0) goto 220
if (j - iq >= jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 continue

if (k > 0) goto 130

return
end subroutine qsortmp
