!*****************************************************************************

!  program testqxfun

!  Revision date:  16 Mar 2023

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired)
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  DESCRIPTION OF THIS PROGRAM:
!   This briefly tests most individual QXFUN special functions, by comparing each
!   result with benchmark results in the file testqxfun.ref.txt, which must be
!   present in the same directory. This is not an exhaustive test of all possible
!   scenarios, but it often detects bugs and compiler issues.

program testqxfun
use qxmodule
implicit none
integer, parameter:: ibz = 3, nb = 60, nfile = 11, np = 2, nq = 3, nrr = 10
real (qxknd), parameter:: eps = 2.e0_qxknd**(ibz-qxnbt)
real (qxknd) aa(np), bb(nq), ber(1:nb), rr(nrr), err, errmx, pi, log2, t3, t4

!   End of declaration

open (nfile, file = 'testqxfun.ref.txt')
rewind nfile
write (6, '(a)') 'QXFUN: Quick check of functions'

errmx = 0.e0_qxknd
pi = qxpi ()
log2 = qxlog2 ()
aa(1) = 0.75e0_qxknd
aa(2) = 1.25e0_qxknd
bb(1) = 1.e0_qxknd
bb(2) = 1.5e0_qxknd
bb(3) = 2.e0_qxknd

call qxberne (nb, ber)
call polylog_ini (-nrr, rr)

write (6, '(a)') 'bessel_i(pi,log2) ='
t3 = bessel_i (pi, log2)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 2, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_in(3,pi) ='
t3 = bessel_in (3, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_in(3,50*pi) ='
t3 = bessel_in (3, 50.e0_qxknd * pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_j(pi,log2) ='
t3 = bessel_j (pi, log2)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_jn(3,pi) ='
t3 = bessel_jn (3, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_jn(3,50*pi) ='
t3 = bessel_jn (3, 50.e0_qxknd * pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_k(pi,log2) ='
t3 = bessel_k (pi, log2)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_kn(3,pi) ='
t3 = bessel_kn (3, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_kn(3,50*pi) ='
t3 = bessel_kn (3, 50.e0_qxknd * pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_y(pi,log2) ='
t3 = bessel_y (pi, log2)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_yn(3,pi) ='
t3 = bessel_yn (3, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_yn(3,50*pi) ='
t3 = bessel_yn (3, 50.e0_qxknd * pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'digamma_be(nb,ber,pi) ='
t3 = digamma_be (nb, ber, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'erf(pi) ='
t3 = erf (pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'erfc(pi) ='
t3 = erfc (pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'expint(pi) ='
t3 = expint (pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'gamma(pi) ='
t3 = gamma(pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'incgamma(pi,2log2) ='
t3 = incgamma (pi, 2.e0_qxknd - log2)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'incgamma(-3,pi) ='
t4 = -3.e0_qxknd
t3 = incgamma (t4, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'hurwitz_zetan(3,1/pi) ='
t4 = 1.e0_qxknd / pi
t3 = hurwitz_zetan (3, t4)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'hurwitz_zetan_be(nb,ber,5,pi) ='
t3 = hurwitz_zetan_be (nb, ber, 5, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'hypergeom_pfq(np,nq,aa,bb,pi) ='
t3 = hypergeom_pfq (np, nq, aa, bb, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'polygamma(3,1/pi) ='
t4 = 1.e0_qxknd / pi
t3 = polygamma (3, 1.e0_qxknd / pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'polygamma_be(nb,ber,5,pi) ='
t3 = polygamma_be (nb, ber, 5, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'polylog_neg (-10, rr, -pi) ='
t3 = polylog_neg (-10, rr, -pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'polylog_pos (10, 1/pi) ='
t4 = 1.e0_qxknd / pi
t3 = polylog_pos (10, 1.e0_qxknd / pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'struve_hn(10,pi) ='
t3 = struve_hn (10, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'zeta(pi) ='
t3 = zeta (pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'zeta_int(10) ='
t3 = zeta_int (10)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, '(a)') 'zeta_be(nb, ber, pi) ='
t3 = zeta_be (nb, ber, pi)
write (6, '(1p,d45.33)') t3
call checkqx (nfile, 1, t3, err)
errmx = max (err, errmx)

write (6, 9) errmx
9 format (/'Max relative error =',1p,d15.6)

if (abs (errmx) < eps) then
  write (6, '(a)') 'ALL TESTS PASSED'
else
  write (6, '(a)') 'ONE OR MORE TESTS FAILED'
endif

stop
end program testqxfun

subroutine checkqx (nfile, i1, t1, err)
use qxmodule
implicit none
integer, intent(in):: nfile, i1
real (qxknd), intent(in):: t1
real (qxknd), intent(out):: err
integer, parameter:: ibz = 3
real (qxknd), parameter:: eps = 2.e0_qxknd**(ibz-qxnbt)
real (qxknd) t2
integer i
character(64) c1

do i = 1, i1
  read (nfile, '(a)') c1
enddo

read (nfile, '(f50.0)') t2
err = abs ((t1 - t2) / t2)
if (abs (err) > eps) write (6, 1) err
1 format ('ERROR:',1pd15.6)

return
end subroutine checkqx

