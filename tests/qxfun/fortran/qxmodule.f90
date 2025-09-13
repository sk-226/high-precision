
!  QXFUN: A quad precision package with special functions

!  High-level language interface module (QXMODULE).

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

!  DESCRIPTION OF THIS MODULE (QXMODULE):
!    This module contains all high-level Fortran language interfaces.

module qxmodule
use qxfune
implicit none

!   These are subroutine names in module QXFUNE (qxfune.f90):

private &
  qxagmr, qxberner, qxpolyadd, qxpolysub, qxpolymul, &
  qxbesselir, qxbesselinr, qxbesseljr, qxbesseljnr, qxbesselkr, &
  qxbesselknr,  qxbesselyr, qxbesselynr, qxdigammabe, qxerfr, &
  qxerfcr, qxexpint, qxgammar, qxhurwitzzetan, qxhurwitzzetanbe, &
  qxhypergeompfq, qxincgammar, qxpolygamma, qxpolygammabe, &
  qxpolylogini, qxpolylogneg, qxpolylogpos, qxstruvehn, qxzetar, &
  qxzetaintr, qxzetabe

!   High-level interfaces for the functions and subroutines below.

interface agm
  module procedure qx_agm
end interface

interface qxberne
  module procedure qx_berne
end interface

interface bessel_i
  module procedure qx_bessel_i
end interface

interface bessel_in
  module procedure qx_bessel_in
end interface

interface bessel_j
  module procedure qx_bessel_j
end interface

interface bessel_jn
  module procedure qx_bessel_jn
end interface

interface bessel_j0
  module procedure qx_bessel_j0
end interface

interface bessel_j1
  module procedure qx_bessel_j1
end interface

interface bessel_k
  module procedure qx_bessel_k
end interface

interface bessel_kn
  module procedure qx_bessel_kn
end interface

interface bessel_y
  module procedure qx_bessel_y
end interface

interface bessel_yn
  module procedure qx_bessel_yn
end interface

interface bessel_y0
  module procedure qx_bessel_y0
end interface

interface bessel_y1
  module procedure qx_bessel_y1
end interface

interface digamma_be
  module procedure qx_digamma_be
end interface

interface qxcssh
  module procedure qx_cssh
end interface

interface qxcssn
  module procedure qx_cssn
end interface

interface qxlog2
  module procedure qx_log2
end interface

interface qxpi
  module procedure qx_pi
end interface

interface erf
  module procedure qx_erf
end interface

interface erfc
  module procedure qx_erfc
end interface

interface expint
  module procedure qx_expint
end interface

interface gamma
  module procedure qx_gamma
end interface

interface hurwitz_zetan
  module procedure qx_hurwitz_zetan
end interface

interface hurwitz_zetan_be
  module procedure qx_hurwitz_zetan_be
end interface

interface hypergeom_pfq
  module procedure qx_hypergeom_pfq
end interface

interface incgamma
  module procedure qx_incgamma
end interface

interface polygamma
  module procedure qx_polygamma
end interface

interface polygamma_be
  module procedure qx_polygamma_be
end interface

interface polylog_ini
  module procedure qx_polylog_ini
end interface

interface polylog_neg
  module procedure qx_polylog_neg
end interface

interface polylog_pos
  module procedure qx_polylog_pos
end interface

interface struve_hn
  module procedure qx_struve_hn
end interface

interface zeta
  module procedure qx_zeta
end interface

interface zeta_be
  module procedure qx_zeta_be
end interface

interface zeta_int
  module procedure qx_zeta_int
end interface

contains

!  QXFUN subroutines and functions, in alphabetical order.

  function qx_agm (qa, qb)
    implicit none
    real (qxknd):: qx_agm
    real (qxknd), intent (in):: qa, qb
    call qxagmr (qa, qb, qx_agm)
    return
  end function

  subroutine qx_berne (nb, rb)
    implicit none
    integer, intent (in):: nb
    real (qxknd), intent (out):: rb(nb)
    call qxberner (nb, rb(1))
    return
  end subroutine

  function qx_bessel_i (qa, ra)
    implicit none
    real (qxknd):: qx_bessel_i
    real (qxknd), intent(in):: qa, ra
    call qxbesselir (qa, ra, qx_bessel_i)
    return
  end function

  function qx_bessel_in (nu, ra)
    implicit none
    real (qxknd):: qx_bessel_in
    integer, intent(in):: nu
    real (qxknd), intent(in):: ra
    call qxbesselinr (nu, ra, qx_bessel_in)
    return
  end function

  function qx_bessel_j (qa, ra)
    implicit none
    real (qxknd):: qx_bessel_j
    real (qxknd), intent(in):: qa, ra
    call qxbesseljr (qa, ra, qx_bessel_j)
    return
  end function

  function qx_bessel_jn (nu, ra)
    implicit none
    real (qxknd):: qx_bessel_jn
    integer, intent(in):: nu
    real (qxknd), intent(in):: ra
    call qxbesseljnr (nu, ra, qx_bessel_jn)
    return
  end function

  function qx_bessel_j0 (ra)
    implicit none
    real (qxknd):: qx_bessel_j0
    integer:: nu
    real (qxknd), intent(in):: ra
    nu = 0
    call qxbesseljnr (nu, ra, qx_bessel_j0)
    return
  end function

  function qx_bessel_j1 (ra)
    implicit none
    real (qxknd):: qx_bessel_j1
    integer:: nu
    real (qxknd), intent(in):: ra
    nu = 1
    call qxbesseljnr (nu, ra, qx_bessel_j1)
    return
  end function

  function qx_bessel_k (qa, ra)
    implicit none
    real (qxknd):: qx_bessel_k
    real (qxknd), intent(in):: qa, ra
    call qxbesselkr (qa, ra, qx_bessel_k)
    return
  end function

  function qx_bessel_kn (nu, ra)
    implicit none
    real (qxknd):: qx_bessel_kn
    integer, intent(in):: nu 
    real (qxknd), intent(in):: ra
    call qxbesselknr (nu, ra, qx_bessel_kn)
    return
  end function

  function qx_bessel_y (qa, ra)
    implicit none
    real (qxknd):: qx_bessel_y
    real (qxknd), intent(in):: qa, ra
    call qxbesselyr (qa, ra, qx_bessel_y)
    return
  end function

  function qx_bessel_yn (nu, ra)
    implicit none
    real (qxknd):: qx_bessel_yn
    integer, intent(in):: nu 
    real (qxknd), intent(in):: ra
    call qxbesselynr (nu, ra, qx_bessel_yn)
    return
  end function

  function qx_bessel_y0 (ra)
    implicit none
    real (qxknd):: qx_bessel_y0
    integer:: nu
    real (qxknd), intent(in):: ra
    nu = 0
    call qxbesselynr (nu, ra, qx_bessel_y0)
    return
  end function

  function qx_bessel_y1 (ra)
    implicit none
    real (qxknd):: qx_bessel_y1
    integer:: nu
    real (qxknd), intent(in):: ra
    nu = 1
    call qxbesselynr (nu, ra, qx_bessel_y1)
    return
  end function

  subroutine qx_cssh (qa, qb, qc)
    implicit none
    real (qxknd), intent (in):: qa
    real (qxknd), intent (out):: qb, qc
    qb = cosh (qa)
    qc = sinh (qa)
    return
  end subroutine

  subroutine qx_cssn (qa, qb, qc)
    implicit none
    real (qxknd), intent (in):: qa
    real (qxknd), intent (out):: qb, qc
    qb = cos (qa)
    qc = sin (qa)
    return
  end subroutine

  function qx_digamma_be (nb, rb, rc)
    implicit none
    integer, intent (in):: nb
    real (qxknd):: qx_digamma_be
    real (qxknd), intent (in):: rb(nb), rc
    call qxdigammabe (nb, rb(1), rc, qx_digamma_be)
    return
  end function

  function qx_egamma ()
    implicit none
    real (qxknd):: qx_egamma
    call qxegamc (qx_egamma)
    return
  end function

  function qx_erf (qa)
    implicit none
    real (qxknd):: qx_erf
    real (qxknd), intent (in):: qa
    call qxerfr (qa, qx_erf)
    return
  end function

  function qx_erfc (qa)
    implicit none
    real (qxknd):: qx_erfc
    real (qxknd), intent (in):: qa
    call qxerfcr (qa, qx_erfc)
    return
  end function

  function qx_expint (qa)
    implicit none
    real (qxknd):: qx_expint
    real (qxknd), intent (in):: qa
    call qxexpint (qa, qx_expint)
    return
  end function

  function qx_gamma (qa)
    implicit none
    real (qxknd):: qx_gamma
    real (qxknd), intent (in):: qa
    call qxgammar (qa, qx_gamma)
    return
  end function

  function qx_hurwitz_zetan (ia, rb)
    implicit none
    real (qxknd):: qx_hurwitz_zetan
    integer, intent (in):: ia
    real (qxknd), intent (in):: rb
    call qxhurwitzzetan (ia, rb, qx_hurwitz_zetan)
    return
  end function

  function qx_hurwitz_zetan_be (nb, rb, is, aa)
    implicit none
    real (qxknd):: qx_hurwitz_zetan_be
    integer, intent (in):: nb, is
    real (qxknd), intent (in):: rb(nb), aa
    call qxhurwitzzetanbe (nb, rb(1), is, aa, &
      qx_hurwitz_zetan_be)
    return
  end function

  function qx_hypergeom_pfq (np, nq, aa, bb, xx)
    implicit none
    real (qxknd):: qx_hypergeom_pfq
    integer, intent (in):: np, nq
    real (qxknd), intent (in):: aa(np), bb(nq), xx
    call qxhypergeompfq (np, nq, aa(1), bb(1), &
      xx, qx_hypergeom_pfq)
    return
  end function

  function qx_incgamma (ra, rb)
    implicit none
    real (qxknd):: qx_incgamma
    real (qxknd), intent (in):: ra, rb
    call qxincgammar (ra, rb, qx_incgamma)
    return
  end function

  function qx_log2 ()
    implicit none
    real (qxknd):: qx_log2
    call qxlog2c (qx_log2)
    return
  end function

  function qx_pi ()
    implicit none
    real (qxknd):: qx_pi
    call qxpic (qx_pi)
    return
  end function    

  function qx_polygamma (nn, ra)
    implicit none
    integer, intent (in):: nn
    real (qxknd), intent (in):: ra
    real (qxknd) qx_polygamma
    call qxpolygamma (nn, ra, qx_polygamma)
    return
  end function

  function qx_polygamma_be (nb, rb, nn, ra)
    implicit none
    integer, intent (in):: nb, nn
    real (qxknd), intent (in):: ra, rb(nb)
    real (qxknd) qx_polygamma_be
    call qxpolygammabe (nb, rb(1), nn, ra, qx_polygamma_be)
    return
  end function

  subroutine qx_polylog_ini (nn, arr)
    implicit none
    integer, intent (in):: nn
    real (qxknd), intent (out):: arr(abs(nn))
    call qxpolylogini (nn, arr(1))
    return
  end subroutine

  function qx_polylog_neg (nn, arr, ra)
    implicit none
    integer, intent (in):: nn
    real (qxknd), intent (in):: arr(abs(nn))
    real (qxknd), intent (in):: ra
    real (qxknd) qx_polylog_neg
    call qxpolylogneg (nn, arr(1), ra, qx_polylog_neg)
    return
  end function

  function qx_polylog_pos (nn, ra)
    implicit none
    integer, intent (in):: nn
    real (qxknd), intent (in):: ra
    real (qxknd) qx_polylog_pos
    call qxpolylogpos (nn, ra, qx_polylog_pos)
    return
  end function

  function qx_struve_hn (nu, ra)
    implicit none
    integer, intent (in):: nu
    real (qxknd):: qx_struve_hn
    real (qxknd), intent (in):: ra
    call qxstruvehn (nu, ra, qx_struve_hn)
    return
  end function

  function qx_zeta (ra)
    implicit none
    real (qxknd):: qx_zeta
    real (qxknd), intent (in):: ra
    call qxzetar (ra, qx_zeta)
    return
  end function

  function qx_zeta_be (nb, rb, rc)
    implicit none
    integer, intent (in):: nb
    real (qxknd):: qx_zeta_be
    real (qxknd), intent (in):: rb(nb), rc
    call qxzetabe (nb, rb(1), rc, qx_zeta_be)
    return
  end function

  function qx_zeta_int (ia)
    implicit none
    real (qxknd):: qx_zeta_int
    integer, intent (in):: ia
    call qxzetaintr (ia, qx_zeta_int)
    return
  end function

end module
