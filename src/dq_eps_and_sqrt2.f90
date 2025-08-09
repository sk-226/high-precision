program dq_eps_and_sqrt2
  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
  use dqmodule
  use dqfun_cwrap, only: dq_tostr
  implicit none
  type(dq_real) :: one, two, eps, prev, val
  character(kind=c_char,len=320) :: cbuf
  character(len=320) :: out

  one = dqreal('1')
  two = dqreal('2')

  ! ! 反復法（unit roundoff の DQ 版）
  ! eps = one
  ! do
  !   prev = eps
  !   eps  = eps / two
  !   if (one + eps == one) exit
  ! end do
  ! call dq_tostr(prev%dqr, 120_c_int, cbuf, len(cbuf, kind=c_int))
  ! call cstr_to_fstr(cbuf, out)
  ! print '(A,1X,A)', 'eps(DQ):', trim(out)

  ! (sqrt(2))^2 - 2 の DQ 版
  val = (sqrt(two)**2) - two
  call dq_tostr(val%dqr, 120_c_int, cbuf, len(cbuf, kind=c_int))
  call cstr_to_fstr(cbuf, out)
  print '(A,1X,A)', 'res(sqrt2^2-2):', trim(out)

contains

  subroutine cstr_to_fstr(cbuf, fstr)
    character(kind=c_char,len=*), intent(in) :: cbuf
    character(len=*), intent(out) :: fstr
    integer :: i, n
    n = len(fstr)
    fstr = ''
    do i = 1, len(cbuf)
      if (cbuf(i:i) == c_null_char) exit
      if (i <= n) fstr(i:i) = achar(iachar(cbuf(i:i)))
    end do
  end subroutine cstr_to_fstr

end program
