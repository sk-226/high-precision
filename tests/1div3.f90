program onediv3
  ! 1/3をdouble, quad, double-double (DD), double-quad (DQ) で計算
  use ddmodule
  use dqmodule
  implicit none

  !=========================================================================
  ! 表示桁数の設定（ここを変更すれば各精度の出力桁数を調整可能）
  !=========================================================================
  integer, parameter :: d_digits  = 69   ! double: 約16桁精度
  integer, parameter :: q_digits  = 69   ! quad: 約34桁精度
  integer, parameter :: dd_digits = 69   ! double-double: 約32桁精度
  integer, parameter :: dq_digits = 69   ! double-quad: 約67桁精度

  ! 出力幅 (digits + 余白)
  integer, parameter :: d_width  = d_digits + 10
  integer, parameter :: q_width  = q_digits + 10
  integer, parameter :: dd_width = dd_digits + 12
  integer, parameter :: dq_width = dq_digits + 12

  ! double precision (約16桁)
  real(kind=8) :: d_one, d_three, d_result
  character(len=d_width) :: d_str

  ! quad precision (約34桁) - IEEE 754 四倍精度
  real(kind=16) :: q_one, q_three, q_result
  character(len=q_width) :: q_str

  ! double-double (DD) - 約32桁精度
  type(dd_real) :: dd_one, dd_three, dd_result

  ! double-quad (DQ) - 約67桁精度
  type(dq_real) :: dq_one, dq_three, dq_result

  ! 出力用バッファ
  character(1) :: dd_buf(dd_width), dq_buf(dq_width)
  character(len=dd_width) :: dd_str
  character(len=dq_width) :: dq_str
  character(len=16) :: fmt_d, fmt_q

  ! ヘッダ
  print '(A)', '1/3 ='

  ! フォーマット文字列を動的に生成
  write(fmt_d, '(A,I0,A,I0,A)') '(F', d_width, '.', d_digits, ')'
  write(fmt_q, '(A,I0,A,I0,A)') '(F', q_width, '.', q_digits, ')'

  ! ========== double ===========
  d_one   = 1.0d0
  d_three = 3.0d0
  d_result = d_one / d_three
  write(*,'(A)', advance='no') '     double        : '
  write(d_str, fmt_d) d_result
  d_str = adjustl(d_str)
  write(*,'(A)') trim(d_str)

  ! ========== quad (long double) ===========
  q_one   = 1.0_16
  q_three = 3.0_16
  q_result = q_one / q_three
  write(*,'(A)', advance='no') '     quad          : '
  write(q_str, fmt_q) q_result
  q_str = adjustl(q_str)
  write(*,'(A)') trim(q_str)

  ! ========== double-double (DD) ===========
  dd_one   = ddreal('1')
  dd_three = ddreal('3')
  dd_result = dd_one / dd_three
  call ddfform(dd_result, dd_width, dd_digits, dd_buf)
  call chars_to_str(dd_buf, dd_str)
  dd_str = adjustl(dd_str)
  write(*,'(A)', advance='no') '     double-double : '
  write(*,'(A)') trim(dd_str)

  ! ========== double-quad (DQ) ===========
  dq_one   = 1.0_dqknd
  dq_three = 3.0_dqknd
  dq_result = dq_one / dq_three
  call dqfform(dq_result, dq_width, dq_digits, dq_buf)
  call chars_to_str(dq_buf, dq_str)
  dq_str = adjustl(dq_str)
  write(*,'(A)', advance='no') '     double-quad   : '
  write(*,'(A)') trim(dq_str)

  print '(A)', ''

contains

  subroutine chars_to_str(buf, str)
    character(1), intent(in) :: buf(:)
    character(len=*), intent(out) :: str
    integer :: i, n
    n = min(size(buf), len(str))
    str = ' '
    do i = 1, n
      str(i:i) = buf(i)
    end do
  end subroutine chars_to_str

end program onediv3
