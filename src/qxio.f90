module qxio
  use qxmodule
  implicit none
  private

  ! Public API: minimal I/O wrappers for QXFUN quad kind
  public :: qxprint, qxprintz, qxsprint, qxsprintz
  public :: qxwrite, qxread, qxeform, qxfform

  character(len=*), parameter :: QX_REAL_FMT    = '(ES80.33E4)'
  character(len=*), parameter :: QX_COMPLEX_FMT = '(ES80.33E4,1X,ES80.33E4)'

  ! Generic interfaces must be declared in the specification part
  interface qxwrite
     module procedure qxwrite_real
     module procedure qxwrite_cplx
  end interface

  interface qxread
     module procedure qxread_real
     module procedure qxread_cplx
  end interface

contains

  ! -------------------- Fixed-format convenience --------------------

  subroutine qxprint(unit, x, iostat)
    integer, intent(in) :: unit
    real(qxknd), intent(in) :: x
    integer, intent(out), optional :: iostat
    integer :: ios
    write(unit, QX_REAL_FMT, iostat=ios) x
    if (present(iostat)) iostat = ios
  end subroutine qxprint

  subroutine qxprintz(unit, z, iostat)
    integer, intent(in) :: unit
    complex(qxknd), intent(in) :: z
    integer, intent(out), optional :: iostat
    integer :: ios
    write(unit, QX_COMPLEX_FMT, iostat=ios) real(z, kind=qxknd), aimag(z)
    if (present(iostat)) iostat = ios
  end subroutine qxprintz

  subroutine qxsprint(buf, x, status)
    character(len=*), intent(out) :: buf
    real(qxknd), intent(in) :: x
    integer, intent(out), optional :: status
    integer :: ios
    write(buf, QX_REAL_FMT, iostat=ios) x
    if (present(status)) status = ios
  end subroutine qxsprint

  subroutine qxsprintz(buf, z, status)
    character(len=*), intent(out) :: buf
    complex(qxknd), intent(in) :: z
    integer, intent(out), optional :: status
    integer :: ios
    write(buf, QX_COMPLEX_FMT, iostat=ios) real(z, kind=qxknd), aimag(z)
    if (present(status)) status = ios
  end subroutine qxsprintz

  ! -------------------- DQ/DD-style API equivalents -----------------

  subroutine qxwrite_real(unit, width, digits, x, iostat)
    integer, intent(in) :: unit, width, digits
    real(qxknd), intent(in) :: x
    integer, intent(out), optional :: iostat
    character(len=64) :: fmt
    integer :: ios, d_eff, w_eff
    call normalize_wd(width, digits, w_eff, d_eff)
    fmt = build_es_format(w_eff, d_eff)
    write(unit, fmt, iostat=ios) x
    if (present(iostat)) iostat = ios
  end subroutine qxwrite_real

  subroutine qxwrite_cplx(unit, width, digits, z, iostat)
    integer, intent(in) :: unit, width, digits
    complex(qxknd), intent(in) :: z
    integer, intent(out), optional :: iostat
    character(len=96) :: fmt
    integer :: ios, d_eff, w_eff
    call normalize_wd(width, digits, w_eff, d_eff)
    fmt = '(' // build_es_component(w_eff, d_eff) // ',1X,' // build_es_component(w_eff, d_eff) // ')'
    write(unit, fmt, iostat=ios) real(z,kind=qxknd), aimag(z)
    if (present(iostat)) iostat = ios
  end subroutine qxwrite_cplx

  subroutine qxread_real(unit, x, iostat)
    integer, intent(in) :: unit
    real(qxknd), intent(out) :: x
    integer, intent(out), optional :: iostat
    integer :: ios
    read(unit, *, iostat=ios) x
    if (present(iostat)) iostat = ios
  end subroutine qxread_real

  subroutine qxread_cplx(unit, z, iostat)
    integer, intent(in) :: unit
    complex(qxknd), intent(out) :: z
    integer, intent(out), optional :: iostat
    integer :: ios
    read(unit, *, iostat=ios) z
    if (present(iostat)) iostat = ios
  end subroutine qxread_cplx

  subroutine qxeform(x, width, digits, s, iostat)
    real(qxknd), intent(in) :: x
    integer, intent(in) :: width, digits
    character(len=*), intent(out) :: s
    integer, intent(out), optional :: iostat
    character(len=64) :: fmt
    integer :: ios, d_eff, w_eff
    call normalize_wd(width, digits, w_eff, d_eff)
    fmt = build_es_format(w_eff, d_eff)
    write(s, fmt, iostat=ios) x
    if (present(iostat)) iostat = ios
  end subroutine qxeform

  subroutine qxfform(x, width, digits, s, iostat)
    real(qxknd), intent(in) :: x
    integer, intent(in) :: width, digits
    character(len=*), intent(out) :: s
    integer, intent(out), optional :: iostat
    character(len=64) :: fmt
    integer :: ios, d_eff, w_eff
    call normalize_wd(width, digits, w_eff, d_eff)
    fmt = build_f_format(w_eff, d_eff)
    write(s, fmt, iostat=ios) x
    if (present(iostat)) iostat = ios
  end subroutine qxfform

  ! -------------------- Helpers --------------------

  pure subroutine normalize_wd(w_in, d_in, w_out, d_out)
    integer, intent(in)  :: w_in, d_in
    integer, intent(out) :: w_out, d_out
    integer, parameter :: QX_MAX_DIGITS = 33
    d_out = max(1, min(d_in, QX_MAX_DIGITS))
    w_out = max(w_in, d_out + 10)
  end subroutine normalize_wd

  pure function build_es_format(w, d) result(fmt)
    integer, intent(in) :: w, d
    character(len=64) :: fmt
    character(len=16) :: sw, sd
    write(sw,'(I0)') w
    write(sd,'(I0)') d
    fmt = '(ES'//trim(sw)//'.'//trim(sd)//'E4)'
  end function build_es_format

  pure function build_es_component(w, d) result(comp)
    integer, intent(in) :: w, d
    character(len=64) :: comp
    character(len=16) :: sw, sd
    write(sw,'(I0)') w
    write(sd,'(I0)') d
    comp = 'ES'//trim(sw)//'.'//trim(sd)//'E4'
  end function build_es_component

  pure function build_f_format(w, d) result(fmt)
    integer, intent(in) :: w, d
    character(len=64) :: fmt
    character(len=16) :: sw, sd
    write(sw,'(I0)') w
    write(sd,'(I0)') d
    fmt = '(F'//trim(sw)//'.'//trim(sd)//')'
  end function build_f_format

end module qxio
