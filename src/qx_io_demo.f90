program qx_io_demo
  use qxmodule
  use qxio
  implicit none

  real(qxknd)    :: x
  complex(qxknd) :: z
  character(len=128) :: buf

  x = 3.1415926535897932384626433832795028841971_qxknd
  z = cmplx(x, -x, kind=qxknd)

  ! Fixed-format helpers
  call qxprint(6, x)
  call qxprintz(6, z)

  call qxsprint(buf, x)
  print *, 'buf =', trim(buf)

  call qxsprintz(buf, z)
  print *, 'bufz=', trim(buf)

  ! DQ/DD-style API equivalents
  call qxwrite(6, 80, 33, x)
  call qxwrite(6, 80, 33, z)
  call qxeform(x, 80, 33, buf)
  print *, 'eform=', trim(buf)
  call qxfform(x, 80, 33, buf)
  print *, 'fform=', trim(buf)

end program qx_io_demo
