program sw_oned
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: N = 100
  real(DP), parameter :: L = 100.0d0
  real(DP), parameter :: XL = 0.0d0
  real(DP), parameter :: GG = 9.81d0
  real(DP), parameter :: X_DAM = 0.4d0 * L
  real(DP), parameter :: DX = L / DFLOAT(N)
  real(DP), parameter :: DEPTH_LEFT = 1.2d0
  real(DP), parameter :: DEPTH_RIGHT = 1.0d0
  real(DP), parameter :: VELOCITY_LEFT = 1.0d0
  real(DP), parameter :: VELOCITY_RIGHT = 1.0d0
  real(DP), parameter :: BC_U_LEFT = -1.0d0
  real(DP), parameter :: BC_U_RIGHT = -1.0d0
  character(len=20), parameter :: OUTPUT_FILE = 'res.dat'

  real(kind=DP), allocatable :: x(:), h(:), u(:)
  allocate(x(0:N+1), h(0:N+1), u(0:N+1))

  call set_mesh(x)
  call set_ic_dambreak(x, h, u)
  call output_solution(OUTPUT_FILE, x, h, u)
  call set_bc(h, u)

  deallocate(x, h, u)

contains

subroutine set_mesh(x)
  real(DP), intent(out) :: x(0:N+1)
  integer :: i

  do i=0,N+1
    x(i) = XL+0.5d0*dx+DFLOAT(i-1)*DX
  enddo
end

subroutine set_ic_dambreak(x, h, u)
  real(DP), intent(in) :: x(0:N+1)
  real(DP), intent(out) :: h(0:N+1), u(0:N+1)
  integer :: i

  do i=1,N
    if (x(i).le.X_DAM) then
      h(i)=DEPTH_LEFT
      u(i)=VELOCITY_LEFT
    else
      h(i)=DEPTH_RIGHT
      u(i)=VELOCITY_RIGHT
    endif
  enddo
end

subroutine set_bc(h, u)
  real(DP), intent(inout) :: h(0: N+1), u(0: N+1)

  h(0) = h(1)
  u(0) = BC_U_LEFT * u(1)
  h(N+1) = h(N)
  u(N+1) = BC_U_RIGHT * u (N)
end

subroutine output_solution(filename, x, h, u)
  character(len=20), intent(in) :: filename
  real(DP), dimension(0:N+1), intent(in) :: x, h, u
  integer :: i

  open(unit=10,file=filename)
  do i=1,N
    write(10,*) x(i), h(i), u(i)
  enddo
  close(10)
end

end program sw_oned
