program sw_oned
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: ITER_MAX = 100000000
  integer, parameter :: N = 100
  real(DP), parameter :: TIME_FIN = 1.0d0
  real(DP), parameter :: CFL = 0.9d0
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

  integer :: it
  real(DP) :: dt=1.d-8, time=0.0d0
  real(DP) :: t_start, t_end
  real(kind=DP), allocatable :: x(:), h(:), u(:)

  allocate(x(0:N+1), h(0:N+1), u(0:N+1))

  call cpu_time(t_start)

  call set_mesh(x)
  call set_ic_dambreak(x, h, u)
  call output_solution(OUTPUT_FILE, x, h, u)
  call set_bc(h, u)
  do while(time<TIME_FIN)
    if (it.ge.ITER_MAX) exit
    dt = timestep(h,u)
    if (dt>TIME_FIN-time) dt = TIME_FIN-time

    it=it+1
    time=time+dt
  enddo

  call cpu_time(t_end)
  deallocate(x, h, u)

  print*, 'cpu time:  ', t_end - t_start
  print*, 'time:      ', time
  print*, 'iterations:', it

contains

! INITIALISATION PROCEDURES

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

real(DP) function timestep(h,u)
  implicit none
  real(kind=DP), dimension(0:N+1), intent(in) :: h, u
  real(kind=DP) :: cmax
  integer :: i

  cmax = 0.0d0
  do i=1,N
    cmax = max( cmax, dabs(u(i)) + dsqrt(GG*h(i)) )
  enddo
  timestep = CFL*dx/cmax
end

! AUX PROCEDURES

subroutine prim_to_cons(h, u, h_cons, u_cons)
  real(DP), intent(in) :: h(0:N+1), u(0:N+1)
  real(DP), intent(out) :: h_cons(0:N+1), u_cons(0:N+1)
  integer :: i

  do i = 1, N
    h_cons(i) = h(i)
    u_cons(i) = h(i) * u(i)
  enddo

end

subroutine cons_to_prim(h_cons, u_cons, h, u)
  real(DP), intent(in) :: h_cons(0:N+1), u_cons(0:N+1)
  real(DP), intent(out) :: h(0:N+1), u(0:N+1)
  integer :: i

  do i = 1, N
    h(i) = h_cons(i)
    u(i) = u_cons(i) / h(i)
  enddo

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
