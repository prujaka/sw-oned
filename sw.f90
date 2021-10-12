program sw_oned
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: ITER_MAX = 100
  integer, parameter :: N = 100
  real(DP), parameter :: TIME_FIN = 8.0d0
  real(DP), parameter :: CFL = 0.45d0
  real(DP), parameter :: XL = 0.0d0
  real(DP), parameter :: L = 100.0d0
  real(DP), parameter :: X_DAM = 40.0d0
  real(DP), parameter :: GG = 9.81d0
  real(DP), parameter :: DX = L / DFLOAT(N)
  real(DP), parameter :: DEPTH_LEFT = 1.8d0
  real(DP), parameter :: DEPTH_RIGHT = 1.0d0
  real(DP), parameter :: VELOCITY_LEFT = 0.0d0
  real(DP), parameter :: VELOCITY_RIGHT = 0.0d0
  real(DP), parameter :: BC_U_LEFT = 1.0d0
  real(DP), parameter :: BC_U_RIGHT = 1.0d0
  character(len=20), parameter :: OUTPUT_FILE = 'res.dat'

  integer :: it=0
  real(DP) :: dt=1.d-8, time=0.0d0
  real(DP) :: t_start, t_end
  real(DP), allocatable :: x(:), h(:), u(:), h_cons(:), u_cons(:)
  real(DP), allocatable :: h_flux(:), u_flux(:)

  allocate(x(0:N+1), h(0:N+1), u(0:N+1), h_cons(0:N+1), u_cons(0:N+1))
  allocate(h_flux(0:N+1), u_flux(0:N+1))

  call cpu_time(t_start)

  call set_mesh(x)
  call set_ic_dambreak(x, h, u)
  do while(time<TIME_FIN)
    if (it.ge.ITER_MAX) exit
    dt = dt_from(h, u)
    if (dt>TIME_FIN-time) dt = TIME_FIN-time
    call set_bc(h, u)
    call riemann_fluxes(h, u, h_flux, u_flux)
    call prim_to_cons(h, u, h_cons, u_cons)
    call godunov(h_cons, u_cons, h_flux, u_flux)
    call cons_to_prim(h_cons, u_cons, h, u)

    it=it+1
    time=time+dt
  enddo
  call output_solution(OUTPUT_FILE, x, h, u)

  call cpu_time(t_end)
  deallocate(x, h, u, h_cons, u_cons, h_flux, u_flux)

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
  u(N+1) = BC_U_RIGHT * u(N)
end

real(DP) function dt_from(h,u)
  implicit none
  real(DP), dimension(0:N+1), intent(in) :: h, u
  real(DP) :: cmax
  integer :: i

  cmax = 0.0d0
  do i=1,N
    cmax = max( cmax, dabs(u(i)) + dsqrt(GG*h(i)) )
  enddo
  dt_from = CFL*dx/cmax
end

subroutine riemann_fluxes(h, u, h_flux, u_flux)
  implicit none
  real(DP), dimension(0:N+1), intent(in) :: h,u
  real(DP), dimension(0:N+1), intent(out) :: h_flux, u_flux
  real(DP) :: hl,ul,hr,ur,f1,f2
  integer :: i

  do i=0,N
    hl=h(i)
    ul=u(i)

    hr=h(i+1)
    ur=u(i+1)

    call rusanov(hl, ul, hr, ur, f1, f2)
    h_flux(i) = f1
    u_flux(i) = f2
  enddo
end

subroutine rusanov(hl,ul,hr,ur,f1,f2)
  implicit none
  real(DP), intent(in) :: hl,ul,hr,ur
  real(DP), intent(out) :: f1,f2
  real(DP) :: pl,al,pr,ar
  real(DP) :: h_cons_l, h_cons_r, u_cons_l, u_cons_r
  real(DP) :: h_flux_l, h_flux_r, u_flux_l, u_flux_r
  real(DP) :: sl,sr,sl1,sr1,sl2,sr2,smid

  pl = 0.5d0*GG*hl*hl
  al = DSQRT(GG*hl)
  pr = 0.5d0*GG*hr*hr
  ar = DSQRT(GG*hr)

  ! Wave speed estimation
  ! Davis S.F.(1988), 'Simplified second order Godunov type methods',
  ! SIAM J. Sci. and Stat. Comput., N3, 445-473.
  sl1=ul-al; sl2=ur-ar
  sr1=ul+al; sr2=ur+ar
  sl=DMIN1(sl1,sl2); sr=DMAX1(sr1,sr2)
  smid=dmax1(dabs(sl),dabs(sr))

  h_cons_l=hl   ;	  h_flux_l=hl*ul
  u_cons_l=hl*ul;	  u_flux_l=hl*ul*ul+pl

  h_cons_r=hr   ;   h_flux_r=hr*ur
  u_cons_r=hr*ur;   u_flux_r=hr*ur*ur+pr

  f1 = 0.5d0*(h_flux_l + h_flux_r) - 0.5d0*smid*(h_cons_r - h_cons_l)
  f2 = 0.5d0*(u_flux_l + u_flux_r) - 0.5d0*smid*(u_cons_r - u_cons_l)
end

subroutine godunov(h_cons, u_cons, h_flux, u_flux)
  implicit none
  real(DP), dimension(0:N+1), intent(in) :: h_flux, u_flux
  real(DP), dimension(0:N+1), intent(inout) :: h_cons, u_cons
  integer :: i
  do i=1,N
    h_cons(i) = h_cons(i) - (h_flux(i) - h_flux(i-1)) * dt / dx
    u_cons(i) = u_cons(i) - (u_flux(i) - u_flux(i-1)) * dt / dx
  enddo
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
