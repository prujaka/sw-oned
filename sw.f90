program sw_oned
  implicit none
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: N = 100
  real(DP), parameter :: gg = 9.81d0
  real(DP), parameter :: L = 100.0d0
  real(DP), parameter :: XL = 0.0d0
  real(DP), parameter :: DX = L / DFLOAT(N)
  character(len=20), parameter :: OUTPUT_FILE = 'res.dat'

  real(kind=DP), allocatable :: x(:)
  allocate(x(0:N+1))

  call set_mesh(x)
  call output_solution(OUTPUT_FILE, x)

  deallocate(x)

contains

subroutine set_mesh(x)
  real(kind=DP), intent(out) :: x(0:N+1)
  integer :: i

  do i=0,N+1
    x(i) = XL+0.5d0*dx+DFLOAT(i-1)*DX
  enddo
end

subroutine output_solution(filename, x)
  character(len=20), intent(in) :: filename
  real(kind=DP), dimension(0:N+1), intent(in) :: x
  integer :: i

  open(unit=10,file=filename)
  do i=1,N
    write(10,*) x(i)
  enddo
  close(10)
end

end program sw_oned
