module FTCS_BTCS_solver
implicit none
!define global variables
integer, parameter  :: dp=selected_real_kind(15,300) !used in most subroutines/functions so defined gloablly here

integer :: istat, N   ! istat global as used in all module subroutines and N global as use outside functions/subroutines
character (len=80) :: error_msg




 contains
!setup
 subroutine open_files(istat, error_msg)
 integer :: istat
 character (len=80) :: error_msg
  
   !open file
   open(unit=10, file="Solver.dat", status="unknown", action="write", &
   & iostat=istat, iomsg=error_msg)
   if (istat /= 0) then
   print *,trim(error_msg)
   stop "Error in open"
   end if

 end subroutine open_files

!Multiplying matrix
 subroutine A_(lambda, matrix,N,delta_t, delta_x,alpha)
  real (kind=dp) :: lambda, delta_t, delta_x,alpha
  real (kind=dp), dimension(:,:) :: matrix
  integer :: N,i

  lambda = lambda_(alpha, delta_x, delta_t)

  !populating
  matrix(1,N)= lambda
  matrix(N,1)= lambda
  do i = 1,N
     matrix(i,i)= 1-2*lambda
     if (1<i.AND.i<N) then
         matrix(i,(i-1))= -lambda
         matrix((i-1),i)= -lambda
     end if
  end do
 end subroutine

 subroutine invert_matrix_B(matrix, N, lambda,delta_t, delta_x,alpha)
    ! Invert the supplied matrix using LAPACK
    implicit none
    real(kind=dp), dimension(:,:), intent(inout) :: matrix
    integer    :: N, LWORK, IERR
    integer,    dimension(:), allocatable  :: IPIV
    real(kind=dp), dimension(:), allocatable      :: WORK
    real (kind=dp) :: lambda, delta_t, delta_x, alpha
    integer :: i
    
    lambda = lambda_(alpha, delta_x, delta_t)
    !populating
    matrix(1,N)= -lambda
    matrix(N,1)= -lambda
    do i = 1,N
       matrix(i,i)= 1+2*lambda
       if (1<i.AND.i<N) then
         matrix(i,(i-1))= -lambda
         matrix((i-1),i)= -lambda
       end if
    end do

    !Matrix inversion
    if (size(matrix,1) /= size(matrix,2)) STOP "Matrix is not square"
    N = size(matrix,1)
    allocate(IPIV(N),stat=IERR)
    if (IERR/=0) STOP "Failed to allocate IPIV"
    LWORK = N**2
    allocate(WORK(LWORK),stat=IERR)
    if (IERR/=0) STOP "Failed to allocate WORK"
    call dgetrf(N,N,matrix,N,IPIV,IERR)
    if (IERR/=0) STOP "Error in dgetrf: Matrix is singular"
    call dgetri(N,matrix,N,IPIV,WORK,LWORK,IERR)
    if (IERR/=0) STOP "Error in dgetri: Matrix is singular"
 end subroutine

! can combine FTCS and BTCS for PBCs into one subroutine as the only area where they differ is in the matrix they use, the form of the equation stays same 
 subroutine PBCs(u, matrix, u_x_t, delta_t, delta_x, delta_t_i, delta_x_i, alpha, lambda)
 real (kind = dp), dimension (:) :: u
 real (kind=dp), dimension (:,:) :: matrix
 real (kind=dp) :: u_x_t, delta_t,  delta_x, delta_t_i, delta_x_i, alpha, lambda
 integer :: solver,i
 write(unit=10, fmt=*, iostat=istat) u_x_t
    do i= 2,N
     if (solver ==1) then
      call A_(lambda, matrix,N,delta_t, delta_x,alpha)
     else
      call invert_matrix_B (matrix, N, lambda,delta_t, delta_x,alpha)
     end if
     u = MATMUL(matrix,u)
     u_x_t = sum(u)/real(N)
     write(unit=10, fmt=*, iostat=istat) u_x_t
     delta_t = delta_t + delta_t_i
     delta_x = delta_x + delta_x_i
    end do
 end subroutine
     
 subroutine FTCS_Dirichlet( u, c, matrix, lambda, u_x_t, delta_t, delta_x, delta_t_i, delta_x_i, alpha)
 real (kind = dp), dimension (:) :: u, c
 real (kind=dp), dimension (:,:) :: matrix
 real (kind=dp) :: lambda, u_x_t, delta_t, delta_x, delta_t_i, delta_x_i, alpha
 integer :: i
    do i= 1,N
      lambda=lambda_(alpha, delta_x, delta_t)
      u = MATMUL(matrix,u) + lambda*c
      u_x_t = sum(u)*delta_t
      !write(unit=10, fmt=*, iostat=istat) u_x_t
      delta_t = delta_t + delta_t_i
      delta_x = delta_x + delta_x_i
    end do
    do i= 1,N    !writing values for ftcs
      write(unit=10, fmt=*, iostat=istat) u(i)
    end do
 end subroutine 

 subroutine BTCS_Dirichlet( u, c, matrix, lambda, u_x_t, delta_t, delta_x, delta_t_i, delta_x_i, alpha)
 real (kind = dp), dimension (:) :: u, c
 real (kind=dp), dimension (:,:) :: matrix
 real (kind=dp) :: lambda, u_x_t, delta_t, delta_x, delta_t_i, delta_x_i, alpha
 integer :: i
    do i= 1,N    
      lambda=lambda_(alpha, delta_x, delta_t)
      u = MATMUL(matrix,u) + MATMUL(matrix, (lambda*c)) 
      u_x_t = sum(u)*delta_t
     ! write(unit=10, fmt=*, iostat=istat) u_x_t
      delta_t = delta_t + delta_t_i
      delta_x = delta_x + delta_x_i
    end do
    do i= 1,N    !writing values for btcs
      write(unit=10, fmt=*, iostat=istat) u(i)
    end do
 end subroutine


 !lambda function
 function lambda_(alpha, delta_x, delta_t) result (lambda)
 real (kind=dp), intent(in) :: alpha, delta_x, delta_t
 real (kind=dp) :: lambda
   !function for lambda
   lambda = (alpha*delta_t)/(delta_x**2)
 end function
  
end module
