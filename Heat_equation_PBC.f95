Program Heat_equation
  use FTCS_BTCS_solver
  implicit none
  !define variables
  !parameter(s)

!heat equation variables
  real (kind=dp) :: alpha, u_x_t, lambda  ! var of eqn, may only use in function
  real (kind=dp) :: delta_x, delta_t, delta_x_i, delta_t_i   ! timestep + iteration 
  real (kind=dp) :: time_duration,  object_length_x, u_0, u_N
  character :: foolishness

!solver arrays
  real (kind=dp), dimension(:), allocatable :: u, c
  real (kind=dp), dimension(:,:), allocatable :: matrix
  
!case variables
  integer :: solver! boundary conditions
  integer :: object_shape
  
  !Declare minor stuff
  integer, parameter :: stdout=6
  integer, parameter :: stdin=5


  ! Building user input, based on code from the previous lab
!User interface:
  !Welcome
  write (stdout,'(/A)') 'Welcome to this Heat Equation Solver Program'  
 
  print*, 'Case to analyze? 1, 1D ring with Periodic boundry conditions (PBCs) or 2, 1D bar with Dirichlet boundary conditions &
  &(specified start and end value). enter integer to choose option.'
  read*, object_shape

  print*, 'Numerical Solver used? 1, Forward-Time Centric solver (FTCS) or 2, Backward-Time Centric solver (BTCS). enter integer&
  & to choose option.'
  read*, solver
  
  print*, 'initial temperature, u_0 (initial heat of object at a point)? will ask again for bar example'
  read*, u_0
  print*, 'Value for x (length change in heat is considered over)?'
  read*, object_length_x
  print*, 'Value for t (time)?'
  read*, time_duration

  print*, 'steps for the FTCS and BTCS solvers are taken using the variable lambda, dictated by the function:&
  & lambda = (alpha*delta_t)/(delta_x**2) . ' !explaining lambda eqn
  print*, 'Value for alpha (thermal diffusivity)? Reccommended lambda is  under 1/2,&
  & for FTCS, BTCS should always be stable.'
  read*, alpha
  print*, 'Numerical solver dx (distance step size)? Reccommended lambda is  under 1/2,&
  & for FTCS, BTCS should always be stable.'
  read*, delta_x ! will test to see if step sizes used create a valid lamba
  delta_x_i= delta_x ! to have initial value for space-steps
 
  print*, 'Numerical solver dt (time size)? Recommended lambda is  under 1/2,&
  & for FTCS, BTCS should always be stable.'
  read*, delta_t ! will test to see if step sizes used create a valid lamba
  delta_t_i= delta_t ! to have initial value for time-steps

!defining dependant variables
  N= int(time_duration/delta_t)
  allocate(u(1:N))

  !allocating matrix
  allocate(matrix(N,N))
  matrix=0

! opening relevant files
  call open_files(istat, error_msg)  
  
  ! create lamba test to see if ratio is stable. If not print warning but provide option for using it unstabally anyway.

  lambda = lambda_(alpha, delta_x, delta_t) 
  print *, 'lambda =', lambda


  
select case (object_shape)
case (1)!ring
   ! setup:
   ! definine initial temp array here
     u = u_0
     u_x_t = u_0
     print *, 'Array of initial heat distribution:', u
     print *, 'initial total heat:', u_0
     print *
   
   select case (solver)
    case (1)!FTCS
     if (lambda >= 0.0_dp.AND.lambda <= 0.5_dp) then
      print*, 'according to von Neumann analysis, chosen lambda is conditionally stable'
     else 
      print*, 'according to von Neumann analysis, chosen lambda is unstable. Do you wish to continue, Y/N?'
      read *, foolishness
      if (foolishness /= 'Y') stop ! Exit program, user has to rerun program to enter valid lambda
     end if
     !solver
     print*, 'calling A_'
     call A_(lambda, matrix,N,delta_t, delta_x,alpha)
           ! as noted in module, PBC cqlculqtion sam form for BTCS and for FTCS just with different matrix
     print*, 'calling PBCs'   
     call PBCs(u, matrix, u_x_t, delta_t, delta_x, delta_t_i, delta_x_i, alpha, lambda)
  case (2)!BTCS
       call invert_matrix_B(matrix, N, lambda,delta_t, delta_x,alpha)
       call PBCs(u, matrix, u_x_t, delta_t, delta_x, delta_t_i, delta_x_i, alpha, lambda)
  case default
      print*, 'hey, something went wrong with the solver case select'

  end select
    ! (The total heat in the ring is the sum of the heat at each grid point multiplied by the step size)
  print*, 'Final total heat:', u_x_t
   
case (2)!bar
   
   ! definine initial temp arrays here

   print*, 'initial temperature, u_0, on one side of bar?'
   read*, u_0
   
   u(1:(N/2)) = u_0

   print*, 'initial temperature, u_N, on other side of bar?'
   read*, u_N
   u(((N/2)+1):N) = u_N

   print *, 'Array of initial heat distribution:', u
   print *
   
   
   !seting up c array
   allocate(c(1:N))
   c=0
   c(1)=0
   c(N)= object_length_x

   !solvers
   select case (solver)
    case (1)!FTCS
     if (lambda >= 0.0_dp.AND.lambda <= 0.5_dp) then
      print*, 'according to von Neumann analysis, chosen lambda is conditionally stable'
     else 
      print*, 'according to von Neumann analysis, chosen lambda is unstable. Do you wish to continue, Y/N?'
      read *, foolishness
      if (foolishness /= 'Y') stop ! Exit program, user has to rerun program to enter valid lambda
     end if
     !solver
     call A_(lambda, matrix,N,delta_t, delta_x,alpha)
     call FTCS_Dirichlet( u, c, matrix, lambda, u_x_t, delta_t, delta_x, delta_t_i, delta_x_i, alpha)
    
     
    case (2)!BTCS
       call invert_matrix_B(matrix, N, lambda,delta_t, delta_x,alpha)
       call BTCS_Dirichlet( u, c, matrix, lambda, u_x_t, delta_t, delta_x, delta_t_i, delta_x_i, alpha)
    case default
      print*, 'hey, something went wrong with the solver case select'
   end select
   
case default
  print*, 'hey, something went wrong with the object_shape case select'  
end select
 

End program  
