program signallab
    implicit none
    !maths variables
    integer, parameter :: dp=selected_real_kind(15,300) 
    real (kind=dp) :: pi, period, dt, t , amp, f_c, L, sigma, mu, exponent
    integer ::  N, k, i, window=0, signal=0
    complex (kind=dp) :: W0
    
    !arrays
    real (kind=dp), dimension(:), allocatable ::    h_n   !vectors
    complex (kind=dp), dimension(:), allocatable::  t_data, h_k ,h_other        !complex vectors
    real (kind=dp), dimension(:,:), allocatable ::  W           !matricies
    complex (kind=dp), dimension(:,:), allocatable ::  win 

    !Declare minor stuff
    integer, parameter :: stdout=6
    integer, parameter :: stdin=5
    integer :: istat

    character (len=80) :: answer, error_msg 
    
    pi=4.d0*datan(1.d0) ! const


    ! User interface
    print*, 'Welcome to this signal processing program! please select what medium you want to&
    & process and fill out the relevant information below.'
    print*, 'N?'
    read*, N
    print*, 'period/length of pulse?'
    read*, period
    L=period !for trianglular and rectangular pulse
    print*, 'amplitude/height of pulse?'
    read*, amp
    !assigning value to dt
    dt=(5.0_dp*period/real(N,kind=dp))

    
    ! signal selection and generation
    print*, 'sine wave(1), rectangular pulse(2), triangular pulse(3), or gaussian pulse(4)? please enter integer.'  ! case select
    read*, signal
    
    !Allocating t_data array
    allocate(t_data(0:N-1))

    select case (signal)
     case(1) ! sine wave
      !populating t_data array (discretised sine function data)
      do i=0,N-1
       t=real(i,kind=dp)*dt
       t_data(i)=cmplx(amp*sin(2.d0*pi*t/period),0.0_dp,kind=dp)
      end do 
    case(2) ! rectangle pulse (func?)
      ! do piecewise
      do i=0,N-1
        if ((i<(int((real(N)-L)/2))) .OR. (i>(int(N-((real(N)-L)/2))))) then
         t_data(i) = 0
        else
         t_data(i)=cmplx(amp, kind=dp)
        end if
      end do 
    case(3) !triangle pulse (func)
      ! do piecewise
      L=L/2.0_dp
      do i=0,N-1
        if ((i<(int((real(N, kind=dp)-L)/2.0_dp))) .OR. (i>(int(N-((real(N, kind=dp)-L)/2.0_dp))))) then
          t_data(i) = 0
        else if ((i>=(int((real(N, kind=dp)-L)/2.0_dp))) .AND. (i<=(int(N-(real(N, kind=dp)/2.0_dp))))) then     
          t_data(i)=cmplx(amp*(1.0_dp-(-((real(i, kind=dp))-N/2.0_dp)/((L+2.0_dp)/2.0_dp))), kind=dp)
          print*, 'running'
        else 
          t_data(i)=cmplx(amp*(1.0_dp-(((real(i, kind=dp))-N/2.0_dp)/((L+2.0_dp)/2.0_dp))), kind=dp)
          print*, 'still'
        end if
      end do 
    case(4) !gauss pulse (func)
      print*, 'standard deviation?' !variables
      read*, sigma
      mu= real(N, kind=dp)/2.0_dp 
      do i=0,N-1 !function
       exponent=-0.5_dp*(((real(i, kind=dp)-mu)/sigma)**2)
       t_data(i)=cmplx(((1.0_dp/(sigma*sqrt(2*pi)))*exp(exponent))*amp, kind=dp)
      end do 

    case default
    print*,  'issue with signal case select.'
    end select

    ! ask about windowing functions and pulse (use case select choosing between them (or not using))
    print*, 'use windowing functions? if so enter integer, triangle(1), cosine(2) or gaussian(3). Else write 0.'  ! case select
    read*, window

    select case (window)
    case (1)
      ! triangle (define operation in functions), can the same function be used for windowing function ond pulse?&
      !& probably? like with different inputs maybe?
      allocate(win(0:N-1,0:N-1))
      allocate(h_other(0:N-1))
      win=cmplx(0, kind=dp)
      h_other=1
      do i=0,N-1
        if (i<(real(N)/2.0_dp)) then
          win(i,i)=cmplx(real(i, kind=dp)/(real(N, kind=dp)/2.0_dp), kind=dp)
          h_other(i)=cmplx(real(i, kind=dp)/(real(N, kind=dp)/2.0_dp), kind=dp)
        else 
          win(i,i)=cmplx(-(real(i, kind=dp)/(real(N, kind=dp)/2.0_dp))+2.0_dp, kind=dp)
          h_other(i)= cmplx(-(real(i, kind=dp)/(real(N, kind=dp)/2.0_dp))+2.0_dp, kind=dp)
        end if
      end do 
      ! mutiply by t_data
      t_data= matmul(t_data, win)
    case (2)
      ! cosine (define operation in functions)
      allocate(win(0:N-1,0:N-1))
      win=cmplx(0, kind=dp)
      do i=0,N-1
        win(i,i)=cmplx(cos(real(i, kind=dp))*(real(N, kind=dp)*2.0_dp), kind=dp)
      end do 
      ! mutiply by t_data
      t_data= matmul(t_data, win)
    case (3)
      ! gaussian (define operation in functions). how to define gaussian function in fortran?
      allocate(win(0:N-1,0:N-1))
      win=cmplx(0, kind=dp)
      sigma = 6.84_dp ! from trail and error on desmos 
      mu= real(N, kind=dp)/2.0_dp
      do i=0,N-1
       exponent=-0.5_dp*(((real(i, kind=dp)-mu)/sigma)**2)
       win(i,i)=cmplx(exp(exponent), kind=dp)
      end do 
      ! mutiply by t_data
      t_data= matmul(t_data, win)
    case default
      print*, 'no windowing function used.' 
    end select

  ! DFT execution  
  !allocating h_n + allocating & populating h_k vectors
    allocate(h_k(0:N-1))
    h_k = t_data  ! the stuff from the lab here
    allocate(h_n(0:N-1))
  
    !allocating & populating W matrix
    allocate(W(0:N-1, 0:N-1))
    W0 = exp(cmplx(0.0_dp, (2.0_dp)/real(N,kind=dp),kind=dp))
    do i = 0, N-1
      do k = 0,N-1
        W(k,i) = W0**(i*k)
      end do
    end do
  
  !apply DFT to data
   ! claculating h_n
      h_n = matmul(W, h_k)
  
  
  ! Nyquist limit:
      f_c = 1/(2*dt) 
      print*, 'Nyquist critical frequency, f_c =', f_c
      print*, 'Sampling frequency, 1/dt =', (1/dt)


! returning answers/ results
    ! open file
    open(unit=12, file="original_sig.csv", status="unknown", action="write", &
    & iostat=istat, iomsg=error_msg)
    if (istat /= 0) then
      print *,trim(error_msg)
      stop "Error in open"
    end if

    open(unit=10, file="DFT_result.csv", status="unknown", action="write", &
    & iostat=istat, iomsg=error_msg)
    if (istat /= 0) then
      print *,trim(error_msg)
      stop "Error in open"
    end if

    !write to file
    write(unit=12, fmt=*, iostat=istat) 'N', ',','real', ',','imag'
    do i=0,N-1 ! line by line
     write(unit=12, fmt=*, iostat=istat) i, ',', real(t_data(i)),',', aimag(t_data(i))  ! original signal
    end do
    !closing t_data file
    close(unit=12, iostat=istat)
    if (istat /= 0) stop "Error closing original_sig.csv"  
    
    write(unit=10, fmt=*, iostat=istat) 'N_dft', ',','real_dft'
    do i=0,N-1 ! line by line
     write(unit=10, fmt=*, iostat=istat) i, ',', h_n (i)  ! DFT
    end do
    !closing DFT file
    close(unit=10, iostat=istat)
    if (istat /= 0) stop "Error closing DFT_result.csv"    
    
    ! print resulting vector
    print*, dt
    print*, 'winshape =', shape(win), 't_datashape =', shape(t_data)
!note: will graph resultant file using python

end program signallab