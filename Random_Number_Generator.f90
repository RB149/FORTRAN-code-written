Program test_random_numbers
 implicit none
 integer :: A, B, M, n=1, x=1, count, istat, I
 integer, parameter :: dp = selected_real_kind(15,300)
 real(kind=dp) :: N_mean, N_vari, sigma, N_mean_err, N_vari_err, x_, pi=3.14159265359, x_min, x_max, dx
 character(len=20)::  error_msg
 real(kind=dp), dimension(:), allocatable :: N_, N_norm, N_gauss
 real(kind=dp), dimension(2):: gaussian_dist_
 integer, dimension(:,:), allocatable :: N_bin_n, N_bin_g
 
!Making the generator (A) 
!requesting values to be inputed
print *, 'A='
read*, A

print *, 'B='
read *, B

print *, 'M='
read*, M

print *, 'count='
read*, count

print *, 'amount of bins='
read*, I

!creats document for the random values to be read into
open(unit=10, file="randomnumbers.txt", status="new", action="write", &
& iostat=istat, iomsg=error_msg)
if (istat /= 0) then
  print *,trim(error_msg)
  stop "Error in open"
end if

!do while loop generates the random numbers
do while (n <= count)
   x=rand_number(A,B,M,x)
   write(unit=10, fmt=*, iostat=istat) x,','
   n=n+1
end do

!-------------------------------------------Normalized RNG------------------------------------------------------

!Normalization and testing (B)
!redefining variables
n=1
x=1

!creats document for the normalized random values to be read into
open(unit=11, file="normrandomnumbers.txt", status="new", action="write", &
& iostat=istat, iomsg=error_msg)
if (istat /= 0) then
  print *,trim(error_msg)
  stop "Error in open"
end if

!Array introduced to facilitate testing the normalized generator's accuracy
allocate(N_(count))
allocate(N_norm(count))

!do while loop generates the random numbers
do while (n <= count)
   x_= norm_random_number(A,B,M)
   write(unit=11, fmt=*, iostat=istat) x_,','
   N_(n)=x_
   x=rand_number(A,B,M,x)
   n=n+1
end do

N_norm = N_

!redefining variables
n=1


!-------------------------------------------Gauss RNG------------------------------------------------------
!Array introduced to facilitate testing the Gaussian generator's accuracy
allocate(N_gauss(count))

!do while loop to gaussianly distribute the random normal numbers
do while (n < count)
   gaussian_dist_= gaussian_dist(n,pi,N_)
 !  write(unit=11, fmt=*, iostat=istat) x_,','
   N_(n)=gaussian_dist_(1)
   N_(n+1)=gaussian_dist_(2)
   n=n+2
end do

N_gauss = N_

!-------------------------------------------Testing----------------------------------------------

!------------------------------------------------STATISTICS DATA NOT PRINTED, ONLY HERE FOR RETRIEVING ANY NECESSARY VARIABLES-------------------

!Finding mean and variance
N_mean= Sum(N_)/count

n = 1
do while (n<= count)
  N_(n)= (N_(n)-N_mean)**2.0
  n=n+1
end do
N_vari= Sum(N_)/(real(count-1))

!Testing the mean and variance error
n = 1
do while (n<= count)
  N_(n)= (N_(n)-N_mean)**2.0
  n=n+1
end do
sigma = (Sum(N_)/real(count))**0.5

N_mean_err = sigma/((real(count))**0.5)
N_vari_err = (sigma**2.0)*((2.0/(real(count)-1))**0.5)

!------------------------------------------------STATISTICS DATA NOT PRINTED-------------------


!---------------------------------------------------Distribution------------------------------
!reset variables
n=1

!allocating bin
allocate(N_bin_n(2,I))

do while (n <= I)
   N_bin_n(1,n)=n
   n=n+1
end do

n=1
allocate(N_bin_g(2,I))

do while (n <= I)
   N_bin_g(1,n)=n
   n=n+1
end do

!creats document for the 1st RNG bin to be read into
open(unit=12, file="normrandombinns.dat", status="new", action="write", &
& iostat=istat, iomsg=error_msg)
if (istat /= 0) then
  print *,trim(error_msg)
  stop "Error in open"
end if

!Distribution for Normalized RNG
x_min= 0.0
x_max= 1.0
dx = (x_max - x_min)/real(I)

!reset variables
n=1
x=1

!binning for Normalized RNG data
do while (n <= count)
  x=int((((N_norm(n))-x_min)/dx)+1)
  x=min(x,100)   !stops overflow
  x=max(x,1)
  N_bin_n(2,x) = N_bin_n(2,x)+1
  n=n+1
end do

!writing to file
n=1
do while (n <= I)
  x = N_bin_g(2,n)
  write(unit=12, fmt=*, iostat=istat) n,',',x
  n=n+1
end do

!creats document for the 2nd RNG bin to be read into
open(unit=13, file="gaussbinns.dat", status="new", action="write", &
& iostat=istat, iomsg=error_msg)
if (istat /= 0) then
  print *,trim(error_msg)
  stop "Error in open"
end if

!Distribution for Gaussian RNG
x_min= -3.0
x_max= 3.0
dx = (x_max - x_min)/real(I)

!reset variables
n=1
x=1
!binning for Gaussian RNG
do while (n <= count)
  x=int((((N_gauss(n))-x_min)/dx)+1)
  x=min(x,100)   !stops overflow
  x=max(x,1)
  N_bin_g(2,x) = N_bin_g(2,x)+1
  n=n+1
end do

!writing to file
n=1
do while (n <= I)
  x = N_bin_g(2,n)
  write(unit=13, fmt=*, iostat=istat) n,',',x
  n=n+1
end do
  
!...functions

contains

function rand_number(A,B,M,x)
  implicit none
  integer, intent(in) :: A,B,M,x
  integer :: rand_number
  rand_number= modulo(A*(x) + B, M) 
end function rand_number

function norm_random_number(A,B,M)
  implicit none
  integer, intent(in) :: A,B,M
  real(kind=dp):: norm_random_number
  norm_random_number= real(rand_number(A,B,M,x))/real(M)
end function norm_random_number

function gaussian_dist(n,pi,N_)
  implicit none
  integer, intent(in) :: n
  real(kind=dp), intent(in) :: pi
  real(kind=dp), dimension(:), intent(in) :: N_
  real(kind=dp), dimension(2):: gaussian_dist
  gaussian_dist(1)=sqrt((-2.0_dp)*(log(N_(n))))*cos(2.0_dp*pi*(N_(n+1)))
  gaussian_dist(2)=sqrt((-2.0_dp)*(log(N_(n))))*sin(2.0_dp*pi*(N_(n+1))) 
end function gaussian_dist
  
End program test_random_numbers
