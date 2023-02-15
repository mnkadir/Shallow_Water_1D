!Fortran project: a model for the steady state 1-D shallow water equations
!This script will calculate the Steady state of 1D shallow Water Depth, Bed elevation and velocity based on different ODE scheme (Eular, RK2 and RK4)
!For selecting different scheme please change the value of ODE in the input file, (USE 1 for Euler, 2 for Runga Kutta 2nd order and 4 for Runga Kutta 4th order)
!Also change the value of dx in the input file for correspoding analytical solution
program SW
   implicit none
!Variable declaration for input data  
      real*8 :: q     !Discharge  [m^2/s]
      real*8 :: Zb    !DS bed elevation  [m]
      real*8 :: y0    !DS water depth  [m]
      real*8 :: Sb    !Bed slope 
      real*8 :: n     !Manning coefficient [m^0.5/s]
      real*8 :: g     !Gravity [m^2/s]
      real*8 :: L     !Domain length [m]
      real*8 :: dx    !Distance between two calculating points
!Variable declaration for processing the data
      integer :: t  !iterration 
      integer :: i  !index
      integer :: ODE !ODE method

      real*8 :: f     ! For function
      real*8 :: k1, k2, c   !Parameter for second order Runge Kutta scheme
      real*8 :: f1, f2, f3, f4 !Parameter for Runge Kutta scheme

      real*8 :: L2   !Norm2
      real*8 :: Lmax  !Norm max
      real*8 :: NA  !default parameter
!Array variable
      real*8, allocatable :: x(:)     !Position on the length  [m] 
      real*8, allocatable :: y(:)     !Water depth [m]
      real*8, allocatable :: z(:)     !Bed elevation [m]  
      real*8, allocatable :: u(:)     !Velocity  [m/s]
      real*8, allocatable :: ana_depth(:)  !analytical solution of water depth  [m]
      real*8, allocatable :: abs_error(:)    ! Absolute error
      real*8, allocatable :: rel_error(:)    ! Relative error
!Character declaration
      character (len=25) :: Result 
      character (len=25) :: Norm 

!reading the variable from input file
      open(11, file='input.txt')   
      read (11,*)
      read (11,*) L, q, y0, Zb, Sb, n, g, dx, k1, k2, c, ODE
      close(11)

!Number of cells (data generating points)
     t = int(L/dx)

!Generating arrays of different variables
      allocate(x(t+1)) 
      allocate(y(t+1))
      allocate(z(t+1))  
      allocate(u(t+1)) 
      allocate(ana_depth(t+1)) 
      allocate(abs_error(t+1))
      allocate(rel_error(t+1))

   !Initial Condition declaration
      x(1)=L  
      y(1)=y0 
      z(1)=Zb 
      u(1)=q/y0

!Calculation of Water Depth, Bed elevation and velocity based on different ODE schemes

      if (ODE==1) then           !For Euler method
         do i=1,t
            y(i+1)=y(i)-dx*f(y(i))
            x(i+1)=x(i)-dx
            z(i+1)=z(i)+Sb*dx
            u(i+1)=q/y(i+1)      
         end do
      else if (ODE==2) then    !For Runga Kutta order 2 scheme
         do i=1,t
            f1=f(y(i))
            f2=f(y(i)-c*dx*f1)
            y(i+1) = y(i) - dx*(k1*f1+k2*f2)
            x(i+1)=x(i)-dx
            z(i+1)=z(i)+Sb*dx
            u(i+1)=q/y(i+1)
         end do
      else if (ODE==4) then     !For Runge Kutta order 4 scheme
         do i=1,t
            f1=f(y(i))
            f2=f(y(i)-(0.5d0*dx*f1))
            f3=f(y(i)-(0.5d0*dx*f2))
            f4=f(y(i)-(dx*f3))
            y(i+1) = y(i) - (1.d0/6.d0*dx)*(f1+2*f2+2*f3+f4)
            x(i+1)=x(i)-dx
            z(i+1)=z(i)+Sb*dx
            u(i+1)=q/y(i+1)               
         end do
      end if

!Calculating the Norm L2 and LMax

      open(22, file='ANALYTICAL_1000CELLS.OUT') !Read(open) the value from analytical solution, Operator needs to change this with corresponding "dx" value
      L2=0
      Lmax=0
      do i=1,t+1
         read(22,*) NA, ana_depth(i)
         abs_error(i)=ana_depth(i)-y(i)
         rel_error(i)= (abs_error(i))/ ana_depth(i)   
      end do  
      L2=sqrt(sum(abs_error**2)/i)
      Lmax=maxval(abs(abs_error))
      close(22)

!Creating files for corresponding results
      if (ODE==1) then
         Result='Euler_Result.txt'
         Norm='Norm_euler.txt'
      else if (ODE==2) then
         Result='RK2_Result.txt'
         Norm='Norm_RK2.txt'
     else if (ODE==4) then
         Result='RK4_Result.txt'
         Norm='Norm_RK4.txt'
      end if
     
! Result generation 
      !This part will generate an output file with the value of Distance, Water Depth, Bed elevation and velocity
      open(33, file=Result, STATUS='REPLACE' )
      write(33,'(a8, a30, a30, a30)')  'x(m)' , 'z(m)', 'y(m)', 'u(m/s)'            
      do i=1,t+1
         write(33,*) x(i), z(i), y(i), u(i)
      end do
      close(33) 
      
      !This part will generate an output file with the Value of Norms and Errors
      open(44, file=Norm, STATUS='REPLACE')
      write(44,'(a15, a25)') 'L2 Norm', 'Lmax Norm'
      write(44,*) L2 , Lmax 
      write(44,*)
      write(44,'(a15, a25, a35, a30)')  'Y_exact', 'Y_numerical', 'Absolute Error', 'Relative Error'
      do i=1,t+1
         write(44,*) ana_depth(i), y(i) , abs_error(i), rel_error(i) 
      end do
      close(44)  

      deallocate(x) 
      deallocate(y)
      deallocate(z)  
      deallocate(u)  

end

!Declaration of an external function
real*8 function f(y)
real*8 :: y, yc, yn, q , Sb, n, g, NA

open(11, file='input.txt')   !take the values from the input file only those are required for the function
 read (11,*)
 read (11,*)   NA, q, NA, NA, Sb, n, g, NA, NA, NA, NA, NA
close(11)

yc=((q**2)/g)**(1.d0/3.d0)         !critical depth calculation
yn=((q*n/(Sb**0.5d0))**3)**0.2d0   !normal depth calculation

f=Sb*((1-(yn/y)**(10.d0/3.d0))/(1-(yc/y)**3)) !Function 

end function f