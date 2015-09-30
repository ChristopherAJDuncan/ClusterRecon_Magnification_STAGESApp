module cosmology
  use Param_Types; 
  implicit none
  
  character(50),private::Rad_Dir
  logical::Verbose


  type Cosmology_Parameters
     real(double)::Om_CDM, Om_B, Om_tot, Om_l, w, h_100
  end type Cosmology_Parameters


  INTERFACE angular_diameter_distance_fromRedshift
     MODULE PROCEDURE angular_diameter_distance_fromRedshift_array, angular_diameter_distance_fromRedshift_scalar
  END INTERFACE angular_diameter_distance_fromRedshift

  INTERFACE luminosity_distance
     MODULE PROCEDURE luminosity_distance_array, luminosity_distance_scalar
  END INTERFACE luminosity_distance
  
contains

  function LCDM
    type(Cosmology_Parameters)::LCDM
    
    LCDM%Om_CDM = 0.2544e0_double
    LCDM%Om_B = 0.0456e0_double
    LCDM%Om_l = 0.7e0_double
    LCDM%Om_tot = LCDM%Om_CDM+LCDM%Om_B+LCDM%Om_l
    LCDM%w = -1.e0_double
    LCDM%h_100 = 0.7e0_double

  end function LCDM

  subroutine print_Cosmology(Cos)
    type(Cosmology_Parameters), optional:: Cos
    type(Cosmology_Parameters):: iCos

    if(present(Cos)) then
       iCos = Cos
    else
       iCos = LCDM()
    end if

    print *, 'Cosmology Parameters as follows:'
    print *, iCos%Om_CDM 
    print *, iCos%Om_B
    print *, iCos%Om_l
    print *, iCos%Om_tot
    print *, iCos%w
    print *, iCos%h_100
    print *, '----------------------------------'
    print *, ' '

  end subroutine print_Cosmology


  function luminosity_distance_array(z1,z2,Cos)
        !-Returned in units of Mpc/h                                                                                                                                                                                                         
    type(Cosmology_Parameters),optional::cos
    real(double),intent(in)::z1
    real(double),intent(in)::z2(:)

    real(double),dimension(size(z2))::luminosity_distance_array

    real(double):: angular_diameter_distance(size(z2))

    type(Cosmology_Parameters)::Internal_Cos

    if(present(cos)) then
       Internal_Cos = Cos
    else
       Internal_Cos = LCDM()
    end if

    angular_diameter_distance =angular_diameter_distance_fromRedshift_array(z1,z2, Internal_Cos)
    luminosity_distance_array = ( ((1+z2)/(1+z1))**2.e0_double )*angular_diameter_distance

  end function luminosity_distance_array

  function luminosity_distance_scalar(z1,z2,Cos)
    !-Returned in units of Mpc/h                                                                                                                                                                                                         
    type(Cosmology_Parameters),optional::cos
    real(double),intent(in)::z1
    real(double),intent(in)::z2

    real(double)::luminosity_distance_scalar

    real(double)::lum_distance(1)

    lum_distance = luminosity_distance_array(z1,(/z2/),Cos)
    luminosity_distance_scalar = lum_distance(1)

  end function luminosity_distance_scalar

  function angular_diameter_distance_fromRedshift_array(z1,z2, Cos)
       !-Returned in units of Mpc/h                                         
    type(Cosmology_Parameters),optional::cos
    real(double),intent(in)::z1
    real(double),intent(in)::z2(:)

    real(double),dimension(size(z2))::angular_diameter_distance_fromRedshift_array

    type(Cosmology_Parameters)::Internal_Cos
    real(double),allocatable::Radius(:)

    if(present(cos)) then
       Internal_Cos = Cos
    else
       if(Verbose) print *, 'Using LCDM to get Distances'
       Internal_Cos = LCDM()
    end if

    allocate(Radius(size(z2)))
    Radius =  getrarray(z1, z2, Internal_Cos)
    
    !--This is flat case only--!
    !--Extra h_100 converts from Mpc to Mpc/h--!
    angular_diameter_distance_fromRedshift_array = (Radius*Internal_Cos%h_100)/(1.e0_double+z2)

    deallocate(Radius)


  end function angular_diameter_distance_fromRedshift_array

 function angular_diameter_distance_fromRedshift_scalar(z1, z2, cos)
   !-Returned in units of Mpc/h
   !--This is returned as a proper distance--!
    type(Cosmology_Parameters),optional::cos
    real(double),intent(in)::z1
    real(double),intent(in)::z2

    real(double)::angular_diameter_distance_fromRedshift_scalar

    type(Cosmology_Parameters)::Internal_Cos
    real(double),allocatable::Radius(:)

    if(present(cos)) then
       Internal_Cos = Cos
    else
       if(Verbose) print *, 'Using LCDM to get Distances'
       Internal_Cos = LCDM()
    end if

    allocate(Radius(1))
    Radius = getrarray(z1, (/z2/), Internal_Cos)
!    print *, '--Got angular diameter distance to lens:', z1, z2, Radius
!!$    if(Radius(1) <= 0.e0_double) then
!!$       print *, '---ERROR: Radius is zero  or negative:', Radius(1), z1, z2
!!$       call print_Cosmology(Internal_Cos)
!!$    end if

    !--Extra h_100 converts from Mpc to Mpc/h--!
    angular_diameter_distance_fromRedshift_scalar = (Radius(1)*Internal_Cos%h_100)/(1.e0_double+z2)

    deallocate(Radius)

  end function angular_diameter_distance_fromRedshift_scalar

  real(double) function Normalised_Hubble_Parameter(z, Cos)
    !-Returns Omega = rho(z)/pcrit(z = 0), which enters the Friedmann Equation, equivalent to the Hubble parameter renomrlised to it's presetn day value-!
    real(double), intent(in)::z
    type(cosmology_parameters),intent(in),optional:: Cos

    type(Cosmology_Parameters):: Internal_Cos

    real(double)::Omega_l(1)

    integer::callcount = 0

    callcount = callcount + 1
    if(present(cos)) then
       Internal_Cos = Cos
    else
       if(Verbose .and. callcount == 1) print *, 'Using LCDM to get Normalised Hubble Parameter'
       Internal_Cos = LCDM()
    end if
    
    Omega_l = Omega_lambda((/z/), Internal_Cos)
    Normalised_Hubble_Parameter = Omega_l(1) +(Internal_Cos%Om_CDM+Internal_Cos%Om_B)*((1+z)**3.e0_double) + (1-Internal_Cos%Om_tot)*((1+z)**2.e0_double)

  end function Normalised_Hubble_Parameter

  function zradiusintegrand(z, Cos)
    use nrtype;
    real(double),dimension(:),intent(in)::z
    type(Cosmology_Parameters),intent(in), optional:: Cos
    real(double),dimension(size(z))::zradiusintegrand
    type(Cosmology_Parameters):: iCos
    !needs generalised to account for how neutrinos vary with redshift
    
    if(present(Cos)) then
       iCos = Cos
    else
       iCos = LCDM()
    end if

    iCos%Om_tot = iCos%Om_l + iCos%Om_B + iCos%Om_CDM !+ Global_Cosmology_Parametersn

    if(iCos%Om_tot .eq. 0) then
       print *, 'zradiusintegrand - iCos%Om_tot equal to zero'
    end if
    if (iCos%h_100 .eq. 0) then
       print *, 'zradiusintegrand - h not set'
       STOP
    end if

    if(any(z .lt. 0)) then
       print *, 'zradiusintegrand - redshift arg less than zero'
       STOP
    END if
    
    !does this need generalised to account for omega_neutrinos
    zradiusintegrand = 3000.e0_DOUBLE/(iCos%h_100*dsqrt(Omega_lambda(z, iCos) +(iCos%Om_CDM+iCos%Om_B)*((1+z)**3.e0_DOUBLE) + (1-iCos%Om_tot)*((1+z)**2.e0_DOUBLE)))

    if(any(zradiusintegrand .lt. 0)) then
       print *, 'zradiusintegrand - less than zero'
       STOP
    end if
    

  end function zradiusintegrand 


  function Omega_lambda(z, Cos)
    use nrtype;
    real(DOUBLE),dimension(:),intent(in)::z
    type(Cosmology_parameters), intent(in):: Cos
    real(DOUBLE),dimension(size(z))::Omega_lambda
    real(DOUBLE), dimension(size(z))::w
    real(DOUBLE)::int
    integer::i
    integer::count = 0
    real(double)::wa = 0.e0_double

    count = count + 1

    Omega_lambda = Cos%Om_l*((1.e0_DOUBLE+z)**(3.e0_DOUBLE*(1+Cos%w)))
    if(wa .ne. 0.e0_double) then !this section needs tested
       STOP 'USING WA /= 0!!!'
       do i = 1, size(z)
          Omega_lambda(i) = Omega_lambda(i)*((1.e0_double+z(i))**(3.e0_double*wa))*dexp((z(i)*wa)/(1.e0_double+z(i)))
       end do
    end if

  end function Omega_lambda

  function Omega_lambda_int(z, Cos)
    !This isn't used yet, what is the integral wrt? Probably a or redshift - would not be used unless wa \= 0?
    use nrtype;
    real(DOUBLE),dimension(:),intent(in)::z
    type(Cosmology_Parameters), intent(in),optional:: Cos
    real(DOUBLE),dimension(size(z))::Omega_lambda_int

    real(DOUBLE), dimension(size(z))::w
    real(DOUBLE)::wa = 0.e0_DOUBLE !in generalisation will be made a global parameter
    type(Cosmology_Parameters):: iCos

    if(present(Cos)) then
       iCos = Cos
    else
       iCos = LCDM()
    end if

    w = iCos%w + (z/(1.e0_DOUBLE+z))*wa
    Omega_lambda_int = (-1.e0_DOUBLE)*(1.e0_DOUBLE + w)/(1 + z)

  end function Omega_lambda_int

  function getrarray(z1,z2,cos) RESULT(r)
    !calculates radius values for z values in zarray, and places them in appropriate slot of rarray
    !- Returned in units of Mpc-!
    !-Calculates CO-MOVING RADIUS [dChi or R_0dr]--!
    use nr; use Integration, only:Integrate_1D; use gridintervals, only:equalscale
    type(cosmology_parameters),optional::cos
    real(double),intent(in)::z1
    real(double)::z2(:)
    real(double), dimension(size(z2)):: r

    type(Cosmology_Parameters):: iCos
    integer::i 
    real(double)::ztemp
    character(len=200)::filename

    !--Intergration Decalrations--!
    real(double),allocatable:: zGrid(:), Integrand(:)

!!$    if(allocated(r)) then
!!$       deallocate(r)
!!$    end if
!!$    allocate(r(size(z2))); r = -1.e0_double



    !set cosmology params that integrand will use. If coss present, these will be set, otherwise we assume global cosmology is set appropriately
    if(present(cos)) then
       iCos = Cos
    else
       iCos = LCDM()
    end if

    !get radius for all z values in z, construct r


    r = 0.e0_double
    do i=1, size(z2)
       if(allocated(zGrid)) deallocate(zGrid)
       call equalscale(z1, z2(i), 300, zGrid)
       allocate(Integrand(size(zGrid))); Integrand = zradiusintegrand(zGrid, iCos)


       if(z2(i) >= z1) then
          r(i) = Simple_TrapInt(zGrid, Integrand)!Integrate_1D(zGrid, Integrand, 2)
          !r(i) = qromo(zradiusintegrand,z1, z2(i), midpnt, 1.e-6_double) !! Removed 20th Feb due to errors in multi-thread run
       else
          r(i) = 0.e0_double
       end if
    end do
    
    if(any(r<=0.e0_double)) then
       print *, 'Error determining distances from redshift - negatives still exist'
       print *, 'Redshifts:', z1, minval(z2), maxval(z2), size(z2)
       do i = 1, size(r)
          if(r(i) <= 0) print *, 'r-', r(i), z1, z2(i)
       end do
       STOP
    END if


  end function getrarray

  function Simple_TrapInt(x,y) RESULT(Area)
    real(double), intent(in):: x(:), y(:)
    real(double):: Area

    integer:: i

    if(size(x) /= size(y)) then
       print *, 'Simple_TrapInt - x and y not conformal'
    end if

    if(size(x) == 1) then
       Area = 0.e0_double
       return
    end if

    Area = 0
    do i = 1, size(y)-1
       Area = Area + 0.5e0_double*(y(i) + y(i+1))*(x(i+1)-x(i))
    end do

  end function Simple_TrapInt

  function convergent_TrapInt(func, limits, Tol) RESULT(Area)
    real(double):: limits(2), Area, Tol
    interface
       function func(x)
         use Param_Types
         real(double), intent(in):: x(:)
         real(double):: func
       end function func
    end interface

    real(double), dimension(1000):: Areas
    real(double), allocatable:: grid(:), integrand(:)
    integer:: nGrid
    integer:: loop, i, R

    integer,parameter:: convergence_limit = 4
    real(double), dimension(convergence_limit):: Diff

    Area = 0.
    Diff = 0.
    Areas = 0.

    nGrid = 2
    do loop = 1, size(Areas)
       nGrid = nGrid*2

       allocate(grid(nGrid))
       Grid = (/ (limits(1) + (i-1)*((limits(2)-limits(1))/(nGrid)), i = 1, nGrid) /)

       allocate(integrand(size(grid))); integrand = func(grid)
       
       Area = Simple_TrapInt(grid,integrand)

       Areas(loop) = Area

       deallocate(grid, integrand)

       if(loop >= convergence_limit+1) then
          Diff = 0.
          do R = loop-convergence_limit, loop
             Diff(R) = Areas(R)-Areas(R-1)
          end do
          if(all(dabs(Diff) <= Tol)) exit
       end if
    end do

    print *, 'Integral Converged on:', loop
    print *, Areas(1:loop)
    read(*,*)

  end function convergent_TrapInt
    

end module cosmology
