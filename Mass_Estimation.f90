
module Mass_Estimation
  !--Contains Subroutines that return mass estimates--!
  use Param_Types
  implicit none

  real(double),private::Default_Source_Redshift = 1.4e0_double !-Used to Calcualte Sigma_Critical when no redhift information is present-!

contains

!!$  subroutine draw_Circular_Aperture( Aperture_Positions, Aperture_Radiuses)
!!$    real(double), intent(in)::Aperture_Positions(:,:), Aperture_Radiuses(:)
!!$
!!$    character(20)::filename_header = 'draw_Aperture'
!!$    character::Cluster_String
!!$
!!$    integer::nSample = 1000
!!$    real(double)
!!$
!!$
!!$    do l = 1, size( Aperture_Positions,1 )
!!$       write(Cluster_String,'(I1)') l
!!$       open(48, file = trim(adjustl(filename_header))//'_'//trim(Cluster_String)//'.dat')
!!$       do n = 1, nSample
!!$          write(48,
!!$       end do
!!$       close(48)
!!$    end do
!!$
!!$
!!$  end subroutine draw_Circular_Aperture

  subroutine Recombine_Mass_Estimates(Mass, Error, Recombined_Mass, Recombined_Mass_Error)
    real(double),intent(in),dimension(:,:,:)::Mass, Error
    real(double),intent(out),dimension(:),allocatable::Recombined_Mass, Recombined_Mass_Error

    integer::Method = 2 !-1: Average, 2: Inverse Variance-!

    integer::i,j,m

    !--Method 2 Declarations--!
    real(double),allocatable::Weighting(:,:)
    real(double),allocatable::Sum_Weighting

    if(allocated(Recombined_Mass)) deallocate(Recombined_Mass)
    allocate(Recombined_Mass(size(Mass,3))); Recombined_Mass = 0.e0_double

    if(allocated(Recombined_Mass_Error)) deallocate(Recombined_Mass_Error)
    allocate(Recombined_Mass_Error(size(Error,3))); Recombined_Mass_Error = 0.e0_double

    if(size(Mass) /= size(Error)) STOP 'Recombine_Mass_Estimates - Mass and Error not the same size, stopping...'
    
    select case(Method)
    case(1)!-Average-!
       STOP 'Recombine_Mass_Estimates - Method 1 not used yet'
    case(2)!-Inverse Variance-!
       allocate(Weighting(size(Mass,1), size(Mass,2)))
       do m = 1,  size(Mass,3)    
          Recombined_Mass(m) = 0.e0_double
          Weighting = 0.e0_double
          do i= 1, size(Mass,1)
             do j = 1, size(Mass,2)
                Weighting(i,j) = 1.e0_double/(Error(i,j,m)*Error(i,j,m))
                Recombined_Mass(m) = Recombined_Mass(m) + Weighting(i,j)*Mass(i,j,m)
             end do
          end do
          Recombined_Mass(m) = Recombined_Mass(m)/sum(Weighting)
          Recombined_Mass_Error(m) = dsqrt(1.e0_double/sum(Weighting))
       end do
       deallocate(Weighting)
    case default
       STOP 'Recombine_Mass_Estimates - Method not valid'
    end select


  end subroutine Recombine_Mass_Estimates

!Projected_Mass, Error_Projected_Mass,
 subroutine Mass_Estimate_Circular_Aperture_byShift(Cat, Aperture_Positions, Aperture_Radiuses, Projected_Mass, Error_Projected_Mass,  global_average_size)
!  subroutine Mass_Estimate_Circular_Aperture_byShift(Cat, Aperture_Positions, Aperture_Radiuses,  global_average_size)
     use Catalogues; use Statistics, only: mean_discrete, variance_discrete; use Cosmology, only:angular_diameter_distance_fromRedshift; use Convergence_Estimation, only: get_Convergence

    !--Aperture Positions must be in (/RA, Dec/), and Aperture Radiuses in DEGREES--!
    type(Catalogue),intent(in)::Cat
    real(double), intent(in)::Aperture_Positions(:,:), Aperture_Radiuses(:)
    real(double),intent(in),optional::Global_Average_Size

    integer,allocatable::Expected_Number_in_Aperture(:)
    integer::i,j
    type(Catalogue),allocatable:: Ap_Cats(:)
    integer,dimension(size(Aperture_Positions,1))::Ap_Counter
    real(double),allocatable::Ap_Galaxy_Convergence(:,:)

    real(double)::Projected_Mass(size(Aperture_Positions,1)), Error_Projected_Mass(size(Aperture_Positions,1))

    real(double),allocatable::iAperture_Radiuses(:)
    real(double)::field_mean_size
    real(double):: Ap_Mean_Size, Ap_Convergence, Sigma_Crit, Ap_Area

    real(double)::D_l, D_ls , D_s, Lens_Redshift = 0.165e0_double, Source_Redshift
    real(double)::Redshift_Tolerance = 1.e-1_double

    character(10):: Mass_String, Error_Mass_String

    integer::Convergence_Estimator = 1 !-1:r/<r>-1, 2: ln(r/<r>)-!

    real(double),dimension(size(Aperture_Positions,1)):: Ap_Convergence_Error

    !--Precursors, decalre internals--!
    allocate(iAperture_Radiuses(size(Aperture_Positions,1)))
    if(size(Aperture_Radiuses) == 1) then
       iAperture_Radiuses = Aperture_Radiuses(1)
    elseif(size(Aperture_Radiuses) /= size(iAperture_Radiuses)) then
       STOP 'Mass_Estimate_Circular_Aperture_Catalogue - FATAL ERROR - Error assigning Aperture_Radiuses internal'
    else
       iAperture_Radiuses = Aperture_Radiuses
    end if

    if(present(global_average_size)) then
       field_mean_size = global_average_size
    else
       print *, 'Calculating the mass in circular aperture using mean physical size of the full catalogue'
       field_mean_size = mean_discrete(Cat%Physical_Sizes)
    end if

    !--Set up Aperture Catalogues--!
    allocate(Expected_Number_in_Aperture(size(Aperture_Positions,1))); Expected_Number_in_Aperture = 0

    allocate(Ap_Cats(size(Aperture_Positions,1)))
    do i = 1, size(Ap_Cats)
       Expected_Number_in_Aperture = count( dsqrt( (Cat%RA-Aperture_Positions(i,1))**2.e0_double +(Cat%Dec-Aperture_Positions(i,2))**2.e0_double ) <= iAperture_Radiuses(i))
       call Catalogue_Construct(Ap_Cats(i), Expected_Number_in_Aperture(i))
    end do

    !--Loop through galaxies in Catalogue and Assign to Relevent Ap_Catalogue--!
    Ap_Counter = 0
    do i = 1, size(Cat%Sizes)
       do j = 1, size(Ap_Cats)
          if( dsqrt( (Cat%RA(i)-Aperture_Positions(j,1))**2.e0_double +(Cat%Dec(i)-Aperture_Positions(j,2))**2.e0_double ) <= iAperture_Radiuses(j) ) then
             Ap_Counter(j) = Ap_Counter(j) + 1
             call Catalogue_Assign_byGalaxy_byCatalogue(Ap_Cats(j), Ap_Counter(j), Cat, i)
          end if
       end do
    end do
    !-Check by expected size-!
    do j =1, size(Ap_Cats)
       if(Ap_Counter(j) /= size(Ap_Cats(j)%RA)) print *, 'WARNING - Possible Error in assigning galaxies to Aperture:', j, ':Expected;Assigned:', size(Ap_Cats(j)%RA), Ap_Counter(j)
    end do
    !--Calculated Convergence Error for each Ap--!
    do j =1, size(Ap_Cats)
         !--Both Methods detailed below agree--!
         !-(Variance of Convergence per Galaxies in Ap)/(Number of Galaxies in Ap) :OR: Varaince of Mean of Convergence of Galaxies in Ap-!
       Ap_Convergence_Error(j) =  dsqrt(variance_discrete(get_Convergence(Ap_Cats(j)%Physical_Sizes,field_mean_size,Convergence_Estimator), get_Convergence(Ap_Cats(j)%Physical_Sizes,field_mean_size,Convergence_Estimator))/size(Ap_Cats(j)%RA))!dsqrt(variance_discrete((Ap_Cats(j)%Physical_Sizes/field_mean_size)-1.e0_double, (Ap_Cats(j)%Physical_Sizes/field_mean_size)-1.e0_double)/size(Ap_Cats(j)%RA))
         !-Variance of: [(Mean of Size in Ap)/(Mean of Size in field) - 1]-!
       !Ap_Convergence_Error(j) = dsqrt(variance_discrete(Ap_Cats(j)%Physical_Sizes, Ap_Cats(j)%Physical_Sizes))/(dsqrt(1.e0_double*size(Ap_Cats(j)%RA))*field_mean_size)
    end do


    !--Galaxies have now been split into Apertures, now statistics for each aperture--!
    Projected_Mass = 0.e0_double
    D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
    do j = 1, size(Ap_Cats)
       Ap_Mean_Size = mean_discrete(Ap_Cats(j)%Physical_Sizes)
       Ap_Convergence = get_Convergence(Ap_Mean_Size,field_mean_size, Convergence_Estimator)!(Ap_Mean_Size/field_mean_size) - 1.e0_double

       !--Convert from Convergence to Mass--!
       Source_Redshift =  mean_discrete(Ap_Cats(j)%Redshift)
       do i = 1, size(Ap_Cats(j)%Redshift)
          if(Ap_Cats(j)%Redshift(i) >= 0.e0_double) then
             Source_Redshift = Source_Redshift + Ap_Cats(j)%Redshift(i)
          else
             Source_Redshift = Source_Redshift + Default_Source_Redshift
          end if
       end do
       Source_Redshift = Source_Redshift/size(Ap_Cats(j)%Redshift)

       if( (Source_Redshift <= 0.e0_double) .or. (dabs(Source_Redshift-Lens_Redshift) <= Redshift_Tolerance)) STOP 'Mass_Estimate_Circular_Aperture_byShift - ERROR IN ASSIGNING SOURCE REDSHIFT'

       D_s = angular_diameter_distance_fromRedshift(0.e0_double, Source_Redshift)
       D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Source_Redshift)
       Sigma_Crit = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-! 
       
       Ap_Area = 3.142e0_double*( (D_l*(iAperture_Radiuses(j)*(3.142e0_double/180.e0_double)))**2.e0_double )

       Projected_Mass(j) = Ap_Convergence*Sigma_Crit*Ap_Area
       Error_Projected_Mass(j) = Sigma_Crit*Ap_Area*Ap_Convergence_Error(j)
    end do

    Projected_Mass = Projected_Mass*1.e18_double !--in MSun/h--!
    Error_Projected_Mass = Error_Projected_Mass*1.e18_double
!!$
!!$    print *, 'Projected Masses by Shift:'
!!$    do j = 1, size(Ap_Cats)
!!$       write(Mass_String, '(e8.2)') Projected_Mass(j); write(Error_Mass_String, '(e8.2)') Error_Projected_Mass(j)
!!$       print *, 'Cluster:', j, ' has mass: ', trim(Mass_String), ' +- ', trim(Error_Mass_String)
!!$    end do

    !--DESTROY--!
    do j = 1, size(Ap_Cats)
       call Catalogue_Destruct(AP_Cats(j))
    end do
    deallocate(Ap_Cats)

  end subroutine Mass_Estimate_Circular_Aperture_byShift


  subroutine Mass_Estimate_Circular_Aperture_byGalaxy(Cat, Profile, Aperture_Positions, Aperture_Radiuses, Projected_Mass, Error_Projected_Mass,global_average_size)
    use Catalogues; use Statistics, only: mean_discrete, variance_discrete; use Cosmology, only:angular_diameter_distance_fromRedshift; use Convergence_Estimation, only: get_Convergence
    !--Aperture Positions must be in (/RA, Dec/), and Aperture Radiuses in DEGREES--!
    !--Each Galaxy must have a redshift assigned, and physical size associated with it--!
    type(Catalogue),intent(in)::Cat
    integer::Profile !-0:Flat, 1:SIS, 2:NFW-!
    real(double), intent(in)::Aperture_Positions(:,:), Aperture_Radiuses(:)
    real(double),intent(in),optional::Global_Average_Size

    !--Internal Declarations--!
    real(double),allocatable::iAperture_Radiuses(:)
    real(double)::field_mean_size
    integer:: Ap, Gal
    real(double),allocatable::Galaxy_Convergence(:)
    real(double)::Convergence_Error, Mean_Convergence

    real(double)::Projected_Mass(size(Aperture_Positions,1)), Error_Projected_Mass(size(Aperture_Positions,1))

    real(double),allocatable::Weight(:,:), Sigma_Crit(:), Area(:), Renormalisation(:,:)
    real(double)::D_l, D_ls , D_s, Lens_Redshift = 0.165e0_double, Source_Redshift

    integer,allocatable:: Aperture_Galaxy_Indexs(:,:)
    integer,allocatable:: Aperture_Counter(:)

    character(10)::Mass_String, Error_Mass_String

    integer::Convergence_Estimator = 1 !-1:r/<r>-1, 2: ln(r/<r>)-!

    !---TESTING DECLARATIONS--!
    integer::nNoRedshift
    integer, dimension(size(Aperture_Positions,1)):: nGal_inAperture

    allocate(iAperture_Radiuses(size(Aperture_Positions,1)))
    if(size(Aperture_Radiuses) == 1) then
       iAperture_Radiuses = Aperture_Radiuses(1)
    elseif(size(Aperture_Radiuses) /= size(iAperture_Radiuses)) then
       STOP 'Mass_Estimate_Circular_Aperture_Catalogue - FATAL ERROR - Error assigning Aperture_Radiuses internal'
    else
       iAperture_Radiuses = Aperture_Radiuses
    end if

    if(any( (/0,1,2/) == Profile)==.false.) STOP 'Mass_Estimate_Circular_Aperture_Catalogue - FATAL ERROR - Supported profiles are 0:Flat, 1:SIS, 2:NFW.'
    
    if(present(global_average_size)) then
       field_mean_size = global_average_size
    else
       print *, 'Calculating the mass in circular aperture using mean physical size of the full catalogue'
       field_mean_size = mean_discrete(Cat%Physical_Sizes)
    end if
    
    !--Get Convergence for each galaxy--!
    allocate(Galaxy_Convergence(size(Cat%Physical_Sizes))); Galaxy_Convergence = 0.e0_double
    do Gal = 1, size(Cat%Physical_Sizes)
       Galaxy_Convergence(Gal) = get_Convergence(Cat%Physical_Sizes(Gal),field_mean_size, Convergence_Estimator)!(Cat%Physical_Sizes(Gal)/field_mean_size) - 1.e0_double
    end do

    Mean_Convergence = mean_discrete(Galaxy_Convergence)
    Convergence_Error = dsqrt(variance_discrete(Galaxy_Convergence, Galaxy_Convergence, Mean_Convergence, Mean_Convergence))

    print *, 'Mean/Variance [Convergence]:', Mean_Convergence, dsqrt(Convergence_Error)
    

    !--Calculate Sigma_Critical for each galaxy, and area in each aperture--!
    nNoRedshift = 0
    allocate(Sigma_Crit(size(Cat%Redshift)))
    D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift) !-In units of Mpc/h                                                                                                                                                         
    do Gal = 1, size(Cat%Redshift) 
       Source_Redshift = Cat%Redshift(Gal)
       if(Cat%Redshift(Gal) < 0.e0_double) then
          nNoRedshift = nNoRedshift + 1
          Source_Redshift = Default_Source_Redshift
       end if
          
       D_s = angular_diameter_distance_fromRedshift(0.e0_double, Source_Redshift)

       D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Source_Redshift)
       Sigma_Crit(Gal) = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!      
    end do
    allocate(Area(size(Aperture_Positions,1))); Area = 0.e0_double
    Area = 3.142e0_double*( (D_l*1.746e-2_double*iAperture_Radiuses)**2.e0_double )
    
    !--Create Indexes of Galaxies in Each Aperture--!
    allocate(Aperture_Galaxy_Indexs(size(Aperture_Positions,1), size(Cat%Physical_Sizes))) ; Aperture_Galaxy_Indexs = -100
    allocate(Aperture_Counter(size(Aperture_Positions,1))); Aperture_Counter = 0
    do Gal = 1, size(Cat%Physical_Sizes)
       do Ap = 1, size(Aperture_Positions,1)
          if( distance_between_points( (/Cat%RA(Gal), Cat%Dec(Gal)/), Aperture_Positions(Ap,:) ) < iAperture_Radiuses(Ap) ) then
             !--Galaxy Lies within Aperture--!
             Aperture_Counter(Ap) = Aperture_Counter(Ap) + 1
             Aperture_Galaxy_Indexs(Ap, Aperture_Counter(Ap)) = Gal
          end if
       end do
    end do

    if(nNoRedshift>0) print *, 'Mass_Estimate_Circular_Aperture_Catalogue - ', nNoRedshift, ' galaxies were assigned the default redshift information of ', Default_Source_Redshift, ' as no redshift information was input'

    nGal_inAperture = 0
    allocate(Weight(size(Aperture_Positions,1),size(Cat%Physical_Sizes))); Weight = 0.e0_double
    allocate(Renormalisation(size(Aperture_Positions,1),size(Cat%Physical_Sizes))); Renormalisation = 0.e0_double
    Projected_Mass = 0.e0_double
    do Ap = 1, size(Aperture_Positions,1)
       select case(Profile)
       case(0)
          !-Sum Convergence-!
          do Gal = 1, count(Aperture_Galaxy_Indexs(Ap,:) > 0)
             Projected_Mass(Ap) = Projected_Mass(Ap) + Galaxy_Convergence(Aperture_Galaxy_Indexs(Ap,Gal))
             Renormalisation(Ap,Gal) = 1.e0_double/(Sigma_Crit(Gal))
          end do
          Projected_Mass(Ap) = Projected_Mass(Ap)/sum(Renormalisation(Ap,:))
          Error_Projected_Mass(Ap) = (dsqrt(1.e0_double*count(Aperture_Galaxy_Indexs(Ap,:) > 0))*Convergence_Error)/sum(Renormalisation(Ap,:))
       end select
       print *, 'Aperture:', Ap, ' had ', count(Aperture_Galaxy_Indexs(Ap,:) > 0), ' galaxies within the Aperture'
    end do
    
!!$    nGal_inAperture = 0
!!$    allocate(Weight(size(Aperture_Positions,1),size(Cat%Physical_Sizes))); Weight = 0.e0_double
!!$    allocate(Renormalisation(size(Aperture_Positions,1),size(Cat%Physical_Sizes))); Renormalisation = 0.e0_double
!!$    Projected_Mass = 0.e0_double
!!$    do Ap = 1, size(Aperture_Positions,1)
!!$          do Gal = 1, size(Cat%Physical_Sizes)
!!$             if( distance_between_points( (/Cat%RA(Gal), Cat%Dec(Gal)/), Aperture_Positions(Ap,:) ) > iAperture_Radiuses(Ap) ) cycle
!!$             nGal_inAperture(Ap) = nGal_inAperture(Ap) + 1
!!$
!!$             select case(Profile)
!!$             case(0) !--Flat Profile--!
!!$                !--Sum Convergence Estimates--!
!!$                Projected_Mass(Ap) = Projected_Mass(Ap)
!!$                Renormalisation(Ap,Gal) = 
!!$
!!$!!!$       do Gal = 1, size(Cat%Physical_Sizes)
!!$!!!$          if( distance_between_points( (/Cat%RA(Gal), Cat%Dec(Gal)/), Aperture_Positions(Ap,:) ) > iAperture_Radiuses(Ap) ) cycle
!!$!!!$          nGal_inAperture(Ap) = nGal_inAperture(Ap) + 1
!!$!!!$          select case(Profile)
!!$!!!$          case(0)
!!$!!!$             
!!$!!!$
!!$!!!$!             Projected_Mass(Ap) = Projected_Mass(Ap) + (Galaxy_Convergence(Gal)/Sigma_Crit(Gal))
!!$!!!$!             Renormalisation(Ap,Gal) = 1.e0_double/(Sigma_Crit(Gal)*Sigma_Crit(Gal))
!!$!!!$!             Weight(Ap,Gal) = 1.e0_double/(Sigma_Crit(Gal))
!!$!!!$          case default
!!$!!!$             STOP 'Mass_Estimate_Circular_Aperture_Catalogue - FATAL ERROR - Supported profiles are 0:Flat'
!!$!!!$          end select
!!$!!!$       end do
!!$
!!$       print *, 'Aperture:', Ap, ' had ', nGal_inAperture(Ap), ' galaxies within the Aperture'
!!$    end do
    
    !--Renormalise--!
    select case(Profile)
    case(0)
!!$       do Ap = 1, size(Projected_Mass)
!!$          if(sum(Renormalisation(Ap,:)) /= 0.e0_double) then
!!$             Projected_Mass(Ap) =  Projected_Mass(Ap)/sum(Renormalisation(Ap,:))
!!$             Error_Projected_Mass(Ap) = Area(Ap)*dsqrt(Convergence_Error)/sum(Renormalisation(Ap,:))
!!$          else
!!$             print *, 'Aperture:', Ap, ' nGal in Aperture:', nGal_inAperture(Ap)
!!$             print*, 'Mass_Estimate_Circular_Aperture_Catalogue - Renormalisation (and Error) for Aperture: Renormalisation is zero. Enter to continue, Ctrl C to stop'
!!$             read(*,*)
!!$          end if
!!$       end do
    end select

    !--Convert from projected Surface Mass Density to Mass--!
    Projected_Mass = Area*Projected_Mass
    Error_Projected_Mass = Area*Error_Projected_Mass

    !--Convert to Units MSun/h--!
    Projected_Mass = Projected_Mass*1.e18_double
    Error_Projected_Mass = Error_Projected_Mass*1.e18_double

!!$    print *, '------------------------------------------------------------------------------'
!!$    print *, 'Masses by galaxy:'
!!$    do Ap = 1, size(Projected_Mass)
!!$       write(Mass_String, '(e8.2)') Projected_Mass(Ap); write(Error_Mass_String, '(e8.2)') Error_Projected_Mass(Ap)   
!!$       print *, 'Mass for cluster:', Ap, ' is: ', trim(Mass_String), ' +- ', trim(Error_Mass_String)
!!$    end do
!!$    print *, '------------------------------------------------------------------------------'

    deallocate(Aperture_Counter, Aperture_Galaxy_Indexs)

  end subroutine Mass_Estimate_Circular_Aperture_byGalaxy
    

  subroutine Mass_Estimate_CircularAperture_byPixel(Convergence_Map, Error_Convergence, x, y, Aperture_Positions, Aperture_Radiuses, Masses, Error_Masses, Source_Redshift, Convergence_Map_Occupation)
    use cosmology
    !-Returns the mass estimate with a circular Aperture--!
    !-Masses are calculated in a model independant way, as detailed in Heymasn et al 2008. Note that in that anaylsis, radius of 0.75' is used-!
    !-Convergence_Map, x, y are considered to be the same struture as those defined in Convergence_Estimates - that is x(i) defines the lower x values for bin i, x(i+1) the upper x value for bin i. ALSO ASSUMED TO BE EQUALLY BINNED

    !--x and y are in DEGREES. Radius should also be consistent (in DEGRESS).  **** dx = dy not necessary, provided it is accounted for on Apix  ****

    !-NOTE IF APERTUE GOES OUTSIDE X/Y RANGE, RESULT WILL NOT BE TRUSTWORTHY-!
    !-Source and Lens Redshift are hardwired. This would need generalised for code which has a different source redshift for each aperature-!
    !-Error_Convergence is used to measure the error on the mass estimates, and must be entered as an RMS value (i.e. sigma, not sigma^2)-!

    !---NOTE: THERE MAY BE A BUG HERE, AS TAKING THE AREA OUT OF THE SUM AND MULTIPLYING BY THE SUMMED TOTAL AREA GIVES DIFFERENT RESULTS. 15 NOV 2013

    real(double),intent(in)::Convergence_Map(:,:), Error_Convergence(:,:), x(:), y(:) !-x is 1st dimension, y is 2nd dimension-!
    real(double), intent(in)::Aperture_Positions(:,:) !-Labels centers of apertures, D1: Center Label, D2: 1: x, 2:y--!
    real(double), intent(in)::Aperture_Radiuses(:) !-Must have same size of D1 as Positions or 1, if 1 then same radius is applied for all apertures-!
    real(double), intent(out)::Masses(:), Error_Masses(:)
    real(double),intent(in),optional::Source_Redshift
    integer, intent(in),optional::Convergence_Map_Occupation(:,:)

    integer::i, j, x_index, y_index, x_rad_index, y_rad_index, x_index_range, y_index_range
    real(double)::Sigma_Crit = 1.e0_double
    real(double)::iSource_Redshift, Lens_Redshift = 0.165e0_double
    real(double)::D_l, D_s, D_ls
    real(double)::A_pix

    real(double),dimension(size(Aperture_Positions,1))::iAperture_Radiuses
    real(double)::nRandomPoints = 1000, NoRandomPoints_Success
    integer:: nPixel, nFullPixel, nEmptyPixel, nPartialPixel !-Testing - contains information on the number of pixels fully in aperture

    real(double)::Total_Area_Enclosed(size(Aperture_Positions,1))

    !-Random Number declarations-!
    real(double),dimension(:),allocatable::Ran(:,:)
    Integer,allocatable::seed(:)
    integer::Nseed, Clock

    real(double),dimension(size(Aperture_Positions,1)):: Sum_Convergence_in_Aperture, NPix_in_Aperture !-Used for testing only. Sums the convergence within the aperture, and the total fraction of pixels in aperture-!

    logical::Inverse_Variance_Weight = .false.
    real(double)::Sum_Weight
    real(double)::Sum_Convergence

    INTERFACE
       subroutine Mass_Estimate_CircularAperture(Convergence_Map, Error_Convergence, x, y, Aperture_Positions, Aperture_Radiuses, Masses, Error_Masses, Source_Redshift, Convergence_Map_Occupation)
         use Param_Types
         real(double),intent(in)::Convergence_Map(:,:), Error_Convergence(:,:), x(:), y(:) !-x is 1st dimension, y is 2nd dimension-!
         real(double), intent(in)::Aperture_Positions(:,:) !-Labels centers of apertures, D1: Center Label, D2: 1: x, 2:y--!              
         real(double), intent(in)::Aperture_Radiuses(:) !-Must have same size of D1 as Positions or 1, if 1 then same radius is applied for all apertures-!
         real(double), intent(out),allocatable::Masses(:), Error_Masses(:)
         
         real(double),intent(in),optional::Source_Redshift
         integer, intent(in),optional::Convergence_Map_Occupation(:,:)
       end subroutine Mass_Estimate_CircularAperture
    END INTERFACE

    print *, 'Estimating Masses:....'

    if(any(isNaN(Convergence_Map))) STOP 'Mass_Estimate_CircularAperture - FATAL ERROR - Convergence Map contains NaNs'
    
    iSource_Redshift = 1.4e0_double
    if(present(Source_Redshift)) iSource_Redshift = Source_Redshift
    print *, 'Estimating Mass with sources at z =', Source_Redshift, 'and lens at z = ', Lens_Redshift

    !--Input Error Catching--!
    if(size(x) /= size(Convergence_Map,1)) then
       print *, size(x), size(Convergence_Map,1)
       STOP 'Error - Mass_Estimate_CircularAperture - x and Kappa Map are not conformal, stopping..'
    end if
    if(size(y) /= size(Convergence_Map,2)) STOP'Error - Mass_Estimate_CircularAperture - x and Kappa Map are not conformal, stopping..'
    if((size(Aperture_Radiuses) /= 1) .and. (size(Aperture_Radiuses) /= size(Aperture_Positions,1))) then
       print *, 'Error - Mass_Estimate_CircularAperture - Aperture_Radiuses is not of the correct size, stopping..'
       print *, size(Aperture_Radiuses), size(Aperture_Positions,1)
       stop
    end if
    if(size(Aperture_Positions,2) /= 2) STOP 'Error - Mass_Estimate_CircularAperture - Aperture_Positions must contain only 2 positions (x,y) for each aperture, stopping..'


    if(size(Aperture_Radiuses)==1) then
       iAperture_Radiuses = Aperture_Radiuses(1)
    else
       iAperture_Radiuses = Aperture_Radiuses
    end if

    !-Determine Sigma_Critical-!
    D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift) !-In units of Mpc/h
    D_s = angular_diameter_distance_fromRedshift(0.e0_double, iSource_Redshift)
    !-Convert from comoving distances to physical distances-!
    !D_l = D_l/(1.e0_double+Lens_Redshift)
    !D_s = D_s/(1.e0_double+Source_Redshift)

    D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, iSource_Redshift)
    Sigma_Crit = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!

    print *, 'Calculating Masses with D_l = ', D_l
    print *, 'and Sigma_Crit = ', Sigma_Crit

    if(Verbose) print *, 'Sigma_Critical = ', Sigma_Crit
!!$    print *, 'Dl = ', D_l
!!$read(*,*)


!    allocate(Masses(size(Aperture_Positions,1))); Masses = 0.e0_double
!    allocate(Error_Masses(size(Aperture_Positions,1))); Error_Masses = 0.e0_double

    Masses = 0.e0_double; Error_Masses = 0.e0_double

    Sum_Convergence_in_Aperture = 0.e0_double; NPix_in_Aperture = 0.e0_double; Total_Area_Enclosed = 0.e0_double
    do i =1, size(Masses)
       if(size(Convergence_Map) == 1) then
          !-Single convergence estimate-!
          print *, 'Calculating Single Convergence Estimate'
          Masses(i) = 3.142e0_double*( (D_l*1.746e-2_double*Aperture_Radiuses(i))**2.e0_double )* Sigma_Crit * Convergence_Map(1,1)
          cycle
       end if

       Sum_Weight = 0.e0_double
       if( (Aperture_Positions(i,1)-iAperture_Radiuses(i) < x(1)) .or. (Aperture_Positions(i,1)+iAperture_Radiuses(i) > x(size(x))) ) print *, 'Warning - Aperture:', i, ' falls outside x range', iAperture_Radiuses(i), x(1), Aperture_Positions(i,1)-iAperture_Radiuses(i),' :',  x(size(x)), Aperture_Positions(i,1)+iAperture_Radiuses(i)
       if( (Aperture_Positions(i,2)-iAperture_Radiuses(i)< y(1)) .or. (Aperture_Positions(i,2)+iAperture_Radiuses(i) > y(size(y))) ) print *, 'Warning - Aperture:', i, ' falls outside y range'

       !--Find x and y positions of the center of the aperture--!
       do x_index = 1, size(x)-1
          if( (x(x_index) <= Aperture_Positions(i,1)) .and. (x(x_index+1) > Aperture_Positions(i,1)) ) then
             exit
          end if
       end do
       do y_index = 1, size(y)-1
          if( (y(y_index) <= Aperture_Positions(i,2)) .and. (y(y_index+1) > Aperture_Positions(i,2)) ) then
             exit
          end if
       end do

       !--Find the number of index points that x and y must cover (defines the boundary of the aperture)--!
       !--ASSUMES EQUAL BINNING IN X, AND ALSO IN Y--!
       x_index_range = int( (iAperture_Radiuses(i)/(x(2)-x(1)) )+2)
       y_index_range = int( (iAperture_Radiuses(i)/(y(2)-y(1)) )+2) !-+2 takes into account that x_index and y_index may not lie on apex

       !--Loop over all pixels within box defined by x/y_index_range--!
       nPixel = 0; nFullPixel = 0; nPartialPixel = 0; nEmptyPixel = 0
       do x_rad_index = maxval((/x_index-x_index_range,1/)), minval((/x_index+x_index_range,size(x)-1/)),1
          do y_rad_index = maxval((/y_index-y_index_range,1/)),  minval((/y_index+y_index_range,size(y)-1/)),1

             nPixel = nPixel + 1
             !--x/y_rad_index now loop over indexs of points within box defined by radius of aperture--!
             !-Calculate A_pix - MUST BE IN CORRECT UNITS, AND ACCOUNT FOR THE UNITS IN WHICH X AND Y ARE ENTERED. WHEN ==1, IT IS IN A PIXEL SCALE-!
             A_pix = (D_l*D_l*(x(x_rad_index+1)-x(x_rad_index))*(y(y_rad_index+1)-y(y_rad_index))*(6.283185e0_double/(360.e0_double))*(6.283185e0_double/(360.e0_double))) !-in (MPC/h)^2. Uses the fact that x and y are in DEGREES, and the small angle formula to convert into distance using steradians


             !--Throw down a set of random points within the pixel, test how many are in the aperture to  estimate the amount of the pixel to include in the mass estimate
             !----Construct Random Numbers----!
             allocate(Ran(2,nint(NRandomPoints))); Ran = 0.e0_double
             call RANDOM_SEED(size = NSeed)
             allocate(Seed(NSeed))
             call SYSTEM_CLOCK(COUNT = Clock)
             seed = Clock + (/ (i-1,i=1,NSeed) /)
             call RANDOM_SEED(PUT = seed)
             deallocate(Seed); NSeed = 0; Clock = 0
             
             call RANDOM_NUMBER(Ran)

             !-Convert Ran to  random positions within pixel-!
             Ran(1,:) = Ran(1,:)*(x(x_rad_index+1)-x(x_rad_index)) + x(x_rad_index) !-x-!
             Ran(2,:) = Ran(2,:)*(y(y_rad_index+1)-y(y_rad_index)) + y(y_rad_index) !-y-!

             NoRandomPoints_Success = 0
             do j = 1, nRandomPoints
                if(distance_between_points(Ran(:,j), Aperture_Positions(i,:)) <= iAperture_Radiuses(i)) then
                   NoRandomPoints_Success = NoRandomPoints_Success + 1
                end if
             end do
             if(NoRandomPoints_Success == 0) then 
                nEmptyPixel = nEmptyPixel + 1
             elseif(NoRandomPoints_Success==nRandomPoints) then
                nFullPixel = nFullPixel + 1
             else
                nPartialPixel = nPartialPixel + 1
             end if

             nPixel = nPixel + (1.e0_double*NoRandomPoints_Success/nRandomPoints)
             Total_Area_Enclosed(i) =  Total_Area_Enclosed(i) + A_pix*(1.e0_double*NoRandomPoints_Success/nRandomPoints)
             if(Inverse_Variance_Weight) then
                   !--Using the occupation numbers as a weight
                if(present(Convergence_Map_Occupation)) then
                   if(i == 1 .and. x_rad_index == maxval((/x_index-x_index_range,1/)) .and. y_rad_index == maxval((/y_index-y_index_range,1/)) ) print *, 'Getting masses by summing using occupation of pixel as weight'
                   if(Convergence_Map_Occupation(x_rad_index,y_rad_index) > 0) print *, Convergence_Map(x_rad_index,y_rad_index), Convergence_Map_Occupation(x_rad_index,y_rad_index)

                   Masses(i) = Masses(i) + ((1.e0_double*NoRandomPoints_Success/nRandomPoints)*Sigma_Crit)*(Convergence_Map(x_rad_index,y_rad_index)*Convergence_Map_Occupation(x_rad_index,y_rad_index))
                   Sum_Convergence_in_Aperture(i) = Sum_Convergence_in_Aperture(i) + (Convergence_Map(x_rad_index,y_rad_index)*Convergence_Map_Occupation(x_rad_index,y_rad_index))
                   Sum_Weight = Sum_Weight + Convergence_Map_Occupation(x_rad_index,y_rad_index)
                else
                   if(i == 1 .and. x_rad_index == maxval((/x_index-x_index_range,1/)) .and. y_rad_index == maxval((/y_index-y_index_range,1/)) ) print *, 'Getting masses by summing using inverse error as weight'
                   if(Error_Convergence(x_rad_index,y_rad_index) <= 0.e0_double) cycle
                   !--This bit is inverse variance weighting using the error bars. This is unsuccessful (possibly due to error bars being incorrect?).
                   Masses(i) = Masses(i) + ((1.e0_double*NoRandomPoints_Success/nRandomPoints)*Sigma_Crit)*(Convergence_Map(x_rad_index,y_rad_index)/(Error_Convergence(x_rad_index,y_rad_index)*Error_Convergence(x_rad_index,y_rad_index)))
                   Sum_Convergence_in_Aperture(i) = Sum_Convergence_in_Aperture(i) + (Convergence_Map(x_rad_index,y_rad_index)/(Error_Convergence(x_rad_index,y_rad_index)*Error_Convergence(x_rad_index,y_rad_index)))
                   Sum_Weight = Sum_Weight + 1.e0_double/((Error_Convergence(x_rad_index,y_rad_index)*Error_Convergence(x_rad_index,y_rad_index)))
                end if
                !--ERROR??--!
                if(i == 1 .and. x_rad_index == maxval((/x_index-x_index_range,1/)) .and. y_rad_index == maxval((/y_index-y_index_range,1/)) ) print *, 'I am not calculating errors on inverse variance weighting'
             else
                if(i == 1 .and. x_rad_index == maxval((/x_index-x_index_range,1/)) .and. y_rad_index == maxval((/y_index-y_index_range,1/)) ) print *, 'Using Model Independant Mass Estimation of Heymans et al'
                !Masses(i) = Masses(i) + A_pix*(1.e0_double*NoRandomPoints_Success/nRandomPoints)*Sigma_Crit*Convergence_Map(x_rad_index,y_rad_index)
                !if(Convergence_Map_Occupation(x_rad_index,y_rad_index) > 0) 

                !--TESTING--!
!!$                print *, 'For i, x, y:', i, x_rad_index, y_rad_index, ' K is:', Convergence_Map(x_rad_index,y_rad_index), Sigma_Crit
!!$                if(x_rad_index ==  minval((/x_index+x_index_range,size(x)-1/)) .and. y_rad_index == minval((/y_index+y_index_range,size(y)-1/)) ) read(*,*)
                

                Masses(i) = Masses(i) +  Sigma_Crit*Convergence_Map(x_rad_index,y_rad_index)
                Sum_Convergence_in_Aperture(i) = Sum_Convergence_in_Aperture(i) + Convergence_Map(x_rad_index,y_rad_index); NPix_in_Aperture(i) = NPix_in_Aperture(i) + (1.e0_double*NoRandomPoints_Success/nRandomPoints)
                !-Add Errors in Quadrature - Error_Masses is Sig^2 until sqrt taken later-!
                Error_Masses(i) = Error_Masses(i) +Sigma_Crit*(Error_Convergence(x_rad_index,y_rad_index)**2.e0_double)
             end if
             deallocate(Ran)
          END do
       END do
!       Masses(i) = Masses(i) * Total_Area_Enclosed(i)
       print *, 'nFullPixel:', nFullPixel, nEmptyPixel, nPartialPixel, nPixel


       if(Inverse_Variance_Weight) then 
          print *, 'Renormalising as using inverse variance weighting'
          Masses(i) = Masses(i) / Sum_Weight
          Sum_Convergence_in_Aperture(i) = Sum_Convergence_in_Aperture(i)/Sum_Weight
       end if

!!$       print *, 'Done Mass for:', i, ' reading'
!!$       print *, 'Summed Convergence = ', Sum_Convergence_in_Aperture(i)
!!$       print *, 'Expected Mass:', ( 3.142e0_double*(D_l*1.746e-2_double*iAperture_Radiuses(i))**2.e0_double )*Sigma_Crit* Sum_Convergence_in_Aperture(i), ' Measured:',  Masses(i), ' Ratio:', ( 3.142e0_double*(D_l*1.746e-2_double*iAperture_Radiuses(i))**2.e0_double )*Sigma_Crit* Sum_Convergence_in_Aperture(i)/Masses(i)
!!$       read(*,*)
    end do

    print *, 'Size:', size(Masses)
    print *, 'Masses and Area:', Total_Area_Enclosed, '::', Masses
    read(*,*)

    Masses = Total_Area_Enclosed*Masses
    Error_Masses = Total_Area_Enclosed*Error_Masses

    !--Convert to Mean--!
    Masses = Masses/nPixel
    Error_Masses = Error_Masses/nPixel

    if(verbose) then
       print *, 'Total Area enclosed in Apertures is:', Total_Area_Enclosed
       print *, 'Which is total fraction of the expected:', Total_Area_Enclosed/( 3.142e0_double*(D_l*1.746e-2_double*iAperture_Radiuses)**2.e0_double )
    end if


!!$    print *, 'Sum of Convergence within ap is:', Sum_Convergence_in_Aperture
!!$    print *, 'With NPix in each aperture of:', NPix_in_Aperture
!!$    read(*,*)

    !--Convert Error_Masses into RMS (sigma)--!
    if(any(Error_Masses < 0.e0_double)) STOP 'Mass_Estimate_CircularAperture - Error on Masses is imaginary, stopping..'
    Error_Masses = dsqrt(Error_Masses)

    print *, 'Finished.'

  end subroutine Mass_Estimate_CircularAperture_byPixel

  function distance_between_points( C1, C2 )
    !-Calculates the straight line distance between grid points, were C1 and C2 are (x,y) pairs-!
    real(double),intent(in)::C1(2), C2(2)

    real(double)::distance_between_points

    distance_between_points = dsqrt( ((C1(1)-C2(1))**2.) + ((C1(2)-C2(2))**2.) )

  end function distance_between_points


end module Mass_Estimation
