module Convergence_Estimation
  use Catalogues; use Param_Types; use FourierTransforms, only:smooth_gaussiankernal
  implicit none

  !--FUNCTION OVERLOADING ON BETA ROUTINES--!
  INTERFACE calculate_beta
     module procedure calculate_beta_fromCatalogue, calculate_beta_fromCatalogue_MagnitudeRedshiftBins, calculate_beta_fromAverageSize_MagnitudeRedshiftBins, calculate_beta_InterpolatefromGrid_File
  END INTERFACE calculate_beta
  INTERFACE correct_for_beta
     module procedure Correct_for_Beta_Scalar, Correct_for_Beta_TwoBins
  END INTERFACE correct_for_beta
  INTERFACE get_Convergence
     module procedure get_Convergence_scalar, get_Convergence_Vector
  END INTERFACE get_Convergence


  integer, private::Default_Convergence_Estimator = 1


contains 

  function get_Convergence_Vector(GSize, GlobalAvSize, Estimator)
    real(double),intent(in)::GSize(:)
    real(double),intent(in)::GlobalAvSize
    integer,intent(in)::Estimator

    real(double),dimension(size(GSize))::get_Convergence_Vector

    integer::i

    do i = 1, size(GSize)
       get_Convergence_Vector(i) = get_Convergence_scalar(GSize(i), GlobalAvSize, Estimator)
    end do

  end function get_Convergence_Vector

  real(double) function get_Convergence_scalar(GSize, GlobalAvSize, Estimator)
    real(double),intent(in)::GSize, GlobalAvSize
    integer,intent(in)::Estimator

    integer,save::callcount = 0
    integer,save::Previous_Method = -1
    logical::Talk_to_Me

    Talk_to_Me = .false.
    if( (Estimator/=Previous_Method) .or. (callcount == 1) ) Talk_to_Me = .true.

    select case (Estimator)
    case(1) !-r/<r>  - 1-!
       if(Talk_To_Me) print *, 'Using r/<r> - 1 as Estimator for convergence'
       get_Convergence_Scalar = (GSize/GlobalAvSize) - 1.e0_double
    case(2) !-ln(r/<r>)-!
       if(Talk_To_Me) print *, 'Using ln(r/<r>) as Estiamtor for convergence'
       get_Convergence_Scalar = dlog(GSize/GlobalAvSize)
    case default

    end select


    Previous_Method = Estimator
    callcount = callcount + 1 

  end function get_Convergence_Scalar

  subroutine construct_RA_Dec_Gridding(Cat,nRA, nDec, RAGrid, DecGrid)
    type(Catalogue),intent(in)::Cat
    integer,intent(in)::nRA, nDec
    real(double), intent(out),allocatable::RAGrid(:), DecGrid(:)

    real(double)::dn_RA, dn_Dec
    integer::i

    if(allocated(RAGrid)) deallocate(RAGrid)
    if(allocated(DecGrid)) deallocate(DecGrid)

    allocate(RAGrid(nRA+1)); RAGrid = 0.e0_double
    allocate(DecGrid(nDec+1)); DecGrid = 0.e0_double


    !--Define RA and Dec Grids such that nRA,nDec label number of bins, and bin i is between limits RA(i) - RA(i+1)--!
    !-Assign maximum and minimum by hand. Make last bin 10% larger to enclose the galaxies with teh largest RA and dec values
    RAGrid(1) = minval(Cat%RA)
    DecGrid(1) = minval(Cat%Dec)
    RAGrid(size(RAGrid)) = maxval(Cat%RA)+0.1e0_double*((maxval(Cat%RA)-minval(Cat%RA))/nRA)
    DecGrid(size(DecGrid)) = maxval(Cat%Dec)+0.1e0_double*((maxval(Cat%Dec)-minval(Cat%Dec))/nDec)
    dn_RA = ((RAGrid(size(RAGrid))-RAGrid(1))/(nRA))
    dn_Dec = ((DecGrid(size(DecGrid))-DecGrid(1))/(nDec))
    do i = 2, maxval((/nRA,nDec/))
       if(i<=nRA) RAGrid(i) = RAGrid(1) + (i-1)*dn_RA !!!
       if(i<=nDec) DecGrid(i) = DecGrid(1) + (i-1)*dn_Dec !!!
    end do

  end subroutine construct_RA_Dec_Gridding


  subroutine Average_Size_in_RA_Dec_Grid(cat, RAGrid, DecGrid, AverageSize, OccupationGrid, KappaEst, Smoothed_OccupationGrid, Size_Type, beta_correction)
    use Matrix_Methods
    type(Catalogue)::Cat
    real(double), intent(in),allocatable::RAGrid(:), DecGrid(:)
    real(double),intent(out),allocatable:: AverageSize(:,:)
    integer,allocatable::OccupationGrid(:,:)
    real(double), allocatable, optional:: Smoothed_OccupationGrid(:,:)
    real(double),intent(out), optional,allocatable::KappaEst(:,:)
    character(*), optional,intent(in)::Size_Type
    real(double), intent(inout),optional::beta_correction !-If present, beta correction is applied. If beta>1000, it is recalculated from the catalogue with a straight line fit to all points. 

    integer::i,j,c
    integer::nFail

    real(double)::global_mean_size
    real(double)::Average_of_Average_Size

    !-Smoothing Scale in Arcminutes-!
    real(double)::Smoothing_Scale = 0.75e0_double

    integer::callcount = 0

    INTERFACE
         subroutine Average_Size_in_RA_Dec_Grid(cat, RAGrid, DecGrid, AverageSize, OccupationGrid, KappaEst, Smoothed_OccupationGrid, Size_Type, beta_correction)
         USE CATALOGUES, ONLY:CATALOGUE; USE PARAM_TYPES
         type(Catalogue)::Cat
         real(double), intent(in),allocatable::RAGrid(:), DecGrid(:)
         real(double),intent(out),allocatable:: AverageSize(:,:)

         integer,allocatable::OccupationGrid(:,:)

         real(double),intent(out), optional,allocatable::KappaEst(:,:)
         real(double),allocatable::Smoothed_OccupationGrid(:,:)
          character(*), optional,intent(in)::Size_Type
          real(double), intent(inout),optional::beta_correction !-If present, beta correction is applied. If beta>1000, it is recalculated from the catalogue with a straight line fit to all points.  
       end subroutine Average_Size_in_RA_Dec_Grid
    END INTERFACE

    callcount = callcount + 1

    if(allocated(RAGrid)==.false. .or. allocated(DecGrid)==.false.) STOP 'FATAL ERROR - Average_Size_in_RA_Dec_Grid - Input RA or Dec grid are not allocated, stopping'

    if(present(Size_Type)) then
       call calculate_Average_Size_in_RA_Dec_Grid(Cat, RAGrid, DecGrid, AverageSize, OccupationGrid, by_size_Type = size_type)
    else
       call calculate_Average_Size_in_RA_Dec_Grid(Cat, RAGrid, DecGrid, AverageSize, OccupationGrid)
    end if

!!$  print *, 'Averages before Smoothing:', sum(AverageSize)/size(AverageSize), sum(Kappa)/size(Kappa)
!!$
!!$  print *, 'Averages after Smoothing:', sum(AverageSize)/size(AverageSize), sum(Kappa)/size(Kappa)

    if(present(KappaEst)) then
       if(allocated(KappaEst)) deallocate(KappaEst)
       allocate(KappaEst(size(AverageSize,1), size(AverageSize,2))); KappaEst = 0.e0_double
       
       call calculate_Convergence_in_RA_Dec_Grid(AverageSize, OccupationGrid, KappaEst)
       if(present(beta_correction)) then
          if(dabs(beta_correction) >= 1.e3_double) then
             print *, 'Calculating beta from catalogue...'
             call calculate_beta_fromCatalogue(Cat, beta_correction)
          end if
          if(callcount ==1) print *, 'I AM NOT CORRECTING FOR BETA'
          !call correct_for_beta(KappaEst, beta_correction)
       end if
       call do_Smoothing(Smoothing_Scale, RAGrid, DecGrid, Av_Size = AverageSize, Kappa = KappaEst, OccGrid = OccupationGrid)
    else
       
       call do_Smoothing(Smoothing_Scale, RAGrid, DecGrid, Av_size = AverageSize, OccGrid = OccupationGrid)
    end if

!    if(callcount == 1) print *, 'I AM NOT SMOOTHING'

    if(present(Smoothed_OccupationGrid)) then
       allocate(Smoothed_OccupationGrid(size(OccupationGrid,1), size(OccupationGrid, 2))); Smoothed_OccupationGrid = 1.e0_double*OccupationGrid
       call do_Smoothing(Smoothing_Scale, RAGrid, DecGrid, SOccGrid = Smoothed_OccupationGrid)
    end if
 
    if(any(isNaN(AverageSize))) STOP 'Average_Size_in_RA_Dec_Grid - FATAL ERROR - average size constains NaNs'
    if(any(isNaN(KappaEst)))STOP 'Average_Size_in_RA_Dec_Grid - FATAL ERROR - KappaEstimate constains NaNs'
   
  end subroutine Average_Size_in_RA_Dec_Grid

  subroutine do_Smoothing(ArcMinScale,RA_In,Dec_In,Av_Size,Kappa, OccGrid, SOccGrid)
    use Smoothing, only: Smooth_GaussianWeight_ConfigurationSpace_2D
    real(double),intent(in)::ArcMinScale
    real(double),intent(in)::RA_In(:),Dec_In(:)
    real(double),dimension(:,:),optional::Av_Size
    real(double),dimension(:,:),optional:: Kappa
    integer, dimension(:,:),optional::OccGrid
    real(double), dimension(:,:),optional::SOccGrid !-Needs to be a real if smoothing-!

    real(double)::Pixel_Scale, Degree_Scale

    real(double),allocatable,dimension(:,:):: S_Av_Size, S_Kappa

    INTERFACE
       subroutine do_Smoothing(ArcMinScale,RA_In,Dec_In,Av_Size,Kappa, OccGrid, SOccGrid)
         use Param_Types
         real(double),intent(in)::ArcMinScale
         real(double),intent(in)::RA_In(:),Dec_In(:)
         
         real(double),dimension(:,:),optional::Av_Size
         real(double),dimension(:,:),optional:: Kappa
         integer, dimension(:,:),optional::OccGrid
         real(double), dimension(:,:),optional::SOccGrid
       end subroutine do_Smoothing
    END INTERFACE

    Degree_Scale = ArcMinScale/60.e0_double

    if(present(Av_Size)) then
!       call smooth_gaussiankernal(Av_Size, Sigma = (/ArcMinScale/60.e0_double, ArcMinScale/60.e0_double/), Box_Length = (/RA_In(size(RA_In,1))-RA_In(1), Dec_In(size(Dec_In,1))-Dec_In(1)/))
       if((present(OccGrid) == .false.)) STOP 'do_Smoothing - Occ Grid not passed'
       call Smooth_GaussianWeight_ConfigurationSpace_2D(Av_Size, RA_In, Dec_In, Degree_Scale, S_Av_Size, OccGrid)
       Av_Size = S_Av_Size
       deallocate(S_Av_Size)
    end if
    if(present(Kappa)) then
!       call smooth_gaussiankernal(Kappa, Sigma = (/ArcMinScale/60.e0_double, ArcMinScale/60.e0_double/), Box_Length = (/RA_In(size(RA_In,1))-RA_In(1), Dec_In(size(Dec_In,1))-Dec_In(1)/))
       if((present(OccGrid) == .false.)) STOP 'do_Smoothing - Occ Grid not passed'
       call Smooth_GaussianWeight_ConfigurationSpace_2D(Kappa, RA_In, Dec_In, Degree_Scale, S_Kappa, OccGrid)
       Kappa = S_Kappa
       deallocate(S_Kappa)
    end if
    if(present(Av_Size) == .false. .and. present(Kappa) == .false. .and. present(SOccGrid)) then
       call smooth_gaussiankernal(SOccGrid, Sigma = (/ArcMinScale/60.e0_double, ArcMinScale/60.e0_double/), Box_Length = (/RA_In(size(RA_In,1))-RA_In(1), Dec_In(size(Dec_In,1))-Dec_In(1)/))
    end if

  end subroutine do_Smoothing

  subroutine calculate_Average_Size_in_RA_Dec_Grid(Cat, RA, Dec, AverageSize, nGrid_In, by_size_Type)
    !-Calculates the averages size in RA and Dec Bins, where RA and Dec has already been constructed
    !--If nGrid_In is present, then nGrid is alos output to file. Therefore when uses to bootstrap randomly assinged galaxy sizes, do not pass in an nGrid
    type(Catalogue),intent(in)::Cat
    real(double),allocatable::AverageSize(:,:)
    real(double),dimension(:),intent(in)::RA,Dec
    integer,allocatable,optionaL::nGrid_In(:,:)
    character(*), optional,intent(in)::by_Size_Type !-'Physical' or 'Pixel'. Default is Pixel-!

    integer,allocatable::nGrid(:,:)
    integer::i,j,c, nFail
    real(double),allocatable::Sizes(:) !-Temporary stores sizes information according to the by_size_Type passed in.-!

    INTERFACE
       subroutine calculate_Average_Size_in_RA_Dec_Grid(Cat, RA, Dec, AverageSize, nGrid_In, by_size_Type)
         use param_types; use Catalogues, only:Catalogue
         type(Catalogue),intent(in)::Cat
         real(double),allocatable::AverageSize(:,:)
         real(double),dimension(:),intent(in)::RA,Dec
         
         integer,allocatable,optionaL::nGrid_In(:,:)
          character(*), optional,intent(in)::by_Size_Type
       end subroutine calculate_Average_Size_in_RA_Dec_Grid
    END INTERFACE

    allocate(Sizes(size(Cat%Sizes))); Sizes = 0.e0_double
    if(present(by_size_Type)) then
       if(by_Size_Type == 'Physical' .or. by_Size_Type == 'physical' .or. by_Size_Type == 'Phys'.or. by_Size_Type == 'phys') then
!          print *, 'Binning by Physical Size..'
          Sizes = Cat%Physical_Sizes
       elseif(by_Size_Type == 'Pixel' .or. by_Size_Type == 'pixel' .or. by_Size_Type == 'Pix'.or. by_Size_Type == 'pix') then
!          print*, 'Binning by Pixel Size..'
          Sizes = Cat%Sizes
       else
          STOP 'calculate_Average_Size_in_RA_Dec_Grid - Input "by_Size_Type" is not valid, stopping...'
       end if
    else
!       print*, 'Binning by Pixel Size..'
       Sizes = Cat%Sizes
    end if

    if(allocated(AverageSize) .and. (size(AverageSize) /= size(RA)*size(Dec))) deallocate(AverageSize)
    if(allocated(AverageSize)==.false.) allocate(AverageSize(size(RA),size(Dec))); 
    AverageSize = 0.e0_double
    
    allocate(nGrid(size(AverageSize,1), size(AverageSize, 2))); nGrid = 0

    !--Assign Galaxies to correct bin, according to RA(i)<=RA<RA(i+1)--!
!!$    do c = 1, size(Cat%Sizes)
!!$!       print *, 'Doing Galaxy number:', c, 'of', size(Cat%Sizes)
!!$       i = int((Cat%RA(c)-RA(1))/dn_RA + 1)
!!$       j = int((Cat%Dec(c)-Dec(1))/dn_Dec + 1)
!!$!       print *, 'i:', i, Cat%RA(c), RA(1)
!!$       AverageSize(i,j) = AverageSize(i,j) + Cat%Sizes(c)
!!$       nGrid(i,j) = nGrid(i,j) + 1
!!$    end do

    do c = 1, size(Cat%Sizes)
       do i = 1, size(RA)-1
          if( (Cat%RA(c)>=RA(i)) .and. (Cat%RA(c)<RA(i+1))) exit
       end do
       do j = 1, size(Dec)-1
          if( (Cat%Dec(c)>=Dec(j)) .and. (Cat%Dec(c)<Dec(j+1))) exit
       end do
       !-Do i and j come out 1 too large?-!
       AverageSize(i,j) = AverageSize(i,j) + Sizes(c)
       nGrid(i,j) = nGrid(i,j) + 1
    end do

    nFail = 0
    do i =1 , size(nGrid,1)
       do j = 1, size(nGrid,2)
          if( (nGrid(i,j) > 0.e0_double) .and. (nGrid(i,j) < 1.e0_double)) then
             nFail = nFail + 1
             print *, 'nFail incremented'
          end if
       end do
    end do
    if(nFail > 0) STOP 'nFail too large'

    if(present(nGrid_In)) then
       allocate(nGrid_In(size(nGrid,1), size(nGrid,2))); nGrid_In = nGrid
    end if

    where(nGrid >= 1)
       AverageSize = AverageSize/(1.e0_double*nGrid)
    end where

    if(any(isnan(AverageSize))) STOP 'calculate_Average_Size_in_RA_Dec_Grid - FATAL ERROR -Average Size contains NaNs'

    deallocate(Sizes)

  end subroutine calculate_Average_Size_in_RA_Dec_Grid

  subroutine calculate_Convergence_in_RA_Dec_Grid(AverageSize, nGrid, Kappa_In)
    !--Returns the convergence, estimated as K = R/<R> - 1, for galaxies binned in RA and Dec, and that averages size taken in each bin
    real(double),intent(in)::AverageSize(:,:)
    integer,intent(in)::nGrid(:,:)
    real(double),allocatable::Kappa_In(:,:)

    real(double)::global_mean_size
    integer::i,j

    
    if(allocated(Kappa_In) .and. size(Kappa_In) /= size(AverageSize)) deallocate(Kappa_In)
    if(allocated(Kappa_In)==.false.) allocate(Kappa_In(size(AverageSize,1), size(AverageSize,2)))
    Kappa_In = 0.e0_double

    global_mean_size = sum(AverageSize*nGrid)/sum(nGrid)!sum(AvSize)/count(nGrid>=1)!sum(AvSize*nGrid)/sum(nGrid)
    !print *, 'Global mean size used in Kappa_In calculation is:', global_mean_size

    do i =1, size(nGrid,1)
       do j = 1, size(nGrid,2)
          if(nGrid(i,j) <= 0) cycle
          Kappa_In(i,j) = get_Convergence(AverageSize(i,j), global_mean_size, Default_Convergence_Estimator)!(AverageSize(i,j)/global_mean_size) - 1.e0_double
       end do
    end do

!!$    where(nGrid >= 1)
!!$       Kappa_In = (AverageSize/global_mean_size) - 1.e0_double
!!$    end where
!!$    print *, 'Number of Kappa_In > 0:', count(Kappa_In>0), 'Number Kappa_In > 1:', count(Kappa_In > 1), '# < 0', count(Kappa_In<0.), ' and # = 0', count(Kappa_In==0.)
!!$    print *, '# empty cells:', count(nGrid < 1)
!!$    print *, 'Kappa_In has means:', sum(Kappa_In)/count(nGrid>=1), sum(Kappa_In)/size(nGrid)

  end subroutine calculate_Convergence_in_RA_Dec_Grid

  subroutine Average_Size_in_RA_Dec_Grid_Errors(Cat_In, global_av_size, RA,Dec, Error_AvSize, Error_Kappa, by_Size_Type)
    use Catalogues, only: Catalogue_Destruct
    !-Calculates the errors per pixel by randomly assinging galaxy sizes to galaxy positions-!
    !-Errors output as sigma, NOT sigma^2
    type(Catalogue),intent(in)::Cat_In
    real(double),allocatable,dimension(:,:)::Error_AvSize, Error_Kappa
    real(double),intent(in),dimension(:),allocatable::RA,Dec
    real(double),intent(in)::global_av_size
    character(*), intent(in),optional:: by_Size_Type !-Phyiscal/Phys, Pixel/Pix-!


    type(Catalogue)::Randomised_Cat
    integer::NBoot
    integer::n,i,j,k
    real(double),allocatable,dimension(:,:,:)::Boot_AvSize, Boot_Kappa
    real(double),allocatable,dimension(:,:)::tBoot_AvSize, tBoot_Kappa
    integer,allocatable:: tBoot_nGrid(:,:)
    real(double),allocatable,dimension(:,:)::Mean_Boot_AvSize, Mean_Boot_Kappa
    integer,allocatable::Boot_nGrid(:,:,:)

    !-Random Number Declarations
    real(double),dimension(:,:),allocatable::Ran
    Integer,allocatable::seed(:)
    integer::Nseed, Clock

    !--Used to test for convergence of BootStrap by Varying nBoot 
    logical::Test_Convergence = .false.
    integer::Convergence_Loop, Convergence_Loop_Maximum
    integer::NBoot_Convergence_Lower
    real(double),allocatable:: Previous_Error_AvSize(:,:), Previous_Error_Kappa(:,:)
    real(double)::Convergence_Tolerance = 1.e-6_double
    logical::Kappa_Error_Convergence, AvSize_Error_Convergence

    INTERFACE
       subroutine Average_Size_in_RA_Dec_Grid_Errors(Cat_In, global_av_size, RA,Dec, Error_AvSize, Error_Kappa, by_Size_Type)
         use Catalogues; use Param_Types
          type(Catalogue),intent(in)::Cat_In
          real(double),allocatable,dimension(:,:)::Error_AvSize, Error_Kappa
          real(double),intent(in),dimension(:),allocatable::RA,Dec
          real(double),intent(in)::global_av_size
          
          character(*), intent(in),optional::by_Size_Type !-Phyiscal/Phys, Pixel/Pix-!     
        end subroutine Average_Size_in_RA_Dec_Grid_Errors
     END INTERFACE

    print *, 'Bootstrapping errors...'

    if(Test_Convergence) then
       Convergence_Loop_Maximum = 6
       allocate(Previous_Error_AvSize(size(RA), size(Dec))); Previous_Error_AvSize = 0.e0_double
       allocate(Previous_Error_Kappa(size(RA), size(Dec))); Previous_Error_Kappa = 0.e0_double
       print *, 'Testing for Convergence:'
       Verbose = .false.
    else
       Convergence_Loop_Maximum = 1
    end if

    allocate(Error_AvSize(size(RA), size(Dec))); Error_AvSize = 0.e0_double
    allocate(Error_Kappa(size(RA), size(Dec))); Error_Kappa = 0.e0_double

    NBoot = 512
    !-Calculate array or random numbers, consisting of NBoot sets of random samplings (with replacement) from the entered galaxy size distribution
    do Convergence_Loop = 1, Convergence_Loop_Maximum
       if(Test_Convergence) NBoot = 64*(2**Convergence_loop)
       
       
       allocate(Ran(NBoot, size(Cat_In%Sizes))); Ran = 0.e0_double
       call RANDOM_SEED(size = NSeed)
       allocate(Seed(NSeed))
       call SYSTEM_CLOCK(COUNT = Clock)
       seed = Clock + (/ (i-1,i=1,NSeed) /)
       call RANDOM_SEED(PUT = seed)
       deallocate(Seed); NSeed = 0; Clock = 0
       
       call RANDOM_NUMBER(Ran)
       !-Convert the random numbers from range 0-1 to 1-Ngal
       Ran = nint(Ran*(size(Cat_In%Sizes)-1) + 1)
       
       !-Construct Average Size grid and Kappa for each bootstrap realisation
       allocate(Boot_AvSize(NBoot, size(RA), size(Dec))); Boot_AvSize = 0.e0_double
       allocate(Boot_Kappa(size(Boot_AvSize,1), size(Boot_AvSize,2), size(Boot_AvSize,3))); Boot_Kappa = 0.e0_double
       allocate(Boot_nGrid(size(Boot_AvSize,1), size(Boot_AvSize,2), size(Boot_AvSize,3))); Boot_nGrid = 0
       do n = 1, NBoot
          !--Set up the Randomised Catalogue--!
          Randomised_Cat = Cat_In; Randomised_Cat%Sizes = 0.e0_double; Randomised_Cat%Physical_Sizes  = 0.e0_double
          do k = 1, size(Cat_In%Sizes)
             Randomised_Cat%Sizes(k) = Cat_In%Sizes(Ran(n,k)); Randomised_Cat%Physical_Sizes(k) = Cat_In%Physical_Sizes(Ran(n,k))
          end do
          if(present(by_Size_Type)) then
             call Average_Size_in_RA_Dec_Grid(Randomised_Cat, RA, Dec, tBoot_AvSize, OccupationGrid = tBoot_nGrid, KappaEst = tBoot_Kappa, Size_Type = by_Size_Type)
          else
             call Average_Size_in_RA_Dec_Grid(Randomised_Cat, RA, Dec, tBoot_AvSize, OccupationGrid = tBoot_nGrid, KappaEst = tBoot_Kappa)
          end if
          
!!$       call calculate_Average_Size_in_RA_Dec_Grid(Randomised_Cat, RA, Dec, tBoot_AvSize, tBoot_nGrid)
!!$       call calculate_Convergence_in_RA_Dec_Grid(tBoot_AvSize, tBoot_nGrid, tBoot_Kappa)
!!$       call do_Smoothing(0.75e0_double,RA,Dec,tBoot_AvSize, tBoot_Kappa)
          Boot_AvSize(n,:,:) = tBoot_AvSize; tBoot_AvSize = 0.e0_double; deallocate(tBoot_AvSize)
          Boot_nGrid(n,:,:) = tBoot_nGrid; tBoot_nGrid = 0; deallocate(tBoot_nGrid)
          Boot_Kappa(n,:,:) = tBoot_Kappa; tBoot_Kappa = 0.e0_double; deallocate(tBoot_Kappa)
          call Catalogue_Destruct(Randomised_Cat)
       end do
       
       !-Bootstrap sampling are not complete - calculate variance of each.
       allocate(Mean_Boot_AvSize(size(Boot_AvSize,2), size(Boot_AvSize,3))); Mean_Boot_AvSize = 0.e0_double
       allocate(Mean_Boot_Kappa(size(Boot_Kappa,2), size(Boot_Kappa,3))); Mean_Boot_Kappa= 0.e0_double

       
       DO n =1, nBoot
          Mean_Boot_AvSize = Mean_Boot_AvSize + Boot_AvSize(n,:,:)
          Mean_Boot_Kappa = Mean_Boot_Kappa + Boot_Kappa(n,:,:)
       end DO
       Mean_Boot_AvSize = Mean_Boot_AvSize/nBoot
       Mean_Boot_Kappa = Mean_Boot_Kappa/nBoot
       
       !print *, 'Mean mean kappa =', sum(Mean_Boot_Kappa)/size(Mean_Boot_Kappa)
       
       !    Mean_Boot_AvSize = sum(Boot_AvSize,1)/size(Boot_AvSize,1)
       !    Mean_Boot_Kappa = sum(Boot_Kappa,1)/size(Boot_Kappa,1)

       Error_AvSize = 0.e0_double; Error_Kappa = 0.e0_double
       do n = 1, nBoot
          !-Sum_i^nBoot (x^i - <x>^i)**2 over 2 dimensions simultaneously
          Error_AvSize = Error_AvSize + (Boot_AvSize(n,:,:) - Mean_Boot_AvSize)**2.
          Error_Kappa = Error_Kappa + (Boot_Kappa(n,:,:) - Mean_Boot_Kappa)**2.
       end do
       
       !-Error output as sigma, NOT sigma^2-!
       Error_AvSize = dsqrt(Error_AvSize/(1.e0_double*(nBoot-1)))
       Error_Kappa = dsqrt(Error_Kappa/(1.e0_double*(nBoot-1)))
       
       !--Set large Error for pixels without information--!
       do i = 1, size(Error_AvSize,1)
          do j = 1, size(Error_Kappa,1)
             if(all(Boot_nGrid(:,i,j) == 0)) then
                !--Set to negative?-!
!                Error_AvSize(i,j) = 1.e30_double; Error_Kappa(i,j) = 1.e30_double
             end if
          end do
       end do
 
       deallocate(Ran, Boot_nGrid, Boot_Kappa, Boot_AvSize, Mean_Boot_AvSize, Mean_Boot_Kappa)

       if(Test_Convergence) then
          if(Convergence_Loop == 1) then
             Previous_Error_AvSize = Error_AvSize
             Previous_Error_Kappa = Error_Kappa
             AvSize_Error_Convergence = .false.; Kappa_Error_Convergence = .false.
             cycle
          end if

          if(all(Error_AvSize-Previous_Error_AvSize <= Convergence_Tolerance) .and. AvSize_Error_Convergence==.false.) then
             !--Success - Converged for all bins--!                                                                                                                                                                                          
             print *, 'Errors in Average Size have converged to tolerance:', Convergence_Tolerance, ' with nBoot = ', nBoot/2
             print *, 'Press Enter to continue:'
!             read(*,*)
             AvSize_Error_Convergence = .true.
             Verbose = .true.
          end if
          if(all(Error_Kappa-Previous_Error_Kappa <= Convergence_Tolerance) .and. Kappa_Error_Convergence==.false.) then
             !--Success - Converged for all bins--!
             print *, 'Errors in Kappa have converged to tolerance:', Convergence_Tolerance, ' with nBoot = ', nBoot/2
             print *, 'Press Enter to continue:'
!             read(*,*)
             Kappa_Error_Convergence = .true.
             Verbose = .true.
          end if
          if(AvSize_Error_Convergence .and. Kappa_Error_Convergence) exit

          if(Convergence_Loop == Convergence_Loop_Maximum) then
             print *, 'Failed to find convergence by nBoot = ', nBoot
          END if

          print *, (count(Error_avSize-Previous_Error_AvSize <= Convergence_Tolerance)*100.e0_double/size(Error_avSize)), '% of AverageSize bins have converged for nBoot = ', nBoot/2,'. Relooping'
          Previous_Error_avSize = Error_AvSize
          print *, (count(Error_Kappa-Previous_Error_Kappa <= Convergence_Tolerance)*100.e0_double/size(Error_Kappa)), '% of  Kappa bins have converged for nBoot = ', nBoot/2,'. Relooping'
          Previous_Error_Kappa = Error_Kappa
          Error_Kappa = 0.e0_double; Error_AvSize = 0.e0_double

       end if

    end do
    if(Test_Convergence)deallocate(Previous_Error_avSize, Previous_Error_Kappa)

    print *, 'Finished Bootstrap'

  end subroutine Average_Size_in_RA_Dec_Grid_Errors

  subroutine Average_Size_in_Luminosity_Redshift_Bins(Cat, Mag_Limits, Red_Limits, SLR)
    !--Returns the average size in bins of luminosity and ABSOLUTE magnitude--!
    type(Catalogue), intent(in)::Cat
    real(double),dimension(:,:)::Red_Limits, Mag_Limits
    real(double),intent(out)::SLR(:,:)

    type(Binned_Catalogue)::BCat_Red
    type(Binned_Catalogue),allocatable::BCat_Lum_Red(:) !-One for each redshift Bin-!
    integer::BinL, BinZ

    character(10)::by_Size_Type = 'Physical'

    if(Verbose) print *, 'Calculating Average Size in Redshift and Luminosity Bins....'

    !-Check Sizes-!
    if(size(SLR,1) /= size(Mag_Limits,1)) stop 'Average_Size_in_Luminosity_Redshift - SLR and Luminosity Bins not conformal...'
    if(size(SLR,2) /= size(Red_Limits,1)) stop 'Average_Size_in_Luminosity_Redshift - SLR and Redshift Bins not conformal...'

    !--Bin first by redshift using pre-defined Redshift--!
    call bin_catalogue_by_redshift(Cat, Red_Limits, BCat_Red)
    print *, 'Average_Size_in_Luminosity_Redshift_Bins - binned in redshift'; read(*,*)

    allocate(BCat_Lum_Red(size(Red_Limits,1)))
    do BinZ = 1, size(BCat_Lum_Red)
       call Bin_Catalogue_by_magnitude(BCat_Red%Cat(BinZ), Mag_Limits, BCat_Lum_Red(BinZ))
    end do

    print *, 'Average_Size_in_Luminosity_Redshift_Bins - binned in magnitude'; read(*,*)

    if(Verbose) print *, 'Finished Binning'

    do BinZ = 1, size(Red_Limits,1)
       do BinL = 1, size(Mag_Limits,1)
          SLR(BinL,BinZ) = global_mean_size(BCat_Lum_Red(BinZ)%Cat(BinL), trim(adjustl(by_Size_Type)))
          if(isNAN(SLR(BinL,BinZ))) SLR(BinL,BinZ) = 0.e0_double
       end do
    end do

     print *, 'Done.'; read(*,*)

    if(Verbose) print *, 'Done.'

  end subroutine Average_Size_in_Luminosity_Redshift_Bins

  subroutine Average_Size_in_Luminosity_Redshift_Errors(Cat, Lum_Limits, Red_Limits, SLR_ERROR)
    !--By Bootstrap--!
    type(Catalogue), intent(in)::Cat
    real(double),dimension(:,:)::Red_Limits, Lum_Limits
    real(double),intent(out)::SLR_ERROR(:,:)
                                                                                            
    !--Random Catalogue Declarations--!
    real(double),dimension(:,:),allocatable::Ran
    Integer,allocatable::seed(:)
    integer::Nseed, Clock
    
    type(Catalogue)::Randomised_Cat
    integer::NBoot
    integer::n,i,j,k
    integer::Boot_Loop
    real(double),allocatable,dimension(:,:,:)::Boot_AvSize
    real(double),allocatable,dimension(:,:)::tBoot_AvSize
    real(double),allocatable::mean_boot_avsize(:,:)

    !--Used to test for convergence of BootStrap by Varying BootStrap
    logical::Test_Convergence = .false.
    integer::Convergence_Loop, Convergence_Loop_Maximum
    integer::NBoot_Convergence_Lower
    real(double),allocatable:: Previous_Error(:,:)
    real(double)::Convergence_Tolerance = 1.e-6_double

    if(Verbose) print *, 'Calculating Bootstrap Error Average Size in Redshift and Luminosity Bins....'

    if(size(SLR_ERROR,1) /= size(Lum_Limits,1)) stop 'Average_Size_in_Luminosity_Redshift - SLR_ERROR and Luminosity Bins not conformal...'
    if(size(SLR_ERROR,2) /= size(Red_Limits,1)) stop 'Average_Size_in_Luminosity_Redshift - SLR_ERROR and Redshift Bins not conformal...'

    if(Test_Convergence) then
       Convergence_Loop_Maximum = 20
       allocate(Previous_Error(size(SLR_Error,1), size(SLR_Error,2))); Previous_Error = 0.e0_double
       print *, 'Testing for Convergence:'
       Verbose = .false.
    else
       Convergence_Loop_Maximum = 1
    end if

    NBoot = 1024;

    do Convergence_Loop = 1, Convergence_Loop_Maximum
       if(Test_Convergence) NBoot = 64*(2**Convergence_loop)

       allocate(Ran(NBoot, size(Cat%Sizes))); Ran = 0.e0_double
       call RANDOM_SEED(size = NSeed)
       allocate(Seed(NSeed))
       call SYSTEM_CLOCK(COUNT = Clock)
       seed = Clock + (/ (i-1,i=1,NSeed) /)
       call RANDOM_SEED(PUT = seed)
       deallocate(Seed); NSeed = 0; Clock = 0
       
       call RANDOM_NUMBER(Ran)
       !-Convert the random numbers from range 0-1 to 1-Ngal                                                                                                                          
       Ran = nint(Ran*(size(Cat%Sizes)-1) + 1)
       
       allocate(tBoot_AvSize(Size(Lum_Limits,1), size(Red_Limits,1)))
       allocate(Boot_AvSize(NBoot, size(tBoot_AvSize,1), size(tBoot_AvSize,2))); Boot_AvSize = 0.e0_double
       do Boot_loop = 1, NBoot
          !--Set up Randomised Catalogue--!
          Randomised_Cat = Cat; Randomised_Cat%Sizes = 0.e0_double
          do k = 1, size(Cat%Sizes)
             Randomised_Cat%Sizes(k) = Cat%Sizes(Ran(Boot_Loop,k))
          end do
          call convert_Size_from_Pixel_to_Physical(Randomised_Cat)
          
          tBoot_AvSize = 0.e0_double
          call Average_Size_in_Luminosity_Redshift_Bins(Randomised_Cat, Lum_Limits, Red_Limits, tBoot_AvSize) 
          Boot_AvSize(Boot_Loop,:,:) = tBoot_AvSize; tBoot_AvSize = 0.e0_double
          
          call Catalogue_Destruct(Randomised_Cat)
       end do
       deallocate(tBoot_AvSize)
       
       allocate(Mean_Boot_AvSize(size(Boot_AvSize,2), size(Boot_AvSize,3))); Mean_Boot_AvSize = 0.e0_double
       do Boot_Loop = 1, NBoot
          Mean_Boot_AvSize = Mean_Boot_AvSize(:,:) + Boot_AvSize(Boot_Loop,:,:)
       end do
       Mean_Boot_AvSize = Mean_Boot_AvSize/nBoot
       
       !--Calculate Variance--!
       do Boot_Loop = 1, NBoot
          SLR_Error = SLR_Error + (Boot_AvSize(Boot_Loop,:,:) - Mean_Boot_AvSize)**2.e0_double
       end do
       SLR_Error = dsqrt(SLR_Error/(1.e0_double*(nBoot-1)))

       !--Allocate large error for bins with zero galaxies--!
       !--Assumes zero average size is applicable only for bins with no information--!
       do i = 1, size(Boot_AvSize,2)
          do j = 1, size(Boot_AvSize,3)
!             if(all(Boot_AvSize(:,i,j) == 0.e0_double)) SLR_Error(i,j) = 1.e30_double
          end do
       end do

       deallocate(Ran, Boot_AvSize, Mean_Boot_AvSize)

       if(Test_Convergence) then
          if(Convergence_Loop == 1) then
             Previous_Error = SLR_Error
             cycle
          end if

          if(all(SLR_Error-Previous_Error <= Convergence_Tolerance)) then
             !--Success - Converged for all bins--!
             print *, 'Errors have converged to tolerance:', Convergence_Tolerance, ' with nBoot = ', nBoot/2
             print *, 'Press Enter to continue:'
             read(*,*)
             Verbose = .true.
             exit
          end if

          if(Convergence_Loop == Convergence_Loop_Maximum) then
             print *, 'Failed to find convergence by nBoot = ', nBoot
             STOP
          END if

          print *, count(SLR_Error-Previous_Error <= Convergence_Tolerance), ' of ', size(SLR_Error), ' bins have converged for nBoot = ', nBoot/2,'. Relooping'
          Previous_Error = SLR_Error
       end if

    end do

    if(Verbose) print *, 'Done.'

  end subroutine Average_Size_in_Luminosity_Redshift_Errors


  !--Correct for Beta--!

  subroutine Correct_for_Beta_Scalar(Conv, beta, Conv_Error)
    real(double), intent(inout)::Conv(:,:)
    real(double), intent(in)::beta
    real(double),intent(inout),optional:: Conv_Error(:,:)

    INTERFACE
        subroutine Correct_for_Beta_Scalar(Conv, beta, Conv_Error)
          use Param_Types
          real(double), intent(inout)::Conv(:,:)
          real(double), intent(in)::beta
          
          real(double),intent(inout),optional:: Conv_Error(:,:)
        END subroutine Correct_for_Beta_Scalar
     END INTERFACE

    if((1.e0_double-2.e0_double*beta)==0.e0_double) STOP 'FATAL ERROR - Correct_for_Beta_Scalar - beta is equal to 0.5, giving infinities'

    Conv = Conv/(1.e0_double-2.e0_double*beta)
    if(present(Conv_Error)) Conv_Error = Conv_Error/(1.e0_double-2.e0_double*beta)

  end subroutine Correct_for_Beta_Scalar

  subroutine Correct_for_Beta_TwoBins(Convergence, beta, Convergence_Error)
    real(double), intent(inout)::Convergence(:,:,:,:)  !First two indeices must be bins
    real(double), intent(in)::beta(:,:)
    real(double),intent(inout),optional:: Convergence_Error(:,:,:,:)

    integer::i,j

    INTERFACE
       subroutine Correct_for_Beta_TwoBins(Convergence, beta, Convergence_Error)
         use Param_Types
         real(double), intent(inout)::Convergence(:,:,:,:)  !First two indeices must be bins
         real(double), intent(in)::beta(:,:)

         real(double),intent(inout),optional:: Convergence_Error(:,:,:,:)
       END subroutine Correct_for_Beta_TwoBins
    END INTERFACE

    do i = 1, size(Convergence,1)
       do j = 1, size(Convergence,2)
          if(present(Convergence_Error)) then
             call Correct_for_Beta_Scalar(Convergence(i,j,:,:), beta(i,j), Convergence_Error(i,j,:,:))
          else
             call Correct_for_Beta_Scalar(Convergence(i,j,:,:), beta(i,j))
          end if
       end do
    end do

  end subroutine Correct_for_Beta_TwoBins

  !---Beta Measurement---!
  subroutine calculate_beta_fromCatalogue(Cat, beta)
    !-Calculates Beta in the usual way, dividing the catalogue into magnitude bins, calculating the average size in each magnitude bin and fitting using the error in average size in each bin--!
    !-Returns a single value of beta, calcaulted for the whole catalogue
    !--This could be easily genrealised to retunr multiple values of beta, going by the size of beta--!

    type(Catalogue),intent(in)::Cat
    real(double), intent(out)::beta

    real(double)::Red_Bin(1,2) !-Always assumes one redshift bin-!

    real(double),allocatable::Mag_Bin(:,:)
    integer::nMag = 5
    real(double),allocatable::SMR(:,:), SMR_Error(:,:)

    print *, 'Calculating Beta from fit across whole catalogue'
    read(*,*)

    print *, 'Cat Abs_Mag:', allocated(Cat%Absolute_Magnitude), size(Cat%Absolute_Magnitude), all(Cat%Absolute_Magnitude==0.e0_double)
    if(allocated(Cat%Absolute_Magnitude)== .false. .or. size(Cat%Absolute_Magnitude) == 0 .or. all(Cat%Absolute_Magnitude==0.e0_double)) then
       print *, 'calculate_beta_fromCatalogue - No magnitude information, returning zero'
       beta = 0.e0_double
       return
    end if

    !-Create Magnitude Bins, roughly as the same number of galxies in each bin. Ideally, nBin shoud iterate towards a value which balances low noise with enough data points-!
    call Calculate_Bin_Limits_by_equalNumber(Cat%Absolute_Magnitude, nMag, Mag_Bin)
    Red_Bin(1,:) = (/ minval(Cat%Redshift), maxval(Cat%Redshift) /)

    print *, 'Got Bin Limits'
    read(*,*)

    !--Calculate beta using all magnitude bins--!
    allocate(SMR(size(Mag_Bin,1), size(Red_Bin,1))); SMR = 0.e0_double
    allocate(SMR_Error(size(Mag_Bin,1), size(Red_Bin,1))); SMR_Error = 0.e0_double

    call Average_Size_in_Luminosity_Redshift_Bins(Cat, Mag_bin, Red_Bin, SMR)
    print *, 'Got SMR'; read(*,*)
    call Average_Size_in_Luminosity_Redshift_Errors(Cat, Mag_Bin, Red_Bin, SMR_Error)
    print *, 'Got SMR Error';read(*,*)
    call calculate_beta_fromAverageSize_AllMagnitudeBins(SMR(:,1), SMR_Error(:,1), Mag_Bin, beta)

    print *, 'Got beta, reading'
    read(*,*)

  end subroutine calculate_beta_fromCatalogue

  subroutine calculate_beta_fromCatalogue_MagnitudeRedshiftBins(Cat, Mag_Bins, Red_Bins, beta)
    !--Wrapper for calculate_beta_fromAverageSize_MagnitudeRedshiftBins, which calculates the average size in the bins provided, and calls the SMR routine
    !--Magnitude bins must be set-!
    type(Catalogue),intent(in)::Cat
    real(double),intent(in),dimension(:,:)::Mag_Bins, Red_Bins
    real(double),intent(out)::beta(:,:)

    real(double),allocatable::SMR(:,:), SMR_Error(:,:)

    print *, 'Calculating beta from Catalogue, via binned Magnitude and Redshift Method'

    if(allocated(Cat%Absolute_Magnitude)== .false. .or. size(Cat%Absolute_Magnitude) == 0 .or. all(Cat%Absolute_Magnitude==0.e0_double)) then
       print *, 'calculate_beta_fromCatalogue - No magnitude information, returning zero'
       beta = 0.e0_double
       return
    end if

    allocate(SMR(size(Mag_Bins,1), size(Red_Bins,1))); SMR = 0.e0_double
    allocate(SMR_Error(size(Mag_Bins,1), size(Red_Bins,1))); SMR_Error = 0.e0_double

    call Average_Size_in_Luminosity_Redshift_Bins(Cat, Mag_bins, Red_Bins, SMR)
    call Average_Size_in_Luminosity_Redshift_Errors(Cat, Mag_Bins, Red_Bins, SMR_Error)
    call calculate_beta_fromAverageSize_MagnitudeRedshiftBins(SMR, SMR_Error, Mag_Bins, beta)

    deallocate(SMR, SMR_Error)

  end subroutine calculate_beta_fromCatalogue_MagnitudeRedshiftBins

  subroutine calculate_beta_fromAverageSize_MagnitudeRedshiftBins(SLR, SLR_Error, Mag_Bins, beta, nFit)
    !--Calulates beta by fitting a straight line to an average size-magnitude-redshift realtion, using mag bins.
    !--Returns a beta for EVERY mag/redshift bin--!
    !--Uses inverse variance weigthing to fit to points. Lowest/Highest magnitude bin fit using only two points--! 
    !--nFit defines the number of points that should be used for each fit, and should be odd-! 
    use nr, only: fit
    real(double), intent(in), dimension(:,:)::SLR, SLR_Error, Mag_Bins
    
    real(double), intent(out), dimension(:,:):: Beta
    integer,intent(in),optional::nFit !-Should be odd-! 

    integer::Method = 1
    real(double)::fit_a, fit_b, fit_siga, fit_sigb, fit_chi2, fit_q
    integer::BinZ, BinL, counter, i, inFit, nFit_EitherSide

    real(double),allocatable::Fit_Mag(:),Fit_Size(:), Fit_Error(:) 
    !--Two Methods 1: Fit using data on eaither side; 2: SubDivide each luminosity bin and fit to subdivisions within bin--!

    if(Verbose) print *, 'Finding Beta'

    inFit = 3
    if(present(nFit)) inFit = nFit

    if(nint(inFit/2.e0_double) /= nint((inFit+1)/2.e0_double)) STOP 'FATAL ERROR - calculate_beta_fromAverageSize_MagnitudeRedshiftBins - inFit is not odd!'

    nFit_EitherSide = nint((nFit - 1)/2.e0_double)

    select case(Method)
    case(1) !-Fit using data on either side-!
       !-Doesn't yet account for grid points with n = 0-!
       do BinZ = 1, size(SLR,2)
          do BinL = 1, size(SLR,1)
             if(BinL ==1) then
                call calculate_beta_fromAverageSize_AllMagnitudeBins(SLR(BinL:BinL+nFit_EitherSide, BinZ), SLR_Error(BinL:BinL+nFit_EitherSide, BinZ), Mag_Bins(BinL:BinL+nFit_EitherSide, :), beta(BinL,BinZ))
             elseif(BinL==size(SLR,1)) then
                call calculate_beta_fromAverageSize_AllMagnitudeBins(SLR(BinL-nFit_EitherSide:BinL, BinZ), SLR_Error(BinL-nFit_EitherSide:BinL, BinZ), Mag_Bins(BinL-nFit_EitherSide:BinL, :), beta(BinL,BinZ))
             else
                call calculate_beta_fromAverageSize_AllMagnitudeBins(SLR(BinL-nFit_EitherSide:BinL+nFit_EitherSide, BinZ), SLR_Error(BinL-nFit_EitherSide:BinL+nFit_EitherSide, BinZ), Mag_Bins(BinL-nFit_EitherSide:BinL+nFit_EitherSide, :), beta(BinL,BinZ))
             end if
          end do
       end do

    case default
       print *, 'calculate_Beta -incorrect method chosen, returning'
       return
    end select

    if(Verbose) print *, 'Found Beta'

  end subroutine calculate_beta_fromAverageSize_MagnitudeRedshiftBins

  subroutine calculate_beta_fromAverageSize_AllMagnitudeBins(SM, SM_Error, Mag_Bins, beta)
    !--THIS CONSTAINS THE MAIN METHOD FOR CALCULATING BETA, GIVEN AVERAGE SIZE ACROSS A SET OF MAGNITUDE BINS--!
    !--Calulates beta by fitting a straight line to an average size-magnitude-redshift realtion, using mag bins.
    !--Returns a singular value of beta by fitting across all magnitude bins
    !--Uses inverse variance weigthing to fit to points. Lowest/Highest magnitude bin fit using only two points--! 
    use nr, only: fit
    real(double), intent(in), dimension(:,:):: Mag_Bins
    real(double), intent(in), dimension(:)::SM, SM_Error
    
    real(double), intent(out):: Beta

    real(double)::fit_a, fit_b, fit_siga, fit_sigb, fit_chi2, fit_q
    integer::counter, i

    real(double),allocatable::Fit_Mag(:),Fit_Size(:), Fit_Error(:) 
    !--Two Methods 1: Fit using data on eaither side; 2: SubDivide each luminosity bin and fit to subdivisions within bin--!

    if(Verbose) print *, 'Finding Beta'

    if(size(SM) /= size(Mag_Bins,1)) STOP 'FATAL ERROR - calculate_beta_fromAverageSize_MagnitudeBins - Size-Magnitude relation is not of the same size as bin information passed in'
    if(size(SM_Error) /= size(SM))  STOP 'FATAL ERROR - calculate_beta_fromAverageSize_MagnitudeBins - Size-Magnitude relation is not of the same size as its Errors' 

       !-Doesn't yet account for grid points with n = 0-!
    allocate(Fit_Size(size(SM))); Fit_Size = 0.e0_double
    allocate(Fit_Mag(Size(SM))); Fit_Mag = 0.e0_double
    allocate(Fit_Error(Size(SM))); Fit_Error = 0.e0_double
             
    counter = 1
    do i = 1, size(SM)
       if(SM(i) == 0.e0_double .or. isNaN(SM_Error(i))) then
          Fit_Size(counter) = 0.e0_double
          Fit_Error(counter) = 1.e30_double
       else
          Fit_Size(counter) = dlog10(SM(i))
          Fit_Error(counter) = 0.434e0_double*SM_Error(i)/SM(i)
       end if
       Fit_Mag(counter) = sum(Mag_Bins(i,:))/2.e0_double
       counter = counter + 1
    end do
     
    call fit(Fit_Mag, Fit_Size, fit_a, fit_b, fit_siga, fit_sigb, fit_chi2, fit_q, fit_Error)
    Beta = -2.5e0_double*fit_b
    
    deallocate(Fit_Size, Fit_Mag, Fit_Error);
 
    print *, 'Found Beta in across all mag bins:', beta

  end subroutine calculate_beta_fromAverageSize_AllMagnitudeBins

  subroutine calculate_beta_InterpolatefromGrid_File(filename, Mag_Bins, Red_Bins, beta_out)
    !--Calculates Beta by reading in the beta_information from a file--!
    !--File format should be 1st row is redshift bin limits, 1st column magnitude bin limits, and 1:,1: is Beta--!
    !--Redshift information, and magnitude information MUST be redshift/ mag LIMITS so that i, j is taken in redshift bin Redshift(i) to Redshift(i+1), and same for magnitude-_!

    character(*),intent(in)::filename
    real(double),intent(out)::beta_out(:,:)
    real(double),intent(in),dimension(:,:)::Mag_Bins, Red_Bins

    STOP 'FATAL ERROR - calculate_beta_InterpolatefromGrid_File - THIS ROUTINE HAS NOT BEEN WRITTEN YET'

    beta_out = 1.e0_8

  end subroutine calculate_beta_InterpolatefromGrid_File

  !--End of Beta Measurement rountines--!

end module Convergence_Estimation


!--DEPRECATED CODE--!


!!$  subroutine calculate_beta_fromAverageSize_MagnitudeRedshiftBins(SLR, SLR_Error, Mag_Bins, beta, nFit)
!!$    !--Calulates beta by fitting a straight line to an average size-magnitude-redshift realtion, using mag bins.
!!$    !--Returns a beta for every mag/redshift bin--!
!!$    !--Uses inverse variance weigthing to fit to points. Lowest/Highest magnitude bin fit using only two points--! 
!!$    use nr, only: fit
!!$    real(double), intent(in), dimension(:,:)::SLR, SLR_Error, Mag_Bins
!!$    
!!$    real(double), intent(out), dimension(:,:):: Beta
!!$    real(double), intent(
!!$
!!$    integer::Method = 1
!!$    real(double)::fit_a, fit_b, fit_siga, fit_sigb, fit_chi2, fit_q
!!$    integer::BinZ, BinL, counter, i
!!$
!!$    real(double),allocatable::Fit_Mag(:),Fit_Size(:), Fit_Error(:) 
!!$    !--Two Methods 1: Fit using data on eaither side; 2: SubDivide each luminosity bin and fit to subdivisions within bin--!
!!$
!!$    if(Verbose) print *, 'Finding Beta'
!!$
!!$    select case(Method)
!!$    case(1) !-Fit using data on either side-!
!!$       !-Doesn't yet account for grid points with n = 0-!
!!$       do BinZ = 1, size(SLR,2)
!!$          do BinL = 1, size(SLR,1)
!!$             if(BinL ==1) then
!!$                allocate(Fit_Size(2)); Fit_Size = 0.e0_double
!!$                allocate(Fit_Mag(2)); Fit_Mag = 0.e0_double
!!$                allocate(Fit_Error(2)); Fit_Error = 0.e0_double
!!$
!!$                counter = 1
!!$                do i = BinL, BinL+1
!!$                   if(SLR(i,BinZ) == 0.e0_double .or. isNaN(SLR_Error(i,BinZ))) then
!!$                      Fit_Size(counter) = 0.e0_double
!!$                      Fit_Error(counter) = 1.e30_double !-Set arbirarily large as no information on this scale-!
!!$                   else
!!$                      Fit_Size(counter) = dlog10(SLR(i, BinZ))
!!$                      Fit_Error(counter) = 0.434e0_double*SLR_Error(i,BinZ)/SLR(i,BinZ)
!!$                   end if
!!$                   Fit_Mag(counter) = sum(Mag_Bins(i,:))/2.e0_double
!!$                   counter = counter + 1
!!$                end do
!!$
!!$
!!$             elseif(BinL==size(SLR,1)) then
!!$                allocate(Fit_Size(2)); Fit_Size = 0.e0_double
!!$                allocate(Fit_Mag(2)); Fit_Mag = 0.e0_double
!!$                allocate(Fit_Error(2)); Fit_Error = 0.e0_double
!!$
!!$                counter = 1
!!$                do i = BinL-1, BinL
!!$                   if(SLR(i,BinZ) == 0.e0_double .or. isNaN(SLR_Error(i,BinZ))) then
!!$                      Fit_Size(counter) = 0.e0_double
!!$                      Fit_Error(counter) = 1.e30_double
!!$                   else
!!$                      Fit_Size(counter) = dlog10(SLR(i, BinZ))
!!$                      Fit_Error(counter) = 0.434e0_double*SLR_Error(i,BinZ)/SLR(i,BinZ)
!!$                   end if
!!$                   Fit_Mag(counter) = sum(Mag_Bins(i,:))/2.e0_double
!!$                   counter = counter + 1
!!$                end do
!!$             
!!$             else
!!$                allocate(Fit_Size(3)); Fit_Size = 0.e0_double
!!$                allocate(Fit_Mag(3)); Fit_Mag = 0.e0_double
!!$                allocate(Fit_Error(3)); Fit_Error = 0.e0_double
!!$
!!$                counter = 1
!!$                do i = BinL-1, BinL+1
!!$                   if(SLR(i,BinZ) == 0.e0_double .or. isNaN(SLR_Error(i,BinZ))) then
!!$                      Fit_Size(counter) = 0.e0_double
!!$                      Fit_Error(counter) = 1.e30_double
!!$                   else
!!$                      Fit_Size(counter) = dlog10(SLR(i, BinZ))
!!$                      Fit_Error(counter) = 0.434e0_double*SLR_Error(i,BinZ)/SLR(i,BinZ)
!!$                   end if
!!$                   Fit_Mag(counter) = sum(Mag_Bins(i,:))/2.e0_double
!!$                   counter = counter + 1
!!$                end do
!!$
!!$             end if
!!$
!!$             call fit(Fit_Mag, Fit_Size, fit_a, fit_b, fit_siga, fit_sigb, fit_chi2, fit_q, fit_Error)
!!$             Beta(BinL, BinZ) = -2.5e0_double*fit_b
!!$
!!$             deallocate(Fit_Size, Fit_Mag, Fit_Error);
!!$          end do
!!$       end do
!!$
!!$    case default
!!$       print *, 'calculate_Beta -incorrect method chosen, returning'
!!$       return
!!$    end select
!!$
!!$    if(Verbose) print *, 'Found Beta'
!!$
!!$  end subroutine calculate_beta_fromAverageSize_MagnitudeRedshiftBins
