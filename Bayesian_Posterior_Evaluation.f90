module Bayesian_Posterior_Evaluation
  !--Contains the methods to evaluate the Bayesian Posterior for an input of free parameters
  use Param_Types

  implicit none

  !-- use lookup renorm does not work, but useful to test run-time dependence of posterior evaluation
  logical, private:: use_lookup = .false., use_lookup_Renorm = .false.
  

  logical, private:: debug = .false.
  character(500), private:: debug_Dir = 'debug/'

  !--Function Overload-------------------------------------
  INTERFACE Likelihood_Evaluation_atVirialRadius_perSource
     module procedure Likelihood_atVirialRadius_SingleCluster_perSource, Likelihood_atVirialRadius_MultipleCluster_perSource
  END INTERFACE Likelihood_Evaluation_atVirialRadius_perSource

  INTERFACE lnLikelihood_Evaluation_atVirialRadius_perSourceSample
     module procedure lnLikelihood_atVirialRadius_SingleCluster_perSourceSample, lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample
  END INTERFACE lnLikelihood_Evaluation_atVirialRadius_perSourceSample

  interface Combine_Posteriors
     module procedure Combine_Posteriors_Scalar, Combine_Posteriors_Vector
  end interface Combine_Posteriors

  interface delens_Source_Size
     module procedure  delens_Source_Size_Array, delens_Source_Size_Scalar
  end interface delens_Source_Size

  interface delens_Source_Magnitude
     module procedure delens_Source_Magnitude_Scalar, delens_Source_Magnitude_Array
  end interface delens_Source_Magnitude

  type RenormLookup
     real(double),allocatable:: muGrid(:), S1_LT(:), S2_LT(:)
  end type RenormLookup
  type(RenormLookup), private:: LT_Renormalisation




contains  

  !---------------------LENSING RELATIONS

  !~~~~~~~~~~~~~Magnitude Relations
  real(double) function delens_Source_Magnitude_Scalar(mu, m)
    !--Returns the DELENSED source magnitude.
    !--mu is magnification factor
    !--m is apparent by default

    real(double), intent(in)::mu, m

    real(double):: dlm(1)
    
    dlm = delens_Source_Magnitude_Array(mu, (/m/))

    delens_Source_Magnitude_Scalar = dlm(1)

  end function delens_Source_Magnitude_Scalar


  function delens_Source_Magnitude_Array(mu, m)
    !--Returns the DELENSED source magnitude.
    !--mu is magnification factor
    !--m is apparent by default

    real(double), intent(in)::mu, m(:)

    real(double), dimension(size(m)):: delens_Source_Magnitude_Array

    delens_Source_Magnitude_Array = m + 2.5e0_double*dlog10(mu)

  end function delens_Source_Magnitude_Array

  !~~~~~~~~~~~~~Size Relations
  function delens_Source_Size_Scalar(mu, T, lnT) RESULT(Res)
    !--Returns the DELENSED source size.
    !--mu is magnification factor
    !--T is size by default, or log-Size if lnT is present and true
    real(double), intent(in):: mu
    real(double), intent(in):: T
    logical, intent(in),optional:: lnT

    real(double):: dlT(1)

    real(double):: Res

    INTERFACE
       function delens_Source_Size_Scalar(mu, T, lnT) RESULT(Res)
         use Param_Types
         real(double), intent(in):: mu
         real(double), intent(in):: T(:)
         logical, intent(in),optional:: lnT

         real(double):: Res
       end function delens_Source_Size_Scalar
    END INTERFACE


    if(present(lnT)) then
       dlT = delens_Source_Size_Array(mu, (/T/), lnT)
    else
       dlT = delens_Source_Size_Array(mu, (/T/))
    end if

    Res = dlT(1)

  end function delens_Source_Size_Scalar

  function delens_Source_Size_Array(mu, T, lnT) RESULT(Res)
    !--Returns the DELENSED source size.
    !--mu is magnification factor
    !--T is size by default, or log-Size if lnT is present and true - uses NATURAL LOGARITHM - must be consistent with defintition og log-Size in catalogue
    real(double), intent(in):: mu
    real(double), intent(in):: T(:)
    logical, intent(in),optional:: lnT

    real(double), dimension(size(T)):: Res


    if(present(lnT) .and. (lnT == .true.)) then
       Res = T - 0.5e0_double*dlog(mu)
    else
       Res = T/dsqrt(mu)
    end if

  end function delens_Source_Size_Array


  !---------------------LIKELIHOOD EVALUATION

  real(double) function lnLikelihood_atVirialRadius_SingleCluster_perSourceSample(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit, use_lnT, Flag, RenormalisationGroup, output_Prefix)
    !--Sigma_Crit must be entered in units (10^18 MSun/h), and the first dimension should hold sigma critical for all lenses considered (not on a lens redshift grid)
    use Cosmology, only: angular_diameter_distance_fromRedshift
    use Mass_Profiles; use Distributions, only: CH08_Redshift_Distribution_Scalar
    use Interpolaters, only: Linear_Interp; use Integration
    real(double),intent(in):: Alpha
    !--Mass Model Declarations
    integer:: Profile
    real(double), intent(in)::  Lens_Redshift, Cluster_Position(:) !- RA/Dec-!
    !--Galaxy Data Input
    real(double), intent(in):: Theta(:), Magnitude(:), Source_Redshift(:), Position(:,:) !-RA,Dec-!
    !--Supplementary
    integer:: Method(:)
    integer, intent(in):: RenormalisationGroup(:)
    real(double), intent(in):: Pr_MagGrid(:), Pr_SizeGrid(:), SM_Prior(:,:), M_Prior(:)
    real(double), intent(in):: Size_Limits(:,:), Magnitude_Limits(:,:)
    real(double), intent(in):: MagnificationGrid(:), M_Pr_Renormalisation(:,:), SM_Pr_Renormalisation(:,:)
    real(double), intent(in):: RedshiftGrid(:), Sigma_Crit(:) !-Source_Redshift-!
    logical:: use_lnT
    character(3):: Flag(:)
    !--Optional
    character(*), optional:: output_Prefix


    real(double),allocatable:: Posterior_perSource(:)

    !--Internal_Declarations--!
    real(double):: tCluster_Position(1,2), tSigma_Crit(1,size(Sigma_Crit))

    tCluster_Position(1,:) = Cluster_Position; tSigma_Crit(1,:) = Sigma_Crit

    if(present(output_Prefix)) then
       lnLikelihood_atVirialRadius_SingleCluster_perSourceSample = lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample((/Alpha/), Method, Profile, tCluster_Position, (/Lens_Redshift/), Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, tSigma_Crit, use_lnT, Flag,RenormalisationGroup, output_Prefix)
    else
       lnLikelihood_atVirialRadius_SingleCluster_perSourceSample = lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample((/Alpha/), Method, Profile, tCluster_Position, (/Lens_Redshift/), Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, tSigma_Crit, use_lnT, Flag, RenormalisationGroup)
    end if

  end function lnLikelihood_atVirialRadius_SingleCluster_perSourceSample

  
  real(double) function lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit, use_lnT, Flag, RenormalisationGroup, output_Prefix)
    !-- Returns the log-Likelihood for a source sample (>1 source) at a single r200 and centroid position, where Alpha (r200) and Cluster_Position (centroid) can contain information for multiple clusters
    !--Sigma_Crit must be entered in units (10^18 MSun/h), and the first dimension should hold sigma critical for all lenses considered (not on a lens redshift grid)
    use Cosmology, only: angular_diameter_distance_fromRedshift
    use Mass_Profiles; use Distributions, only: CH08_Redshift_Distribution_Scalar
    use Interpolaters, only: Linear_Interp; use Integration
    real(double),intent(in):: Alpha(:)
    !--Mass Model Declarations
    integer, intent(in):: Profile
    real(double), intent(in)::  Lens_Redshift(:), Cluster_Position(:,:) !-Cluster, RA/Dec-!
    !--Galaxy Data Input
    real(double), intent(in):: Theta(:), Magnitude(:), Source_Redshift(:), Position(:,:) !-RA,Dec-!
    !--Supplementary
    integer, intent(in):: Method(:), RenormalisationGroup(:)
    real(double), intent(in):: Pr_MagGrid(:), Pr_SizeGrid(:), SM_Prior(:,:), M_Prior(:)
    real(double), intent(in):: Size_Limits(:,:), Magnitude_Limits(:,:)
    real(double), intent(in):: MagnificationGrid(:), M_Pr_Renormalisation(:,:), SM_Pr_Renormalisation(:,:)
    real(double), intent(in):: RedshiftGrid(:), Sigma_Crit(:,:) !-Lens, Source_Redshift-!
    logical:: use_lnT
    character(3):: Flag(:)
    !--Optional
    character(*), optional:: output_Prefix

    real(double),allocatable:: Posterior_perSource(:)

    integer:: nGal_Ignored_MagLimits, nGal_Ignored_NaN

    integer::c

    real(double):: Res

    real(double), allocatable:: Posterior_M(:), Posterior_SM(:), Posterior_S(:), data_S(:,:), data_SM(:,:), data_M(:,:)
    real(double), dimension(3):: combined_Posterior_byMethod
    integer, dimension(3):: methodCount
    logical:: split_byMethod = .true. !--Can be turned off to speed up run-time, at expense of removing possible output of posterior per source split by method

    character(10):: fmt
    character(5):: output_Type_Label
    logical:: here

    allocate(Posterior_perSource(size(Theta))); Posterior_perSource = dsqrt(-1.e0_double)

    if(split_byMethod) then
       allocate(Posterior_S(count(Method == 1))); Posterior_S = 0.e0_double/0.e0_double
       allocate(Posterior_SM(count(Method == 2))); Posterior_SM = 0.e0_double/0.e0_double
       allocate(Posterior_M(count(Method == 3))); Posterior_M = 0.e0_double/0.e0_double !--NaNs as default

       !--Store information on sample as split, used on output
       allocate(data_S(5,size(Posterior_S))); data_S = 0.e0_double
       allocate(data_SM(5,size(Posterior_SM))); data_SM = 0.e0_double
       allocate(data_M(5,size(Posterior_M))); data_M = 0.e0_double
    end if

    !---Set up debugging output
    if(debug) then

       if(present(output_Prefix)) then
          debug_Dir = trim(adjustl(output_Prefix))//'debug/'
       else
          debug_Dir = 'debug/'
       end if

       inquire(directory = debug_Dir, exist = here)
       if(here == .false.) then
          print *, '--- Debugging Directory:', trim(adjustl(debug_Dir)), ' does not exist, creating....'
          call system('mkdir '//trim(adjustl(debug_Dir)))
       end if
    end if

    methodCount = 0
    do c = 1, size(Theta)
       Posterior_perSource(c) = 1.e0_double

!!$       !--Skip Evaluation if Galaxy Input falls outside limits -- Removed as should be done at posterior evaluation level: e.g. magnitude-only should not care about whether size is a NaN
!!$       if( (Magnitude(c) > Magnitude_Limits(2)) .or. (Magnitude(c) < Magnitude_Limits(1))) then
!!$          nGal_Ignored_MagLimits = nGal_Ignored_MagLimits + 1
!!$          cycle
!!$       end if
!!$       if(isNaN(Theta(c)) .or. isNaN(Magnitude(c))) then
!!$          nGal_Ignored_NaN = nGal_Ignored_NaN + 1
!!$          cycle
!!$       end if
       
       Posterior_perSource(c) = Likelihood_atVirialRadius_MultipleCluster_perSource(Alpha, Method(c), Profile, Cluster_Position, Lens_Redshift, Theta(c), Magnitude(c), Source_Redshift(c), Position(c,:), Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits(RenormalisationGroup(c),:), Magnitude_Limits(RenormalisationGroup(c),:), MagnificationGrid, M_Pr_Renormalisation(RenormalisationGroup(c),:), SM_Pr_Renormalisation(RenormalisationGroup(c),:), RedshiftGrid, Sigma_Crit, use_lnT, Flag(c))
       
       if(split_byMethod) then
          select case(Method(c))
          case(1)
             methodCount(1) = methodCount(1) + 1
             data_S(:,methodCount(1)) = (/Theta(c), Magnitude(c), Source_Redshift(c), Position(c,:)/)
             Posterior_S(methodCount(1)) = Posterior_perSource(c)
          case(2)
             methodCount(2) = methodCount(2) + 1
             data_SM(:,methodCount(2)) = (/Theta(c), Magnitude(c), Source_Redshift(c), Position(c,:)/)
             Posterior_SM(methodCount(2)) = Posterior_perSource(c)
          case(3)
             methodCount(3) = methodCount(3) + 1
             data_M(:,methodCount(3)) = (/Theta(c), Magnitude(c), Source_Redshift(c), Position(c,:)/)
             Posterior_M(methodCount(3)) = Posterior_perSource(c)
          end select
       end if
    end do

    !--Error Catching
    if(methodCount(1) /= size(Posterior_S)) STOP 'lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample - Error in split counts for method 1'
    if(methodCount(2) /= size(Posterior_SM)) STOP 'lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample - Error in split counts for method 2'
    if(methodCount(3) /= size(Posterior_M)) STOP 'lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample - Error in split counts for method 3'

    !--Recombine
    call Combine_Posteriors_Scalar(Alpha(1), Posterior_perSource(:), .true., Renormalise = .false., Return_lnP = .true., Combined_Posterior = Res)    

    combined_Posterior_byMethod = 0.e0_double
    if(split_byMethod) then
       if(methodCount(1) /= 0) call Combine_Posteriors_Scalar(Alpha(1), Posterior_S(:), .true., Renormalise = .false., Return_lnP = .true., Combined_Posterior = combined_Posterior_byMethod(1))
       if(methodCount(2) /= 0) call Combine_Posteriors_Scalar(Alpha(1), Posterior_SM(:), .true., Renormalise = .false., Return_lnP = .true., Combined_Posterior = combined_Posterior_byMethod(2))
       if(methodCount(3) /= 0) call Combine_Posteriors_Scalar(Alpha(1), Posterior_M(:), .true., Renormalise = .false., Return_lnP = .true., Combined_Posterior = combined_Posterior_byMethod(3))
    end if

    if(present(output_Prefix)) then

       !--Output Posterior per source. Each row corresponds to the single alpha entered here. Therefore, this only works for a single alpha value, and may not necessarily be ordered in alpha 

       !--This may be slow as opend and closed for each alpha
       !--General Combined output
       inquire(file = trim(adjustl(output_Prefix))//'_Posterior_per_Source.dat', exist = here)
       if(here == .false.) then
          !Write header consisting of all sources in sample
          write(fmt, '(I10)') size(Theta)+2
          open(unit = 24, file = trim(adjustl(output_Prefix))//'_Posterior_per_Source.dat', status = 'new')
          write(24, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Sizes:', 0., 0., Theta
          write(24, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Mags:', 0., 0., Magnitude
          close(24)
       end if

       open(unit = 24, file = trim(adjustl(output_Prefix))//'_Posterior_per_Source.dat', Access = 'append')
       write(fmt, '(I10)') size(Alpha)+size(Posterior_perSource)+1
       write(24, '('//trim(fmt)//'(e14.7,x))') Alpha, Res, Posterior_perSource
       close(24)

       if(split_byMethod) then
          if(methodCount(1) /= 0) then

             output_Type_Label = 'S'
             inquire(file = trim(adjustl(output_Prefix))//'_Posterior_per_Source_'//trim(adjustl(output_Type_Label))//'.dat', exist = here)
             if(here == .false.) then
                !Write header consisting of all sources in sample
                write(fmt, '(I10)') size(Alpha)+size(data_S,2)+1
                open(unit = 26, file = trim(adjustl(output_Prefix))//'_Posterior_per_Source_'//trim(adjustl(output_Type_Label))//'.dat', Access = 'append')
                write(26, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Sizes:', 0., 0., data_S(1,:)
                write(26, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Mags:', 0., 0., data_S(2,:)
                write(26, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Z:', 0., 0., data_S(3,:)
                write(26, '(A,x,'//trim(fmt)//'(e14.7,x))') '#RA:', 0., 0., data_S(4,:)
                write(26, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Dec:', 0., 0., data_S(5,:)
                close(26)
             end if


             write(fmt, '(I10)') size(Alpha)+size(Posterior_S)+1
             write(26, '('//trim(fmt)//'(e14.7,x))') Alpha, Res, Posterior_S
             close(26)
          end if
          
          if(methodCount(2) /= 0) then
             output_Type_Label = 'SM'
             inquire(file = trim(adjustl(output_Prefix))//'_Posterior_per_Source_'//trim(adjustl(output_Type_Label))//'.dat', exist = here)
             if(here == .false.) then
                !Write header consisting of all sources in sample
                write(fmt, '(I10)') size(Alpha)+size(data_SM,2)+1
                open(unit = 27, file = trim(adjustl(output_Prefix))//'_Posterior_per_Source_'//trim(adjustl(output_Type_Label))//'.dat', Access = 'append')
                write(27, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Sizes:', 0., 0., data_SM(1,:)
                write(27, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Mags:', 0., 0., data_SM(2,:)
                write(27, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Z:', 0., 0., data_SM(3,:)
                write(27, '(A,x,'//trim(fmt)//'(e14.7,x))') '#RA:', 0., 0., data_SM(4,:)
                write(27, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Dec:', 0., 0., data_SM(5,:)
                close(27)
             end if

             open(unit = 27, file = trim(adjustl(output_Prefix))//'_Posterior_per_Source_'//trim(adjustl(output_Type_Label))//'.dat', Access = 'append')
             write(fmt, '(I10)') size(Alpha)+size(Posterior_SM)+1
             write(27, '('//trim(fmt)//'(e14.7,x))') Alpha, Res, Posterior_SM
             close(27)
          end if
          
          if(methodCount(3) /= 0) then
             output_Type_Label = 'M'
             inquire(file = trim(adjustl(output_Prefix))//'_Posterior_per_Source_'//trim(adjustl(output_Type_Label))//'.dat', exist = here)
             if(here == .false.) then
                !Write header consisting of all sources in sample
                write(fmt, '(I10)') size(Alpha)+size(data_M,2)+1
                open(unit = 28, file = trim(adjustl(output_Prefix))//'_Posterior_per_Source_'//trim(adjustl(output_Type_Label))//'.dat', Access = 'append')
                write(28, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Sizes:', 0., 0., data_M(1,:)
                write(28, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Mags:', 0., 0., data_M(2,:)
                write(28, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Z:', 0., 0., data_M(3,:)
                write(28, '(A,x,'//trim(fmt)//'(e14.7,x))') '#RA:', 0., 0., data_M(4,:)
                write(28, '(A,x,'//trim(fmt)//'(e14.7,x))') '#Dec:', 0., 0., data_M(5,:)
                close(28)
             end if

             open(unit = 28, file = trim(adjustl(output_Prefix))//'_Posterior_per_Source_'//trim(adjustl(output_Type_Label))//'.dat', Access = 'append')
             write(fmt, '(I10)') size(Alpha)+size(Posterior_M)+1
             write(28, '('//trim(fmt)//'(e14.7,x))') Alpha, Res, Posterior_M
             close(28)
          end if

          deallocate(Posterior_M, Posterior_SM, Posterior_S, data_S, data_M, data_SM)

          !--Output Combined
          open(unit = 29, file = trim(adjustl(output_Prefix))//'_CobminedPosterior_byMethod.dat', Access = 'append')
          write(fmt, '(I10)') size(Alpha)+size(combined_Posterior_byMethod)
          write(29, '('//trim(fmt)//'(e14.7,x))') Alpha, combined_Posterior_byMethod
          close(29)
       end if

    end if

    lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample = Res

  end function lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample

  !--Evalution of the posterior on a single virial radius grid point, per source
  real(double) function Likelihood_atVirialRadius_SingleCluster_perSource(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit, use_lnT, Flag)
    !~~~~~DESCRIPTION: Wrapper routine which returns the posterior evulated at a single point, //where only a single cluster is modelled//. Uses a call to the standard multi-cluster routine to get it's result.
    !--Sigma_Crit must be entered in units (10^18 MSun/h), and the first dimension should hold sigma critical for all lenses considered (not on a lens redshift grid)
    real(double),intent(in):: Alpha
    !--Mass Model Declarations
    integer:: Profile
    real(double), intent(in)::  Lens_Redshift, Cluster_Position(:) !-Cluster, RA/Dec-!
    !--Galaxy Data Input
    real(double), intent(in):: Theta, Magnitude, Source_Redshift, Position(2) !-RA,Dec-!
    !--Supplementary
    integer:: Method
    real(double), intent(in):: Pr_MagGrid(:), Pr_SizeGrid(:), SM_Prior(:,:), M_Prior(:)
    real(double), intent(in):: Size_Limits(:), Magnitude_Limits(:)
    real(double), intent(in):: MagnificationGrid(:), M_Pr_Renormalisation(:), SM_Pr_Renormalisation(:)
    real(double), intent(in):: RedshiftGrid(:), Sigma_Crit(:)
    logical:: use_lnT
    character(3):: Flag

    !--Internal_Declarations--!
    real(double):: tCluster_Position(1,2), tSigma_Crit(1,size(Sigma_Crit))

    tCluster_Position(1,:) = Cluster_Position; tSigma_Crit(1,:) = Sigma_Crit

    Likelihood_atVirialRadius_SingleCluster_perSource = Likelihood_atVirialRadius_MultipleCluster_perSource((/Alpha/), Method, Profile, tCluster_Position, (/Lens_Redshift/), Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, tSigma_Crit, use_lnT, Flag)

  end function Likelihood_atVirialRadius_SingleCluster_perSource

  real(double) function Likelihood_atVirialRadius_MultipleCluster_perSource(Alpha, Method, Profile, Cluster_Position, Lens_Redshift, Theta, Magnitude, Source_Redshift, Position, Pr_MagGrid, Pr_SizeGrid, SM_Prior, M_Prior, Size_Limits, Magnitude_Limits, MagnificationGrid, M_Pr_Renormalisation, SM_Pr_Renormalisation, RedshiftGrid, Sigma_Crit, use_lnT, Posterior_Flag)
    !--Produces likelihood for a single source
    !-This is the workhorse of the routine as it stands. This constructs the likelihood of the data (theta and magnitude) given a set of cluster mass profile free parameters (alpha). The input priors MUST be normalised to unity in the size and magnitude range of the survey/source data set. Nearly everything reuiwred here is discernable fromt eh catalogue itself, or can be obtained using a call to get_Likelihood_Evaluation_Precursors.
    !--SM_Prior is understood to be evaluated over the prior size and magnitude limits. M_Prior is evaluated as the SM prior integrated over the entered *PRIOR* size limits. NOTE that survey limits are not taken into consideration in this method
    !--Sigma_Crit must be entered in units (10^18 MSun/h), and the first dimension should hold sigma critical for all lenses considered (not on a lens redshift grid)
    !--If use_lnSize is present and true, then prior is assumed to be ln size prior, and Theta assumed to be ln-Theta 
    use Cosmology, only: angular_diameter_distance_fromRedshift
    use Mass_Profiles; use Distributions, only: CH08_Redshift_Distribution_Scalar, redshift_Distribution_byLookup
    use Interpolaters, only: Linear_Interp; use Integration
    real(double),intent(in):: Alpha(:)
    !--Mass Model Declarations
    integer:: Profile
    real(double), intent(in)::  Lens_Redshift(:), Cluster_Position(:,:) !-Cluster, RA/Dec-!
    !--Galaxy Data Input
    real(double), intent(in):: Theta, Magnitude, Source_Redshift, Position(2) !-RA,Dec-!
    !--Supplementary
    integer:: Method
    real(double), intent(in):: Pr_MagGrid(:), Pr_SizeGrid(:), SM_Prior(:,:), M_Prior(:)
    real(double), intent(in):: Size_Limits(:), Magnitude_Limits(:) !-First dimension allows byCluster!!
    real(double), intent(in):: MagnificationGrid(:), M_Pr_Renormalisation(:), SM_Pr_Renormalisation(:)
    real(double), intent(in):: RedshiftGrid(:), Sigma_Crit(:,:) !-Lens, Source_Redshift-!
    logical:: use_lnT
    character(3), intent(out):: Posterior_Flag

    
    !--Internal Declarations
    real(double), allocatable:: Posterior_perGalaxy_Redshift(:)
    real(double):: Effective_Magnification
    logical:: Known_Redshift
    integer:: z, m, j, l
    real(double):: Galaxy_Sigma_Critical(size(Sigma_Crit,1))
    real(double):: RedshiftPDF, Renorm
    real(double), dimension(size(Alpha)):: Distance_From_Mass_Center
    real(double):: D_l
    real(double):: Theta_0, Magnitude_0 !--Delensed Versions
    integer:: m0Index, T0Index

    real(double), dimension(2):: Renormalisation_Size_Limits, Renormalisation_Magnitude_Limits
    real(double), dimension(2):: iSize_Limits, iMag_Limits
    real(double), dimension(size(SM_Prior,1), size(SM_Prior,2)):: Kappa_Renormalised_Prior, Size_Only_Mag_Prior
    real(double), dimension(size(M_Prior)):: Kappa_Renormalised_MagPrior
    real(double),allocatable:: Size_Only_Prior(:), Mag_Only_Prior(:) !-Used in mag-only and size-only analyses

    integer:: Galaxy_Posterior_Method !-Converts from inoput posterior method to the individual method, based on the Method input
    integer:: nMagPosterior, nSizePosterior, nSizeMagPosterior
    integer:: nGal_Ignored_SizeLimits
    integer:: Unlensed_Mag_Index !--Used to store position of unlensed measurement on mag grid (used in Mag-Only)


    integer::nZ

    !--Internal Decalrations that could eventually be set/passed in
    logical:: Enforce_Weak_Lensing = .false., Cuts_Renormalise_Likelihood = .true.

    !---Testing Declarations
    real(double):: Time00, Time01, Time0, Time1, Time2, Time3, Time4, Time5, Time6, Time7
    integer:: tt, ttu = 1

    Posterior_Flag = '0'

    !--Error Catching on Input---!
    if(size(Alpha) /= size(Distance_From_Mass_Center) .or. (size(Alpha) /= size(Lens_Redshift))) then
       print *,  'Likelihood_Evaluation_atVirialRadius - Error on input of Mass Profile variables:'
       print *, 'Sizes (Alpha, Distance, Lens_Redshift):', size(Alpha), size(Distance_From_Mass_Center), size(lens_Redshift)
       STOP
    END if
    
    !--Redshift Knowledge--!
    if(Source_Redshift > 0.e0_double) then
       Known_Redshift = .true.
       nZ = 1
    else
       Known_Redshift = .false.
       nZ = size(RedshiftGrid)
    end if
    
    !--Get Distance of source from each Cluster
    do l = 1, size(Alpha)
       D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift(l))
       Distance_From_Mass_Center(l) = dsqrt( (dabs(Position(1)-Cluster_Position(l, 1))**2.e0_double) + (dabs(Position(2)-Cluster_Position(l, 2))**2.e0_double) ) !-in Degrees-!
       Distance_From_Mass_Center(l) = D_l*Distance_From_Mass_Center(l)*(3.142e0_double/180.e0_double) !-In Mpc/h-!
    end do

    Likelihood_atVirialRadius_MultipleCluster_perSource = 1.e-100_double

    !--Set which method will be used, including effective survey limits
    select case(Method)
    case(0)
       !--Do Nothing, return a constant
       write(*,'(A)') 'WARNING - Galaxy posterior method is zero - this should have been removed!!'
       Likelihood_atVirialRadius_MultipleCluster_perSource = 1.0e0_double
    case(1) !--SizeOnly--!                                                                                                                                                                              
       !iSize_Limits = Size_Limits(1,:); iMag_Limits = (/minval(Pr_MagGrid), maxval(Pr_MagGrid)/) !--Mag limits should encompass the whole data set
       nSizePosterior = nSizePosterior + 1
       Galaxy_Posterior_Method = 1
    case(2)!--SizeMag--!                                                                                                                                                                                    
       !iSize_Limits = Size_Limits(1,:); iMag_Limits = Magnitude_Limits(1,:)
       nSizeMagPosterior = nSizeMagPosterior + 1
       Galaxy_Posterior_Method = 2
    case(3) !-Magnitude Only--!                                                                                                                                                          
       !iMag_Limits = Magnitude_Limits(1,:); iSize_Limits = (/minval(Pr_SizeGrid), maxval(Pr_SizeGrid)/) !--Should encompase the whole data set                    
       nMagPosterior = nMagPosterior + 1
       
       Galaxy_Posterior_Method = 3
    case(4) !-Size mag above size data limit, Mag-Only below size data limit-!                                                                                                                              
       STOP 'Galaxy Posterior Method 4 has been disabled at the level where the posterior is evaluated - Instead, the sample should be split before call, and the correct posterior method and prior fed in for each galaxy. Fatal.'

       !!$STOP 'Posterior Method 4 has been disabled as it needs further thought into its construction. In particular, with respect to the construction of the prior from the size-mag distribution, the effect of prior cuts, and correct renormalisation. The skeleton code has been added for this, but disabled for now. See commented code labeled MagRenorm'
       !--Note, in this case the data-renormalisation should take inot account that the mag only is constructed from a sample which takes small, faint galaxies                                             
       !if(Theta < Size_Limits(1)) then
          !--Use Magnitude information only--!
          !iSize_Limits = (/0.e0_double, Size_Limits(2)/); iMag_Limits = Magnitude_Limits(1,:) !!-I believe this is wrong, but I am testing it anyway
       !   Galaxy_Posterior_Method = 3
       !   nMagPosterior = nMagPosterior + 1
       !else
          !--Use full magnitude-size information--!
          !iSize_Limits = Size_Limits(1,:); iMag_Limits = Magnitude_Limits(1,:)
       !   Galaxy_Posterior_Method = 2
       !   nSizeMagPosterior = nSizeMagPosterior + 1
       !end if
    case default
       STOP 'Likelihood_atVirialRadius_MultipleCluster_perSource - Entered Posterior Method invalid, FATAL'
    end select

    !---Error Catching: Return Likelihood == 0 for all alpha if size or magnitude falls outside set analysis limits - This should not occur as this body should be removed from the sample
    !if( (Magnitude < iMag_Limits(1)) .or. (Magnitude > iMag_Limits(2)) .or. isNaN(Magnitude)) then
    !   Posterior_Flag = '1.1'
    !   Likelihood_atVirialRadius_MultipleCluster_perSource = 1.e-100_double
    !   print *, 'Source Magnitude, Limits:', Magnitude, iMag_Limits
    !   STOP 'Magnitude of entered source falls outside limits or NaN, this should not happen - Check source selection'
    !   return
    !end if
    if((Galaxy_Posterior_Method == 1) .or. (Galaxy_Posterior_Method == 2)) then
       !if( (Theta < iSize_Limits(1)) .or. (Theta > iSize_Limits(2)) .or. isNaN(Theta)) then
       !   Posterior_Flag = '1.2'
       !   Likelihood_atVirialRadius_MultipleCluster_perSource = 1.e-100_double
       !   STOP 'Size of entered source falls outside limits or NaN, this should not happen - Check source selection'
       !   return
       !end if
    end if

    !---Error Catching: Check for cases where the source lies outside evaluated prior limits - This should not occur as this body should be removed from the sample
    if( (Magnitude < Pr_MagGrid(1)) .or. (Magnitude > Pr_MagGrid(size(Pr_MagGrid)-1)) .or. isNaN(Magnitude)) then
       Posterior_Flag = '1.1'
       Likelihood_atVirialRadius_MultipleCluster_perSource = 1.e-100_double
       print *, 'Source Magnitude, Prior Grid Limits:', Magnitude, Pr_MagGrid(1), Pr_MagGrid(size(Pr_MagGrid)-1)
       STOP 'Magnitude of entered source falls outside prior limits or NaN, this should not happen - Check source selection'
       return
    end if
    if((Galaxy_Posterior_Method == 1) .or. (Galaxy_Posterior_Method == 2)) then
       if( (Theta < Pr_SizeGrid(1)) .or. (Theta > Pr_SizeGrid(size(Pr_SizeGrid)-1)) .or. isNaN(Theta)) then
          Posterior_Flag = '1.2'
          Likelihood_atVirialRadius_MultipleCluster_perSource = 1.e-100_double
          print *, 'Source Size, Prior Grid Limits:', Theta, Pr_SizeGrid(1), Pr_SizeGrid(size(Pr_SizeGrid)-1)
          STOP 'Size of entered source falls outside limits or NaN, this should not happen - Check source selection'
          return
       end if
    end if
 

    allocate(Posterior_perGalaxy_Redshift(nZ)); Posterior_perGalaxy_Redshift = 0.e0_double

    !---Set up debugging output
    if(debug) then
       open(56, file = trim(adjustl(debug_Dir))//'Magnification_Factor_List.dat', access = 'APPEND')
    end if
    
    do z = 1, nZ !-Vary source redshift
       !--Set Sigma Critical for that galaxy, across all lenses


       if(Known_Redshift) then
          Do l = 1, size(Sigma_Crit,1) !--Get Sigma_Critical for all lenses by interpolation
             Galaxy_Sigma_Critical(l) = Linear_Interp(Source_Redshift, RedshiftGrid, Sigma_Crit(l,:))
          end Do
       else
          !--As RedshiftPDF and Sigma_Crit are on the same grid, take value directly
          Galaxy_Sigma_Critical = Sigma_Crit(:,z)
       END if

       if(any(isNaN(Galaxy_Sigma_Critical)) .or. any((Galaxy_Sigma_Critical > huge(1.e0_double)))) STOP 'DM_Profile_Variable_Posterior - FATAL - Galaxy Sigma Critical not set correctly'
       
       !~~~~Ths section will need editing to account for multiple clusters in flat and SIS case
       select case(Profile)
       case(1) !-Flat-!                                                                                                                                                                                      
          STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
          Effective_Magnification = 1.e0_double+2.e0_double*(Alpha(1)/Galaxy_Sigma_Critical(1))
       case(2) !-SIS-!                                                                                                                                                                                       
          STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
          Effective_Magnification = 1.e0_double+2.e0_double*(SMD_SIS(Alpha(1), Distance_From_Mass_Center(1))/(Galaxy_Sigma_Critical(1)*1.e18_double))
       case(3) !-NFW-!                                                                                                                                                                                       
          if(Enforce_Weak_Lensing) then
             STOP 'Weak Lensing disabled for NFW for now, until the equivalent for multiple clusters is coded [dont worry, its actually easy]'
             Effective_Magnification = 1.e0_double + 2.e0_double*(SMD_NFW(Distance_From_Mass_Center(1), Lens_Redshift(1), Alpha(1))/(Galaxy_Sigma_Critical(1)*1.e18_double))
          else
             Effective_Magnification = Total_MagnificationFactor_MultipleClusters(3, Position, Cluster_Position, Alpha, Lens_Redshift, Sigma_Crit = Galaxy_Sigma_Critical*1.e18_double)
          end if
       case default
          STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
       end select

       !--Take absolute value of mu, to account for sources whihc lie within Einstein radius (this shouldn't really happen all that often, but may be an issue where the centorid is allowed to vary, or the cluster mass is large)
       Effective_Magnification = dabs(Effective_Magnification)

       if(debug) then
          write(56, '(e14.7,x)') Effective_Magnification
       end if

       if(isNaN(Effective_Magnification)) then
          print *, 'Magnification Factor is a NaN:', Distance_From_Mass_Center, Alpha, Lens_Redshift, Galaxy_Sigma_Critical*1.e18_double
          STOP
       end if
       
       if(Effective_Magnification < minval(MagnificationGrid) .or. Effective_Magnification > maxval(MagnificationGrid)) then
          !---Use this for debugging only
!!$          print *, 'Magnification outwith limits detected for galaxy:', Theta, Magnitude, Effective_Magnification
!!$          print *, 'And for lensing parameters:'
!!$          print *, Position, Cluster_Position, Alpha, Lens_Redshift, Sigma_Crit*1.e18_double
!!$          read(*,*)

          !--Skipping as outside limits on which magnification was evaluated--!                                                                                                                              
          Posterior_Flag = '02' 
          Posterior_perGalaxy_Redshift(z) = 1.e-100_double
          cycle
       end if
       
       if(delens_Source_Magnitude(Effective_Magnification, Magnitude) > maxval(Pr_MagGrid)) then
          !--Skipping as outside limits on which magnitude prior was evaluated (In effect, this should be close to the limit of the input catalogue magnitude)--! 
          !--NOTE: THis could be a source of bias, where a single lensed galaxy is close to the magnitude limit of the survey: delensing it takes it out of the limits in which the magnitude was evaluated--!
          Posterior_Flag = '03'
          Posterior_perGalaxy_Redshift(z) = 1.e-100_double
          cycle
       end if

       !--Get de-lensed version of measured quantities
       Theta_0 = delens_Source_Size(Effective_Magnification, Theta, use_lnT)
       Magnitude_0 = delens_Source_Magnitude(Effective_Magnification, Magnitude)
          
       !---Set RedshiftPDF, to unity if redshift known, and to a point in a particular distribution if it is not

       if(Known_Redshift) then
          RedshiftPDF = 1.e0_double
       else

          if(use_lookup) then
             RedshiftPDF = redshift_Distribution_byLookup(Magnitude_0, RedshiftGrid(z))
          else
             RedshiftPDF = CH08_Redshift_Distribution_Scalar(Magnitude_0, RedshiftGrid(z))
          end if

          !redshift_Distribution_byLookup(Magnitude_0, RedshiftGrid(z))
          !--Use this to evaluate at every level CH08_Redshift_Distribution_Scalar(Magnitude + 2.5e0_double*dlog10(Effective_Magnification), RedshiftGrid(z))
       end if
       
       !--Construct the correctly renormalised prior
       !Renormalisation_Size_Limits =  delens_Source_Size(Effective_Magnification, iSize_Limits, use_lnT)!!DEPRECATED iSize_Limits/dsqrt(Effective_Magnification)
       !Renormalisation_Magnitude_Limits = delens_Source_Magnitude(Effective_Magnification, iMag_Limits) !!DEPRECATED iMag_Limits + 2.5e0_double*dlog10(Effective_Magnification)
       

       
       if(Cuts_Renormalise_Likelihood) then             

!          if(z ==1) print *, 'KAPPA RENORMALISATION TURNED ON'
          if(Galaxy_Posterior_Method == 1 .or. Galaxy_Posterior_Method == 2) then
             if(use_lookup_Renorm) then
                Renorm = Renormalisation_Lookup(Effective_Magnification, 1)
             else
                Renorm = Linear_Interp(Effective_Magnification, MagnificationGrid, SM_Pr_Renormalisation)
             end if
             if(Renorm == 0) then
                Kappa_Renormalised_Prior = 0.e0_double
             else
                Kappa_Renormalised_Prior = SM_Prior/Renorm
             end if
             Renorm = 0.e0_double
          elseif(Galaxy_Posterior_Method == 3) then
             if(use_lookup_Renorm) then
                Renorm = Renormalisation_Lookup(Effective_Magnification, 1)
             else
                Renorm = Linear_Interp(Effective_Magnification, MagnificationGrid, M_Pr_Renormalisation) ! DO NOT ALLOW TO EXTRAPOLATE
             end if
             if(Renorm == 0) then
                Kappa_Renormalised_MagPrior = 0.e0_double
             else
                Kappa_Renormalised_MagPrior = M_Prior/Renorm
             end if
             Renorm = 0.e0_double
          end if

          if((use_lookup_Renorm == .false.) .and. (Effective_Magnification < minval(MagnificationGrid) .or. (Effective_Magnification > maxval(MagnificationGrid)))) then
             print *, 'Magnification outwith limits detected for galaxy:', Theta, Magnitude
             print *, 'And for lensing parameters:'
             print *, Position, Cluster_Position, Alpha, Lens_Redshift, Sigma_Crit*1.e18_double
             print *, 'Giving:', Renorm
             read(*,*)
          end if

       else
 !         if(i == 1 .and. z ==1 .and. c == 1) print *, 'KAPPA RENORMALISATION TURNED OFF'
          Kappa_Renormalised_Prior =  SM_Prior
          Kappa_Renormalised_MagPrior = M_Prior
       end if

       !--Find where m_0 and T_0 lie on the prior grid - speed-up for interpolation
       m0Index = 1 + nint((Magnitude_0-Pr_MagGrid(1))/(Pr_MagGrid(2)-Pr_MagGrid(1)))
       T0Index = 1 + nint((Theta_0-Pr_SizeGrid(1))/(Pr_SizeGrid(2)-Pr_SizeGrid(1)))

       select case(Galaxy_Posterior_Method)
       case(1) !--Size-Only--!
          !--Evaluate p_{theta_0, m_0}*p_{z|m_0} for the whole magnitude grid                                                                                                                                
          do m =  1, size(Pr_MagGrid)
             !--m_0, theta_0--!                                                                                                                                                                              
             if((Pr_MagGrid(m) < Renormalisation_Magnitude_Limits(1)) .or. (Pr_MagGrid(m) > Renormalisation_Magnitude_Limits(2))) then
                !-Since Integrand does not extend over this region anyway-!                                                                                                                                  
                Posterior_Flag = '11'
                Size_Only_Mag_Prior(m,:) = 1.e-100_double
             else
                if(Known_Redshift) then
                   Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)
                else
                   if(use_lookup) then
                      Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)*redshift_Distribution_byLookup(Pr_MagGrid(m), RedshiftGrid(z))
                   else
                      Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)*CH08_Redshift_Distribution_Scalar(Pr_MagGrid(m), RedshiftGrid(z))
                   end if
                   !--Use to evaulate at every position *CH08_Redshift_Distribution_Scalar(Pr_MagGrid(m), RedshiftGrid(z))
                end if
             end if
          end do
          
          if((Theta_0 > maxval(PR_SizeGrid)) .or. (Theta_0 < minval(PR_SizeGrid))) then
             !--Extrapolation--!                      
             Posterior_Flag = '11'
             Posterior_perGalaxy_Redshift(z) = 1.e-100_double
             cycle
          else
             
             !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-time by minimising the number of integrations required                                       
             allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
             do j = 1, size(PR_SizeGrid)-1
                if( (PR_SizeGrid(j)<= Theta_0) .and. ( PR_SizeGrid(j+1) > Theta_0) ) then
                   print *, 'Calling Integrate'
                   Size_Only_Prior(1) = Integrate(Pr_MagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
                   Size_Only_Prior(2) = Integrate(Pr_MagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
                   exit
                end if
             end do
             
             !--EVALUATE SIZE-ONLY LIKELIHOOD
             if(use_lnT) then
                Posterior_perGalaxy_Redshift(z) = Linear_Interp(Theta_0, PR_SizeGrid(j:j+1), Size_Only_Prior(:), ExValue =  1.e-100_double)
             else
                Posterior_perGalaxy_Redshift(z) = Linear_Interp(Theta_0, PR_SizeGrid(j:j+1), Size_Only_Prior(:), ExValue =  1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification))
             end if
          end if
          deallocate(Size_Only_Prior)
          
       CASE(2) !--Joint Size-Magnitude Analysis--!
          !--Error Catching
          if(Magnitude_0 < 21.e0_double) then
             print *, 'Possible problem with de-lensed magnitude - falls outwith bright limit of Scrabback Fit (21)'
             print *, Magnitude_0
          end if
          

          
          !--EVALUATE SIZE-MAGNITUDE LIKELIHOOD
!!$          if(use_lnT) then
!!$             Posterior_perGalaxy_Redshift(z) = Linear_Interp(Magnitude_0, Theta_0, Pr_MagGrid(m0Index-1:m0Index+1), PR_SizeGrid(T0Index-1:T0Index+1), Kappa_Renormalised_Prior(m0Index-1:m0Index+1,T0Index-1:T0Index+1), ExValue =  1.e-100_double)*RedshiftPDF
!!$          else
!!$             Posterior_perGalaxy_Redshift(z) = Linear_Interp(Magnitude_0, Theta_0, Pr_MagGrid(m0Index-1:m0Index+1), PR_SizeGrid(T0Index-1:T0Index+1), Kappa_Renormalised_Prior(m0Index-1:m0Index+1,T0Index-1:T0Index+1), ExValue =  1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification))*RedshiftPDF
!!$          end if

!!$          DEPREACTED FOR ISOALTION ROUTINE
!!$          !--EVALUATE SIZE-MAGNITUDE LIKELIHOOD
          if(use_lnT) then
             Posterior_perGalaxy_Redshift(z) = Linear_Interp(Magnitude_0, Theta_0, Pr_MagGrid, PR_SizeGrid, Kappa_Renormalised_Prior, ExValue = 1.e-100_double)*RedshiftPDF
          else
             Posterior_perGalaxy_Redshift(z) = Linear_Interp(Magnitude_0, Theta_0, Pr_MagGrid, PR_SizeGrid, Kappa_Renormalised_Prior, ExValue = 1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification))*RedshiftPDF
          end if
          
       case(3) !--Magnitude-only
          !allocate(Mag_Only_Prior(size(Pr_MagGrid))); Mag_Only_Prior = 0.e0_double                                                                                     
        
          !---DEPRACTED AFTER SOURCE SAMPLE SEPERATION CHANGES (9Apr2015)
!!$          if(Method == 4) then !-if Size_Mag+Mag and here, then this lies below the lower size limit and *must* therefore be correctly renormalised. This requires p_[m_0] to be calculated over the whole grid [Note: this may then take longer]
!!$
!!$             print *, 'RECONSTRUCTING MAG PRIOR DISTRIBUTION AS PART OF SIZEMAG+MAG ANALYSIS'
!!$
!!$             do j = 1, size(Pr_MagGrid)          
!!$                !--Construct p_[m_0] = int p_[theta_0, m_0], where int is understood to integrate over magnifcation-dependant survey size limits.
!!$                if( ( Pr_MagGrid(j)<=  Magnitude_0 ) .and. (  Pr_MagGrid(j+1) > Magnitude_0) ) then
!!$                   Unlensed_Mag_Index = j
!!$                end if
!!$                Mag_Only_Prior(j) = Integrate(Pr_SizeGrid, Kappa_Renormalised_Prior(j,:), 2, lim = Renormalisation_Size_Limits)        
!!$             end do
!!$             !--Renormalise
!!$             Mag_Only_Prior = Mag_Only_Prior/Integrate(Pr_MagGrid, Mag_Only_Prior, 2, lim = Renormalisation_Magnitude_Limits)
!!$             
!!$             
!!$          else !--This is the case of Magnitude only analysis
!!$
!!$             !--Construct p_[m_0] as the magnitude distribution for a sample of galaxies between cuts-corrected survey size limits                                                
!!$
!!$             !--Find point where de-lensed magnitude lie on intrinsic distribution - could also use locate (NR).                           
!!$             do j = 1, size(Pr_MagGrid)-1          
!!$                !--Construct p_[m_0] = int p_[theta_0, m_0], where int is understood to integrate over magnifcation-dependant survey size limits.
!!$                !--Renormalisation here considered to be set by the overall renormalisation of the size-magnitude distribution. This is not true if considering the mag-only case of Size-Mag+Mag
!!$                if( ( Pr_MagGrid(j)<=  Magnitude_0 ) .and. (  Pr_MagGrid(j+1) > Magnitude_0) ) then
!!$                   Unlensed_Mag_Index =j
!!$                   Mag_Only_Prior(j) = Integrate(Pr_SizeGrid, Kappa_Renormalised_Prior(j,:), 2, lim = Renormalisation_Size_Limits)                                                                      
!!$                   Mag_Only_Prior(j+1) = Integrate(Pr_SizeGrid, Kappa_Renormalised_Prior(j+1,:), 2, lim = Renormalisation_Size_Limits)                                                                    
!!$                   exit                                                                                                                                                                                  
!!$                end if
!!$             end do
!!$            
!!$          end if
!!$
!!$
!!$          !---EVALUATE MAGNITUDE-ONLY LIKELIHOOD                   
!!$          Posterior_perGalaxy_Redshift(z) = Linear_Interp(Magnitude_0, Pr_MagGrid(Unlensed_Mag_Index:Unlensed_Mag_Index+1), Mag_Only_Prior(Unlensed_Mag_Index:Unlensed_Mag_Index+1), ExValue =  1.e-100_double)*RedshiftPDF
!!$
!!$          !--Could also reasonably pass in the whole Mag_Only_Prior, but interpoaltion may take longer
!!$          deallocate(Mag_Only_Prior)                                                                                                                                                                   

          !---Added after source sample split (9th April 2015)
          Posterior_perGalaxy_Redshift(z) = Linear_Interp(Magnitude_0, Pr_MagGrid(:), Kappa_Renormalised_MagPrior(:), ExValue =  1.e-100_double)*RedshiftPDF


          !---------------------------------------------------------------------------------------------------------------------                                                                             
          
          !--The following gives unbiased results if no size cuts are used, however is known to be biased in the presence of size cuts, I believe that this is due to the fact that the prior itself should depend on the size cuts in the presence of a size-magnitude correlation. THIS SHOULD BE DELETED IF CODE IS VERIFIED AS UNBIASED AFTER 23RD JAN 2015
          
          !--Renormalise Magnitude Prior Distribution--!                                                                                                                                                     
!!$          if(Cuts_Renormalise_Likelihood) then
!!$             !--This method assumes not cuts on galaxy size, only cuts on magnitude
!!$             Renorm = Linear_Interp(Effective_Magnification, MagnificationGrid,  M_Pr_Renormalisation, ExValue = 1.e30_double)
!!$             if(Renorm == 0) then
!!$                Kappa_Renormalised_MagPrior = 0.e0_double
!!$             else
!!$                Kappa_Renormalised_MagPrior = M_Prior/Renorm
!!$             end if
!!$          end if
!!$          
!!$          Posterior_perGalaxy_Redshift(z) = Linear_Interp(Magnitude + 2.5e0_double*dlog10(Effective_Magnification), Pr_MagGrid, Kappa_Renormalised_MagPrior, ExValue =  1.e-100_double)*RedshiftPDF
          
       case default
          STOP 'DM_Profile_Variable_Posterior - Error in choosing method of finding posterior'
       END select
       
 
    END do !--End of Redshift Loop
    
    !-Integrate over nuisance parameters (redshift)
    if(size(Posterior_perGalaxy_Redshift,1) == 1) then
       !--Redshift was taken to be exact                                                                                                                                                                     
       Likelihood_atVirialRadius_MultipleCluster_perSource = Posterior_perGalaxy_Redshift(1)
    else
       !--Integrate over the redshift Information--!                                                                                                                                                         
       Likelihood_atVirialRadius_MultipleCluster_perSource = TrapInt(RedshiftGrid, Posterior_perGalaxy_Redshift(:))
    end if

    !--Deallocate Internal Declarations
    deallocate(Posterior_perGalaxy_Redshift)
    
    !--Close open files
    if (debug) then
       close(56)
    end if
    
  end function Likelihood_atVirialRadius_MultipleCluster_perSource

  subroutine get_Likelihood_Evaluation_Precursors(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid,Renormalisation_by_Magnification, MagOnly_Renormalisation_by_Magnification, PriorSizeGrid, PriorMagGrid, MagPrior, TwoDMagPrior, Prior, Survey_Magnitude_Limits, Survey_Size_Limits, MagSample_Magnitude_Limits, MagSample_Size_Limits, Lens_Redshift, Redshift_Limit, Output_Prefix, use_lnSize)
    use Integration; use Cosmology, only: angular_diameter_distance_fromRedshift
    use gridintervals
    use Distributions, only: create_redshiftdistribution_lookuptable
    !--Routine that:
    !---Normalises the Prior Distributions to integrate to one within the survey limits
    !---Sets up Sigma_Critical evalaution for all lenses considered
    !---Sets up a redshift grid on which the redshiftPDF (marginalisation) and Sigma_Critical will be evaluated
    !---Determines the normalisation of the prior as a function of magnification factor in the presence of cuts.
    !---

    !---MAG PRIOR NOTES:
    !-- MagPrior and TwoDMagPrior refer to the sample selection where the measurement will be taken on magnitudes only. Where non-trivial source seperation is used (e.g. by SNR or size cuts)
    !--Part of the code, marked with **, has been disabled (10Jul2015), contains the method by which the mag-only sample is renormalised accounting for the application of a size cut. This was verified unbiased in teh case where all sources had sizes, however is biased when some sources do not have a bias attached (as lnT = NaN is not taken into account in prior distribution construction). It may be possible to edit ot accoutn for this
    !---This should *always* be called before any of the evaluation routines above

    !--Input
    real(double):: PriorSizeGrid(:), PriorMagGrid(:), TwoDMagPrior(:,:), Prior(:,:) !-Magnitude, Size-!
    real(double), intent(inout),allocatable::  MagPrior(:) 
    real(double):: Survey_Magnitude_Limits(:,:), Survey_Size_Limits(:,:), Redshift_Limit, MagSample_Magnitude_Limits(:,:), MagSample_Size_Limits(:,:) !--First dimension allows byCluster
    real(double):: Lens_Redshift(:)
    character(*):: Output_Prefix
    
    logical,optional:: use_lnSize

    !--Output
    real(double),allocatable:: Survey_Renormalised_Prior(:,:), Survey_Renormalised_MagPrior(:)
    real(double),allocatable:: RedshiftGrid(:)
    real(double),allocatable:: Sigma_Crit(:,:), MagnificationGrid(:),Renormalisation_by_Magnification(:,:), MagOnly_Renormalisation_by_Magnification(:,:) !--First dimension allows byCluster
    real(double):: Renorm !Discardable

    !--Internal
    real(double), dimension(size(TwoDMagPrior,1), size(TwoDMagPrior,2)):: Survey_Renormalised_2DMagPrior
    real(double),dimension(2):: Renormalisation_Magnitude_Limits, Renormalisation_Size_Limits
    integer:: z, C, i
    real(double):: D_l, D_ls, D_s
    logical:: use_lnT

    integer,parameter:: nRedshift_Sampling = 75
    real(double)::Redshift_Lower, Redshift_Higher = 6.e0_double

    integer:: nMagnificationGrid = 2000
    real(double):: MagFactorGridLower = 0.9e0_double, MagFactorGridHigher = 65.e0_double

    character(5):: fmt

    print *, ' '
    write(*,'(A)') '________ Setting up Precusors to Bayesian Analysis_______________'

    use_lnT = .false.
    if(present(use_lnSize)) use_lnT = use_lnSize

    print *, ' '
    print *, 'NOTE: Renormalisation defined with respect to renormalisation for cluster 1 across the whole survey'
    print *, ' '

    allocate(Survey_Renormalised_Prior(size(Prior,1),size(Prior,2))); Survey_Renormalised_Prior = 0.e0_double
    !--Renormalise the prior within these size and magnitude limits--!                                                                                                           
    Renorm = Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits(1,:), lim2 = Survey_Size_Limits(1,:))
    print *, 'Renormalisation of the intrinsic distribution:', Renorm
    print *, sum(Survey_Renormalised_Prior), Prior(1:5,1:5)
    Survey_Renormalised_Prior = Prior/Renorm
    !--Allow for seperate renormalisation of the magnitude prior. This is required if the magnitude prior is constructed from galaxies which are excluded from the joint size-magnitude analysis, e.g. due to size cuts--!                                    

    allocate(Survey_Renormalised_MagPrior(size(PriorMagGrid))); Survey_Renormalised_MagPrior = MagPrior
    Renorm = Integrate(PriorMagGrid, Survey_Renormalised_MagPrior, 2, lim = (/MagSample_Magnitude_Limits(1,:)/))
    Survey_Renormalised_MagPrior = Survey_Renormalised_MagPrior/Renorm
    print *, 'Magnitude Distribution integrated between:', MagSample_Magnitude_Limits, ' with renormalisation:', Renorm
    print *, ' '
    
    !---** ---------------------------------------------------------------------------------------------------------------------------------**------------!
!!$    Renorm = Integrate(PriorMagGrid, PriorSizeGrid, TwoDMagPrior, 2, lim1 = MagSample_Magnitude_Limits, lim2 = MagSample_Size_Limits)
!!$    print *, 'Renormalisation of the intrinsic magnitude distribution:', Renorm
!!$    Survey_Renormalised_2DMagPrior = TwoDMagPrior/Renorm
    !---** ---------------------------------------------------------------------------------------------------------------------------------**------------!
    

    !--Set up the redshift grid--!                                                                                                                                                                                
    Redshift_Lower = Redshift_Limit
    allocate(RedshiftGrid(nRedshift_Sampling)); RedshiftGrid = 0.e0_double
    do z = 1, nRedshift_Sampling
       RedshiftGrid(z) = Redshift_Lower + (z-1)*((Redshift_Higher-Redshift_Lower)/(nRedshift_Sampling-1))
    end do

    !--Get Sigma_Critical for each point on the Redshift PDF grid--!                                        

    allocate(Sigma_Crit(size(Lens_Redshift),size(RedshiftGrid))); Sigma_Crit = 0.e0_double
    do C = 1, size(Lens_Redshift)
       do z = 1, size(RedshiftGrid)
          if(RedshiftGrid(z) <= Lens_Redshift(C)) then
             Sigma_Crit(C,z) = 1.e100_double !as multiplied by 10^18 before passing to mass profile code
             cycle
          end if
          D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift(C))
          D_s = angular_diameter_distance_fromRedshift(0.e0_double, RedshiftGrid(z))
          D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift(C), RedshiftGrid(z))
          Sigma_Crit(C,z) = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!
       end do
    end do

    allocate(MagnificationGrid(nMagnificationGrid)); MagnificationGrid = 0.e0_double
    allocate(Renormalisation_by_Magnification(size(Survey_Size_Limits,1),size(MagnificationGrid))); Renormalisation_by_Magnification = 0.e0_double
    allocate(MagOnly_Renormalisation_by_Magnification(size(magSample_Size_limits,1),size(MagnificationGrid))); MagOnly_Renormalisation_by_Magnification = 0.e0_double

    !---Construct the Renormalisation of the prior distribution as a function of magnification. Prior is allowed to be non-zeros outside the survey limits, since these limits should only apply to lensed quantities [not that survey limits here are most often used to detail user-imposed limits, such as magnitude or size cuts]. Only cuts on the prior distribution are necessary
    MagFactorGridHigher = min(1.0e0_double*(10.e0_double**((maxval(PriorMagGrid)-minval(Survey_Magnitude_Limits(:,1)))/2.5e0_double)), 1000.) !Allow maximum to be 100 only

    if(use_lookup_Renorm) then
       STOP 'MAGNIFICATION GRID LOOKUP HAS BEEN DISABLED UPON USE OF BYCLUSTER RENORMALISATION'
       call create_Renormalisation_Lookup((/MagFactorGridLower, MagFactorGridHigher, 1.e0_double*size(MagnificationGrid)/), Survey_Size_Limits, Survey_Magnitude_Limits, MagSample_Magnitude_Limits,  PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, use_lnT)
    else
       
       call logscale_min(MagFactorGridLower, MagFactorGridHigher, size(MagnificationGrid), MagnificationGrid)
       print *, 'Getting Renormalisation'
       DO C = 1, size(Survey_Size_Limits,1) 
          do i = 1, size(MagnificationGrid)
             !MagnificationGrid(i) = MagFactorGridLower + (i-1)*((MagFactorGridHigher- MagFactorGridLower)/(size(MagnificationGrid)-1))
             Renormalisation_Size_Limits = delens_Source_Size(MagnificationGrid(i), Survey_Size_Limits(c,:), use_lnT)
             Renormalisation_Magnitude_Limits = delens_Source_Magnitude(MagnificationGrid(i),Survey_Magnitude_Limits(c,:))
             
             !----TESTING----!
!!$       print *, 'Survey Size Limits:', Survey_Size_Limits, ' Mu:', MagnificationGrid(i), Renormalisation_Size_Limits
!!$       print *, 'Survey Magnitude Limits:', Survey_Magnitude_Limits, ' Mu:', MagnificationGrid(i), Renormalisation_Magnitude_Limits
!!$       print *, 'Grids:', minval(PriorMagGrid), maxval(PriorMagGrid), minval(PriorSizeGrid), maxval(PriorSizeGrid), sum(Survey_Renormalised_Prior)
             !---------------!
             Renormalisation_by_Magnification(c,i) = Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
             
             
             !--Get Magnitufaction -dependent renomralisation for mag-prior, usually constructed from sample for which size-information is not presnt
             Renormalisation_Magnitude_Limits = delens_Source_Magnitude(MagnificationGrid(i), MagSample_Magnitude_Limits(c,:))
             MagOnly_Renormalisation_by_Magnification(c,i) = Integrate(PriorMagGrid, Survey_Renormalised_MagPrior, 2, lim = Renormalisation_Magnitude_Limits)
             
             !---- ** Get magnitude-only signal renormalisation where the sample has been chosen from the main sample using size cuts - In this case, the effect fo size limits needs to be taken into account to produce unbiased results ** ----!
!!$       Renormalisation_Size_Limits = delens_Source_Size(MagnificationGrid(i), MagSample_Size_Limits, use_lnT)
!!$       Renormalisation_Magnitude_Limits = delens_Source_Magnitude(MagnificationGrid(i), MagSample_Magnitude_Limits)
!!$
!!$       MagOnly_Renormalisation_by_Magnification(i) = Integrate(PriorMagGrid,  PriorSizeGrid, Survey_Renormalised_2DMagPrior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
             !----**------------------------------------------------------------------------------------------------------------------------------------------------------------------------------**--------------!
          end do
       end DO
       print *, 'Got Renormalisation'
    end if

    write(fmt,'(I5)') size(Renormalisation_by_Magnification,1)
    open(unit = 53, file = trim(adjustl(Output_Prefix))//'Renormalisation_by_Magnification.dat')
    do i = 1, size(Renormalisation_by_Magnification,2)
       write(53, '('//trim(adjustl(fmt))//'(e14.7,x))') MagnificationGrid(i), Renormalisation_by_Magnification(:,i)
    end do
    write(*,'(2(A))') 'File output: ', trim(adjustl(Output_Prefix))//'Renormalisation_by_Magnification.dat'

    open(unit = 53, file = trim(adjustl(Output_Prefix))//'Magnitude_Renormalisation_by_Magnification.dat')
    do i = 1, size(MagOnly_Renormalisation_by_Magnification,2)
       write(53, '('//trim(adjustl(fmt))//'(e14.7,x))') MagnificationGrid(i), MagOnly_Renormalisation_by_Magnification(:,i)
    end do
    write(*,'(2(A))') 'File output: ', trim(adjustl(Output_Prefix))//'Magnitude_Renormalisation_by_Magnification.dat'


    !---Set up integrated MagPrior as the integrand of the 2D Mag Prior
    !-- ** Use this where the mag-only case comes as the result of a size cut
!!$    print *, 'Getting MagPrior'
!!$    if(allocated(MagPrior) == .false.) allocate(MagPrior(size(PriorMagGrid)))
!!$    do i= 1, size(PriorSizeGrid)
!!$       MagPrior(i) = Integrate(PriorSizeGrid, Survey_Renormalised_2DMagPrior(i,:), 2, lim = MagSample_Size_Limits)
!!$    end do
!!$    print *, 'Got MagPrior'
!!$    allocate(Survey_Renormalised_MagPrior(size(MagPrior))); Survey_Renormalised_MagPrior = MagPrior/Integrate(PriorMagGrid, MagPrior, 2, lim = MagSample_magnitude_Limits)
!!$    print *, 'Got Renormalised MagPrior'
    !----**------------------------------------------------------------------------------------------------------------------------------------------------------------------------------**--------------!

    !---- Set up Lookup Table to p(z|m) is used - on the same grid as that used to evaluate posterior
    if(use_lookup) call create_redshiftDistribution_LookupTable((/minval(PriorMagGrid), 1.2*maxval(PriorMagGrid), 0.01e0_double/), (/minval(RedshiftGrid), maxval(RedshiftGrid), (RedshiftGrid(2)-RedshiftGrid(1))/), OutputDirectory = trim(adjustl(Output_Prefix)))

!!$    !---- Set up Lookup Table to p(z|m) is used - on a finer same grid as that used to evaluate posterior (could be useful for interpolation, but probably uneccessary)
!!$    if(use_lookup) call create_redshiftDistribution_LookupTable((/minval(PriorMagGrid), 1.2*maxval(PriorMagGrid), 0.01e0_double/), (/minval(RedshiftGrid), maxval(RedshiftGrid), 0.5e0_double*(RedshiftGrid(2)-RedshiftGrid(1))/), OutputDirectory = trim(adjustl(Output_Prefix)))

    write(*,'(A)') '________ Successfully got Precusors_______________'
    print *, ' '

  end subroutine get_Likelihood_Evaluation_Precursors

  subroutine create_Renormalisation_Lookup(GridLimits, Survey_T_Limits, Survey_M_Limits, S2_M_Limits,  PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, use_lnSize)
    use integration, only: integrate

    real(double), intent(in):: GridLimits(3) !--Lower, Upper, nGrid
    real(double), intent(in):: Survey_T_Limits(2), Survey_M_Limits(2), S2_M_Limits(2) 
    real(double), intent(in):: PriorMagGrid(:), PriorSizeGrid(:), Survey_Renormalised_Prior(:,:), Survey_Renormalised_magPrior(:)
    logical, intent(in):: use_lnSize

    !--Internal
    integer:: nM, i
    real(double):: dM

    real(double):: Renormalisation_Size_Limits(2), Renormalisation_Magnitude_Limits(2)

    if(allocated(LT_Renormalisation%muGrid)) deallocate(LT_Renormalisation%muGrid)
    if(allocated(LT_Renormalisation%S1_LT)) deallocate(LT_Renormalisation%S1_LT)
    if(allocated(LT_Renormalisation%S2_LT)) deallocate(LT_Renormalisation%S2_LT)

    nM = nint(GridLimits(3))

    allocate(LT_Renormalisation%muGrid(nM)); LT_Renormalisation%muGrid = 0.e0_double
    allocate(LT_Renormalisation%S1_LT(nM)); LT_Renormalisation%S1_LT = 0.e0_double
    allocate(LT_Renormalisation%S2_LT(nM)); LT_Renormalisation%S2_LT = 0.e0_double

    !-- Set up muGrid - This needs to be replicated in the lookup part
    dM = log(GridLimits(2)/GridLimits(1))/(nM-1)
    do i =0, nM-1
       LT_Renormalisation%muGrid(i+1) = GridLimits(1)*dexp(i*dM)
    end do

    !--Construct Renormalisation Part
    print *, 'Getting Renormalisation'
    do i = 1, size(LT_Renormalisation%muGrid)
       Renormalisation_Size_Limits = delens_Source_Size(LT_Renormalisation%muGrid(i), Survey_T_Limits, use_lnSize)
       Renormalisation_Magnitude_Limits = delens_Source_Magnitude(LT_Renormalisation%muGrid(i),Survey_M_Limits)

       LT_Renormalisation%S1_LT(i) = Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)

       !--Get Magnification -dependent renomralisation for mag-prior, usually constructed from sample for which size-information is not presnt
       Renormalisation_Magnitude_Limits = delens_Source_Magnitude(LT_Renormalisation%muGrid(i), S2_M_Limits)
       LT_Renormalisation%S2_LT(i) = Integrate(PriorMagGrid, Survey_Renormalised_MagPrior, 2, lim = Renormalisation_Magnitude_Limits)

    end do
    print *, 'Got Renormalisation'


  end subroutine create_Renormalisation_Lookup


  function Renormalisation_Lookup(mu, sample)
    use Interpolaters, only: Linear_Interp
    real(double), intent(in):: mu
    integer, intent(in):: sample

    real(double):: Renormalisation_Lookup

    !--Internal
    real(double):: dM
    integer:: nM, index
    real(double), dimension(size(LT_Renormalisation%muGrid)):: Renorm

    logical:: doInterpolate = .false.

    nM = size(LT_Renormalisation%muGrid)
    dM = log(LT_Renormalisation%muGrid(nM)/LT_Renormalisation%muGrid(1))/(nM-1)

    index = 1 + int(dlog(mu/LT_Renormalisation%muGrid(1))/dM)

    if(sample == 1) then
       Renorm = LT_Renormalisation%S1_LT
    elseif(sample == 2) then
       Renorm =LT_Renormalisation%S2_LT
    else
       STOP 'Renormalisation_Lookup - Invalid Sample entered'
    end if
       
    if((index > size(Renorm)) .or. (index == 0)) then
       !-- This could be dealt with by setting to a large number
       print *, 'Extrapolation error with lookup:', mu, minval(LT_Renormalisation%muGrid), maxval(LT_Renormalisation%muGrid)
       STOP
    END if

    if(doInterpolate) then
       Renormalisation_Lookup = Linear_Interp(mu, LT_Renormalisation%muGrid(index-1:index+1), Renorm(index-1:index+1))
    else
       Renormalisation_Lookup = Renorm(index)
    end if

  end function Renormalisation_Lookup

 subroutine Combine_Posteriors_Scalar(GridValue, Posteriors,  Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
    !--Combines the posterior on a single grid value (free parameter alpha), using a call to the "normal" (vector) combined posterior subroutine                                                                  
     real(double), intent(in):: GridValue,Posteriors(:) !-Galaxy-!                                                                                                                                                
     real(double), intent(out):: Combined_Posterior
     logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP

     !--Internal Declarations                                                                                                                                                                                     
     real(double),dimension(1):: tGrid, tCombinedPosterior
     real(double), dimension(size(Posteriors),1):: tPosteriors

     !--Set up internals                                                                                                                                                                                          
     tGrid = GridValue
     tPosteriors(:,1) = Posteriors

     call Combine_Posteriors(tGrid, tPosteriors, Combine_by_ln, .false., Return_lnP, tCombinedPosterior)
     Combined_Posterior = tCombinedPosterior(1)

   end subroutine Combine_Posteriors_Scalar

  subroutine Combine_Posteriors_Vector(PosteriorGrid, Posteriors, Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
    !--Combines posteriors by looping over the first dimension--!                                                                                                                                                 
    real(double), intent(in):: PosteriorGrid(:),Posteriors(:,:) !-Galaxy, Grid/Value-!                                                                                                                            
    real(double), intent(out):: Combined_Posterior(:)
    logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP

    integer::c, j
    real(double)::Renorm, Combination_Normalisation
    logical::iDoRenormalise

    integer::nPosteriorsSkipped

    INTERFACE
       subroutine Combine_Posteriors_Vector(PosteriorGrid, Posteriors, Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
         use Param_Types
         real(double), intent(in):: PosteriorGrid(:), Posteriors(:,:)
         real(double), intent(out):: Combined_Posterior(:)
         logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP
       END subroutine Combine_Posteriors_Vector
    END INTERFACE
    
    !--Set renormalisation by input. If a single alpha value entered, then do not renormalise                                                                                                                     
    iDoRenormalise = Renormalise
    if(size(PosteriorGrid) == 1) iDoRenormalise = .false.
    
    nPosteriorsSkipped = 0
    if(Combine_by_ln == .false.) then
       Combination_Normalisation = 1.e0_double/maxval(Posteriors)!or 1.e0_double/(0.5e0_double*maxval(Posteriors(c,:))) within loop                                                                               
       Combined_Posterior = 1.e0_double
       do c = 1, size(Posteriors,1) !-Loop over galaxies-!                                                                                                                                                        
          !--Skip if zero (lnP not defined then) or NaN                                                                                                                                                           
          if(all(Posteriors(c,:) == 0.e0_double) .or. all(isNaN(Posteriors(c,:)))) then
             nPosteriorsSkipped = nPosteriorsSkipped + 1
             cycle
          end if

          !--Skip when `renormalised' posteriors are invalid (zero/negative)                                                                                                                                      
          if(all(Posteriors(c,:)*Combination_Normalisation == 0.e0_double)) then
             print *, 'Invalid Posterior for galaxy:', c, ' (==0) press [ENTER] to output an stop..'; READ(*,*)
             print *, Posteriors(c,:)
             STOP
          end if
          if(any(Posteriors(c,:)*Combination_Normalisation < 0.e0_double)) then
             print *, 'Invalid Posterior for galaxy:', c, ' (<0) presS [ENTER] to output an stop..'; READ(*,*)
             print *, Posteriors(c,:)
             STOP
          end if

          Combined_Posterior(:) = Combined_Posterior(:)*(Posteriors(c,:)*Combination_Normalisation)

          if(all(Combined_Posterior == 0.e0_double)) then
             print *, 'Invalid CPosterior for galaxy:', c, ' press (==0) [ENTER] to output an stop..'; READ(*,*)
             print *, Combination_Normalisation
             read(*,*)
             print *, Combined_Posterior(:)
             read(*,*)
             print *, Posteriors(c,:)
             STOP
          end if
          if(anY(Combined_Posterior < 0.e0_double)) then
             print *, 'Invalid CPosterior for galaxy:', c, ' press (<0) [ENTER] to output an stop..'; READ(*,*)
             print *, Posteriors(c,:)
             STOP
          end if

       end do
    else
       if(return_lnP) iDoRenormalise = .false.
       !       print *, 'Combining using logs'                                                                                                                                                                   
       !-Set lnP to a large negative value as default, equivalent to P ~ 0                                                                                                                                        
       Combined_Posterior = -1000.0_double
       Combination_Normalisation = size(Posteriors,1)
       do c = 1, size(Posteriors,1) !-Loop over galaxies-!                                                                                                                                                         
          !-Error Catching--!                        
          if(all(Posteriors(c,:) == 1.e-75_double) .or. all(isNaN(Posteriors(c,:)))) then
             nPosteriorsSkipped = nPosteriorsSkipped + 1
             cycle
          end if
          !--Sum log posteriors--!                                                                                                                                                                                
!!$          print *, 'CombinePosterior, loop:', c                                                                                                                                                                
!!$          print *, Combined_Posterior, dlog(Posteriors(c,:))                                                                                                                                                   
          where(Posteriors(c,:) <= 1.e-100_double)
             Combined_Posterior = Combined_Posterior - 1000.e0_double
          elsewhere
             Combined_Posterior = Combined_Posterior + dlog(Posteriors(c,:)) + 1.e0_double
          end where


          if(any(isNAN(Combined_Posterior(:)))) then
             print *, 'Any NaNs in Combined Posterior?, galaxy:',c, any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
             STOP
          end if
       end do
       !--Convert to PDF, not ln(PDF)--!                                                                                                                                                                          
       if(Return_lnP) then
          return
       else
          if(size(Combined_Posterior)/= 1) then
             Combined_Posterior = dexp(Combined_Posterior - maxval(Combined_Posterior))
          else
             Combined_Posterior = dexp(Combined_Posterior)
          end if
       end if
    end if

    if(any(isNAN(Combined_Posterior(:)))) then
       print *, 'Any NaNs in Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
       print *, 'Stopping'
       STOP
    END if

    if(nPosteriorsSkipped > 0) write(*,'(A,I4,A,I4,A)') '################### Combine_Posteriors - ', nPosteriorsSkipped, ' of ', size(Posteriors,1), ' posteriors were skipped as they we invlaid - NaNs or 0 ##########################'
    if((1.e0_double*nPosteriorsSkipped)/size(Posteriors,1) > 0.1) STOP 'Combine_Posteriors - number of skipped posteriors too large, stopping!'

    !--Renormalise--!                                                                                                                                                                                             
    if(iDoRenormalise) then
       Renorm = 0.e0_double
       do j = 1, size(Combined_Posterior)-1
          Renorm = Renorm + 0.5e0_double*(Combined_Posterior(j) + Combined_Posterior(j+1))*(PosteriorGrid(j+1)-PosteriorGrid(j))
       end do
       if(Renorm <= 0.e0_double) then
          print *, 'Renormalisation:', Renorm
          STOP 'Combine_Posteriors - Invalid Renormalisation for combined Posterior'
       end if
       Combined_Posterior(:) = Combined_Posterior(:)/Renorm
    end if

    if(any(isNAN(Combined_Posterior(:)))) then
       print *, 'Any NaNs in Renormalised Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
       print *, 'Stopping'
       STOP
    END if

  end subroutine Combine_Posteriors_Vector


end module Bayesian_Posterior_Evaluation
