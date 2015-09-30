module Distributions
  use Param_Types
implicit none

character(500),private::Reference_Catalogue = 'Catalogues/STAGES_COMBO17_gsMag_size_matched.pzcat'
integer::Reference_Catalogue_Columns(13) = (/-1,10,11,-1,-1,8,-1,-1,-1,14,16,17,5/)

INTERFACE CH08_redshift_distribution
   module procedure CH08_redshift_distribution_Scalar, CH08_redshift_distribution_Array
END INTERFACE CH08_redshift_distribution

INTERFACE redshift_Distribution_byLookup
   module procedure redshift_Distribution_byLookup_Array, redshift_Distribution_byLookup_Scalar
END INTERFACE redshift_Distribution_byLookup

type LookupTable
   real(double),allocatable:: xGrid(:), yGrid(:)
   real(double), allocatable:: LT(:,:)
end type LookupTable
type(LookupTable), private:: LT_RedshiftDistribution


contains


  subroutine randomly_Sample_From_Distribution(Grid, PDF, Res)
    use Common_Functions, only: return_Random_Set
    !--This could easily be modifed to use library routines, but this has not been done (19 Dec 2014)
    !--Populates the array Res by Randomly Sampling from the PDF described by Grid and PDF--!
    !--Assumes that the grid is equally binned, and the PDF is renormalised--!
    real(double), intent(in)::Grid(:), PDF(:)
    real(double), intent(out):: Res(:)

    real(double)::dGrid

    !--Random Number Generator Declarations--!
    real(double),dimension(:),allocatable::Ran
    Integer,allocatable::seed(:)
    integer::Nseed, Clock
    integer::NRandom 

    integer::i, F

    !--Intrinsic PDF declarations--!
    real(double), allocatable::CumulativePDF(:)
    real(double), allocatable::Renormalised_PDF(:)
    real(double)::AreaUnderCurve

    if(size(Grid) /= size(PDF)) STOP 'randomly_Sample_From_Distribution - Grid and PDF entered not conformal'

    NRandom = size(Res) !-1 random number for each galaxy-!
    allocate(Ran(NRandom)); Ran = 0.e0_double
    call RANDOM_SEED(size = NSeed)
    allocate(Seed(NSeed))
    call SYSTEM_CLOCK(COUNT = Clock)
    seed = Clock + (/ (i-1,i=1,NSeed) /)
    call RANDOM_SEED(PUT = seed)
    deallocate(Seed); NSeed = 0; Clock = 0
    
    call RANDOM_NUMBER(Ran)

!!$    allocate(Ran(nRandom)); Ran = 0.0e0_ldp
!!$    Ran = return_Random_Set(nRandom) !#Method = 0

    if(size(Grid) < 4) STOP 'randomly_Sample_From_Distribution - Grid is too small, something is possibly up'
    dGrid =  Grid(2)-Grid(1)!0.5e0_double*(Grid(4)-Grid(2))

    !--Renormalise PDF such that CumulativePDF always reaches 1--!
    AreaUnderCurve = dGrid*PDF(1)
    do i =2, size(PDF)
       dGrid = (Grid(i)-Grid(i-1)) !-0.5e0_double*(Grid(i+1)-Grid(i-1))-!
       AreaUnderCurve = AreaUnderCurve + dGrid*PDF(i)
    end do
    allocate(Renormalised_PDF(size(PDF))); Renormalised_PDF = 0.e0_double
    Renormalised_PDF = PDF/AreaUnderCurve

    !-Construct Cumulative PDF-!
    dGrid =Grid(2)-Grid(1)
    allocate(CumulativePDF(size(PDF))); CumulativePDF = 0.e0_double
    CumulativePDF(1) = Renormalised_PDF(1)*dGrid
    do i = 2, size(Renormalised_PDF)
       dGrid = (Grid(i)-Grid(i-1))
       CumulativePDF(i) = CumulativePDF(i-1)+(Renormalised_PDF(i)*dGrid)
       if(CumulativePDF(i)-1.e0_double > 1.e-6_double) then
          print *, CumulativePDF(i)-1.e0_double
          STOP 'randomly_Sample_From_Distribution - Cumulative PDF is getting too large, larger than rounding tolerance of 1'
       end if
    end do
    !--If by integration error it is within a rounding tolerance of 1, then the cumualtive PDF should be set to one at it's max value

    if(dabs(CumulativePDF(size(CUmulativePDF))-1.e0_double) <= 1.e-6_double) then
       CumulativePDF(size(CUmulativePDF)) = 1.e0_double
    end if

    if(any(Ran < 0.e0_double) .or. any(Ran > 1.e0_double)) STOP 'randomly_Sample_From_Distribution - FATAL ERROR - Random array values lies outside bounds (0,1)'
    if(any(CumulativePDF < 0.e0_double) .or. any(CumulativePDF > 1.e0_double+1.e-3_double)) then
       print *, 'Cumulative PDF:', CumulativePDF
       print *, any(CumulativePDF < 0.e0_double), any(CumulativePDF > 1.e0_double)
       STOP 'randomly_Sample_From_Distribution - FATAL ERROR - Cumulative PDF array values lies outside bounds (0,1)'
    end if

    if(maxval(CumulativePDF) < 1.e0_double) then
       print *, maxval(CumulativePDF)
       STOP 'randomly_Sample_From_Distribution - Cumulatie PDF has not reached appropriate precision, and has not reached one, stopping'
    end if

    !--Randomly Sample from PDF-!
    !--Assumes that the PDF is NORMALISED--!
    do i= 1, size(Ran)
       do F = 1, size(CumulativePDF)-1
          !--First two cases account for fintie size of Cumulative PDF and Rounding Errors--!
          if(Ran(i) < CumulativePDF(1)) then
             Res(i) = Grid(1)
          elseif(Ran(i) > maxval(CumulativePDF)) then
             print *, Ran(i), maxval(CumulativePDF)
             STOP 'rANDOM NUMBER IS LARGER THAN CUMULATIVE PDF'
             Res(i) = Grid(size(Grid))
          elseif(CumulativePDF(F) == Ran(i)) then
             Res(i) = Grid(F)
             exit
          elseif( (CumulativePDF(F) < Ran(i)) .and. (CumulativePDF(F+1) > Ran(i)) ) then
             !-Linearly Interpolate-!
             Res(i) = Grid(F) + ((Grid(F+1)-Grid(F))/(CumulativePDF(F+1)-CumulativePDF(F)))*(Ran(i)-CumulativePDF(F))
             exit
          elseif(CumulativePDF(F+1) == Ran(i)) then
             Res(i) = Grid(F+1)
             exit
          elseif(F == size(CumulativePDF)-1) then
             print *, 'Loop:', i
             print *, 'Max Cumulative, Random:', CumulativePDF(size(CUmulativePDF)), Ran(i)
             print *, 'Min Cumulative, Random:', CumulativePDF(1), Ran(i)
             STOP 'Assign_Intrinsic_Sizes - FATAL ERROR - Cannot assign magnitude'
          end if
       end do
    end do 

  end subroutine randomly_Sample_From_Distribution


  !---------------REDSHIFT DISTRIBUTIONS--------------------------------------------------------!

  function redshift_Distribution_byLookup_Scalar(Apparent_Magnitude, Redshift, Renorm_Limits)
    real(double), intent(in)::Apparent_Magnitude, Redshift
    real(double):: redshift_Distribution_byLookup_Scalar
    real(double),intent(in),optional:: Renorm_Limits(2)

    real(double):: PDF(1)

    if(present(Renorm_Limits)) then
       PDF = redshift_Distribution_byLookup_Array(Apparent_Magnitude, (/Redshift/), Renorm_Limits)
    else
       PDF = redshift_Distribution_byLookup_Array(Apparent_Magnitude, (/Redshift/))
    end if

    redshift_Distribution_byLookup_Scalar = PDF(1)

  end function redshift_Distribution_byLookup_Scalar
    

  function redshift_Distribution_byLookup_Array(Apparent_Magnitude, Grid, Renorm_Limits) RESULT(PDF)
    !--Wrapper Routine for use of redshift distribution lookup. If lookup does not exisit/ has not been created, then this will create it
    use Interpolaters, only: Linear_Interp; use Integration, only: Integrate
    real(double), intent(in)::Apparent_Magnitude
    real(double), intent(in)::Grid(:)
    real(double), intent(in),optional::Renorm_Limits(2)
    real(double)::PDF(size(Grid))

    logical:: doInterpolate = .false.

    !--Internal Declarations
    integer:: nM, nZ, Z, zz
    real(double):: dz

    !--Lookup Declarations
    integer:: Index, zIndex
    real(double):: Width

    !--Temporary Declarations
    real(double), allocatable:: tPDF(:), LT_PDF(:)

    !---Error Catching Declarations
    integer, save:: nOutsideMagRange = 0

    !--Error Catching
    if(allocated(LT_RedshiftDistribution%LT) == .false.) then
       STOP 'redshift_Distribution_byLookup_Array: FATAL: Lookup Table not initialised'
    end if

    if(allocated(LT_RedshiftDistribution%xGrid) == .false. .or. allocated(LT_RedshiftDistribution%yGrid) == .false.) then
       STOP 'redshift_Distribution_byLookup_Array - Error in setting lookup table use: Grids not allocated'
    end if

    if(size(LT_RedshiftDistribution%xGrid) /= size(LT_RedshiftDistribution%LT,1) .or. size(LT_RedshiftDistribution%yGrid) /= size(LT_RedshiftDistribution%LT,2)) STOP 'redshift_Distribution_byLookup_Array - Error in setting lookup table use: Grids not conformal with Table'


    !---Getting Distribution
    if(Apparent_Magnitude > maxval(LT_RedshiftDistribution%xGrid) .or. Apparent_Magnitude < minval(LT_RedshiftDistribution%xGrid)) then
       nOutsideMagRange = nOutsideMagRange + 1
       print *, 'redshift_Distribution_byLookup_Array - WARNING - Apparent Magnitude falls outside lookup range. Constructing PDF individually for this case. This is the ', nOutsideMagRange, ' time'
       print *, 'Input, Range:', Apparent_Magnitude, minval(LT_RedshiftDistribution%xGrid), maxval(LT_RedshiftDistribution%xGrid)

       call CH08_redshift_distribution_Array(Apparent_Magnitude, Grid, tPDF)
       PDF = tPDF

       deallocate(tPDF)
    end if

    !--Find Magnitude on Grid
    Width =  LT_RedshiftDistribution%xGrid(2)- LT_RedshiftDistribution%xGrid(1)
    Index = nint((Apparent_Magnitude-LT_RedshiftDistribution%xGrid(1))/Width) + 1

    if(Index < 1 .or. Index > size(LT_RedshiftDistribution%xGrid)) then
       STOP 'LookupTable: Error getting apparent magnitude of grid'
    end if

    !--Store LT value for that magnitude
    allocate(LT_PDF(size(LT_RedshiftDistribution%LT,2))); LT_PDF = LT_RedshiftDistribution%LT(Index,:)

    if(present(Renorm_Limits)) then
       !--Renormalise PDF between entered limits. As this will be evaluated for each source, this is likely to correspond to a noticable increase in run-time if used. Alternative is, if a single renormalisation is used for the whole system is to do this once when the table is constructed
       LT_PDF = LT_PDF/Integrate(LT_RedshiftDistribution%yGrid, LT_PDF, 2, lim = Renorm_Limits)
    end if

    if(doInterpolate) then
       !--This could be reduced down if size of grid == 1 to reduce searching needed to interpolate
       PDF = Linear_Interp(Grid, LT_RedshiftDistribution%yGrid, LT_PDF, ExValue = 0.e0_double)
    else
       !--zGrid constructed as: (/ (zLimits(1) + (M-1)*zLimits(3), M = 1, size(LT_RedshiftDistribution%yGrid)) /)
       dz = (LT_RedshiftDistribution%yGrid(size(LT_RedshiftDistribution%yGrid))-LT_RedshiftDistribution%yGrid(1))/(size(LT_RedshiftDistribution%yGrid)-1)

       do zz = 1, size(Grid)
          zIndex = nint(1.e0_double+((Grid(zz)-LT_RedshiftDistribution%yGrid(1))/dz))

          if(zIndex > size(LT_PDF)) then
             print *, 'zGrid check (exceeding):', zindex, LT_RedshiftDistribution%yGrid(size(LT_RedshiftDistribution%yGrid)), grid(zz), dz
          end if

          !--Evaluate on coarse grid
          PDF(zz) = LT_PDF(zIndex)
       end do
    end if

    !--2D interpolation - More accurate (but different only < 1% when dm = 0.01 over index approach), but marginally slower
!!$    do z = 1, size(Grid)
!!$       PDF(z) = Linear_Interp(Apparent_Magnitude, Grid(z), LT_RedshiftDistribution%xGrid, LT_RedshiftDistribution%yGrid, LT_RedshiftDistribution%LT, ExValue = 0.e0_double)
!!$    end do

  end function redshift_Distribution_byLookup_Array


  subroutine create_redshiftDistribution_LookupTable(MagLimits, zLimits, OutputDirectory)
    !___--- Construct the redshift distribution lookup table (stored as a private global parameter)
    real(double), intent(in):: MagLimits(3), zLimits(3) !--Lower, Upper, Width
    character(*), intent(in),optional:: OutputDirectory

    real(double), allocatable:: tPDF(:)

    integer:: nZ, nM, M
    character(500):: fmt, filename

    INTERFACE 
       subroutine create_redshiftDistribution_LookupTable(MagLimits, zLimits, OutputDirectory)
         use Param_Types
         real(double), intent(in):: MagLimits(3), zLimits(3) !--Lower, Upper, Width
         
         character(*), intent(in),optional:: OutputDirectory
       end subroutine create_redshiftDistribution_LookupTable
    END INTERFACE

    print *, ' '
    print *, 'Constructing Redshift Distribution Look-Up Table............................'
    print *, ' '


    if(allocated(LT_RedshiftDistribution%xGrid)) deallocate(LT_RedshiftDistribution%xGrid)
    if(allocated(LT_RedshiftDistribution%yGrid)) deallocate(LT_RedshiftDistribution%yGrid)
    if(allocated(LT_RedshiftDistribution%LT)) deallocate(LT_RedshiftDistribution%LT)

    nM = 1 + nint((MagLimits(2)-MagLimits(1))/MagLimits(3))
    nZ = 1 + nint((zLimits(2)-zLimits(1))/zLimits(3))
    
    allocate(LT_RedshiftDistribution%xGrid(nM)); LT_RedshiftDistribution%xGrid = 0.
    allocate(LT_RedshiftDistribution%yGrid(nZ)); LT_RedshiftDistribution%yGrid = 0.
    allocate(LT_RedshiftDistribution%LT(nM, nZ)); LT_RedshiftDistribution%LT = 0.

    LT_RedshiftDistribution%xGrid = (/ (MagLimits(1) + (M-1)*MagLimits(3), M = 1, size(LT_RedshiftDistribution%xGrid)) /)
    LT_RedshiftDistribution%yGrid = (/ (zLimits(1) + (M-1)*zLimits(3), M = 1, size(LT_RedshiftDistribution%yGrid)) /)

    do M = 1, size(LT_RedshiftDistribution%xGrid)
       call CH08_redshift_distribution_Array(LT_RedshiftDistribution%xGrid(M), LT_RedshiftDistribution%yGrid, tPDF)

       if(size(tPDF)/= size( LT_RedshiftDistribution%LT,2)) STOP 'create_redshiftDistribution_LookupTable: Error in setting Redshift distribution for given magnitude: Returned PDF not conformal'

        LT_RedshiftDistribution%LT(M,:) = tPDF

       deallocate(tPDF)
    end do

    !--Output
    if(present(OutputDirectory)) then
       Filename = trim(adjustl(OutputDirectory))//'Redshift_Distribution_Lookup.dat'
       open(unit = 31, file = Filename)

       !-_Write Header_-!
       write(31, '(A)') '# First Row is Redshift Grid, First Column is magnitude, Other is p(z|m)'
       
       write(fmt,*) size(LT_RedshiftDistribution%yGrid)+1
       fmt = '('//trim(adjustl(fmt))//'(e14.7,x))'
       !-_Write First Row_-!
       write(31, fmt) 0., LT_RedshiftDistribution%yGrid
       do M = 1, size(LT_RedshiftDistribution%xGrid)
          write(31, fmt) LT_RedshiftDistribution%xGrid(M), LT_RedshiftDistribution%LT(M,:)
       end do

       close(31)
       write(*,'(A)') '--- Output Redshift Distribution Lookup to:', trim(Filename)
    end if


    print *, ' '
    print *, 'Redshift Distribution Look-Up Table Constructed'
    print *, ' '


  end subroutine create_redshiftDistribution_LookupTable

  subroutine CH08_redshift_distribution_Array(Apparent_Magnitude, Grid, PDF)
    !--Returns a PDF for the redshift distribution, used in CH08, which is a Smail et al. distribution with apparent-magnitude-dependent median redshift--!
    real(double), intent(in)::Apparent_Magnitude
    real(double), intent(in)::Grid(:)
    real(double), intent(out),allocatable::PDF(:)
    
    real(double):: zmed, alpha = 2.e0_double, beta = 1.5e0_double

    integer:: j

    integer,save:: callcount = 0

    if(allocated(PDF)) deallocate(PDF)
    allocate(PDF(size(Grid))); PDF = 0.e0_double

!    if(callcount == 0) PRINT *, 'MAGNITUDE DEPENDANCE OF P(Z|M) TURNED OFF'
    callcount = callcount + 1


    zmed = 0.29e0_double*(Apparent_Magnitude - 22.e0_double) + 0.31e0_double
!!$    if(zmed < 0.e0_double) then
!!$       print *, 'CH08_redshift_distributions - Median Redshift returned is negative, suggesting that galaxies which are too bright have been entered'
!!$       print *, 'zmed, m:', zmed, Apparent_Magnitude
!!$       STOP
!!$    end if

    call Analytic_Source_Redshift_Distribution(alpha, beta, zmed, Grid, PDF)

  end subroutine CH08_redshift_distribution_Array

  real(double) function CH08_redshift_distribution_Scalar(Apparent_Magnitude, Redshift)
    use nr, only: gammln
    
    real(double), intent(in)::Apparent_Magnitude, Redshift
    real(double):: zmed, alpha = 2.e0_double, beta = 1.5e0_double
    real(double)::PDF
    real(double):: z_0, Norm
    integer,save:: callcount  = 0

!    if(callcount == 0) PRINT *, 'MAGNITUDE DEPENDANCE OF P(Z|M) TURNED OFF'
    callcount = callcount + 1

    zmed = 0.29e0_double*(Apparent_Magnitude - 22.e0_double) + 0.31e0_double
    if(zmed < 0.e0_double) then
       print *, 'CH08_redshift_distributions - Median Redshift returned is negative, suggesting that galaxies which are too bright have been entered'
       print *, 'zmed, m:', zmed, Apparent_Magnitude
       STOP
    end if

    z_0 = zmed/1.412e0_double
    Norm = (beta/(z_0*dexp(gammln( (alpha+1.e0_double)/(beta) ))))

    PDF = Norm*((Redshift/z_0)**alpha)*dexp(-(Redshift/z_0)**beta)

!!    call Analytic_Source_Redshift_Distribution(alpha, beta, zmed, (/Redshift/), PDF)
    CH08_redshift_distribution_Scalar = PDF

  end function CH08_redshift_distribution_Scalar


  subroutine Analytic_Source_Redshift_Distribution(alpha, beta, z_med, Redshift, pdf, Output_Directory)
    use nr, only: gammln
    real(double), intent(in)::Redshift(:)
    real(double), intent(out),allocatable::pdf(:)
    real(double),intent(in):: alpha, beta, z_med
    character(*),optional::Output_Directory

    integer::Method = 1 !-1:Smail-!                                                                                                                                                                                                           
    integer::i

    real(double)::AreaUnderCurve, dRedshift

    !--Smail-!      
    real(double):: z_0, Norm
    
    if(allocated(PDF) == .false.) then
       allocate(PDF(size(Redshift)));
    end if
    PDF = -1.e0_double !--Set to negative to test for unassigned elements later

    AreaUnderCurve = 0.e0_double
    select case(Method)
    case(1) !-Smail-!
       !print *, 'Producing a Smail et al source redshift distribution with alpha = ', alpha, '; beta = ', beta, '; z_med  = ', z_med 
       if(z_med <= 0.) then
          print *, 'WARNING: Analytic Redshift Distribution: Likelihood returned as zero as median redshift is zero or negative'
          PDF = 1.e-100_double
          return
       end if

       z_0 = z_med/1.412e0_double
       Norm = (beta/(z_0*dexp(gammln( (alpha+1.e0_double)/(beta) ))))

       if(Norm < 0) STOP 'Analytic_Source_Redshift_Distribution - Error in renormalising the Smail distribution'

       do i = 1, size(Redshift)
          PDF(i) = Norm*((Redshift(i)/z_0)**alpha)*dexp(-(Redshift(i)/z_0)**beta)
!          if(i>1) AreaUnderCurve = AreaUnderCurve + 0.5e0_double*(Redshift(i)-Redshift(i-1))*(PDF(i)+PDF(i-1))
       end do
       !-Renormalise-!
!       PDF = PDF/AreaUnderCurve
    end select

    if(any(PDF < 0.e0_double)) STOP 'Source_Redshift_Distribution - FALTAL ERROR - PDF contains negatives'

   !-Output-!                                                                                                                                                                                                                                
    if(present(Output_Directory)) then
       open(32, file = trim(adjustl(Output_Directory))//'Source_Redshift_Distribution_Smail.dat')
       !print *, 'Source Redshift Distribution output to:' ,trim(adjustl(Output_Directory))//'Source_Redshift_Distribution_Smail.dat'                                                                                                         
    else
       open(32, file = 'Distributions/Source_Redshift_Distribution_Smail.dat')
       !print *, 'Source Redshift Distribution output to: Distributions/Source_Redshift_Distribution_Smail.dat'   
    end if
    do i = 1, size(PDF)
       write(32, *) Redshift(i), PDF(i)
    end do
    close(32)

  end subroutine Analytic_Source_Redshift_Distribution



  subroutine get_Redshift_Distribution_MagnitudeBinning_byCatalogue(MagBins, Redshifts, PDF, RefCat, Output_Dir)
    !--Constructs the redhsift distribution p(z,m) [or p(z|m)?] from a reference catalogue
    use Catalogues; use Statistics, only: variance_discrete, median
    real(double), intent(out),allocatable:: Redshifts(:)
    real(double), intent(out),allocatable::pdf(:,:) !-Redshift Bin, GridValue-!
    real(double),dimension(:,:),allocatable,intent(inout)::MagBins
    character(*), intent(in), optional:: Output_Dir
    type(Catalogue), intent(in)::RefCat

    character(500)::Catalogue_Filename
    integer::Catalogue_Columns(13)

    type(Catalogue)::Cat
    
    integer::nMags 
    real(double)::Renormalisation

    integer,parameter::nZ = 20
    real(double)::Z_lower, Z_Higher, dZ
    real(double),dimension(:,:),allocatable::ZBins
    
    integer::i,j,c
    
    character(5)::fmtstring
    
    type(Binned_Catalogue)::BCat
    real(double),allocatable::Temporary_Z_Array(:)

    !--Fit declarations [ to z_m =  a*m + b]--!
    real(double),allocatable::median_redshift(:)
    real(double),allocatable::mean_Mags(:)
    real(double)::fita, fitb, siga, sigb, chi2, fitq

    Catalogue_Filename = trim(Reference_Catalogue)
    Catalogue_Columns = Reference_Catalogue_Columns
    
    PRINT *, 'Reconstructing p(z,m) from Reference Catalogue:'
    Cat = RefCat
!!$    
!!$    PRINT *, 'Reconstructing p(z,m) distribution from the Catalogue:', trim(Catalogue_Filename),'....'
!!$    !--Read in Reference Catalogue--!
!!$    call Catalogue_ReadIn(Cat, trim(Catalogue_Filename), 'COMBO', Catalogue_Columns)

    !--Test for redshifts for all galaxies in reference catalogue--!

    Z_Higher = maxval(Cat%redshift); Z_Lower = 0.e0_double
    allocate(Redshifts(nZ)); Redshifts = 0.e0_double
    allocate(ZBins(nZ,2)); ZBins = 0.e0_double
    dZ = (( Z_Higher- Z_Lower )/(1.e0_double*(nZ-1)) )
    do i = 1, nZ
       !--Use i-2 so that no galaxies fall into the first bin--!                                                                                                                                                                         
       ZBins(i,1) = Z_Lower + (i-2)*dZ
       ZBins(i,2) = ZBins(i,1) + dZ
       if(i==nZ) ZBins(i,2) = ZBins(i,2) + 1.e-3_double*dZ
       Redshifts(i) = 0.5e0_double*(ZBins(i,1)+ZBins(i,2))
    end do
    Redshifts(1) = ZBins(1,2)

    !----START HERE----!
    if(allocated(MagBins) == .false.) then
       nMags = 5
       print *, 'Magnitude Limits not passed so getting Mag Limitsf ro nBin', nMags
       call Calculate_Bin_Limits_by_equalNumber(Cat%Absolute_Magnitude, nMags, MagBins)
    end if
    nMags = size(MagBins,1)
    !--Bin Catalogue by Magnitude--!
    call bin_catalogue_by_magnitude(Cat,MagBins,BCat)

    do i =1 , nMags
       print *, 'Mag Bin ', i, ' : ', MagBins(i,:)
    end do
    
    allocate(PDF(nMags, nZ)); PDF = 0.e0_double
    do i = 1, nMags
       allocate(Temporary_Z_Array(size(BCat%Cat(i)%Redshift))); Temporary_Z_Array = BCat%Cat(i)%Redshift
       do c = 1, size(Temporary_Z_Array)
          do j = 1, nZ
             if( (Temporary_Z_Array(c) > ZBins(j,1)) .and. (Temporary_Z_Array(c) <= ZBins(j,2)) ) then
                PDF(i,j) = PDF(i,j) + 1
                exit
             end if
          end do
       end do
       deallocate(Temporary_Z_Array)
    end do

    !--Renormalise for each Magnitude Bin--!                                                                                                                                                                                                 
    do i =1, nMags
       Renormalisation = 0.e0_double
       do j = 1, size(PDF,2)
          Renormalisation = Renormalisation + PDF(i,j)*(ZBins(i,2)-ZBins(i,1))
       end do
       if(Renormalisation> 0.e0_double) PDF(i,:) = PDF(i,:)/Renormalisation
    end do

    !--Output PDFS--!                                                                                                                                                                                                                   
    if(present(output_dir)) then
       open(unit = 49, file = trim(Output_Dir)//'Redshift_Distribution_MagnitudeBinning_Catalogue.dat')
           print *, 'Redshift Distribution, by Magnitude Bin, output to '//trim(Output_Dir)//'Redshift_Distribution_MagnitudeBinning_Catalogue.dat'
    else
       open(unit = 49, file = 'Distributions/Redshift_Distribution_MagnitudeBinning_Catalogue.dat')
       print *, 'Size Distribution, by Magnitude Bin, output to Distributions/Redshift_Distribution_MagnitudeBinning_Catalogue.dat'
    end if
    !--Write Header--!
    do j = 1, size(MagBins,1)
       write(49, '(A1, 2(e14.7,x))') '#', MagBins(j,:)
    end do
    write(49,'(A)')
    write(fmtstring, '(I5)') size(PDF,1)+1
    do j = 1, size(PDF,2)
       write(49, '('//trim(fmtstring)//'(e14.7,x))') Redshifts(j), PDF(:,j)
    end do
    close(49)

!!$    print *, 'Getting fit of median redshift to *absolute magnitude*:'
!!$    allocate(median_redshift(nMags)); allocate(mean_Mags(nMags));
!!$    do i = 1, nMags
!!$       mean_Mags(i) = 0.5e0_double*(sum(MagBins(i,:)))
!!$       median_redshift(i) = median(pdf(i,:), redshifts)
!!$    end do
!!$    call fit(mean_Mags, median_redshift, fita, fitb, siga, sigb, chi2, fitq) !--Doesn't use error for median redshift--!
!!$    print *, 'Redshift distribution fit to: z_m = ',fita, '*m +', fitb

    deallocate(median_redshift, mean_Mags)

  end subroutine get_Redshift_Distribution_MagnitudeBinning_byCatalogue

  subroutine angular_Diameter_Distance_Distribution(Cat, Magnitude_Bins, Output_Dir)
    use Catalogues; use Cosmology
    type(Catalogue),intent(in)::Cat
    real(double),optional:: MAgnitude_Bins(:,:)
    character(*), intent(in),optional::Output_Dir

    type(Binned_Catalogue)::BCat
    real(double),allocatable::MagBins(:,:)
    integer::nBin
    integer::MAgnitude_Type  = 1

    real(double)::AD_Lower, AD_Higher, dAD
    integer::nAD = 20
    real(double),allocatable::ADBins(:,:)
    real(double),allocatable::AD_Grid(:)
    real(double),allocatable::PDF(:,:) !-Bin, PDF-!
    
    real(double),allocatable::AD_PerGalaxy(:)
    real(double)::Renormalisation

    integer::b, c, ad_loop

    logical::here
    character(5)::fmtstring

    if(present(Magnitude_Bins)) then
       nBin = size(Magnitude_Bins,1)
       allocate(MagBins(nBin, size(Magnitude_Bins,2))); MagBIns = 0.e0_double
       MagBins = Magnitude_Bins
    else
       STOP 'not sorted yet'
    end if
    
    call bin_catalogue_by_magnitude(Cat,MagBins,BCat, Magnitude_Type)

    AD_Higher = angular_diameter_distance_fromRedshift_scalar(0.e0_double, 2.0e0_double*maxval(Cat%redshift))
    AD_lower = 1.e10_double
    do c = 1, size(Cat%Redshift)
       if(Cat%Redshift(c) < 0.e0_double) cycle
       if(Cat%Redshift(c) < AD_lower) AD_lower = Cat%Redshift(c)
    end do
    AD_lower = angular_diameter_distance_fromRedshift_scalar(0.e0_double, AD_Lower)

    allocate(AD_Grid(nAD)); AD_Grid = 0.e0_double
    allocate(ADBins(nAD,2)); ADBins = 0.e0_double
    dAD = (( AD_Higher- AD_Lower )/(1.e0_double*(nAD-1)) )
    if(dAD == 0.e0_double) STOP 'get_AD_Distribution_MagnitudeBinning_byCatalogue - Error setting up AD Grid, dAD = 0'
    do b = 1, nAD
       !--Use i-2 so that no galaxies fall into the first bin--!                
       ADBins(b,1) = AD_Lower + (b-2)*dAD
       ADBins(b,2) = ADBins(b,1) + dAD
       if(b==nAD) ADBins(b,2) = ADBins(b,2) + 1.e-3_double*dAD
       AD_Grid(b) = 0.5e0_double*(ADBins(b,1)+ADBins(b,2))
    end do
    AD_Grid(1) = ADBins(1,2)       
    
    
    allocate(PDF(nBin, size(AD_Grid))); PDF = 0.e0_double
    do b = 1, nBin
       !-Calculate Angular Diamater Distance per Galaxy in Catalogue-!
       allocate(AD_PerGalaxy(size(BCat%Cat(b)%RA))); AD_PerGalaxy = -1.e0_double
       do c = 1, size(AD_PerGalaxy)
          !-_Test for negative?--!
          AD_perGalaxy(c) = angular_diameter_distance_fromRedshift_scalar(0.e0_double, BCat%Cat(b)%Redshift(c))
          do ad_loop = 1, nAD
             if( (AD_PerGalaxy(c) > ADBins(ad_loop,1)) .and. (AD_PerGalaxy(c) <= ADBins(ad_loop,2)) ) then
                PDF(b,ad_loop) = PDF(b,ad_loop) + 1
                exit
             end if
          end do
       end do
       deallocate(AD_PerGalaxy)
    end do
                                           
    do b =1, nBin
       Renormalisation = 0.e0_double
       do ad_loop = 1, size(PDF,2)
          Renormalisation = Renormalisation + PDF(b,ad_loop)*(ADBins(b,2)-ADBins(b,1))
       end do
       if(Renormalisation> 0.e0_double) PDF(b,:) = PDF(b,:)/Renormalisation
    end do
                                           
    if(present(output_dir)) then

       !--Check for existence--!
       inquire(directory =  trim(Output_Dir), exist = here)
       if(here == .false.) call system('mkdir '//trim(Output_Dir))

       open(unit = 49, file = trim(Output_Dir)//'AngularDiameter_Distribution.dat')
           print *, 'Angular Diameter Distribution, by Magnitude Bin, output to '//trim(Output_Dir)//'AngularDiameter_Distribution.dat'
    else
       open(unit = 49, file = 'Distributions/AngularDiameter_Distribution.dat')
       print *, 'Angular Diameter Distribution, by Magnitude Bin, output to Distributions/AngularDiameter_Distribution.dat'
    end if
    !--Write Header--!
    do b = 1, size(MagBins,1)
       write(49, '(A1, 2(e14.7,x))') '#', MagBins(b,:)
    end do
    write(49,'(A)')
    write(fmtstring, '(I5)') size(PDF,1)+1
    do ad_loop = 1, size(PDF,2)
       write(49, '('//trim(fmtstring)//'(e14.7,x))') AD_Grid(ad_loop), PDF(:,ad_loop)
    end do
    close(49)


  end subroutine angular_Diameter_Distance_Distribution

  !----------------Magnitude Distributions------------------------------------------------------!
  subroutine produce_Magnitude_Distribution(MagGrid, PDF, RefCat, KDE_Smooth)
    use Catalogues; use Statistics, only: Discrete_Covariance; use Smoothing, only:KDE_Univariate_Gaussian; use Integration, only:Integrate
    !--Produces a magnitude distribution for the data--!
    !-Uses covariance for size-magnitude to set smoothing length. This could be edited to a simpler routine, but kept as is for speed's sake--!
    real(double), intent(inout),allocatable,dimension(:):: MagGrid
    real(double), intent(out),allocatable::pdf(:)
    type(Catalogue), intent(in):: RefCat
    logical, intent(in):: KDE_Smooth

    integer:: nSmoothed_Sampling_Mag = 110, nHist = 60
    real(double):: Higher, Lower, Mag_Limit_Convergence_Buffer = 0.3e0_double
    real(double), allocatable:: Bins(:)

    integer::i, c

    !--KDE Declarations--!
    real(double),allocatable:: KDE_Gaussian_Covariance(:,:), Data_Vectors(:,:)
    real(double):: KDE_Gaussian_Covariance_Reduction = 0.01e0_double

    if(KDE_Smooth == .false.) STOP "produce_Magnitude_Distribution - I haven't coded up a non-KDE version of this yet... Stopping"

    if(allocated(MagGrid) == .false.) then
       Higher =maxval(RefCat%MF606W)+ 2.17e0_double*Mag_Limit_Convergence_Buffer; Lower = minval(RefCat%MF606W)
       
       if(KDE_Smooth) then
          !--KDE Declarations--!
          allocate(MagGrid(nSmoothed_Sampling_Mag)); MagGrid = 0.e0_double
          do i =1, nSmoothed_Sampling_Mag
             MagGrid(i) = Lower + (i-1)*((Higher-Lower)/(nSmoothed_Sampling_Mag-1))
          end do
       else
          !--Histogram Decalrations
          allocate(Bins(nHist+1)); Bins = 0.e0_double
          allocate(MagGrid(nHist)); MagGrid = 0.e0_double
          Bins(1) = Lower
          do i =2, nHist + 1
             Bins(i) = Bins(i-1) + ((Higher-Lower)/(nHist))
             MagGrid(i-1) = 0.5e0_double*(Bins(i-1) + Bins(i))
          end do
       end if
    end if
       
    print *, 'Producing Magnitude Distribution:'
    print *, 'Catalogue Limits:', minval(RefCat%MF606W), maxval(RefCat%MF606W)
    print *, 'Grid Limits:', minval(MagGrid), maxval(MagGrid)

    if(KDE_Smooth) then
       allocate(Data_Vectors(2,size(RefCat%Sizes))); Data_Vectors(1,:) = RefCat%MF606W; Data_Vectors(2,:) = RefCat%Sizes
       call Discrete_Covariance(Data_Vectors, KDE_Gaussian_Covariance)
       KDE_Gaussian_Covariance = KDE_Gaussian_Covariance_Reduction*KDE_Gaussian_Covariance
       
       allocate(PDF(size(MagGrid))); PDF = 0.e0_double
       call KDE_Univariate_Gaussian(RefCat%MF606W, dsqrt(KDE_Gaussian_Covariance(1,1)), MagGrid, PDF)
       
       deallocate(KDE_Gaussian_Covariance, Data_Vectors)
    else
       allocate(PDF(size(MagGrid))); PDF = 0.e0_double
       do c = 1, size(RefCat%MF606W)
          do i = 1, nHist
             if( (RefCat%MF606W(c) >= Bins(i)) .and. (RefCat%MF606W(c) < Bins(i+1)) ) then
                PDF(i) = PDF(i) + 1.e0_double
                exit
             end if
          end do
       end do
       if(sum(PDF) /= size(RefCat%MF606W)) STOP 'produce_Magnitude_Distribution - FATAL ERROR - Not all galaxies accounted for in histogram'

       deallocate(Bins)
    end if

    !--Renormalise Distribution
    PDF = PDF/Integrate(MagGrid, PDF, 2, lim = (/minval(MagGrid), maxval(MagGrid)/))

  end subroutine produce_Magnitude_Distribution


  !----------------SIZE DISTRIBUTIONS-----------------------------------------------------------!
  subroutine produce_Joint_Size_Magnitude_Distribution(SizeGrid, MagGrid, PDF, RefCat, use_Physical_sizes, Magnitude_Type, Output_Dir, ln_size_Distribution, KDE_Smooth, SizeLimits, MagLimits)
    use Catalogues; use Statistics, only: variance_discrete, Discrete_Covariance; use Smoothing, only:KDE_Bivariate_Gaussian; use Integration, only: Integrate, TrapInt, RectangularIntegration
    !--Essentially Wrapper routine - Does the same as get_Size_Distribution_MagnitudeBinning_byCatalogue, but binning type is set here, and output is slightly varied--!
    !--Returns p(R,m)
    !~~~~ Q: How does the result vary with Magnitude Binnning and Size Binning? I.e. if a magnitude bin has few galaxies, then the distribution will be noisy when taking 60 size bins!
    !~~~To Do: 
    !~~ Normalise across magnitudes
    !~~ Smoothing - Inverse Variance?
    real(double), intent(inout),allocatable,dimension(:)::SizeGrid, MagGrid
    real(double), intent(out),allocatable::pdf(:,:) !-Magnitude, Size-!
    type(Catalogue),intent(in)::RefCat
    logical,intent(in)::use_Physical_Sizes
    integer, intent(in)::Magnitude_Type
    character(*), intent(in), optional:: Output_Dir
    logical,intent(in),optional::ln_size_Distribution, KDE_Smooth
    real(double), intent(in), optional:: SizeLimits(2), MagLimits(2)
    
    !----Binning Decalarations-----!
    integer,parameter::nSizes = 60, nMags = 60
    real(double):: Lower, Higher, dParam
    real(double):: Mag_Limit_Convergence_Buffer = 0.3e0_double
    real(double),allocatable,dimension(:,:):: SizeBins, MagBins

    integer:: I,J
    character(20)::fmtstring

    character(500):: iOutput_Dir
    logical:: iln_size_Distribution

    real(double)::Renormalisation

    !--Smoothed Version--!
    real(double),allocatable::Smoothed_Grid_Size(:), Smoothed_Grid_Mag(:), Smoothed_PDF(:,:)
    real(double):: Sig = 0.5e0_double !-Want a different width for magnitude and size?
    integer:: nSmoothed_Sampling = 100, nSmoothed_Sampling_Mag = 110
    real(double),allocatable::TCatMags(:), TCatSizes(:) !-Temprary storage for Magnitude and sizes-!
    real(double),allocatable:: Smoothed_SizeOnly_PDF(:), Smoothed_MagOnly_PDF(:) !--Used for Testing-!
    real(double),allocatable:: KDE_Gaussian_Covariance(:,:), Data_Vectors(:,:)
    real(double):: KDE_Gaussian_Covariance_Reduction = 0.01e0_double !-How much is sig^2 which give KDE width reduced from measured covariance?--!
    

    real(double),allocatable:: Histogram_PDF(:,:)
    real(double),allocatable:: Histogram_SizeOnly_PDF(:), Histogram_MagOnly_PDF(:)
    real(double),allocatable:: Histogram_MagGrid(:), Histogram_SizeGrid(:)
    
    logical:: here

    INTERFACE
       subroutine produce_Joint_Size_Magnitude_Distribution(SizeGrid, MagGrid, PDF, RefCat, use_Physical_sizes, Magnitude_Type, Output_Dir, ln_size_Distribution, KDE_Smooth, SizeLimits, MagLimits)
         use Param_Types         
         use Catalogues
         real(double), intent(inout),allocatable,dimension(:)::SizeGrid, MagGrid
         real(double), intent(out),allocatable::pdf(:,:) !-Magnitude, Size-!                                                                                                                                            
         type(Catalogue),intent(in)::RefCat
         logical,intent(in)::use_Physical_Sizes
         integer, intent(in)::Magnitude_Type

         character(*), intent(in), optional:: Output_Dir
         logical,intent(in),optional::ln_size_Distribution, KDE_Smooth
         real(double), intent(in), optional:: SizeLimits(2), MagLimits(2)
       END subroutine produce_Joint_Size_Magnitude_Distribution
    END INTERFACE

    print *, 'Inside Joint Size.. Distribution'

    if(present(ln_size_Distribution) .and. ln_size_Distribution .and. RefCat%log_sizes) STOP 'produce_Joint_Size_Magnitude_Distribution - log-Size specified and catalogue already in log-Size. This can be rectified by setting TCat%Sizes = dexp(RefCat%Sizes) in this routine, but has been disabled for safety'

    iOutput_Dir = 'Distributions/'
    if(present(Output_Dir)) iOutput_Dir = Output_Dir

    inquire(Directory = trim(iOutput_Dir), exist = here)
    if(here == .false.) call system('mkdir '//trim(iOutput_Dir))

    iln_size_Distribution = .false.
    if(present(ln_size_Distribution)) iln_size_Distribution = ln_size_Distribution


    print *, 'PRoducing joint size-magnitude distribution: Limits:', minval(RefCat%Sizes), maxval(RefCat%Sizes)

    allocate(TCatMags(size(RefCat%RA))); TCatMags = 0.e0_double
    allocate(TCatSizes(size(RefCat%RA))); TCatSizes = 0.e0_double

    print *, 'Getting Limits'
    if(present(SizeLimits)) then
       Lower = SizeLimits(1); Higher = SizeLimits(2)
    else
       !--Set up SizeGrid--!
       if(use_Physical_Sizes) then
          if(present(ln_size_Distribution) .and. ln_size_Distribution) then
             Higher = dlog(maxval(RefCat%Physical_Sizes)); Lower = dlog(minval(RefCat%Physical_Sizes))
             TCatSizes = dlog(RefCat%Physical_Sizes)
          else
             Higher = maxval(RefCat%Physical_Sizes); Lower = 0.e0_double
             TCatSizes = RefCat%Physical_Sizes
          end if
       else
          if(present(ln_size_Distribution) .and. ln_size_Distribution) then
             Higher = dlog(maxval(1.1e0_double*RefCat%Sizes)); Lower = dlog(minval(RefCat%Sizes/(1.e0_double+ 2.e0_double*Mag_Limit_Convergence_Buffer)))
             print *, 'lnSize Dist: Limits:', Lower, Higher
             TCatSizes = dlog(RefCat%Sizes)
          else
             !--Min 0 Allows user to pass in lnT even if ln_size_Distribution == F
             Higher = maxval(1.1e0_double*RefCat%Sizes); Lower = min(0.e0_double, minval(RefCat%Sizes)-0.5e0_double*dlog(1.e0_double+ 2.e0_double*Mag_Limit_Convergence_Buffer))
             print *, 'Size Dist: Limits:', Lower, Higher
             TCatSizes = RefCat%Sizes
          end if
       end if
    end if
    if(Higher > 100.e0_double) STOP 'Upper limit on size distribution too large, check catalogue'
    
    !--Produce KDE Smoothed Version--!
    if(allocated(SizeGrid)) then
       allocate(Smoothed_Grid_Size(size(SizeGrid))); Smoothed_Grid_Size = SizeGrid
    else
       allocate(Smoothed_Grid_Size(nSmoothed_Sampling)); Smoothed_Grid_Size = 0.e0_double
       do i =1, nSmoothed_Sampling
          Smoothed_Grid_Size(i) = Lower + (i-1)*((Higher-Lower)/(nSmoothed_Sampling-1))
       end do
    end if

!    allocate(Histogram_SizeGrid(nSizes)); Sizes = 0.e0_double
    allocate(SizeBins(nSizes,2)); SizeBins = 0.e0_double
    dParam = (( Higher- Lower )/(1.e0_double*(nSizes-1)) )
    if(dParam == 0.e0_double) STOP 'get_Size_Distribution_MagnitudeBinning_byCatalogue - Error setting up Size Grid, dSize = 0'
    do i = 1, nSizes
       !--Use i-2 so that no galaxies fall into the first bin--!                                                                                                                                             
       SizeBins(i,1) = Lower + (i-2)*dParam
       SizeBins(i,2) = SizeBins(i,1) + dParam
       if(i==nSizes) SizeBins(i,2) = SizeBins(i,2) + 1.e-3_double*dParam
    end do

        
    !--Set up Magnitude Binning--!
    !-Set Higher by + 2.17e0_double*Mag_Limit_Convergence_Buffer as will be going along a +2.17Kappa de-lensing line - Ensures at least kappa = 0.3 is acheivable for *every* galaxy in sample
    if(present(MagLimits)) then
       Lower = MagLimits(1); Higher = MagLimits(2)
       if(Magnitude_Type == 1) then !-Absolute_Magnitude-!
          TCatMags = RefCat%Absolute_Magnitude
       elseif(Magnitude_Type == 2) then
          TCatMags = RefCat%MF606W
       else
          STOP 'get_Size_Distribution_MagnitudeBinning_byCatalogue - Invalid Magnitude Type entered'
       end if

    else
       if(Magnitude_Type == 1) then !-Absolute_Magnitude-!
          !       call Calculate_Bin_Limits_by_equalNumber(Cat%Absolute_Magnitude, nMags, MagBins)
          Higher = maxval(RefCat%Absolute_Magnitude)+ 2.17e0_double*Mag_Limit_Convergence_Buffer; Lower = minval(RefCat%Absolute_Magnitude)
          TCatMags = RefCat%Absolute_Magnitude
       elseif(Magnitude_Type == 2) then
          !       call Calculate_Bin_Limits_by_equalNumber(Cat%MF606W, nMags, MagBins)
          Higher =maxval(RefCat%MF606W)+ 2.17e0_double*Mag_Limit_Convergence_Buffer; Lower = minval(RefCat%MF606W)
          TCatMags = RefCat%MF606W
       else
          STOP 'get_Size_Distribution_MagnitudeBinning_byCatalogue - Invalid Magnitude Type entered'
       end if
    end if

    !--Set up magnitude histogram binning
    if(allocated(MagGrid)) then
       !--Is Grid passed in, then set this as priority
       allocate(Histogram_MagGrid(size(MagGrid))); Histogram_MagGrid = MagGrid
       !--Convert to mag binning
       allocate(MagBins(size(Histogram_MagGrid),2)); MagBins = 0.0e0_double
       MagBins(1,:) = (/Histogram_MagGrid(1) - 0.5e0_double*(Histogram_MagGrid(2)-Histogram_MagGrid(1)), Histogram_MagGrid(2) + 0.5e0_double*(Histogram_MagGrid(2)-Histogram_MagGrid(1))/)
       do i = 2, size(MagBins,1)-1
          MagBins(i,:) = (/MagBins(i-1,2),  Histogram_MagGrid(i) + 0.5e0_double*(Histogram_MagGrid(i+1)-Histogram_MagGrid(i))/)
       end do
       MagBins(size(MagBins,1),:) = (/MagBins(size(MagBins,1)-1,2), Histogram_MagGrid(size(MagBins,1)) + 0.5e0_double*(Histogram_MagGrid(size(MagBins,1))-Histogram_MagGrid(size(MagBins,1)-1))/)
    else
       dParam = (( Higher- Lower )/(1.e0_double*(nMags-1)) )
       allocate(MagBins(nMags,2)); MagBins = 0.e0_double
       allocate(Histogram_MagGrid(nMags)); Histogram_MagGrid = 0.e0_double
       do i = 1, nMags
          !--Use i-2 so that no galaxies fall into the first bin--!
          !       MagBins(i,1) = Lower + (i-1)*dParam
          MagBins(i,1) = Lower + (i-2)*dParam
          MagBins(i,2) = MagBins(i,1) + dParam
          if(i==nSizes) MagBins(i,2) = MagBins(i,2) + 1.e-3_double*dParam
          Histogram_MagGrid(i) = 0.5e0_double*(sum(MagBins(i,:)))
       end do
    end if
    

    !--KDE Declarations--!
    if(allocated(MagGrid)) then
       allocate(Smoothed_Grid_Mag(size(MagGrid))); Smoothed_Grid_Mag = MagGrid
    else
       allocate(Smoothed_Grid_Mag(nSmoothed_Sampling_Mag)); Smoothed_Grid_Mag = 0.e0_double
       do i =1, nSmoothed_Sampling_Mag
          Smoothed_Grid_Mag(i) = Lower + (i-1)*((Higher-Lower)/(nSmoothed_Sampling_Mag-1))
       end do
    end if

    print *, 'Magnitude Limits (BinMin, BinMax, SmoothedGrid Min/Max, Catalogue Min/Max):', MagBins(1,1), MagBins(size(MagBIns,1), 2), Smoothed_Grid_Mag(1), maxval(Smoothed_Grid_Mag), minval(RefCat%MF606W), maxval(RefCat%MF606W)

    !--Get KDE Smoothed Version--!
    if(present(KDE_Smooth)) then
       if(KDE_Smooth) then
          allocate(Data_Vectors(2,size(TCatMags))); Data_Vectors(1,:) = TCatMags; Data_Vectors(2,:) = TCatSizes

          if(any(isNaN(TCatMags))) STOP 'NaNs found in Magnitude'
          if(any(isNaN(TCatSizes))) STOP 'NaNs found in Sizes'

          call Discrete_Covariance(Data_Vectors, KDE_Gaussian_Covariance)
          KDE_Gaussian_Covariance = KDE_Gaussian_Covariance_Reduction*KDE_Gaussian_Covariance
          deallocate(Data_Vectors)

          allocate(Smoothed_PDF(size(Smoothed_Grid_Mag), size(Smoothed_Grid_Size))); Smoothed_PDF = 0.e0_double
          call KDE_Bivariate_Gaussian(TCatMags, TCatSizes, Smoothed_Grid_Mag, Smoothed_Grid_Size, Smoothed_PDF, Covariance = KDE_Gaussian_Covariance)
          deallocate(KDE_Gaussian_Covariance)

          !--Renormalise--!
          allocate(Smoothed_sizeOnly_PDF(size(Smoothed_Grid_Size))); Smoothed_sizeOnly_PDF = 0.e0_double
          do i =1, size(Smoothed_sizeOnly_PDF)
             Smoothed_sizeOnly_PDF(i) = TrapInt(Smoothed_Grid_Mag, Smoothed_PDF(:,i))
          end do
          Renormalisation = TrapInt(Smoothed_Grid_Size, Smoothed_SizeOnly_PDF)
          if(Renormalisation == 0.e0_double) STOP 'Joint_Size_Magnitude_Distribution - Renormalisation of smoothed PDF is zero, stopping...'
          Smoothed_PDF = Smoothed_PDF/Renormalisation
          deallocate(Smoothed_sizeOnly_PDF)


          !--Output--!
          open(unit = 45, file  = trim(iOUtput_Dir)//'Smoothed_Size_Mag_Dist.dat')
          write(fmtstring, '(I5)') size(Smoothed_Grid_Mag) + 1
          write(45, '('//trim(fmtstring)//'(e14.7,x))') 0.0e0_double, Smoothed_Grid_Mag
          do i= 1, size(Smoothed_Grid_Size)
             write(45, '('//trim(fmtstring)//'(e14.7,x))') Smoothed_Grid_Size(i), Smoothed_PDF(:,i)
          end do
          close(45)
          print *, 'Output to: ', trim(iOUtput_Dir), ' Smoothed_Size_Mag_Dist.dat'
          
          
          allocate(Smoothed_sizeOnly_PDF(size(Smoothed_Grid_Size))); Smoothed_sizeOnly_PDF = 0.e0_double
          do i =1, size(Smoothed_sizeOnly_PDF)
             Smoothed_sizeOnly_PDF(i) = TrapInt(Smoothed_Grid_Mag, Smoothed_PDF(:,i))
          end do
          Smoothed_sizeOnly_PDF = Smoothed_sizeOnly_PDF/TrapInt(Smoothed_Grid_Size, Smoothed_sizeOnly_PDF)
          !--Output--!
          open(unit=4, file = trim(iOUtput_Dir)//'Smoothed_SizeDist_From_BivariateKDE.dat')
          do i = 1, size(Smoothed_sizeOnly_PDF)
             write(4, *) Smoothed_Grid_Size(i), Smoothed_sizeOnly_PDF(i)
          end do
          print *, 'Sucessfully output to: ', trim(iOUtput_Dir), ' Smoothed_SizeDist_From_BivariateKDE.dat'
          close(4)
          deallocate(Smoothed_sizeOnly_PDF)

          allocate(Smoothed_MagOnly_PDF(size(Smoothed_Grid_Mag))); Smoothed_MagOnly_PDF = 0.e0_double
          do i =1, size(Smoothed_MagOnly_PDF)
             Smoothed_MagOnly_PDF(i) = TrapInt(Smoothed_Grid_Size, Smoothed_PDF(i,:))
          end do
          Smoothed_MagOnly_PDF = Smoothed_MagOnly_PDF/TrapInt(Smoothed_Grid_Mag, Smoothed_MagOnly_PDF)
          !--Output--!
          open(unit=4, file = trim(iOUtput_Dir)//'Smoothed_MagDist_From_BivariateKDE.dat')
          do i = 1, size(Smoothed_MagOnly_PDF)
             write(4, *) Smoothed_Grid_Mag(i), Smoothed_MagOnly_PDF(i)
          end do
          print *, 'Sucessfully output to: ', trim(iOUtput_Dir), ' Smoothed_MagDist_From_BivariateKDE.dat'
          close(4)
          deallocate(Smoothed_MagOnly_PDF)
       end if
    end if


    !--Get the Histogram--!
    call get_Size_Distribution_MagnitudeBinning_byCatalogue(MagBins, Histogram_SizeGrid, Histogram_PDF, RefCat, use_Physical_sizes,  Magnitude_Type, ln_size_Distribution = iln_size_Distribution, SizeBins = SizeBins, Renormalise = .false., KDE_Smooth = .false.)

    !--Output--!
    allocate(Histogram_sizeOnly_PDF(size(Histogram_SizeGrid))); Histogram_sizeOnly_PDF = 0.e0_double
    do i =1, size(Histogram_sizeOnly_PDF)
       Histogram_sizeOnly_PDF(i) = RectangularIntegration(Histogram_MagGrid, HISTOGRAM_PDF(:,i))!sum(PDF(:,i))!TrapInt(Histogram_MagGrid, PDF(:,i))
    end do
    Histogram_sizeOnly_PDF = Histogram_sizeOnly_PDF/RectangularIntegration(Histogram_SizeGrid, Histogram_sizeOnly_PDF)
    !--Output--!
    open(unit=4, file = trim(iOUtput_Dir)//'Histogram_SizeDist.dat')
    do i = 1, size(Histogram_sizeOnly_PDF)
       write(4, *) Histogram_SizeGrid(i), Histogram_sizeOnly_PDF(i)
    end do
    print *, 'Sucessfully output to: ', trim(iOUtput_Dir),'Histogram_SizeDist.dat'
    close(4)
    deallocate(Histogram_sizeOnly_PDF)

    allocate(Histogram_MagOnly_PDF(size(Histogram_MagGrid))); Histogram_MagOnly_PDF = 0.e0_double       
    do i =1, size(Histogram_MagOnly_PDF)
       Histogram_MagOnly_PDF(i) = RectangularIntegration(Histogram_SizeGrid, Histogram_PDF(i,:))
    end do
    Histogram_MagOnly_PDF = Histogram_MagOnly_PDF/TrapInt(Histogram_MagGrid, Histogram_MagOnly_PDF)
    !--Output--!
    open(unit=4, file = trim(iOUtput_Dir)//'Histogram_MagDist.dat')
    do i = 1, size(Histogram_MagOnly_PDF)
       write(4, *) Histogram_MagGrid(i), Histogram_MagOnly_PDF(i)
    end do
    print *, 'Sucessfully output to: ',trim(iOutput_Dir),' Histogram_MagDist.dat'
    close(4)
    deallocate(Histogram_MagOnly_PDF)

    !--Renormalise across magnitudes??-!
    Renormalisation = 0.e0_double
    do i= 1, size(HISTOGRAM_PDF,1) !--Mags--!
       do j =1, size(HISTOGRAM_PDF,2) !-- Size--!
          Renormalisation = Renormalisation + (MagBins(i,2)-MagBins(i,1))*(SizeBins(j,2)-SizeBins(j,1))*HISTOGRAM_PDF(i,j)
       end do
    end do
    if(Renormalisation <= 0.e0_double) STOP 'get_Joint_Size_Magnitude_Distribution - Renormalisation is invalid, <=0. Stopping'
    HISTOGRAM_PDF = HISTOGRAM_PDF/Renormalisation

    !--Output--!
    if(iln_size_Distribution) then
       open(unit = 37, file = trim(iOutput_Dir)//'Joint_lnSize_Magnitude_Distribution_Histogram.dat')
    else
       open(unit = 37, file = trim(iOutput_Dir)//'Joint_Size_Magnitude_Distribution_Histogram.dat')
    end if
    !--Header--!
    write(37, '(A)', advance = 'no') '#'
    do i =1, size(SizeBins,1)-1
       write(37, '(e14.7,x)', advance = 'no') SizeBins(i,1)
    end do
    write(37, '(e14.7,x)') SizeBins(size(SizeBins,1),2)
    write(37, '(A)', advance = 'no') '#'
    do i =1, size(MagBins,1)-1
       write(37, '(e14.7,x)', advance = 'no') MagBins(i,1)
    end do
    write(37, '(e14.7,x)') MagBins(size(MagBins,1),2)
    write(37, *)    
    !--First line has Magnitude Binning Information--!
    write(fmtstring, '(I3)') size(MagGrid) + 1
    write(37, '('//trim(fmtstring)//'(e14.7,x))') 0.0e0_double, MagGrid
    do i= 1, size(Histogram_SizeGrid)
       write(37, '('//trim(fmtstring)//'(e14.7,x))') Histogram_SizeGrid(i), HISTOGRAM_PDF(:,i)
    end do
    deallocate(SizeBins, MagBins)
    close(37)

    if(allocated(SizeGrid)) deallocate(SizeGrid)
    if(allocated(MagGrid)) deallocate(MagGrid)
    !--Set return values--!Histogram_SizeGrid
    if(present(KDE_Smooth) .and. KDE_Smooth) then
       !--Set to Smoothed Version--!
       allocate(SizeGrid(size(Smoothed_Grid_Size))); SizeGrid = Smoothed_Grid_Size
       allocate(MagGrid(size(Smoothed_Grid_Mag))); MagGrid = Smoothed_Grid_Mag
       allocate(PDF(size(Smoothed_PDF,1), size(Smoothed_PDF,2))); PDF = Smoothed_PDF

       deallocate(Smoothed_PDF, Smoothed_Grid_Mag, Smoothed_Grid_Size)
    else
       !--Set to Histogram--!   
       allocate(SizeGrid(size(Histogram_SizeGrid))); SizeGrid = Histogram_SizeGrid
       allocate(MagGrid(size(Histogram_MagGrid))); MagGrid = Histogram_MagGrid
       allocate(PDF(size(Histogram_PDF,1), size(Histogram_PDF,2))); PDF = Histogram_PDF
    end if

    deallocate(Histogram_SizeGrid, Histogram_MagGrid, Histogram_PDF)
    deallocate(TCatMags, TCatSizes)


  end subroutine produce_Joint_Size_Magnitude_Distribution

  subroutine get_Size_Distribution_MagnitudeBinning_byCatalogue(MagBins, Sizes, PDF, RefCat, use_Physical_sizes,  Magnitude_Type, Output_Dir, ln_size_Distribution, SizeBins, Renormalise, KDE_Smooth)
    use Smoothing, only: KDE_Univariate_Gaussian
    !--Returns the size distribution according to the magnitude binning of MagBins, of Magnitude_Type (1:Abs, 2:App)
    !---If ln_size_Distribution present and true, then a ln size distribution is returned, in which case sizes is really ln(Sizes)
    !-- If magbins allocated, then distribution produced for the magnitude bins set out in magbins.
    !-- If SizeBin entered, then histrogram produced for the Size Binning used in SizeBins
    !-- Distribution in sizes assumed unless ln_Size_Distribution is present and true
    !-- Renormalisation is done by default, and is only not done if Renormalise is present and false. If Renormalised, returns p(R|m)
    !~~Edits:
    !~~~~~ Takes in a Size Bin and returns a SizeGrid - This is unneccessary - instead the Sizes part should be deleted and *only* a size bins and Mag bins returned.
    use Catalogues; use Statistics, only: variance_discrete
       real(double), intent(inout),allocatable:: Sizes(:)
       real(double), intent(out),allocatable::pdf(:,:) !-Mag Bin, GridValue-!
       real(double),dimension(:,:),allocatable,intent(inout)::MagBins
       character(*), intent(in), optional:: Output_Dir
       type(Catalogue),intent(in)::RefCat
       logical,intent(in)::use_Physical_Sizes
       integer, intent(in)::Magnitude_Type
       logical,intent(in),optional::ln_size_Distribution
       real(double), intent(in),optional::SizeBins(:,:)
       logical,intent(in),optional::Renormalise, KDE_Smooth

       character(500)::Catalogue_Filename
       integer::Catalogue_Columns(13)

       type(Catalogue)::Cat
       
       integer::nMags
       logical:: iRenormalise
       real(double)::Renormalisation
       
       !--Cuts--!  May not be implemented as want the result to resemble the catalogue as much as possible                                                                                                                                  
       integer::nSizes 
       logical::Apply_Size_Cuts = .false.
       real(double):: Size_Cuts(2) = (/0.,20./)
       real(double)::Size_lower, Size_Higher
       !    integer:: Catalogue_Size_Column = 14 !!!!!                                                                                                                                                                                               
       real(double)::dSize
       real(double),dimension(:,:),allocatable::iSizeBins
       
       integer::i,j,c
       
       character(5)::fmtstring
       
       type(Binned_Catalogue)::BCat
       real(double),allocatable::Temporary_Sizes_Array(:)

       !---Smoothed Distributions---!
       real(double),allocatable:: Smoothed_Grid(:)
       integer:: nSampling_Smooth = 500
       real(double):: Sig
       real(double),allocatable:: Smoothed_PDF(:,:)

       logical::here

       INTERFACE
          subroutine get_Size_Distribution_MagnitudeBinning_byCatalogue(MagBins, Sizes, PDF, RefCat, use_Physical_sizes, Magnitude_Type, Output_Dir, ln_size_Distribution, SizeBins, Renormalisation, KDE_Smooth)
            use Param_Types; use Catalogues
            real(double), intent(inout),allocatable:: Sizes(:)
            real(double), intent(out),allocatable::pdf(:,:) !-Mag Bin, GridValue-!      
            real(double),dimension(:,:),allocatable,intent(inout)::MagBins
            type(Catalogue),intent(in)::RefCat
            logical,intent(in)::use_Physical_Sizes
            integer, intent(in)::Magnitude_Type

            character(*), intent(in), optional:: Output_Dir
            logical,intent(in),optional::ln_size_Distribution
            real(double), intent(in),optional::SizeBins(:,:)
            logical, intent(in),optional:: Renormalisation, KDE_Smooth
          end subroutine get_Size_Distribution_MagnitudeBinning_byCatalogue
       END INTERFACE

!!$       Catalogue_Filename = trim(Reference_Catalogue)
!!$       Catalogue_Columns = Reference_Catalogue_Columns
!!$       
!!$       PRINT *, 'Reconstructing magnitude distribution from the Catalogue:', trim(Catalogue_Filename), ' to get size distributions'
!!$       !--Read in Catalogue--!                                                                                                                                                                                                                   
!!$       call Catalogue_ReadIn(Cat, trim(Catalogue_Filename), 'COMBO', Catalogue_Columns)
       print *, 'Reconstructing size-magnitude distribution from the reference catalogue'
       Cat = RefCat

       if(allocated(Sizes)) then
          nSizes = size(Sizes)
          !--Construct iSizeBins
          allocate(iSizeBins(size(Sizes),2)); iSizeBins = 0.0e0_double
          iSizeBins(1,:) = (/Sizes(1) - 0.5e0_double*(Sizes(2)-Sizes(1)), Sizes(2) + 0.5e0_double*(Sizes(2)-Sizes(1))/)
          do i = 2, size(iSizeBins,1)-1
             iSizeBins(i,:) = (/iSizeBins(i-1,2),  Sizes(i) + 0.5e0_double*(Sizes(i+1)-Sizes(i))/)
          end do
          iSizeBins(size(iSizeBins,1),:) = (/iSizeBins(size(iSizeBins,1)-1,2), Sizes(size(iSizeBins,1)) + 0.5e0_double*(Sizes(size(iSizeBins,1))-Sizes(size(iSizeBins,1)-1))/)
          
          Size_Higher = maxval(iSizeBins); Size_Lower = minval(iSizeBins)
       elseif(present(SizeBins)) then
          !--Assumed that SizeBins is set up to be the right range etc.. i.e. if using ln sizes then size bins should already be set up deal with this
          !--Testing? i.e shold only be negative if using lnSizes etc?
          nSizes = size(SizeBins,1)
          Size_Higher = maxval(SizeBins); Size_Lower = minval(SizeBins)
          allocate(Sizes(nSizes)); Sizes = 0.e0_double
          allocate(iSizeBins(size(SizeBins,1), size(SizeBins,2))); iSizeBins = SizeBins
          
          Sizes(1) = iSizeBins(1,2)
          do i = 2, nSizes
             Sizes(i) = 0.5e0_double*(iSizeBins(i,1)+iSizeBins(i,2))
          end do
       else
          nSizes = 120
          if(use_Physical_Sizes) then
             if(present(ln_size_Distribution) .and. ln_size_Distribution) then
                Size_Higher = dlog(maxval(Cat%Physical_Sizes)); Size_Lower = dlog(minval(Cat%Physical_Sizes))
             else
                Size_Higher = maxval(Cat%Physical_Sizes); Size_Lower = 0.e0_double!minval(Cat%Physical_Sizes)
             end if
          else
             if(present(ln_size_Distribution) .and. ln_size_Distribution) then
                Size_Higher = dlog(maxval(Cat%Sizes)); Size_Lower = dlog(minval(Cat%Sizes))
                print *, 'lnSize Dist: Limits:', Size_Lower, Size_Higher
             else
                Size_Higher = maxval(Cat%Sizes); Size_Lower = 0.e0_double!minval(Cat%Sizes);
                print *, 'Size Dist: Limits:', Size_Lower, Size_Higher
             end if
          end if
          
          allocate(Sizes(nSizes)); Sizes = 0.e0_double
          allocate(ISizeBins(nSizes,2)); ISizeBins = 0.e0_double
          dSize = (( Size_Higher- Size_Lower )/(1.e0_double*(nSizes-1)) )
          if(dSize == 0.e0_double) STOP 'get_Size_Distribution_MagnitudeBinning_byCatalogue - Error setting up Size Grid, dSize = 0'
          do i = 1, nSizes
             !--Use i-2 so that no galaxies fall into the first bin--!  
             iSizeBins(i,1) = Size_Lower + (i-2)*dSize
             iSizeBins(i,2) = iSizeBins(i,1) + dSize
             if(i==nSizes) iSizeBins(i,2) = iSizeBins(i,2) + 1.e-3_double*dSize
             Sizes(i) = 0.5e0_double*(iSizeBins(i,1)+iSizeBins(i,2))
          end do
          Sizes(1) = iSizeBins(1,2)       
       end if
       !--Bin by equal number density in magnitude bins--!                                                                                                                                                                
       if(allocated(MagBins) == .false.) then
          nMags = 5
          print *, 'Size Distribution, Getting Mag Limits', nMags
          if(all(Cat%Absolute_Magnitude >= 0.e0_double)) STOP 'size_Disribution by magnitude:, Absolute magnitude not set'
          if(Magnitude_Type == 1) then !-Absolute_Magnitude-!
             call Calculate_Bin_Limits_by_equalNumber(Cat%Absolute_Magnitude, nMags, MagBins)
          elseif(Magnitude_Type == 2) then
             call Calculate_Bin_Limits_by_equalNumber(Cat%MF606W, nMags, MagBins)
          else
             STOP 'get_Size_Distribution_MagnitudeBinning_byCatalogue - Invalid Magnitude Type entered'
          end if
       end if
       !--Bin Catalogue by Magnitude--!
       nMags = size(MagBins,1)
       call bin_catalogue_by_magnitude(Cat,MagBins,BCat, Magnitude_Type)

       print *, 'Contructing from total of:', size(Cat%RA), ' galaxies, binned according to:', BCat%Occupation

       print *, 'In producing Histogram of Size-Magnitude Distribution:'
       do i =1 , nMags
          if(present(ln_size_Distribution) .and. ln_size_Distribution) then
!             if(any(BCat%Cat(i)%Physical_Sizes <= 0.e0_double)) STOP 'get_Size_Distribution_MagnitudeBinning_byCatalogue - Physical Size contains zeros or negatives, stopping'
             print *, 'Mag Bin ', i, ' : ', MagBins(i,:)!, ' Variance in log Physical Sizes:', dsqrt(variance_discrete(dlog(BCat%Cat(i)%Physical_Sizes), dlog(BCat%Cat(i)%Physical_Sizes)))
          else
             print *, 'Mag Bin ', i, ' : ', MagBins(i,:)!, ' Variance in Physical Sizes:', dsqrt(variance_discrete(BCat%Cat(i)%Physical_Sizes, BCat%Cat(i)%Physical_Sizes))
          end if
       end do

       !--Set up Smoothing Result--!
       allocate(Smoothed_Grid(nSampling_Smooth))
       do i = 1, nSampling_Smooth
          Smoothed_Grid(i) = minval(Sizes) + (i-1)*( (maxval(Sizes)-minval(Sizes))/(nSampling_Smooth -1))
       end do
       !--This is editable. 0.5 Worrks well for STAGES--!
       Sig = 0.5e0_double
       allocate(Smoothed_PDF(nMags,size(Smoothed_Grid))); Smoothed_PDF = 0.e0_double
       !--------------------------!

       allocate(PDF(nMags, nSizes)); PDF = 0.e0_double
           
!!$       if(Cat%log_sizes == .false.) then
!!$          if(any(Cat%Sizes <= 0.e0_double)) then
!!$             print *, 'ERROR - Some of the reference catalogue sizes are negative or zero:', count(Cat%Sizes<= 0.e0_double)
!!$             STOP
!!$          END if
!!$          if(use_Physical_Sizes .and. any(Cat%Physical_Sizes <= 0.e0_double)) then
!!$             print *, 'ERROR - Some of the reference catalogue physical sizes are negative or zero:', count(Cat%Physical_Sizes<= 0.e0_double)
!!$             STOP
!!$          END if
!!$       end if
       
       
       do i = 1, nMags
          allocate(Temporary_Sizes_Array(size(BCat%Cat(i)%Sizes)));
          if(use_Physical_Sizes) then
             Temporary_Sizes_Array = BCat%Cat(i)%Physical_Sizes
          else
             Temporary_Sizes_Array= BCat%Cat(i)%Sizes
          end if
          
!!$          if(any(Temporary_Sizes_Array <= 0.e0_double)) then
!!$             print *, count(Temporary_Sizes_Array <= 0.e0_double)
!!$             STOP 'get_Size_Distribution_MagnitudeBinning_byCatalogue - Negative Sizes in Temporary array, stopping'
!!$          end if
          if(present(ln_size_Distribution) .and. ln_size_Distribution) then
             print *, 'produce_joint_size.....Producing log-Sizes Array, and log-Size distriubtion!'
             Temporary_Sizes_Array = dlog(Temporary_Sizes_Array)
          end if
          
          !-Get KDE Smoothed PDF--!
          if(present(KDE_Smooth)) then
             if(KDE_Smooth) then
                call KDE_Univariate_Gaussian(Temporary_Sizes_Array, Sig, Smoothed_Grid, Smoothed_PDF(i,:))
             end if
          end if
          
          do c = 1, size(Temporary_Sizes_Array)
             do j = 1, nSizes
                if( (Temporary_Sizes_Array(c) > ISizeBins(j,1)) .and. (Temporary_Sizes_Array(c) <= ISizeBins(j,2)) ) then
                   PDF(i,j) = PDF(i,j) + 1
                   exit
                end if
             end do
          end do
          deallocate(Temporary_Sizes_Array)
       end do
       
       !--Renormalise for each Magnitude Bin--!
       if(present(Renormalise)) then
          iRenormalise = Renormalise
       else
          iRenormalise = .True.
       end if
       if(iRenormalise) then
          do i =1, nMags
             !--Renormalise Histogram--!
             Renormalisation = 0.e0_double
             do j = 1, size(PDF,2)
                Renormalisation = Renormalisation + PDF(i,j)*(ISizeBins(i,2)-ISizeBins(i,1))
             end do
             if(Renormalisation> 0.e0_double) PDF(i,:) = PDF(i,:)/Renormalisation
             !--Renormalise KDE Smoothed PDF--!
             if(present(KDE_Smooth)) then
                if(KDE_Smooth) then
                   Renormalisation = 0.e0_double
                   do j = 1, size(Smoothed_PDF,2)-1
                      Renormalisation = Renormalisation + Smoothed_PDF(i,j)*(Smoothed_Grid(i+1)-Smoothed_Grid(i))
                   end do
                   if(Renormalisation> 0.e0_double) Smoothed_PDF(i,:) = Smoothed_PDF(i,:)/Renormalisation
                end if
             end if
          end do
       end if
       
       !--Output Size PDFS--!                                                                                                                                                                                        
       if(present(output_dir)) then
          
          !--Check for existence--!
          inquire(directory =  trim(Output_Dir), exist = here)
          if(here == .false.) call system('mkdir '//trim(Output_Dir))
          
          if(present(ln_size_Distribution) .and. ln_size_Distribution) then
             open(unit = 49, file = trim(Output_Dir)//'lnSize_Distribution_MagnitudeBinning_Catalogue.dat')
             print *, 'Size Distribution, by Magnitude Bin, output to '//trim(Output_Dir)//'lnSize_Distribution_MagnitudeBinning_Catalogue.dat'
             
             !--Smoothed--!
             open(71, file = trim(Output_Dir)//'Smoothed_lnSize_Distribution_KDE.dat')
             print *, 'Smoothed Size Distribution output to: ', trim(Output_Dir)//'Smoothed_lnSize_Distribution_KDE.dat'
          else
             open(unit = 49, file = trim(Output_Dir)//'Size_Distribution_MagnitudeBinning_Catalogue.dat')
             print *, 'Size Distribution, by Magnitude Bin, output to '//trim(Output_Dir)//'Size_Distribution_MagnitudeBinning_Catalogue.dat'
             
             !--Smoothed--!
             open(71, file = trim(Output_Dir)//'Smoothed_Size_Distribution_KDE.dat')
             print *, 'Smoothed Size Distribution output to: ', trim(Output_Dir)//'Smoothed_Size_Distribution_KDE.dat'
          end if
       else
          if(present(ln_size_Distribution) .and. ln_size_Distribution) then
             open(unit = 49, file = 'Distributions/lnSize_Distribution_MagnitudeBinning_Catalogue.dat')
             print *, 'Size Distribution, by Magnitude Bin, output to Distributions/lnSize_Distribution_MagnitudeBinning_Catalogue.dat'
             
             !--Smoothed--!
             open(71, file = 'Distributions/Smoothed_lnSize_Distribution_KDE.dat')
             print *, 'Smoothed Size Distribution output to: ', 'Distributions/Smoothed_lnSize_Distribution_KDE.dat'
             
          else
             open(unit = 49, file = 'Distributions/Size_Distribution_MagnitudeBinning_Catalogue.dat')
             print *, 'Size Distribution, by Magnitude Bin, output to Distributions/Size_Distribution_MagnitudeBinning_Catalogue.dat'
             
             open(71, file = 'Distributions/Smoothed_Size_Distribution_KDE.dat')
             print *, 'Smoothed Size Distribution output to: ', 'Distributions/Smoothed_Size_Distribution_KDE.dat'
             
          end if
          
          
       end if
       !--Write Header--!
       do j = 1, size(MagBins,1)
          write(49, '(A1, 2(e14.7,x))') '#', MagBins(j,:)
       end do
       write(49,'(A)')
       write(fmtstring, '(I5)') size(PDF,1)+1
       do j = 1, size(PDF,2)
          write(49, '('//trim(fmtstring)//'(e14.7,x))') Sizes(j), PDF(:,j)
       end do
       close(49)
       !--Output Smoothed Version--!
       do j = 1, size(MagBins,1)
          write(71, '(A1, 2(e14.7,x))') '#', MagBins(j,:)
       end do
       write(fmtstring, '(I5)') size(PDF,1)+1
       do j = 1, size(Smoothed_PDF,2)
          write(71, '('//trim(fmtstring)//'(e14.7,x))') Smoothed_Grid(j), Smoothed_PDF(:,j)
       end do
       close(71)
       
       if(present(KDE_Smooth)) then
          if(KDE_Smooth) then
             deallocate(Sizes); allocate(Sizes(size(Smoothed_Grid))); Sizes = Smoothed_Grid
             deallocate(PDF); allocate(PDF(size(Smoothed_PDF,1),size(Smoothed_PDF,2))); PDF = Smoothed_PDF
          end if
       end if
       deallocate(Smoothed_Grid, Smoothed_PDF)
       
       deallocate(iSizeBins)
       
     end subroutine get_Size_Distribution_MagnitudeBinning_byCatalogue
     
     
   end module Distributions
