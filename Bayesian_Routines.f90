module  Bayesian_Routines
  !--Uses Bayesian Means to Calculate Posterior for the Dark Matter Profile variables--!
  !--Method must have redshifts for every galaxy (20Jan2014)--!
  use Param_Types; use Catalogues
  implicit none

  character(500):: Bayesian_Routines_Output_Directory
  
  real(double)::Default_Source_Redshift = 1.4e0_double
  logical:: Analyse_with_Physical_Sizes = .false. !#Default Value#

  logical::Combine_log_Posteriors = .true. !-If False, then combined by multiplication-!
  logical:: use_lnSize = .true.

  logical,private::Debug_Mode = .true.

  integer::Surface_Mass_Profile = 3 !-1:Flat, 2:SIS, 3:NFW-!

  !--Method: 1: Size-Only, 2: Size-Magnitude, 3:Magnitude Only, 4: SizeMag - MagOnly Combination (SizeLimits)--!
  integer:: Posterior_Method = 2

  logical:: Enforce_Weak_Lensing = .false.
  real(double),allocatable:: Core_Cut_Position(:,:), Core_Cut_Radius(:)

  logical:: use_KDE_Smoothed_Distributions = .true., KDE_onTheFly = .false., allow_KDE_Extrapolation = .false.
  logical:: Cuts_Renormalise_Likelihood = .true.
  !--Survey Limits are applied only to the source sample
  real(double),dimension(2):: Survey_Magnitude_Limits = (/23.e0_double, 27.5e0_double/), Survey_Size_Limits = (/0.e0_double, 100.e0_double/), Survey_SNR_Limits = (/0.e0_double, 10000.e0_double/)
  real(double), dimension(2):: magSample_Magnitude_Limits, magSample_Size_Limits
  real(double), dimension(:,:),allocatable:: S1_Magnitude_Limits_byCluster, S2_Magnitude_Limits_byCluster, S1_Size_Limits_byCluster, S2_Size_Limits_byCluster, Size_Limits_byCluster, Magnitude_Limits_byCluster
  !--Global Limits are applied to both the source sample and field sample
  real(double), dimension(2):: global_SNR_Limits(2) = (/0., 100000./)
  !--Prior Limits are applied to field only - Deprecated, not applied, but not deleted to enable use of old input files
  real(double),dimension(2):: Prior_Magnitude_Limits = (/23.e0_double, 27.5e0_double/), Prior_Size_Limits = (/0.0e0_double, 100.e0_double/) !-3.3  

  real(double),parameter:: Lower_Redshift_Cut = 0.21

  !--Parallisation Decalaration (To Be Passed In)
  integer:: nOMPThread = 1

  !----MCMC Declarations-----!
  integer:: nBurnin = 50
  integer:: nChains = 6, nChainOut = 6, nMinChain = 1050
  logical:: tune_MCMC = .false.
  logical:: allow_Automatic_Exit_MCMC = .false.
  logical:: output_Burnin = .true.
  integer:: fit_Parameter(3) = (/1,0,0/) !-r_200, concentration, position
  real(double), allocatable:: centroid_Prior_Width(:)
  real(double), dimension(4):: MCMC_Proposal_Width = (/0.02e0_double, 0.0e0_double, 0.0005e0_double, 0.0005e0_double/) !-r200, c, RA, Dec-!
  !----Proposal Width:
  !--- 4 Cluster: (/0.04e0_double, 0.0e0_double, 0.0015e0_double, 0.0015e0_double/)
  !--- 6 Cluster: (/0.02e0_double, 0.0e0_double, 0.0005e0_double, 0.0005e0_double/)

  !--Overload function for Combine Posteriors
!!$  interface Combine_Posteriors
!!$     module procedure Combine_Posteriors_Scalar, Combine_Posteriors_Vector
!!$  end interface Combine_Posteriors


contains

!!$  subroutine Combine_Posteriors_Scalar(GridValue, Posteriors,  Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
!!$    !--Combines the posterior on a single grid value (free parameter alpha), using a call to the "normal" (vector) combined posterior subroutine
!!$     real(double), intent(in):: GridValue,Posteriors(:) !-Galaxy-!  
!!$     real(double), intent(out):: Combined_Posterior
!!$     logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP
!!$
!!$     !--Internal Declarations
!!$     real(double),dimension(1):: tGrid, tCombinedPosterior
!!$     real(double), dimension(size(Posteriors),1):: tPosteriors
!!$
!!$     !--Set up internals
!!$     tGrid = GridValue
!!$     tPosteriors(:,1) = Posteriors
!!$
!!$     call Combine_Posteriors(tGrid, tPosteriors, Combine_by_ln, .false., Return_lnP, tCombinedPosterior)
!!$     Combined_Posterior = tCombinedPosterior(1)
!!$
!!$   end subroutine Combine_Posteriors_Scalar
!!$
!!$  subroutine Combine_Posteriors_Vector(PosteriorGrid, Posteriors, Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
!!$    !--Combines posteriors by looping over the first dimension--!
!!$    real(double), intent(in):: PosteriorGrid(:),Posteriors(:,:) !-Galaxy, Grid/Value-!
!!$    real(double), intent(out):: Combined_Posterior(:)
!!$    logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP
!!$
!!$    integer::c, j
!!$    real(double)::Renorm, Combination_Normalisation
!!$    logical::iDoRenormalise
!!$
!!$    integer::nPosteriorsSkipped
!!$
!!$    INTERFACE
!!$       subroutine Combine_Posteriors_Vector(PosteriorGrid, Posteriors, Combine_by_ln, Renormalise, Return_lnP, Combined_Posterior)
!!$         use Param_Types
!!$         real(double), intent(in):: PosteriorGrid(:), Posteriors(:,:)
!!$         real(double), intent(out):: Combined_Posterior(:)
!!$         logical, intent(in):: Combine_by_ln, Renormalise, Return_lnP
!!$       END subroutine Combine_Posteriors_Vector
!!$    END INTERFACE
!!$
!!$    !--Set renormalisation by input. If a single alpha value entered, then do not renormalise
!!$    iDoRenormalise = Renormalise
!!$    if(size(PosteriorGrid) == 1) iDoRenormalise = .false.
!!$
!!$    nPosteriorsSkipped = 0
!!$    if(Combine_by_ln == .false.) then
!!$       Combination_Normalisation = 1.e0_double/maxval(Posteriors)!or 1.e0_double/(0.5e0_double*maxval(Posteriors(c,:))) within loop
!!$       Combined_Posterior = 1.e0_double
!!$       do c = 1, size(Posteriors,1) !-Loop over galaxies-!   
!!$          !--Skip if zero (lnP not defined then) or NaN
!!$          if(all(Posteriors(c,:) == 0.e0_double) .or. all(isNaN(Posteriors(c,:)))) then
!!$             nPosteriorsSkipped = nPosteriorsSkipped + 1
!!$             cycle
!!$          end if
!!$
!!$          !--Skip when `renormalised' posteriors are invalid (zero/negative)
!!$          if(all(Posteriors(c,:)*Combination_Normalisation == 0.e0_double)) then
!!$             print *, 'Invalid Posterior for galaxy:', c, ' (==0) press [ENTER] to output an stop..'; READ(*,*)
!!$             print *, Posteriors(c,:)
!!$             STOP
!!$          end if
!!$          if(any(Posteriors(c,:)*Combination_Normalisation < 0.e0_double)) then
!!$             print *, 'Invalid Posterior for galaxy:', c, ' (<0) presS [ENTER] to output an stop..'; READ(*,*)
!!$             print *, Posteriors(c,:)
!!$             STOP
!!$          end if
!!$            
!!$          Combined_Posterior(:) = Combined_Posterior(:)*(Posteriors(c,:)*Combination_Normalisation)
!!$
!!$          if(all(Combined_Posterior == 0.e0_double)) then
!!$             print *, 'Invalid CPosterior for galaxy:', c, ' press (==0) [ENTER] to output an stop..'; READ(*,*)
!!$             print *, Combination_Normalisation
!!$             read(*,*)
!!$             print *, Combined_Posterior(:)
!!$             read(*,*)
!!$             print *, Posteriors(c,:)
!!$             STOP
!!$          end if
!!$          if(anY(Combined_Posterior < 0.e0_double)) then
!!$             print *, 'Invalid CPosterior for galaxy:', c, ' press (<0) [ENTER] to output an stop..'; READ(*,*)
!!$             print *, Posteriors(c,:)
!!$             STOP
!!$          end if
!!$          
!!$       end do
!!$    else
!!$       if(return_lnP) iDoRenormalise = .false.
!!$!       print *, 'Combining using logs'
!!$       !-Set lnP to a large negative value as default, equivalent to P ~ 0
!!$       Combined_Posterior = -100_double
!!$       Combination_Normalisation = size(Posteriors,1)
!!$       do c = 1, size(Posteriors,1) !-Loop over galaxies-!
!!$          !-Error Catching--!
!!$          if(all(Posteriors(c,:) == 1.e-75_double) .or. all(isNaN(Posteriors(c,:)))) then
!!$             nPosteriorsSkipped = nPosteriorsSkipped + 1
!!$             cycle
!!$          end if
!!$          !--Sum log posteriors--!
!!$!!!$          print *, 'CombinePosterior, loop:', c
!!$!!!$          print *, Combined_Posterior, dlog(Posteriors(c,:))
!!$
!!$          where(Posteriors(c,:) == 0.e0_double)
!!$             Combined_Posterior = Combined_Posterior - 100.e0_double
!!$          elsewhere
!!$             Combined_Posterior = Combined_Posterior + dlog(Posteriors(c,:)) + 1.e0_double
!!$          end where
!!$
!!$
!!$          if(any(isNAN(Combined_Posterior(:)))) then
!!$             print *, 'Any NaNs in Combined Posterior?, galaxy:',c, any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
!!$             STOP
!!$          end if
!!$       end do
!!$       !--Convert to PDF, not ln(PDF)--!
!!$       if(Return_lnP) then
!!$          return
!!$       else
!!$          if(size(Combined_Posterior)/= 1) then
!!$             Combined_Posterior = dexp(Combined_Posterior - maxval(Combined_Posterior))
!!$          else
!!$             Combined_Posterior = dexp(Combined_Posterior)
!!$          end if
!!$       end if
!!$    end if
!!$ 
!!$    
!!$
!!$    if(any(isNAN(Combined_Posterior(:)))) then
!!$       print *, 'Any NaNs in Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
!!$       print *, 'Stopping'
!!$       STOP
!!$    END if
!!$
!!$    if(nPosteriorsSkipped > 0) write(*,'(A,I4,A,I4,A)') '################### Combine_Posteriors - ', nPosteriorsSkipped, ' of ', size(Posteriors,1), ' posteriors were skipped as they we invlaid - NaNs or 0 ###########################'
!!$    if((1.e0_double*nPosteriorsSkipped)/size(Posteriors,1) > 0.1) STOP 'Combine_Posteriors - number of skipped posteriors too large, stopping!'
!!$
!!$    !--Renormalise--!
!!$    if(iDoRenormalise) then
!!$       Renorm = 0.e0_double
!!$       do j = 1, size(Combined_Posterior)-1
!!$          Renorm = Renorm + 0.5e0_double*(Combined_Posterior(j) + Combined_Posterior(j+1))*(PosteriorGrid(j+1)-PosteriorGrid(j))
!!$       end do
!!$       if(Renorm <= 0.e0_double) then
!!$          print *, 'Renormalisation:', Renorm
!!$          STOP 'Combine_Posteriors - Invalid Renormalisation for combined Posterior'
!!$       end if
!!$       Combined_Posterior(:) = Combined_Posterior(:)/Renorm
!!$    end if
!!$
!!$    if(any(isNAN(Combined_Posterior(:)))) then
!!$       print *, 'Any NaNs in Renormalised Combined Posterior?', any(isNAN(Combined_Posterior)), count(isNAN(Combined_Posterior) == .true.)
!!$       print *, 'Stopping'
!!$       STOP
!!$    END if
!!$
!!$  end subroutine Combine_Posteriors_Vector

  subroutine cut_Sample_aroundCluster(Cat, sampleFlag, Cluster_Pos, Ap_Rad)
    !-- Removes a full sample (Method = 1,2,3) around a particualr cluster
    !-- Designed to cut the faint sources without sizes around A901b
    type(Catalogue), intent(inout):: Cat
    integer:: sampleFlag
    real(double), intent(in):: Cluster_Pos(2)
    real(double), intent(in):: Ap_Rad

    type(Catalogue):: tCat

    real(double), allocatable:: sep(:)
    integer:: c, counter, nCut

    print *, '\n Sources with Posterior Method == ', sampleFlag
    print *, '  with ', Ap_Rad ,' of ', cluster_pos, ' will be cut from sample'

    allocate(sep(size(Cat%RA))); sep = dsqrt( (Cat%RA-Cluster_Pos(1))**2.e0_double +(Cat%Dec-Cluster_Pos(2))**2.e0_double )


    print *, Cat%Posterior_Method
    print *, count((sep <= Ap_Rad)), count(Cat%Posterior_Method == sampleFlag), count((sep <= Ap_Rad) .and. (Cat%Posterior_Method == sampleFlag))
    nCut = count((sep <= Ap_Rad) .and. (Cat%Posterior_Method == sampleFlag))

    tCat = Cat
   
    call Catalogue_Destruct(Cat)
    call Catalogue_Construct(Cat, size(tCat%RA)-nCut)

    counter = 0
    do c = 1, size(sep)
      if(((sep(c) <= Ap_Rad) .and. (tCat%Posterior_Method(c) == sampleFlag)) == .false.) then
         counter = counter + 1
         call Catalogue_Assign_byGalaxy_byCatalogue(Cat, counter, tCat, c)
      end if
   end do

   if(counter /= size(Cat%RA)) then
      print *, 'Fatal Error - cut_Sample_aroundCluster - Error setting Catalogue size:', counter, size(Cat%RA)
      STOP
   END if

   call Catalogue_Destruct(tCat)

   print *, '\n Finished'

  end subroutine cut_Sample_aroundCluster
    


  subroutine split_Sample(iCat, oCat, asPrior, S2_Magnitude_Limits, S2_Size_Limits)
    !---Splits the sample into a magnitude-only catalogue, and a size-magnitude catalogue, according to the entered size limits, and SNR limits - anything outside size and SNR limits is appended to the magnitude-only catalogue

    !---- if asPrior == True, then the retruned catalogue array (oCat) will have as it's first element the case *without* size or magnitude cuts. This is because for the size-magnitude analysis, the size-mag catalogue from which the prior is constructed must not be cut by size or magnitude. The second elemet will stil have cuts applied, as in the SizeMag+Mag case it is defined by these cuts (in size)

    type(Catalogue), intent(in):: iCat
    type(Catalogue), intent(out):: oCat(2)
    logical, intent(in), optional:: asPrior
    real(double), intent(out), dimension(2), optional:: S2_Magnitude_Limits, S2_Size_Limits

    type(Catalogue):: tCat, inputCat
    logical:: iasPrior

    INTERFACE
       subroutine split_Sample(iCat, oCat, asPrior, S2_Mag_Limits, S2_Size_Limits)
         use Catalogues
         type(Catalogue), intent(in):: iCat
         type(Catalogue), intent(out):: oCat(2)
         logical, intent(in), optional:: asPrior

         real(double), intent(in), dimension(2):: S2_Mag_Limits, S2_Size_Limits
       end subroutine split_Sample
    END INTERFACE


    write(*,'(A)') '_______________________________________________________SPLITTING CATALOGUE________________________________________________'

    if(present(asPrior)) then
       iasPrior = asPrior
    else
       iasPrior = .false.
    end if

    inputCat = iCat

    PRINT *, 'Input Catalogue Check:'
    print *, 'Sizes:', minval(inputCat%Sizes), maxval(inputCat%Sizes), count(isNaN(inputCat%Sizes)), count(inputCat%Sizes <= 0.)
    print *, 'Magnitude:', minval(inputCat%MF606W), maxval(inputCat%MF606W)


    if((Catalogue_Constructed(inputCat) == .false.) .or. (size(inputCat%RA) == 0)) STOP 'FATAL ERROR - split_Sample - Input Catalogue is not allocated, or zero'

    !-----------------------------GENERAL CUTS - APPLY TO PRIOR ALSO

    !---Cut negative magnitudes on MF606W -- This as completely removed, and do not enther either sample (could reasonably be a size-only sample but can't be bothered. They are automatically removed anyway due to bright cut!!
    print *, '--- \n General cuts - Applied to both Source and Field sample ---'

    print *, 'Cutting by magnitude (removing invalid magnitudes for all samples):'
    call Cut_By_Magnitude(inputCat, 0.e0_double, 100000.e0_double)


    !--Split by Size Limits
    print *, 'Cutting by Pixel Size'

    !--Remove NaNs from Size array in SM case (i.e. this should always take badly measured szes into account, even in prior determination)
    if(use_lnSize) then
       print *, 'Invalid Size Check:', count(isNaN(inputCat%Sizes)), count(inputCat%Sizes <= -1000.0_double)
       call Cut_By_PixelSize(inputCat, -10000.e0_double, 10000.e0_double, tCat)
    else
       print *, 'Invalid Size Check:', count(inputCat%Sizes <= 0.e0_double)
       call Cut_By_PixelSize(inputCat, 0.e0_double, 10000.e0_double, tCat)
    end if
    !--Place NaNs in oCat(2)
    call Concatonate_Catalogues(oCat(2), tCat); call Catalogue_Destruct(tCat)

    !--Split by SNR Limits
    call Cut_By_SNR(inputCat, Survey_SNR_Limits(1), Survey_SNR_limits(2), tCat)
    call Concatonate_Catalogues(oCat(2), tCat)
    call Catalogue_Destruct(tCat)

    print *, '--- \n End of General cuts - Applied to both Source and Field sample ---'
    
    !----------------------------END OF GENERAL CUTS----------------------------------------------------------!

    print *, ' ' 
    print *, 'Split By SNR (mag limits):'
    print *, Survey_SNR_Limits
    print *, 'oCat(1):', minval(inputCat%MF606W), maxval(inputCat%MF606W)
    print *, 'oCat(2):', minval(oCat(2)%MF606W), maxval(oCat(2)%MF606W)
    !-----------------------------END GENERAL CUTS 

    !--At this stage, inputCat contains all sizes, but is cut by SNR
    !-- oCat(2) contains those invalid sizes as well as those failing SNR cut

    !--

    if(iasPrior) then
       !--Do not allow cut by size when evaluating prior. This should be taken into account when integrating 2D distribtion to get magnitude distriubtion, but size information is necessary for correct renormalisation even in this case. - This si true for the case where a size cut is used, and therefore the magnitde distribution for sample 2 can be obtained my integrating over the full sample between size limits. Where sample 2 also contains invalid sizes, the prior sample must be chosen as the full sample will be.
       !-oCat(2) = inputCat !-- oCat(2) contains size and SNR cuts

       !--In this case, set oCat(1) before the cut, so that no size cuts are applied to the Size sample
       !--oCat(2) will have sources which fail the size cut applied, to ensure that when sample 2 is contructed from both sources failing size cut and invalid sizes, that the prior is constructed from a representative sample of both
       oCat(1) = inputCat       
       call Cut_By_PixelSize(inputCat, Survey_Size_Limits(1), Survey_Size_Limits(2), tCat)
       call Concatonate_Catalogues(oCat(2), tCat)
       call Catalogue_Destruct(tCat)
    else
       !-- In this case, set oCat(1) after the cut (i.e. size cut applied)
       call Cut_By_PixelSize(inputCat, Survey_Size_Limits(1), Survey_Size_Limits(2), tCat)
       print *, 'Pixel Size Cut, Mag Check:', minval(tCat%MF606W), maxval(tCat%MF606W), count(tCat%MF606W <= 5.), count(tCat%MF606W <= 0.1), count(tCat%MF606W <= 0.)
       call Concatonate_Catalogues(oCat(2), tCat)
       call Catalogue_Destruct(tCat)

       oCat(1) = inputCat
    end if

    print *, ' '
    print *, 'Split By Size:'
    print *, Survey_Size_Limits
    print *, 'oCat(1):', minval(inputCat%Sizes), maxval(inputCat%Sizes)
    print *, 'oCat(2):', minval(oCat(2)%Sizes), maxval(oCat(2)%Sizes)


    call Catalogue_Destruct(inputCat)

    
    !--Assign Galaxy Posterior Methods to each
    if(Posterior_Method == 1) then !-Size-Only
       oCat(1)%Posterior_Method = 1
       oCat(2)%Posterior_Method  = 0
    elseif(Posterior_Method == 2) then !-Size-Mag
       oCat(1)%Posterior_Method = 2
       oCat(2)%Posterior_Method  = 0
    elseif(Posterior_Method == 3) then !-Mag-Only
       !--Edit Mag-Only catalogue to contain every catalogue (this whole routine could really be by-passed in this case)
       !--Use the following to produce magnitude-only on the whole catalogue
       print *, 'Sample 2 uses the full catalogue, as considering mag-only analysis'
       oCat(2) = iCat

       oCat(1)%Posterior_Method = 0
       oCat(2)%Posterior_Method  = 3
    elseif(Posterior_Method == 4) then !- Size-Mag where appropriate, M-Only thereafter
       oCat(1)%Posterior_Method = 2
       oCat(2)%Posterior_Method  = 3
    end if

    print *, ' '
    print *, 'oCat(1) [Valid Sizes]:'
    print *, 'nGal:', size(oCat(1)%RA)
    print *, 'Mag:',  minval(oCat(1)%MF606W), maxval(oCat(1)%MF606W)
    print *, 'Size:',  minval(oCat(1)%Sizes), maxval(oCat(1)%Sizes)
    print *, 'Posterior Method:', count(oCat(1)%Posterior_Method == 0), count(oCat(1)%Posterior_Method == 1), count(oCat(1)%Posterior_Method == 2), count(oCat(1)%Posterior_Method == 3), count(oCat(1)%Posterior_Method == 4)

    print *, ' '
    write(*,'(A)') 'oCat(2) [Invalid Sizes, Valid Magnitudes]:'
    print *, 'nGal:', size(oCat(1)%RA)
    print *, 'Mag:',  minval(oCat(2)%MF606W), maxval(oCat(2)%MF606W)
    print *, 'Size:',  minval(oCat(2)%Sizes), maxval(oCat(2)%Sizes)
    print *, 'Invalid Sizes?:', count(isNaN(oCat(2)%Sizes)), count(oCat(2)%Sizes<=0.e0_double)
    print *, 'Posterior Method:', count(oCat(2)%Posterior_Method == 0), count(oCat(2)%Posterior_Method == 1), count(oCat(2)%Posterior_Method == 2), count(oCat(2)%Posterior_Method == 3), count(oCat(2)%Posterior_Method == 4)

    !--Set Survey Limits - this should only be done on source sample (i.e. when iasPrior == false)
    if(iasPrior == .false.) then
       if(present(S2_Magnitude_Limits)) S2_Magnitude_Limits = (/minval(oCat(2)%MF606W), maxval(oCat(2)%MF606W)/)
       if(present(S2_Size_Limits)) S2_Size_Limits = (/minval(oCat(2)%Sizes), maxval(oCat(2)%Sizes)/)

       print *, 'Sample 1 has limits set as:'
       print *, 'Size Limits:', Survey_Size_Limits
       print *, 'Magnitude Limits:', Survey_Magnitude_Limits
      
       print *, 'Sample 2 [oCat(2)] has limits set as:'
       print *, 'Magnitude: ', S2_Magnitude_Limits
       print *, 'Size: ', S2_Size_Limits
    end if

    print *, ' '
    write( *,'(A)') '__________________________________ Succesfully Split Catalogue _________________________________________________'
    print *, 'Size-Mag Cat has:', size(oCat(1)%RA), ' of ', size(iCat%RA), ' galaxies'
    print *, 'Mag-Only Cat has:', size(oCat(2)%RA), ' of ', size(iCat%RA), ' galaxies'
    write(*,'(A)') '__________________________________ Succesfully Split Catalogue _________________________________________________'
    print *, ' '



  end subroutine split_Sample

  subroutine recombine_Source_Sample(iCat, oCat)
    !--Un-does split sample. Removes galaxies with Posterior_Method == 0
    type(Catalogue), intent(in):: iCat(:)
    type(Catalogue), intent(out):: oCat

    type(Catalogue):: tCat

    integer:: nTotal, i, nPass

    !--Recombine Sample
    call Catalogue_Destruct(oCat)

    do i = 1, size(iCat)
       call Concatonate_Catalogues(tCat, iCat(i))
    end do

    print *, '__________________________________________Recombining Source Sample___________________________________________'
    print *, 'Recombining sample: Concatonated sample has:', size(tCat%RA), ' galaxies'
    print *, '     ', count(tCat%Posterior_Method == 0), ' will be removed (0 Post. Method)'

    !--Remove galaxies without an assigned posterior method
    call Catalogue_Construct(oCat, count(tCat%Posterior_Method /= 0))
    nPass = 0
    do i = 1, size(tCat%RA)
       if(tCat%Posterior_Method(i) /= 0) then
          nPass = nPass + 1
          call Catalogue_Assign_byGalaxy_byCatalogue(oCat, nPass, tCat, i)
       end if
    end do

    print *, 'Recombined Sample has:', size(oCat%RA), ' galaxies'
    print *, '___________________________________________Recombined Source Sample___________________________________________'

  end subroutine recombine_Source_Sample

  subroutine Posterior_Statistics(PosteriorGrid, Posterior, MeanVal, ModeVal, Error, AntiSymm_Error)
    use Statistics, only: mean, mode_distribution, variance_distribution, Antisymmetric_Variance_Distribution
    real(double),intent(in)::Posterior(:), PosteriorGrid(:)
    real(double), intent(out),optional::MeanVal, ModeVal, Error, AntiSymm_Error(2)

    INTERFACE
       subroutine Posterior_Statistics(PosteriorGrid, Posterior, MeanVal, ModeVal, Error, AntiSymm_Error)
         use Param_Types
         real(double),intent(in)::Posterior(:), PosteriorGrid(:)
         
         real(double), intent(out),optional::MeanVal, ModeVal, Error
         real(double), intent(out),optional:: AntiSymm_Error(2)
       end subroutine Posterior_Statistics
    END INTERFACE

    if(any(isNaN(PosteriorGrid)) .or. any(PosteriorGrid > huge(1.e0_double))) STOP 'Posterior_Statistics: NaNs or Inf in Posterior Grid'

    if(present(ModeVal)) ModeVal = mode_distribution(PosteriorGrid, Posterior)
    if(present(MeanVal)) MeanVal = mean(Posterior, PosteriorGrid)
    if(present(Error)) Error = dsqrt(variance_distribution(PosteriorGrid, Posterior))
    if(present(AntiSymm_Error)) then
       if(present(ModeVal)) then
          !--Take error about mode as a preference--!
          AntiSymm_Error = Antisymmetric_Variance_Distribution(PosteriorGrid, Posterior, ModeVal)
       elseif(present(Meanval)) then
          !--Both the following use the mean value, however re-use if already calculated--!
          AntiSymm_Error = Antisymmetric_Variance_Distribution(PosteriorGrid, Posterior, MeanVal)
       else
          AntiSymm_Error = Antisymmetric_Variance_Distribution(PosteriorGrid, Posterior)
       end if
    end if

  end subroutine Posterior_Statistics

  subroutine Maximise_Convergence_byShifts_inAperture(Cat, BFCat, Ap_Pos, Ap_Radius, Core_Radius)
    use Statistics, only: mean_discrete
    type(Catalogue), intent(in)::Cat, BFCat
    real(double),intent(in)::Ap_Pos(:), Ap_Radius
    real(double), intent(in), Optional:: Core_Radius
    
    type(Catalogue)::Ap_Cat
    integer:: i,j
    integer,parameter:: nSide = 2500
    real(double)::Start_Pos(2)
    real(double):: Try_Pos(2), dTheta
    real(double):: Global_Mean_Size, Aperture_Mean_Size, Global_Mean_Mag, Aperture_Mean_Mag
    real(double), dimension(nSide, nSide, 2):: Convergence

    integer:: Maximum_Position(2)

    !--Cuts on Prior--!

    Global_Mean_Size = mean_discrete(BFCat%Sizes)
    Global_Mean_Mag =  mean_discrete(BFCat%MF606W)

    Start_Pos = Ap_Pos - Ap_Radius

    dTheta = Ap_Radius/(nSide-1)

    do i = 1, nSide
       Try_Pos(1) = Start_Pos(1) + (i-1)*dTheta
       do j = 1, nSide
          Try_Pos(2) = Start_Pos(2) + (j-1)*dTheta

          call Identify_Galaxys_in_Circular_Aperture(Cat, Try_Pos, Ap_Radius, Ap_Cat)

          Aperture_Mean_Size = mean_discrete(Ap_Cat%Sizes)
          Aperture_Mean_Mag = mean_discrete(Ap_Cat%MF606W)

          call Catalogue_Destruct(Ap_Cat)

          Convergence(i,j,1) = (Aperture_Mean_Size/Global_Mean_Size) - 1.e0_double !-Size-!
          Convergence(i,j,2) = (Global_Mean_Mag-Aperture_Mean_Mag)/2.17e0_double !-Magnitude-!
       end do
    end do
    !--In below, if not dabs, then maximum positive convergence is returned--!
    print *, 'Start Pos for Aperture:', Ap_Pos
    !--Output Maximum Pos for Mag--!
    Maximum_Position = maxloc(Convergence(:,:,2))
    print *, '** Maximum convergence for magnitude at:', Start_Pos(1)+(Maximum_Position(1)-1)*dTheta, Start_Pos(2)+(Maximum_Position(2)-1)*dTheta
    print *, '   which has value:', Convergence(Maximum_Position(1), Maximum_Position(2), 2)

    !--Output Maximum Position for Size--!
    Maximum_Position = maxloc(Convergence(:,:,1))
    print *, '** Maximum convergence for magnitude at:', Start_Pos(1)+(Maximum_Position(1)-1)*dTheta, Start_Pos(2)+(Maximum_Position(2)-1)*dTheta
    print *, '   which has value:', Convergence(Maximum_Position(1), Maximum_Position(2), 1)

    !--Output Maximum Position for both--!

  end subroutine Maximise_Convergence_byShifts_inAperture



  subroutine Convert_Alpha_Posteriors_to_VirialMass(SMD_Profile, Posteriors, Mass_Posteriors, Lens_Redshift, Output_Label)
    !--Converts from posteriors on DM free parameter (alpha) to posteriors on virial mass using the conservation of probability
    use Mass_Profiles, only: Halo_Mass, virial_Radius_from_ProfileFreeParameter; use Integration, only: Integrate; use gridintervals, only: equalscale
    
    integer,intent(in):: SMD_Profile
    real(double), intent(in):: Posteriors(:,:) !-Posteriors(1,:) assumed to be grid of DM Free parameters, (2,:) the actual posterior
    real(double), intent(out):: Mass_Posteriors(:,:)
    real(double):: Lens_Redshift !-Used only in NFW-!
    character(*), intent(in), optional:: Output_Label

    integer:: i
    real(double):: Discardable(1)
    character(20)::fmt
    real(double):: Redshift

    !-Posteriors entered must be in the form of the output mass posterior: that is, Dimension 1 labels grid/cluster, so that element A of Dim. 1 is Grid if A==1, or Cluster A-1 if A>1

    if(size(Posteriors,1) < 2) STOP 'Convert_VirialRadius_Posteriors_to_VirialMass - Error in input posterior - no grid/value'
    if(size(Mass_Posteriors) /= size(Posteriors)) STOP 'Convert_VirialRadius_Posteriors_to_VirialMass - Error in size of mass posterior, not equal to Original posterior' !Can be deleted in the case of interpolation

    print *, 'Converting alpha posterior to Mass Posterior'

    Redshift = Lens_Redshift

    !--Assumes that all the clusters entered are at teh same redshift--!
    do i = 1, size(Posteriors,2)
       call Halo_Mass(SMD_Profile, Posteriors(1,i), (/0.e0_double/), Mass_Posteriors(1,i), Discardable, Redshift)
       
       !--p(M)dM = p(r)dr -> p(M) \propto p(r)/r^2 -- Proportionality just affects renormalisation
!       Mass_Posteriors(2:,i) = Posteriors(2:,i)/(virial_Radius_from_ProfileFreeParameter(SMD_Profile, Posteriors(1,i))**2.e0_double)
       Mass_Posteriors(2:,i) = Posteriors(2:,i)/(Mass_Posteriors(1,i)**(2.e0_double/3.e0_double))
    end do

    print *, '**Mass Posterior output on grid of 10^14 Msun/h'
    Mass_Posteriors(1,:) = Mass_Posteriors(1,:)/1.e14_double

    !--Renormalise
    do i = 2, size(Posteriors,1)
       Mass_Posteriors(i,:) = Mass_Posteriors(i,:)/Integrate(Mass_Posteriors(1,:), Mass_Posteriors(i,:), 2, lim = (/minval(Mass_Posteriors(1,:)), maxval(Mass_Posteriors(1,:))/))
    end do

    if(present(Output_Label)) then
       open(unit = 78, file = trim(Output_Label)//'VirialMass_Posterior.dat')
       write(fmt, *) size(Mass_Posteriors,1)
       do i = 1, size(Mass_Posteriors,2)
          write(78, '('//trim(fmt)//'(e9.3,x))') Mass_Posteriors(:,i)
       end do
       write(*,'(2A)') '* Output file to: ', trim(Output_Label)//'VirialMass_Posterior.dat'
       close(78)
    end if
    print *, '** Virial Mass Posterior output to: ', trim(Output_Label)//'VirialMass_Posterior.dat'

    print *, '---Finished Conversion'

  end subroutine Convert_Alpha_Posteriors_to_VirialMass


  function Select_Source_Sample(Cat, Ap_Pos, Ap_Radius, Group_Index, Mag_Limits, Size_Limits, Redshift_Limit)
    !--Issues:
    !--- 3rd Sept: As well as unresovled seg fault in foreground contamination (removed with addition of print statement...), in this case when size limit cuts were applied, the min size of the source sample was lower than the input sample. This cannot be the case, and is possibly worrying
    !--Ap Radius must be in degrees
    !--Redshift Limit is, for now, a LOWER limit only
    type(Catalogue), intent(in):: Cat
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)
    integer, intent(in), optional:: Group_Index(:)
    real(double), dimension(:,:), intent(in), optional:: Mag_Limits, Size_Limits
    real(double), intent(in), optional:: Redshift_Limit(:)
    !! DEPRECATED real(double), intent(in), optional:: Redshift_Limit(:,:)
    type(Catalogue):: Select_Source_Sample

    integer, allocatable:: iGroup_Index(:)

    integer:: C, i
    type(Catalogue):: tCat, tSource_Catalogue, Group_Cat

    INTERFACE
        function Select_Source_Sample(Cat, Ap_Pos, Ap_Radius, Group_Index, Mag_Limits, Size_Limits, Redshift_Limit)
          use Param_Types; use Catalogues, only: Catalogue
          !--Ap Radius must be in degrees                                                                                                                                                
          type(Catalogue), intent(in):: Cat
          real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:)

          integer, intent(in), optional:: Group_Index(:)
          real(double), intent(in), optional:: Redshift_Limit(:)
          real(double), dimension(:,:), intent(in), optional:: Mag_Limits, Size_Limits
          type(Catalogue):: Select_Source_Sample
        end function Select_Source_Sample
    END INTERFACE


    write(*,'(A)') '__________________________________SELECTING SAMPLE__________________________________________________'

    print *, 'Sample input has size and magnitude limits:', minval(Cat%Sizes), maxval(Cat%Sizes), minval(Cat%MF606W), maxval(Cat%MF606W)

    
    if(present(Redshift_Limit)) then
       print *, 'Lower Redshift_Limit of ', Redshift_Limit, ' will be applied to sample selection'
       if(size(Redshift_Limit) /= size(Ap_Pos,1)) STOP 'Select_Source_Sample - FATAL ERROR - Redshift Limit is not of correct size'
    end if

    if(present(Group_Index)) then
       allocate(iGroup_Index(size(Group_Index))); iGroup_Index = Group_Index
    else
       allocate(iGroup_Index(size(Ap_Pos,1))); iGroup_Index = (/(i,i=1,size(Ap_Pos,1))/)
    end if

    !--Identify the source sample as the combination of sources in each aperture
    tCat = Cat
    do C = 1, size(iGroup_Index)
       !--Get sources for that aperture
       call Identify_Galaxys_in_Circular_Aperture(tCat, Ap_Pos(iGroup_Index(C),:), Ap_Radius(iGroup_Index(C)), TSource_Catalogue)!, Core_Radius = Core_Cut_Radius(Group_Index(C))/60.e0_double)

       !--Set renormalisation group (used in this case as one per cluster)
       TSource_Catalogue%RenormalisationGroup = iGroup_Index(C)

       !--Cut those sources which are undesirable (Usually by redshift)
       if(present(Redshift_Limit)) then
          print *, 'SOURCE SELECTION - CUTTING by 0.2 for all clusters!!!: THIS NEEDS DEALT WITH'
          call Cut_By_PhotometricRedshift(TSource_Catalogue, 0.2e0_double, 10000.e0_double)!DEPRECATED , Redshift_Limit(iGroup_Index(C),2))
          !call Cut_By_PhotometricRedshift(TSource_Catalogue, Redshift_Limit(iGroup_Index(C)), 10000.e0_double)!DEPRECATED , Redshift_Limit(iGroup_Index(C),2))
       end if

       if(present(Mag_Limits)) then
          print *, 'Cutting by input cluster magnitude limits:', Mag_Limits(iGroup_Index(C),:)
          call Cut_by_Magnitude(tSource_Catalogue, Mag_Limits(iGroup_Index(C),1), Mag_Limits(iGroup_Index(C),2))
       end if


!!$ Removed as also removed NaN sizes, which should be present due to Mag-Only sample
       if(present(Size_Limits)) then
          print *, 'Select Sample - Size Limit selection around clusters disabled until edit to ensure invalid sizes remain' 
!!$          print *, 'Cutting by input cluster size limits:', Size_Limits(iGroup_Index(C),:), minval(tSource_Catalogue%Sizes), maxval(tSource_Catalogue%Sizes)
!!$          !call Cut_by_PixelSize(tSource_Catalogue, Size_Limits(iGroup_Index(C),1), Size_Limits(iGroup_Index(C),2))
!!$          call Cut_by_PixelSize(tSource_Catalogue, -1000.e0_double, Size_Limits(iGroup_Index(C),2))
       end if

       call Concatonate_Catalogues(Group_Cat, TSource_Catalogue) 
       call Catalogue_Destruct(TSource_Catalogue)
       !--Mask Source Galaxies in temporary catalogue for that aperture to ensure no double counting
       print *, '--__ Applying mask to temporary Catalogue as part of Sample Selection, to ensure no double counting...'
       call Mask_Circular_Aperture(tCat, Ap_Pos(iGroup_Index(C),:),  Ap_Radius(iGroup_Index(C)))          
    end do
    call Catalogue_Destruct(tCat); call Catalogue_Destruct(tSource_Catalogue)

    Select_Source_Sample = Group_Cat

    print *, 'Got source sample. Sample contains:', size(Group_Cat%RA), ' sources'
    print *, 'and has limits (Mag, Size):', minval(Group_Cat%MF606W), maxval(Group_Cat%MF606W), minval(Group_Cat%Sizes), maxval(Group_Cat%Sizes)
    write(*,'(A)') '_____________________________________________________________________________________________________'

    call Catalogue_Destruct(Group_Cat); deallocate(iGroup_Index)
    
  end function Select_Source_Sample

  !----------------------------------------------------------POSTERIOR PRODUCTION ROUTINES-------------------------------------------------------------------------------------------------------------!

  subroutine DM_Profile_Fitting_Simultaneous_MCMC(iCat, Ap_Pos, Ap_Radius, Lens_Redshift, Marginalised_Posteriors, Distribution_Directory, reproduce_Prior, Fit_Group, Output_Prefix, Blank_Field_Catalogue)
    use Bayesian_Posterior_Evaluation, only: lnLikelihood_Evaluation_atVirialRadius_perSourceSample, get_likelihood_evaluation_precursors, lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample
    use Foreground_Clusters, only: get_distance_between_Clusters
    use Interpolaters, only: Linear_Interp; use MCMC; use Statistics, only: create_Histogram; use Common_Functions, only: Pyth_distance_between_Points; use Matrix_Methods, only: Sort_Array; use Variable_Tests, only: real_Ordered_Array_Decreasing
    !-Routine that produces posteriors on dark matter profile free parameters, by simultaneously fitting to all clusters that belong to the same 'Fit_Group'.
    !--Fit_Group should be in assending order starting from 1, with no numerical gaps, and should hold a value for each aperture considered. The number of free parameters for each fit is therefore determined by the number of clusters considered in each group. 

    type(Catalogue), intent(in):: iCat(:)
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:), Lens_Redshift(:)
    real(double),intent(out),allocatable::Marginalised_Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! 
    character(*), intent(in):: Distribution_Directory
    logical, intent(in):: reproduce_Prior
    type(Catalogue), intent(in),optional::Blank_Field_Catalogue(:)
    integer:: Fit_Group(:)
    character(*),intent(in):: Output_Prefix

    type(Catalogue):: tSource_Catalogue, tCat, Cat
    type(Catalogue):: Group_Cat
    real(double), dimension(Size(Ap_Pos,1))::iAp_Radius
    integer, allocatable:: Group_Index(:)

    integer:: C,G,i,j, counter

    character(500):: Group_Output_Prefix

    real(double),allocatable:: Source_Positions(:,:)
    real(double), dimension(size(Ap_Pos,1), size(Ap_Pos,1)):: angularDistance_betweenClusters

    !--Result Output Decalrations
    integer::nAlpha_Out = 500
    real(double):: AlphaLimit_Out(2) = (/0.05e0_double, 3.e0_double/)

    character(500):: Combined_Chain_Output, Convergence_Test_Output, Acceptance_Rate_Output
    character(25):: fmt, conv_fmt, acc_fmt
    character(100):: Filename
    character(3),allocatable:: Posterior_Flag(:)

    !--Distribution Declarations
    real(double),allocatable:: Joint_Size_Magnitude_Distribution(:,:), Magnitude_Distribution(:), TwoD_MagPrior_Distribution(:,:)
    real(double),allocatable::SizeGrid(:), MagGrid(:)

    !--Temporary Allocations--!
    real(double), allocatable:: tSigma_Crit(:,:), tLens_Redshift(:)

    !--Likelihood Evaluation Precursor Declarations
    real(double),dimension(:,:),allocatable:: Survey_Renormalised_Prior
    real(double),dimension(:),allocatable:: Survey_Renormalised_MagPrior
    real(double),allocatable:: RedshiftGrid(:)
    real(double),allocatable::Sigma_Crit(:,:)    
    real(double),allocatable:: MagnificationGrid(:), Renormalisation_by_Magnification(:,:), MagOnly_Renormalisation_by_Magnification(:,:)

    !--MCMC routines
    type(Array_of_MCMCChains),allocatable::ChainArray(:)
    integer:: M, out, nParameter_per_Cluster
    real(double),allocatable:: ChainWidths(:)
    real(double),allocatable:: Parameter_Start_Limits(:,:)
    character(300),allocatable:: ChainsOutput_Filename(:)
    integer:: chainCount
    real(double),allocatable:: Likelihood(:,:)
    real(double),allocatable:: Acceptance_Rate(:)
    logical:: Chains_Converged
    real(double),allocatable:: tAlpha(:), tAp_Pos(:,:), centroid_Prior_Offset(:)
    logical:: Acceptance
    integer:: nConvergenceTestPoint = 10
    real(double),allocatable:: ConvergenceStatistic(:)
    integer:: nMaxChain = 1000000
    logical, dimension(:),allocatable:: ifree_Parameter !used as a handle to tell convergence routine which parameter to check

    !---Mass Grouping
    type MassGrouping
       integer, allocatable:: ClusterIndex(:)
    end type MassGrouping
    type(MassGrouping),allocatable:: MassOrderingGroup(:)
    integer, allocatable:: cluster_MassOrdering_Groups(:)
    real(double), allocatable:: MO_tr200(:)
    logical:: Enforce_Mass_Ordering = .false.

    !--Temporary Posterior construction
    real(double),allocatable:: tMarginalised_Posterior_Grid(:), tMarginalised_Posterior(:)
    real(double), allocatable:: Combined_Chain(:,:), Combined_Likelihood(:)

    !--Paralleisation Declarations
    real(double), dimension(:), allocatable:: tSizes, tMF606W, tRedshift
    integer, allocatable:: tPosterior_Method(:), tRenormalisationGroup(:)
    integer:: nOpenThreads, OMP_GET_NUM_THREADS

    !--Testing Declarations:
    real(double):: Time1, Time2
    real:: Time_ChainLink_Start, Time_ChainLink_End
    
    if(size(Ap_Radius)==1) then
       iAp_Radius = Ap_Radius(1)
    else
       iAp_Radius = Ap_Radius
    end if

    !--Recombine Source Sample
    call recombine_Source_Sample(iCat, Cat)

    !--Error Catching on Fit Group
    if(size(Fit_Group) /= size(Ap_Pos,1)) STOP 'DM_Profile_Fitting_Simultaneous - Fit_Group entered is not of the correct size.'
    if(minval(Fit_Group) /= 1) STOP 'DM_Profile_Fitting_Simultaneous - Fit_Group entered does not satisfy minimum value conditions (should be 1)'
    
    if(size(Lens_Redshift) /= size(Ap_Pos,1)) STOP 'DM_Profile_Fitting_Simultaneous - Lens Redshift Entered not of correct size:, one redshift per cluster should be supplied'

    !_____________________________________________________________________________Process of Posterior Evaluation_____________________________________________________________________________!
    !--Construct/Read in prior distribution
    !~~~Repeat for all groups in Fit_Group:
    !-Process: For each group identify source sample: Should be all galaxies in aperture radius entered. If apertures do not overlap, then output an error message but continue with evaluation
    !- Set up new 'propose' point on chain for all free parameters
    !- Evaluate Posterior on chain.
    !- Accept/Reject point on chain
    !- Output Joint Posterior for each group to file, taken as the histogram of each column in the chain
    !- Assign Marginalised Posterior output through linear interpolation
    !~~~~ Repeat for next group
    !_____________________________________________________________________________Process of Posterior Evaluation_____________________________________________________________________________!


    !--Get Prior Distribution
    !--If a catalogue for the blank field is passed in, then use this catalogue to get the intrinsic distributions--!
    if(reproduce_Prior) then
       if(present(Blank_Field_Catalogue) == .false.) STOP 'DM_Profile_Variable_Posteriors_CircularAperture - Blankf Field Catalogue must be entered to allow for the production of the prior on a grid'
       write(*,'(A)') 'Producing Distribution from Catalogue'
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, TwoD_MagPrior_Distribution, Distribution_Directory, .true., Blank_Field_Catalogue)
    else
       write(*,'(A)') 'Reading in distribution from:', Distribution_Directory
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, TwoD_MagPrior_Distribution, Distribution_Directory, .false., Blank_Field_Catalogue)
       print *, 'Success:', Distribution_Directory
    end if

    !--Get Precursors
    call get_Likelihood_Evaluation_Precursors(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid, Renormalisation_by_Magnification, MagOnly_Renormalisation_by_Magnification, SizeGrid, MagGrid, Magnitude_Distribution, TwoD_MagPrior_Distribution, Joint_Size_Magnitude_Distribution, S1_Magnitude_Limits_byCluster, S1_Size_Limits_byCluster, S2_Magnitude_limits_byCluster, S2_Size_Limits_byCluster, Lens_Redshift, Lower_Redshift_Cut, Output_Prefix, use_lnSize)

    deallocate(Joint_Size_Magnitude_Distribution, Magnitude_Distribution, TwoD_MagPrior_Distribution)

    !--Set up output grid, on which the marginalised posteriors will be output, as the linear interpolation of the constructed marginalised posterior below
    allocate(Marginalised_Posteriors(size(Ap_Pos,1), 2, nAlpha_Out));
    do i =1, nAlpha_Out
       Marginalised_Posteriors(:,1,i) = (/(AlphaLimit_Out(1) + (i-1)*((AlphaLimit_Out(2)-AlphaLimit_Out(1))/(nAlpha_Out-1)),j=1,size(Marginalised_Posteriors,1))/)
    end do

    !--Get the angular distance between clusters being considered
    angularDistance_betweenClusters = get_distance_between_Clusters(Ap_Pos)

    if(count(Fit_Group == 1) /= size(Fit_Group)) then
       write(*,'(A)') 'Entered fit group does not simulataneously fit all clusters. Are you sure? <Enter> to continue'
       read(*,*)
    end if

    do G = 1, 100 !-Assume no more than 100 groups will be present
       !--Find all clusters with that grouping
       !---Group Index contains the index of the clusters within that group, and can be used to access the correct Aperture Location and Radius. Size of Group_Index is the number of free parameters being used
       allocate(Group_Index(count(Fit_Group == G)));
       !--Exit when all groups exhausted
       if(size(Group_Index) == 0) then
          if(maxval(Fit_Group) > G) then
             deallocate(Group_Index)
             cycle
          else
             exit
          end if
       end if

       print *, ' '
       print *, 'Joint Fitting Cluster Group ', G, ' of ', maxval(Fit_Group),'....'
       
       !--Store the index of the clusters being fit, for easy identification in the following
       i = 0
       do C = 1, size(Fit_Group)
          if(Fit_Group(C) == G) then
             i = i + 1
             Group_Index(i) = C
          end if
       end do

       !-- Mag_Limits, Size_Limits, Redshift_Limit

       Group_Cat = Select_Source_Sample(Cat, Ap_Pos, iAp_Radius, Group_Index, Mag_Limits = Magnitude_Limits_byCluster, Size_Limits = Size_Limits_byCluster, Redshift_Limit =  1.2e0_double*Lens_Redshift)

!!$ !! DEPRECATED
!!$       !--Identify the source sample as the combination of sources in each aperture
!!$       tCat = Cat
!!$       do C = 1, size(Group_Index)
!!$          call Identify_Galaxys_in_Circular_Aperture(tCat, Ap_Pos(Group_Index(C),:), iAp_Radius(Group_Index(C)), TSource_Catalogue)!, Core_Radius = Core_Cut_Radius(Group_Index(C))/60.e0_double)
!!$          call Concatonate_Catalogues(Group_Cat, TSource_Catalogue) 
!!$          call Catalogue_Destruct(TSource_Catalogue)
!!$          !--Mask Source Galaxies for that aperture to ensure no double counting
!!$          print *, '--__ Applying mask to temporary Catalogue as part of MCMC, to ensure no double counting...'
!!$          call Mask_Circular_Aperture(tCat, Ap_Pos(Group_Index(C),:),  iAp_Radius(Group_Index(C)))          
!!$       end do
!!$       call Catalogue_Destruct(tCat)
!!$
!!$       print *, 'Source Catalogue contains:', size(Group_Cat%RA), ' of ', size(Cat%RA), 'galaxies'

       if(allocated(Source_Positions)) deallocate(Source_Positions)
       allocate(Source_Positions(size(Group_Cat%RA),2));
       Source_Positions(:,1) = Group_Cat%RA; Source_Positions(:,2) = Group_Cat%Dec

       !--Set up temporary arrays of Sigma_Critical (usually aperture position)
       allocate(tSigma_Crit(size(Group_Index),size(Sigma_Crit,2))); !tSigma_Crit(1,:) = Sigma_Crit(Group_index(1),:); tSigma_Crit(2,:) = Sigma_Crit(Group_index(2),:) 
       allocate(tLens_Redshift(size(Group_Index)));
       do C = 1, size(Group_index)
          tSigma_Crit(C,:) = Sigma_Crit(Group_Index(C),:)
          tLens_Redshift(C) = Lens_Redshift(Group_Index(C))
       end do

       !--Check for overlap between clusters

       !--Set up original chains (nChain in total, used to calculate convergence)
       allocate(ChainArray(nChains))
       allocate(Acceptance_Rate(size(ChainArray))); Acceptance_Rate = 0

       !--Set output filename
       allocate(ChainsOutput_Filename(size(ChainArray)))
       write(ChainsOutput_Filename(1), '(I3)') G
       Combined_Chain_Output = trim(adjustl(Output_Prefix))//'Group'//trim(adjustl(ChainsOutput_Filename(1)))//'_MCMC_CombinedChain.dat'
       Convergence_Test_Output = trim(adjustl(Output_Prefix))//'Group'//trim(adjustl(ChainsOutput_Filename(1)))//'_MCMC_ConvergenceTest_R.dat'
       Acceptance_Rate_Output = trim(adjustl(Output_Prefix))//'Group'//trim(adjustl(ChainsOutput_Filename(1)))//'_MCMC_AcceptanceRate.dat'
       ChainsOutput_Filename = trim(adjustl(Output_Prefix))//'Group'//trim(adjustl(ChainsOutput_Filename(1)))//'_MCMC_Chain_'
       do M = 1, size(ChainArray)
          write(ChainsOutput_Filename(M),'(A,I1,A)') trim(adjustl(ChainsOutput_Filename(M))),M,'.dat'
       end do


       !---Edit centroid_Prior_Width to avoid overalp between nearby clusters
       if(allocated(centroid_Prior_Width) == .false.) then
          STOP '\n DM_Profile_Fitting_Simultaneous_MCMC - FATAL ERROR - centroid_Prior_Width is not allocated. Exiting.'
       end if

!!$ !!---This section of code edits the centroid prior width entered for the cluster to avoid cluster overlap----!
!!$       do i = 1, size(Group_Index)
!!$          do j = 1, size(Group_Index)
!!$             if (i == j) cycle
!!$             if((angularDistance_betweenClusters(Group_Index(i), Group_Index(j)) < centroid_Prior_Width(Group_Index(i))) .or. (angularDistance_betweenClusters(Group_Index(i), Group_Index(j)) < centroid_Prior_Width(Group_Index(j)))) then !-If overlap then
!!$                !--Set priorWidth for one cluster to half the distance between clusters, only if that is smaller than the input
!!$                centroid_Prior_Width(Group_Index(i)) = min(centroid_Prior_Width(Group_Index(i)), 0.5e0_double*angularDistance_betweenClusters(Group_Index(i),Group_Index(j)))
!!$                !-- set other clusters priorWidth to the remaining distance between the clusters
!!$                centroid_Prior_Width(Group_Index(j)) = min(centroid_Prior_Width(Group_Index(j)),angularDistance_betweenClusters(Group_Index(i),Group_Index(j))-centroid_Prior_Width(Group_Index(i)))
!!$             end if
!!$          end do
!!$       end do
!!$ !!--------------------------------------------------------------------------------------------------------------!

       !!--- Set Size Ordering groups for the clusters ---!!
       if(enforce_Mass_Ordering) then
          !---Find Groups which are close together (within centroid prior limits)
          allocate(cluster_MassOrdering_Groups(size(Group_Index))); cluster_MassOrdering_Groups = 0.
          do i = 1, size(Group_Index)
             do j = i+1, size(Group_Index)
                if(i == j) cycle
                if((angularDistance_betweenClusters(Group_Index(i), Group_Index(j)) < (centroid_Prior_Width(Group_Index(i)) + centroid_Prior_Width(j)) )) then !-If overlap then
                   if(cluster_MassOrdering_Groups(i) == 0) then
                      !--Create New Group
                      cluster_MassOrdering_Groups(i) = maxval(cluster_MassOrdering_Groups) + 1
                   end if
                   !--Append to Group
                   cluster_MassOrdering_Groups(j) = cluster_MassOrdering_Groups(i)
                end if
             end do
          end do
          
          !--Store mass ordering group
          allocate(MassOrderingGroup(maxval(cluster_MassOrdering_Groups)))
          do i = 1, size(MassOrderingGroup)
             allocate(MassOrderingGroup(i)%ClusterIndex(count(cluster_MassOrdering_Groups == i))); MassOrderingGroup(i)%ClusterIndex = 0
             counter = 0
             do j = 1, size(cluster_MassOrdering_Groups)
                if(cluster_MassOrdering_Groups(j) == i) then
                   counter = counter + 1
                   MassOrderingGroup(i)%ClusterIndex(counter) = j
                end if
             end do

             print *, ' '
             print *, 'Mass Grouping:', i, ':', MassOrderingGroup(i)%ClusterIndex
             print *, ' '

          end do


          
          deallocate(cluster_MassOrdering_Groups)
       end if


       !-Set Parameter_Start_Limits, which also sets which parameters are being varied
       !--Chain is set up to include all parameters by default, however if they are not marginalised over then the proposal distribution has width zero in that parameters direction
       !--Chain Width specifies the typical width of the proposal distribution in that axis of parameter space
       !--ifree_Parameter labels whether all parameters in chain are free (true) or fixed (false). When fixed, chain width is set to zero, and further in code, starting point is set to the fixed value
       nParameter_per_Cluster = size(fit_Parameter)+1 !-Increment as position requires two parameters
       allocate(Parameter_Start_Limits(size(Group_Index)*nParameter_per_Cluster, 2)); Parameter_Start_Limits = 0.e0_double
       allocate(ChainWidths(size(Group_Index)*nParameter_per_Cluster)); ChainWidths = 0.e0_double
       allocate(ifree_Parameter(size(Group_Index)*nParameter_per_Cluster)); ifree_Parameter = .true.
       do C = 1, size(Group_Index)

          Parameter_Start_Limits((C-1)*nParameter_per_Cluster+1,:) = (/0.4e0_double, 1.5e0_double/) !-r200
          Parameter_Start_Limits((C-1)*nParameter_per_Cluster+2,:) = (/0.0e0_double, 0.0e0_double/) !-Concentration
          Parameter_Start_Limits((C-1)*nParameter_per_Cluster+3,:) = (/maxval((/Ap_Pos(Group_Index(C),1)-centroid_Prior_Width(Group_Index(C)), 148.7707e0_double/)), minval((/Ap_Pos(Group_Index(C),1)+centroid_Prior_Width(Group_Index(C)), 149.3524e0_double/))/) !-RA
          Parameter_Start_Limits((C-1)*nParameter_per_Cluster+4,:) = (/maxval((/Ap_Pos(Group_Index(C),2)-centroid_Prior_Width(Group_Index(C)), -10.291e0_double/)), minval((/Ap_Pos(Group_Index(C),2)+centroid_Prior_Width(Group_Index(C)), -9.748e0_double/))/) !-Dec

          !-0.1 works well for 6 Cluster (r200 only).
          ChainWidths((C-1)*nParameter_per_Cluster+1) = MCMC_Proposal_Width(1) !-r200
          ChainWidths((C-1)*nParameter_per_Cluster+2) = MCMC_Proposal_Width(2) !-Concentration (not coded up yet)
          ChainWidths((C-1)*nParameter_per_Cluster+3) = MCMC_Proposal_Width(3) !-RA
          ChainWidths((C-1)*nParameter_per_Cluster+4) = MCMC_Proposal_Width(4) !-Dec

          !--Set values when not free
          if(fit_Parameter(1)==0) then !--r200
             ifree_Parameter((C-1)*nParameter_per_Cluster+1) = .false.
          end if
          if(fit_Parameter(2)==0) then !--Concentration
             ifree_Parameter((C-1)*nParameter_per_Cluster+2) = .false.
          end if
          if(fit_Parameter(3)==0) then !--Position
             ifree_Parameter((C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4) = .false.
             ChainWidths((C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4) = 0.e0_double
          end if
       end do
       

       !--Set starting link for chain, by randomly placing within set parameter limits
       do M = 1,size(ChainArray)
          call MCMC_StartChain_Random(Parameter_Start_Limits, ChainArray(M)%Chain)
          !--Start chain postions seperately, since prior is defined on an aperture (Reset over those positions set in previous step)
          print *, 'WARNING, starting position set very small indeed to test acceptance rate, THIS NEEDS FIXED'
          if(fit_Parameter(3) == 1) then !--If fitting centroid position
             do C = 1, size(Group_Index)
                ChainArray(M)%Chain(1,(C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4) = Start_MCMC_Chain_TopHatAperture(Ap_Pos(Group_Index(C),:), minval( (/0.25e0_double,0.25e0_double*centroid_Prior_Width(Group_Index(C))/))) !-maximum width of 0.25 ensures that arbitrarily large  centroid_Prior can be used without placing apeture too far from center, prefactor on centroid prior notes that usually the entered position is well known
                !--Ensure that centroid is not positioned outside survey limits
                !-(Ignored for now, as probably shouldn't affec the result, however this would need relaxed if prior on survey RA/Dec limits imposed)
                if(Pyth_distance_between_Points(ChainArray(M)%Chain(1,(C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4), Ap_Pos(Group_Index(C),:)) > centroid_Prior_Width(Group_Index(C))) STOP 'Error in setting initial starting position for cluster'
             end do
          end if
       end do

       deallocate(Parameter_Start_Limits)


       !--Rearrange Starting r200 if mass ordering is enforced
       if(enforce_Mass_Ordering) then
          do i = 1, size(MassOrderingGroup)
             do M = 1,size(ChainArray)
                allocate(MO_tr200(size(MassOrderingGroup(i)%ClusterIndex))); 
                do j = 1, size(MO_tr200)
                   MO_tr200(j) = ChainArray(M)%Chain(1,(MassOrderingGroup(i)%ClusterIndex(j)-1)*nParameter_per_Cluster+1)
                end do
                   !--Sort Cluster Index by starting r200
                MO_tr200 = Sort_Array(MO_tr200, Decreasing = .true.)

                do j = 1, size(MO_tr200)
                   ChainArray(M)%Chain(1,(MassOrderingGroup(i)%ClusterIndex(j)-1)*nParameter_per_Cluster+1) = MO_tr200(j)
                end do

                deallocate(MO_tr200)
             end do
          end do
       end if


       !--If a certain parameter is not being fit, set to input values (never true for r200, but possibly true for position and concentration)
       if(fit_Parameter(1) == 0) then
          !-r200
          print *, 'Default r200 value has been set by hand, THIS SHOULD BE EDITED'
          do C = 1, size(Group_Index)
             ChainWidths((C-1)*nParameter_per_Cluster+1) = 0.e0_double
             do M = 1, size(ChainArray)
                ChainArray(M)%Chain(1,(C-1)*nParameter_per_Cluster+1) = 1.2e0_double
             END do
          END do
          STOP 'Whilst r200 may not necessarily need to be fit, I cant see why one wouldnt. Therefore, I am stopping'
       end if
       if(fit_Parameter(2) == 0) then
          !-Concentration
       end if
       if(fit_Parameter(3) == 0) then
          !--Position
          do C = 1, size(Group_Index)
             do M = 1, size(ChainArray)
                ChainArray(M)%Chain(1,(C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4) = Ap_pos(Group_Index(C),:)
             end do
          end do
       end if

       !--Testing
       print *, '------------------------------------------------------------------'
       print *, 'Chains starting at:'
       do M = 1, size(ChainArray)
          print *, M, ':', ChainArray(M)%Chain(1,:)
       end do
       print *, '------------------------------------------------------------------'

       print *, 'Chain width used:'
       print *, ChainWidths
       print *, ' '

       !_____________________________________________________Proceed with chain_______________________________________________________________!

       !----Open output files
       !------Ensure no attempt is made to output more chains than exists
       nChainOut = minval((/nChains, nChainOut/))
       do out = 1, nChainOut
          open(unit = 30+out, file = ChainsOutput_Filename(out))
          write(*,'(2(A))') '---Outputting Chain to file: ', ChainsOutput_Filename(out)
          !--Header--!
          write(30+out,'(A,I3,A)') '# Output is: ChainLink (r200, c, RA, DEC)x', size(Group_Index), ' clusters. Final column is ln-likelihood (not renormalised).'
          if(output_Burnin) write(30+out, '(A)') '# Burnin is included.'
       end do
       write(fmt, '(I3)') size(ChainArray(1)%Chain,2) + 1
       fmt = '(I7,x,'//trim(adjustl(fmt))//'(e16.9,x))'
       
       !-Convergence_Test
       open(unit = 28, file = Convergence_Test_Output)
       !--Header--!
       write(28, '(A)') '# Included is R (Gelman-Rubin) across parameters'
       write(conv_fmt, '(I3)') size(ChainArray(1)%Chain,2)
       conv_fmt = '('//trim(adjustl(conv_fmt))//'(e12.5,x))'

       open(unit = 27, file = Acceptance_Rate_Output)
       write(27, '(A)') '# Chain Link; Acceptance Rate by chain'
       write(acc_fmt,'(I4)') size(Acceptance_Rate)
       acc_fmt = '(I5,x,'//trim(adjustl(acc_fmt))//'(e10.3,x))'

       !--Set minimum chain length so that it is at least 20 times the burnin
       nMinChain = maxval((/nMinChain, 20*nBurnin/))

       write(*,'(A,x,I4,x,A)') 'A minimum of:', nMinChain, ' chain links will be evaluated.'

       !--Set Paralellisation definitions
       if(nOMPThread <= 0) then
          nOpenThreads = size(ChainArray)
          !nOpenThreads = minval((/OMP_GET_NUM_THREADS(), size(ChainArray)/))
       else
          nOpenThreads = nOMPThread
       end if
       write(*,'(A,I3,A)') '!---- Running Parallelisation over: ', nOpenThreads, ' threads.'

       !--Initialise counters
       chainCount = 0; allocate(Likelihood(size(ChainArray),nMaxChain)); Acceptance_Rate = 0.e0_double
       chains_Converged = .false.
       do
          !--Exit Conditions
          if(Chains_Converged .and. chainCount >= nMinChain .and. allow_Automatic_Exit_MCMC) exit
          chainCount = chainCount + 1


          print *, 'Considering Chain Link:', chainCount

          if(chainCount > nMaxChain) then
             print *, 'WARNING: Reached the maximum number of allowed chains without convergence, exiting ~~~~~~~~~'
             exit
          end if
          
          !--Reset Acceptance Rate after burnin
          if(chainCount == nBurnin) Acceptance_Rate = 0.e0_double

          !_______________________Propose New Point for each chain_____________________________________________________!
          if(chainCount > 1) then
             do M = 1, size(ChainArray)
                !--Take new point (skip on first run to allow for evaluation of first point
                call MCMC_Propose_Point_TopHat(ChainArray(M)%Chain, ChainWidths)
             end do
          end if

          !______________________EVALUATE LIKELIHOOD AT NEW CHAIN POINT________________________________________________!
          !--Set up Temporary Source Sample array (Parallelisation does not like derived types)
          if(allocated(tMF606W)) deallocate(tMF606W); if(allocated(tRedshift)) deallocate(tRedshift); if(allocated(tSizes)) deallocate(tSizes);
          if(allocated(tPosterior_Method)) deallocate(tPosterior_Method); if(allocated(Posterior_Flag)) deallocate(Posterior_Flag)
          allocate(tMF606W(size(Group_Cat%MF606W))); tMF606W = Group_Cat%MF606W
          allocate(tRedshift(size(Group_Cat%Redshift))); tRedshift = Group_Cat%Redshift
          allocate(tSizes(size(Group_Cat%Sizes))); tSizes = Group_Cat%Sizes
          allocate(tRenormalisationGroup(size(Group_Cat%RenormalisationGroup))); tRenormalisationGroup = Group_Cat%RenormalisationGroup
          allocate(tPosterior_Method(size(Group_Cat%Posterior_Method))); tPosterior_Method = Group_Cat%Posterior_Method
          allocate(Posterior_Flag(size(Group_Cat%Posterior_Method))); Posterior_Flag = '000'

          call CPU_TIME(Time_ChainLink_Start)

          call OMP_SET_NUM_THREADS(nOpenThreads)
          !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(M, tAlpha, tAp_Pos, centroid_Prior_Offset, MO_tr200)
          !$OMP DO
          do M = 1, size(ChainArray)
             if(allocated(tAlpha)) deallocate(tAlpha); if(allocated(tAp_Pos)) deallocate(tAp_Pos)

             !--Set up temporary arrays to pass in free parameters
             allocate(tAlpha(size(Group_Index))); allocate(tAp_Pos(size(Group_Index),2))
             do C = 1, size(Group_Index)
                tAlpha(C) = ChainArray(M)%Chain(chainCount, (C-1)*nParameter_per_Cluster+1)
                tAp_Pos(C,:) = ChainArray(M)%Chain(chainCount, (C-1)*nParameter_per_Cluster+3:(C-1)*nParameter_per_Cluster+4)
             end do

             !__________________________________________Evaluate Posterior
             !--Priors on Parameters
             !---If using position, to avoid degeneracy need alpha1>alpha2>...., and then 1 not longer necessarily labels the 1st cluster, but the largest (identifiable by location)
             if(any(tAlpha <= 0.05)) then
                !--Put a low probability on such low values of alpha so they are never accepted
                Likelihood(M,chainCount) = -1.e0_ldp*huge(1.e0_ldp)!-100000
                cycle
             end if
             if(fit_Parameter(3) == 1) then
                !--If fitting centroid position....
                                
                !--Apply prior on position of centroid. Note: Care must be taken here when two clusters can share the same parameter space, as this will lead to degeneracies which should be accounted for by setting an ordering prior on alpha
                if(allocated(centroid_Prior_Offset)) deallocate(centroid_Prior_Offset)
                allocate(centroid_Prior_Offset(size(Group_Index))); centroid_Prior_Offset = 0.e0_double
                !--Set offset as the distance between the centroid and the centroid position entered as a first guess
                centroid_Prior_Offset = (/ (Pyth_distance_between_Points(tAp_Pos(C,:), Ap_Pos(Group_Index(C),:)), C = 1, size(Group_Index)) /)
                if(any(centroid_Prior_Offset > centroid_Prior_Width)) then
                   !--If any of the centroids fall outside a circular top hat centered on the entered aperture position, with width given by centroid_Prior_Width...
                   Likelihood(M,chainCount) = -1.e0_ldp*huge(1.e0_ldp)
                   cycle
                end if


                if(enforce_Mass_Ordering) then
                   do i = 1, size(MassOrderingGroup) !--For each group
                      if(allocated(MO_tr200)) deallocate(MO_tr200)
                      allocate(MO_tr200(size(MassOrderingGroup(i)%ClusterIndex)));

                      !--Isolate r200 values for that group
                      do j = 1, size(MO_tr200)
                         MO_tr200(j) = tAlpha(MassOrderingGroup(i)%ClusterIndex(j))
                      end do

                      !--Test that r200 is decreasing
                      if(real_Ordered_Array_Decreasing(MO_tr200) == .false.) then
                         !--If fitting position, and Cluster N is larger than N-1 then set probability to zero
                         !----This allows parameters to be fit over the whole field space, where labelling of cluster is not conserved (i.e. returned cluster 1 is the most massive)
                         Likelihood(M,chainCount) = -1.e0_ldp*huge(1.e0_ldp)
                         cycle
                      end if
                      deallocate(MO_tr200)
                   end do
                end if
                
             end if

             !lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample
             !--Likelihood evaluation provided priors have been passed
             Likelihood(M,chainCount) = lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample(tAlpha, tPosterior_Method, Surface_Mass_Profile, tAp_Pos, tLens_Redshift, tSizes, tMF606W, tRedshift, Source_Positions, MagGrid, SizeGrid, Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, Size_Limits_byCluster, Magnitude_Limits_byCluster, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification, RedshiftGrid, tSigma_Crit, use_lnSize, Posterior_Flag, tRenormalisationGroup)

             if(allocated(centroid_Prior_Offset)) deallocate(centroid_Prior_Offset)
             if(allocated(tAlpha)) deallocate(tAlpha); if(allocated(tAp_Pos)) deallocate(tAp_Pos)
          end do
          !$OMP END PARALLEL
          deallocate(tMF606W, tRedshift, tSizes, tRenormalisationGroup, tPosterior_Method, Posterior_Flag)
          !--Check that initial value of the chain has not been placed outside priors
          if(chainCount == 1 .and. any(Likelihood(:,chainCount) <= -(huge(1.e0_double)+100.e0_double))) then
             print *, 'MCMC - Possible error with initial positioning of chain, outside prior limits. FATAL.'
             do M = 1, size(Likelihood,1)
                if(Likelihood(M,1) <= (huge(1.e0_double)+100.e0_double)) then
                   print *, '!--Chain:', M, ChainArray(M)%Chain(1,:)
                   print *, ' '
                end if
             end do
             read(*,*)
          end if

          !___________________________________ACCEPT/REJECT NEW POINT ON EACH CHAIN_________________________________________________!
          do M = 1, size(ChainArray)
             !--Test for acceptance
             if(chainCount > 1) then
                call MCMC_accept_reject_Point_Metropolis(ChainArray(M)%Chain(chainCount-1:chainCount,:), Likelihood(M,chainCount-1:chainCount), acceptance, .true., Acceptance_Rate(M), chainCount-1)

                !--Test Behaviour of priors
                if(any( (/(ChainArray(M)%Chain(chainCount, (C-1)*nParameter_per_Cluster+1), C = 1, size(Group_Index))/) < 0.e0_double) .and. Likelihood(M,chainCount) > -100000) then
                   print *, 'Group, Chain, chainLink:', G, M, chainCount
                   STOP 'MCMC- Fatal- A negative Value of Alpha has been accepted!'
                end if
             end if

             if(chainCount > 1 .and. tune_MCMC) then
                write(*,'(A,x,I2,x,I5,x,e12.5,x,L)') 'Chain, ChainLink, Acceptance Rate: ', M, chainCount, Acceptance_Rate(M), acceptance
                if(M == size(ChainArray)) print *, ' '
             end if

          end do

          !--Output acceptance rate to file
          if(chainCount > 1) write(27, acc_fmt) chainCount, Acceptance_Rate

          !--Output Point to file
          if(output_Burnin .or. chainCount > nBurnin) then
             do out = 1, nChainOut
                write(30+out, fmt) chainCount, ChainArray(out)%Chain(chainCount, :), Likelihood(out,chainCount) 
             end do
          end if

          !--Cycle on first point in chain
          if(chainCount < 2) cycle

          !--Every nConvergenceTestPoint test for convergence across the chains
          if(chainCount > nBurnin .and. mod(chainCount-nBurnin, nConvergenceTestPoint) == 0) then
             Chains_Converged = convergence_Test_GelmanRubin(ChainArray, (/nBurnin/), ConvergenceStatistic, ifree_Parameter, 1.03e0_double)
             write(28,conv_fmt) ConvergenceStatistic
          end if
          if(allocated(ConvergenceStatistic)) deallocate(ConvergenceStatistic)

          call CPU_TIME(Time_ChainLink_End)
          print *, '--ChainLink took:', Time_ChainLink_End-Time_ChainLink_Start, ' seconds'

       end do !--End of chain loop
       print *, 'Finished Chain' !--Could output maximum value using maxloc, maybe average across chains
       do out = 1, nChainOut
          close(unit = 30+out)
       end do
       close(28)
       close(27)

       !--Combine all chains onto on large chain
       allocate(Combined_Chain(size(ChainArray)*size(ChainArray(1)%Chain(nBurnin+1:,1),1), size(ChainArray(1)%Chain,2))); Combined_Chain = dsqrt(-1.e0_double)
       allocate(Combined_Likelihood(size(Combined_Chain,1))); Combined_Likelihood = dsqrt(-1.e0_double)
       do M = 1, size(ChainArray)
          Combined_Chain((M-1)*size(ChainArray(1)%Chain(nBurnin+1:,1),1)+1:M*size(ChainArray(1)%Chain(nBurnin+1:,1),1),:) = ChainArray(M)%Chain(nBurnin+1:,:)
          Combined_Likelihood((M-1)*size(ChainArray(1)%Chain(nBurnin+1:,1),1)+1:M*size(ChainArray(1)%Chain(nBurnin+1:,1),1)) = Likelihood(M,nBurnin+1:)
       end do
       if(any(isNAN(Combined_Chain)) .or. any(isNAN(Combined_Likelihood))) STOP'FATAL: NaNs in combined chain (MCMC)'

       open(unit = 29, file = Combined_Chain_Output)
       write(29,'(A,I3,A)') '# Output is: ChainLink (r200, c, RA, DEC)x', size(Group_Index), ' clusters. Final column is likelihood.'
       write(29, '(A)') '# Burnin is NOT included.'
       do i = 1, size(Combined_Chain,1)
          write(29, fmt) i, Combined_Chain(i, :), Combined_Likelihood(i)
       end do
       close(29)

       !--Construct Marginalised Posteriors from one chain
       print *, 'I am only constructing Marginalised Posterior of Alpha......' !--For other variables, perhaps consider looping over P, and check that parent routine loops over this value
       !--Could use nChainLinks/n for second arguement (nBins), where the ''n'' here sets the number of chain links in each bin
       do C = 1, size(Group_Index)
          call create_Histogram(Combined_Chain(:,(C-1)*nParameter_per_Cluster+1), 40, tMarginalised_Posterior_Grid, tMarginalised_Posterior)
          Marginalised_Posteriors(Group_Index(C), 1+1, :) = Linear_Interp(Marginalised_Posteriors(Group_Index(C),1,:), tMarginalised_Posterior_Grid, tMarginalised_Posterior, ExValue =1.e-100_double)
          deallocate(tMarginalised_Posterior_Grid, tMarginalised_Posterior)
       end do

       deallocate(Combined_Chain, Combined_Likelihood)

       deallocate(ChainsOutput_Filename, Likelihood, Source_Positions, tSigma_Crit, Group_Index, ChainWidths, tLens_Redshift)
       call Catalogue_Destruct(Group_Cat)

       do M = 1, size(ChainArray)
          deallocate(ChainArray(M)%Chain)
       end do
       deallocate(ChainArray)
       deallocate(Acceptance_Rate)

    end do !-End of group loop

    deallocate(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification, MagGrid, SizeGrid, Posterior_Flag)



  end subroutine DM_Profile_Fitting_Simultaneous_MCMC

  function Start_MCMC_Chain_TopHatAperture(Ap_Center, Width)
    use Common_Functions, only: return_Random_Set
    real(double), intent(in):: Ap_Center(:)
    real(double), intent(in):: Width

    real(double), dimension(size(Ap_Center)):: Start_MCMC_Chain_TopHatAperture
    real(double),dimension(size(Ap_Center)):: Ran
    real(double):: Sum2, Pen_Width
    integer::i

    Ran = return_Random_Set(size(Ran))
    !--Shift random onto interval [-1,1]
    Ran = 2.e0_double*(Ran-0.5e0_double)

    !--Randomly set 1st N-1 parameters
    Sum2 = 0.e0_double
    do i = 1, size(Ap_Center)-1
       Start_MCMC_Chain_TopHatAperture(i) = Ap_Center(i) + Ran(i)*Width
       Sum2 = Sum2 + (Start_MCMC_Chain_TopHatAperture(i)-Ap_Center(i))**2.e0_double
    end do

    !--Set width of penultimate co-ordinate as y = sqrt(Width - sum(x^2))
    Pen_Width = dsqrt(Width*Width - Sum2)
    i = size(Start_MCMC_Chain_TopHatAperture)
    Start_MCMC_Chain_TopHatAperture(i) = Ap_Center(i) + Ran(i)*Pen_Width

  end function Start_MCMC_Chain_TopHatAperture

  subroutine DM_Profile_Fitting_Simultaneous_2Cluster(iCat, Ap_Pos, Ap_Radius, Lens_Redshift, Marginalised_Posteriors, Distribution_Directory, reproduce_Prior, Fit_Group, Alpha_Limit, Grid_Tolerance, Alpha_Tolerance, Output_Prefix, Blank_Field_Catalogue)
    use Bayesian_Posterior_Evaluation, only: lnLikelihood_Evaluation_atVirialRadius_perSourceSample, get_likelihood_evaluation_precursors, lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample
    use gridintervals, only: equalscale; use Integration, only: Integrate; use Interpolaters, only: Linear_Interp
    !-Routine that produces posteriors on dark matter profile free parameters, by simultaneously fitting to all clusters that belong to the same 'Fit_Group'.
    !--Fit_Group should be in assending order starting from 1, with no numerical gaps, and should hold a value for each aperture considered. The number of free parameters for each fit is therefore determined by the number of clusters considered in each group.
    !--Alpha_Limit should be a one-dimensional array, size of 2*nCluster. Each pair or elements labels the alpha limits on which the posterior should be evaluated. Outwith these limits, the posterior is set to zero, so these limits CAN BE THOUGHT OF AS THE LIMITS OF A FLAT PRIOR ON ALPHA.
    !--Grid_tolerance is a nCl size 1D array with sets the grid size on which the posterior is evaluated
    !--Alpha Tolerance sets the tolerance to which the maximum of the grid is found. If Alpha_Tolerance < Grid_Tolerance, then a recursive call evaluates the grid on a finer grid around the location of the maximum (either in terms of the marginalised posteriors, or the likelihood)
    !--First implementation will evaluate on a basic grid. Extensions will search for the mode value by: Evaluating on a finer grid around maxloc; Numerical Recipes Simplex method around mode;
    !---Further extension could implement MCMC if more than two free parameters are to be evaluated
    !-- Output is marginalised posterior interpolated for true output. Difficulty in passing out joint posterior results from the unknown size of the joint distribution (rank), and how many need passed out. Instead, they are written to file.

    type(Catalogue), intent(in):: iCat(:)
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:), Lens_Redshift(:)
    real(double),intent(out),allocatable::Marginalised_Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! 
    character(*), intent(in):: Distribution_Directory
    logical, intent(in):: reproduce_Prior
    type(Catalogue), intent(in),optional::Blank_Field_Catalogue(:)
    integer:: Fit_Group(:)
    real(double):: Alpha_Limit(:), Grid_tolerance(:), Alpha_Tolerance
    character(*),intent(in):: Output_Prefix

    
    real(double),dimension(2*size(Ap_Pos,1)):: iAlpha_Limit
    type(Catalogue):: tSource_Catalogue, tCat, Cat
    type(Catalogue):: Group_Cat
    real(double), dimension(Size(Ap_Pos,1))::iAp_Radius
    integer, allocatable:: Group_Index(:)

    integer:: C,G,i,j

    character(500):: Group_Output_Prefix

    real(double),allocatable:: Source_Positions(:,:)

    !--Result Output Decalrations
    integer::nAlpha_Out = 500
    real(double):: AlphaLimit_Out(2) = (/0.05e0_double, 3.e0_double/)

    character(15):: fmt
    character(100):: Filename

    !--Distribution Declarations
    real(double),allocatable:: Joint_Size_Magnitude_Distribution(:,:), Magnitude_Distribution(:), TwoD_MagPrior_Distribution(:,:)
    real(double),allocatable::SizeGrid(:), MagGrid(:)

    !--Posterior Grid Declarations
    integer, allocatable:: nAlpha(:)
    real(double),allocatable:: AlphaGrid1(:), AlphaGrid2(:)

    real(double),allocatable:: Likelihood(:,:)
    character(3),allocatable:: Posterior_Flag(:,:,:)
    real(double),allocatable:: Marginalised_Cl1(:), Marginalised_Cl2(:)

    !--Temporary Allocations--!
    real(double), allocatable:: tSigma_Crit(:,:), tAp_Pos(:,:)
    !----As Part of OPENMP
    real(double), allocatable:: tMF606W(:), tRedshift(:), tSizes(:)
    integer, allocatable:: tPosterior_Method(:), tRenormalisationGroup(:)

    !--Likelihood Evaluation Precursor Declarations
    real(double),dimension(:,:),allocatable:: Survey_Renormalised_Prior
    real(double),dimension(:),allocatable:: Survey_Renormalised_MagPrior
    real(double),allocatable:: RedshiftGrid(:)
    real(double),allocatable::Sigma_Crit(:,:)    
    real(double),allocatable:: MagnificationGrid(:), Renormalisation_by_Magnification(:,:), MagOnly_Renormalisation_by_Magnification(:,:)

    !--Fine Grid declaration
    integer,dimension(2):: Max_Point
    real(double),allocatable:: Max_Marginalised_Posteriors(:,:,:)
    real(double), dimension(2*size(Ap_Pos,1)):: Max_Limit

    !--Single Cluster Fitting Declarations
    real(double),allocatable:: Posterior_Single(:,:)

    !-- Parallelisation
    integer:: OMP_GET_NUM_THREADS

    if(size(Ap_Radius)==1) then
       iAp_Radius = Ap_Radius(1)
    else
       iAp_Radius = Ap_Radius
    end if

    if(size(Alpha_Limit) == 2) then
       do i = 1, size(Ap_Pos,1)
          iAlpha_Limit(2*i-1:2*i) = Alpha_Limit
       end do
    elseif(size(Alpha_Limit) /= size(Ap_Pos,1)) then
       STOP '2Cluster - Alpha_limit Entered is not of the correct size'
    else
       iAlpha_Limit = Alpha_Limit
    end if

    !--Recombine Source Sample
    call recombine_Source_Sample(iCat, Cat)

!DELETE    Output_Prefix = trim(adjustl(Bayesian_Routines_Output_Directory))

    !--Error Catching on Fit Group
    if(size(Fit_Group) /= size(Ap_Pos,1)) STOP 'DM_Profile_Fitting_Simultaneous - Fit_Group entered is not of the correct size.'
    if(minval(Fit_Group) /= 1) STOP 'DM_Profile_Fitting_Simultaneous - Fit_Group entered does not satisfy minimum value conditions (should be 1)'

    !_____________________________________________________________________________Process of Posterior Evaluation_____________________________________________________________________________!
    !--Construct/Read in prior distribution
    !~~~Repeat for all groups in Fit_Group
    !-Process: For each group identify source sample: Should be all galaxies in aperture radius entered. If apertures do not overlap, then output an error message but continue with evaluation
    !- Set up coarse grid for all free parameters
    !- Evaluate Posterior on coarse grid.
    !- Output Joint Posterior for each group to file
    !- Calculate and output marginalised posterior for each group.
    !- Assign Marginalised Posterior output through linear interpolation
    !~~~~ Repeat for next group
    !_____________________________________________________________________________Process of Posterior Evaluation_____________________________________________________________________________!


    !--Get Prior Distribution
    !--If a catalogue for the blank field is passed in, then use this catalogue to get the intrinsic distributions--!
    if(reproduce_Prior) then
       if(present(Blank_Field_Catalogue) == .false.) STOP 'DM_Profile_Variable_Posteriors_CircularAperture - Blank Field Catalogue must be entered to allow for the production of the prior on a grid'
       write(*,'(A)') 'Producing Distribution from Catalogue'
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, TwoD_MagPrior_Distribution, Distribution_Directory, .true., Blank_Field_Catalogue)
    else
       write(*,'(A)') 'Reading in distribution from:', Distribution_Directory
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, TwoD_MagPrior_Distribution, Distribution_Directory, .false., Blank_Field_Catalogue)
       print *, 'Success:', Distribution_Directory
    end if

    !--Get Precursor
    call get_Likelihood_Evaluation_Precursors(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid, Renormalisation_by_Magnification, MagOnly_Renormalisation_by_Magnification, SizeGrid, MagGrid, Magnitude_Distribution, TwoD_MagPrior_Distribution, Joint_Size_Magnitude_Distribution, S1_Magnitude_Limits_byCluster, S1_Size_Limits_byCluster, S2_Magnitude_limits_byCluster, S2_Size_Limits_byCluster, Lens_Redshift, Lower_Redshift_Cut, Output_Prefix, use_lnSize)


    !--Set up output grid, on which the marginalised posteriors will be output, as the linear interpolation of the constructed marginalised posterior below
    allocate(Marginalised_Posteriors(size(Ap_Pos,1), 2, nAlpha_Out));
    do i =1, nAlpha_Out
       Marginalised_Posteriors(:,1,i) = (/(AlphaLimit_Out(1) + (i-1)*((AlphaLimit_Out(2)-AlphaLimit_Out(1))/(nAlpha_Out-1)),j=1,size(Marginalised_Posteriors,1))/)
    end do

    do G = 1, 100 !-Assume no more than 100 group will be present
       !--Find all clusters with that grouping
       !---Group Index contains the index of the clusters within that group, and can be used to access the correct Aperture Location and Radius. Size of Group_Index is the number of free parameters being used
       allocate(Group_Index(count(Fit_Group == G)));
       !--Exit when all groups exhausted
       if(size(Group_Index) == 0) then
          if(maxval(Fit_Group) > G) then
             deallocate(Group_Index)
             cycle
          else
             exit
          end if
       end if
       !--Exit if the grouping is too large. This is limited mainly by the inability to write a simultaneous fitting routine for a generally sized array - BOLLOCKS TO FORTRAN
       if(size(Group_Index) > 2) STOP 'DM_Profile_Fitting_Simultaneous - I have not been written to deal with the simultaneous fitting of more than two clusters at once, stopping.'

       print *, ' '
       print *, 'Joint Fitting Cluster Group ', G, ' of ', maxval(Fit_Group),'....'
       
       i = 0
       do C = 1, size(Fit_Group)
          if(Fit_Group(C) == G) then
             i = i + 1
             Group_Index(i) = C
          end if
       end do

       Group_Cat = Select_Source_Sample(Cat, Ap_Pos, iAp_Radius, Group_Index, Mag_Limits = Magnitude_Limits_byCluster, Size_Limits = Size_Limits_byCluster, Redshift_Limit = 1.2e0_double*Lens_Redshift)

!!$       tCat = Cat
!!$       do C = 1, size(Group_Index)
!!$!          print *, 'Cutting on a core radius of:', Core_Cut_Radius(Group_Index(C)), ' arcminutes for aperture:', Group_Index(C)
!!$          !--Get galaxies in aperture
!!$          call Identify_Galaxys_in_Circular_Aperture(tCat, Ap_Pos(Group_Index(C),:), iAp_Radius(Group_Index(C)), TSource_Catalogue)!, Core_Radius = Core_Cut_Radius(Group_Index(C))/60.e0_double)
!!$          call Concatonate_Catalogues(Group_Cat, TSource_Catalogue) 
!!$          call Catalogue_Destruct(TSource_Catalogue)
!!$          !--Mask Source Galaxies for that aperture to ensure no double counting
!!$          print *, '--__Applying mask as part of 2Cluster, to enusre no double counting...'
!!$          call Mask_Circular_Aperture(tCat, Ap_Pos(Group_Index(C),:),  iAp_Radius(Group_Index(C)))          
!!$       end do
!!$       call Catalogue_Destruct(tCat)
!!$       print *, 'Got source sample'

       if(allocated(Source_Positions)) deallocate(Source_Positions)
       allocate(Source_Positions(size(Group_Cat%RA),2));
       Source_Positions(:,1) = Group_Cat%RA; Source_Positions(:,2) = Group_Cat%Dec

       !--Check for overlap between clusters
       
       !--Set up a coarse grid on which the posterior will be evaluated
       !---For this point on the code is limited to two dimensions
       if(allocated(nAlpha)) deallocate(nAlpha)
       allocate(nAlpha(2)); nAlpha = 1
       do C= 1, size(Group_Index)
          nAlpha(C) = nint((IAlpha_Limit(2*Group_index(C)) - IAlpha_Limit(2*Group_index(C)-1))/Grid_tolerance(Group_index(C))+1)
       end do
       print *, 'Simultaneous fit for Clusters:', Group_Index, ' will take ', nAlpha(1)*nAlpha(2), 'grid points'
       call equalscale(IAlpha_Limit(2*Group_index(1)-1), IAlpha_Limit(2*Group_index(1)), nAlpha(1), AlphaGrid1)
       
       if(size(Group_Index) == 1) then
          !--Do single fitting of cluster
          print *, 'Just attempting a single fit...'
          allocate(Posterior_Single(2,size(AlphaGrid1))); Posterior_Single = dsqrt(-1.e0_double); Posterior_Single(1,:) = AlphaGrid1
          call DM_Profile_Variable_Posterior_SingleFit(Group_Cat, Surface_Mass_Profile, Lens_Redshift(Group_Index(1)), Ap_Pos(Group_Index(1),:), Posterior_Single, Output_Prefix, .false., PriorMagGrid = MagGrid, PriorSizeGrid = SizeGrid, Prior = Joint_Size_Magnitude_Distribution, MagPrior = Magnitude_Distribution)
          Marginalised_Posteriors(Group_Index(1),2,:) = Linear_Interp(Marginalised_Posteriors(Group_Index(1),1,:), Posterior_Single(1,:), Posterior_Single(2,:), ExValue =1.e-100_double)
          deallocate(Posterior_Single, AlphaGrid1, Group_Index, Source_Positions)
          call Catalogue_Destruct(Group_Cat)
          print *, 'Got single fit of cluster'

          cycle
       end if

       call equalscale(IAlpha_Limit(2*Group_index(2)-1), IAlpha_Limit(2*Group_index(2)), nAlpha(2), AlphaGrid2)
       
       print *, ' '
       write(*,'(A,2(e9.2,x),A,I4,A,I1,A)') '--Alpha is being evaluated between:', minval(AlphaGrid1), maxval(AlphaGrid1), ', in ', nAlpha(1), ' grid points (Cluster ', Group_index(1), ')'
       write(*,'(A,2(e9.2,x),A,I4,A,I1,A)') '--Alpha is being evaluated between:', minval(AlphaGrid2), maxval(AlphaGrid2), ', in ', nAlpha(2), ' grid points (Cluster ', Group_index(2), ')'
       print *, ' '

       if(allocated(Likelihood)) deallocate(Likelihood)
       allocate(Likelihood(nAlpha(1), nAlpha(2))); Likelihood = 0.e0_double
       if(allocated(Posterior_Flag)) deallocate(Posterior_Flag)
       allocate(Posterior_Flag(size(Likelihood,1), size(Likelihood,2), size(Group_Cat%Posterior_Method))); Posterior_Flag = '000'

       allocate(tSigma_Crit(2,size(Sigma_Crit,2))); tSigma_Crit(1,:) = Sigma_Crit(Group_index(1),:); tSigma_Crit(2,:) = Sigma_Crit(Group_index(2),:) 
       allocate(tAp_Pos(2, 2)); tAp_pos(1,:) = Ap_Pos(Group_index(1),:); tAp_Pos(2,:) = Ap_Pos(Group_index(2),:)

       !-----Evaluate Posterior - Parallelisation used, however can be switched off by setting the number of threads to 1. All varabiables passed to likelihood evaluation need to be intent(in) only, and be careful of library calls which depend on global declarations.
       allocate(tMF606W(size(Group_Cat%MF606W))); tMF606W = Group_Cat%MF606W
       allocate(tRedshift(size(Group_Cat%Redshift))); tRedshift = Group_Cat%Redshift
       allocate(tSizes(size(Group_Cat%Sizes))); tSizes = Group_Cat%Sizes
       allocate(tPosterior_Method(size(Group_Cat%Posterior_Method))); tPosterior_Method = Group_Cat%Posterior_Method
       allocate(tRenormalisationGroup(size(Group_Cat%RenormalisationGroup))); tRenormalisationGroup = Group_Cat%RenormalisationGroup
       !--Evaluate Posterior

       call OMP_SET_NUM_THREADS(nOMPThread)
       write(*,'(A,I2,A)') '---Paralellising over:', nOMPThread, ' threads.'
       !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(i,j)
       !--OMP DO parallelises the *outer loop only*- is nalpha(1) == nalpha(2), this could be collapsed using collapse(2)
       !$OMP DO
       do i = 1, nAlpha(1)
          do j = 1, nAlpha(2)
             Likelihood(i,j) = lnLikelihood_atVirialRadius_MultipleCluster_perSourceSample((/AlphaGrid1(i),AlphaGrid2(j)/), tPosterior_Method, Surface_Mass_Profile, tAp_Pos, (/(Lens_Redshift(Group_Index(i)), i = 1, 2)/), tSizes, tMF606W, tRedshift, Source_Positions, MagGrid, SizeGrid, Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, Size_Limits_byCluster, Magnitude_Limits_byCluster, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification, RedshiftGrid, tSigma_Crit, use_lnSize, Posterior_Flag(i,j,:), tRenormalisationGroup)
             print *, 'Done:', i, j,  Likelihood(i,j)
          end do
       end do
       !$OMP END PARALLEL
       deallocate(tRedshift, tMF606W, tSizes, tPosterior_Method, tRenormalisationGroup)

       !--Convert from lnP to P
       Likelihood = dexp(Likelihood - maxval(Likelihood))
       
       !--Renormalise Posterior
       !-(Flat prior assumed at this point)
       !--This may not be strictly necessary
       Likelihood = Likelihood/Integrate(AlphaGrid1, AlphaGrid2, Likelihood, 2, lim1 = (/minval(AlphaGrid1), maxval(AlphaGrid1)/), lim2 = (/minval(AlphaGrid2), maxval(AlphaGrid2)/))
       
       deallocate(tSigma_Crit, Source_Positions, tAp_Pos)
       call Catalogue_Destruct(Group_Cat)
       
       !--Output Full Likelihood (or posterior, assuming flat prior bounded by grid)
       write(Filename, '(I1)') G
       open(unit = 21, file = trim(adjustl(Output_Prefix))//'_SimultaneousFitting_Likelihood_Group'//trim(Filename)//'.dat')
       write(21,'(2(A,x,I2,x))') '# Column 1 contains Virial Radius for cluster:', Group_Index(1), ' whilst Row 1 constains cluster:', Group_Index(2)

       !--Set Output Format
       write(fmt, '(I3)') size(Likelihood,2)+1
       fmt = '('//trim(fmt)//'(e12.5,x))'
       !--Write first row containing grid along second cluster
       write(21, fmt) 0.e0_double, AlphaGrid2
       do i = 1, size(Likelihood,1)
          write(21, fmt) AlphaGrid1(i), Likelihood(i,:)
       end do
       write(*,'(A,I3,A,A)') '---Likelihood for group:', G, ' output to ', trim(adjustl(Output_Prefix))//'_SimultaneousFitting_Likelihood_Group'//trim(Filename)//'.dat'
       close(21)

       !--Output Flag for posterior evaluation - REMOVED UNTIL A WAY OF DEALING WITH ALL GALXIES IS CONSIDERED
!!$       open(unit = 22, file = trim(adjustl(Output_Prefix))//'_Likelihood_Evalutation_Flag'//trim(Filename)//'.dat')
!!$       write(fmt, '(I3)') size(Likelihood,2)+1
!!$       fmt = '('//trim(fmt)//'(A3,x))'
!!$       !--Write first row containing grid along second cluster
!!$       write(22, fmt) 0.e0_double, AlphaGrid2
!!$       do i = 1, size(Posterior_Flag,1)
!!$          write(22, fmt) AlphaGrid1(i), Posterior_Flag(i,:)
!!$       end do
!!$       write(*,'(A,I3,A,A)') '---Likelihood Flags for group:', G, ' output to ', trim(adjustl(Output_Prefix))//'_SimultaneousFitting_Likelihood_Group'//trim(Filename)//'.dat'
!!$       close(22)
!!$       deallocate(Posterior_Flag)
       
       !--Get Marginalised Distributions
       !(At this point, a flat prior is assumed)
       allocate(Marginalised_Cl1(nAlpha(1))); Marginalised_Cl1 = 0.e0_double
       do i = 1, size(Marginalised_Cl1)
          Marginalised_Cl1(i) = Integrate(AlphaGrid2, Likelihood(i,:), 2, lim = (/minval(AlphaGrid2), maxval(AlphaGrid2)/))
       end do
       
       write(Filename, '(I1)') Group_Index(1)
       open(unit = 22, file = trim(adjustl(Output_Prefix))//'_SimultaneousFitting_MarginalisedLikelihood_Cluster'//trim(Filename)//'.dat')
       do i= 1, size(AlphaGrid1)
          write(22, '(2(e12.5,x))') AlphaGrid1(i), Marginalised_Cl1(i)
       end do
       print *, '--- Marginalised Likehood for Cluster:', Group_Index(1),' of group:', G, ' output to ', trim(adjustl(Output_Prefix))//'_SimultaneousFitting_MarginalisedLikelihood_Cluster'//trim(Filename)//'.dat'
       close(22)

       allocate(Marginalised_Cl2(nAlpha(2))); Marginalised_Cl2 = 0.e0_double
       do i = 1, size(Marginalised_Cl2)
          Marginalised_Cl2(i) = Integrate(AlphaGrid1, Likelihood(:,i), 2, lim = (/minval(AlphaGrid1), maxval(AlphaGrid1)/))
       end do

       write(Filename, '(I1)') Group_Index(2)
       open(unit = 22, file = trim(adjustl(Output_Prefix))//'_SimultaneousFitting_MarginalisedLikelihood_Cluster'//trim(Filename)//'.dat')
       do i= 1, size(AlphaGrid2)
          write(22, '(2(e12.5,x))') AlphaGrid2(i), Marginalised_Cl2(i)
       end do
       print *, '--- Marginalised Likehood for Cluster:', Group_Index(2),' of group:', G, ' output to ', trim(adjustl(Output_Prefix))//'_SimultaneousFitting_MarginalisedLikelihood_Cluster'//trim(Filename)//'.dat'
       close(22)
       
       !--Find Maximum of Likihood
       Max_Point=  maxloc(likelihood)
       !--Use of Group Tolerance ensures that the true maximum is bounded within the points considered
       Max_Limit(2*Group_Index(1)-1:2*Group_Index(1)) = (/AlphaGrid1(Max_Point(1))-Grid_Tolerance(Group_Index(1)),AlphaGrid1(Max_Point(1))+Grid_Tolerance(Group_Index(1))/)
       Max_Limit(2*Group_Index(2)-1:2*Group_Index(2)) = (/AlphaGrid1(Max_Point(2))-Grid_Tolerance(Group_Index(2)),AlphaGrid1(Max_Point(2))+Grid_Tolerance(Group_Index(2))/)
       
       !--Output Marginalised Posteriors
       !-Cluster, Grid/Value, Value
       Marginalised_Posteriors(Group_Index(1),2,:) = Linear_Interp(Marginalised_Posteriors(Group_Index(1),1,:), AlphaGrid1, Marginalised_Cl1, ExValue = 1.e-100_double)
       Marginalised_Posteriors(Group_Index(2),2,:) = Linear_Interp(Marginalised_Posteriors(Group_Index(2),1,:), AlphaGrid2, Marginalised_Cl2, ExValue = 1.e-100_double)

       !-Reset for next group
       deallocate(Marginalised_Cl1, Marginalised_Cl2, Likelihood, AlphaGrid1, AlphaGrid2)
       deallocate(Group_Index)
    end do

    deallocate(Joint_Size_Magnitude_Distribution, Magnitude_Distribution, SizeGrid, MagGrid)
    deallocate(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification)
    
    !--Recursive call to find maximum of likelihood to a given tolerance
    !--Switched off for now to test ability of code to simultaneously fit two clusters
!!$    if(any(Grid_Tolerance > Alpha_Tolerance)) then
!!$       !--This could be edited so that the prior does not need to be read in again, or indeed that the precursors do not need to be read in
!!$       !--Pass new Alpha tolerance ten times greater to enusre recursive loop does not continue more than once
!!$       print *, 'Max limit is:', Max_Limit
!!$
!!$       call DM_Profile_Fitting_Simultaneous_2Cluster(Cat, Ap_Pos, Ap_Radius, Max_Marginalised_Posteriors, Distribution_Directory, .false., Fit_Group, Max_Limit, (/(Alpha_Tolerance, i =1, size(Ap_Pos,1))/), Alpha_Tolerance*10.e0_double, trim(adjustl(Output_Prefix))//'_FineGrid_AroundMaximum')
!!$       
!!$       deallocate(Max_Marginalised_Posteriors, Max_Limit)
!!$    end if
    
    
  end subroutine DM_Profile_Fitting_Simultaneous_2Cluster

  

  subroutine DM_Profile_Variable_Posteriors_CircularAperture(iCat, Ap_Pos, Ap_Radius, Lens_Redshift, Posteriors, Distribution_Directory, reproduce_Prior, Blank_Field_Catalogue)
    use Statistics, only: mean_discrete, mode_distribution, variance_distribution; use Mass_profiles
    use Distributions; use Cosmology, only: angular_diameter_distance_fromRedshift; use Interpolaters, only: Linear_Interp; use gridintervals, only: equalscale
    !--Main routine that returns the Posteriors on DM free parameter (alpha) over all apertures.
    !--If the posterior construction is such that grid may vary, this interpolates onto a finer grid (linear interpolation). Thus returned posteriors on on the same grid
    !--Ap_Radius in DEGREES
    !--If Blank_Field_Catalogue is entered, then intrinsic distributions are produced form this Catalogue
    !--To Do: 
    !~~Recent edits to use a joint size magnitude have ignored magnitude binning


    type(Catalogue), intent(in)::iCat(:)
    real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:), Lens_Redshift(:)
    real(double),intent(out),allocatable::Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! - Second dimension allows for a different grid for each aperture
    character(*), intent(in):: Distribution_Directory
    logical, intent(in):: reproduce_Prior
    type(Catalogue), intent(in),optional::Blank_Field_Catalogue(:)


    real(double),dimension(Size(Ap_Pos,1))::iAp_Radius

    type(catalogue):: Cat
    type(catalogue),dimension(size(Ap_pos,1))::Ap_Cats
    integer::i, ap, m, j

    !--Size Grid Declarations--!     
    integer::nSizeGrid = 100        
!    real(double)::SizeGrid_Lower = 0.e0_double, SizeGrid_Upper = 70.e0_double, dSizeGrid !-Pixel-! 
    real(double)::SizeGrid_Lower = 0.e0_double, SizeGrid_Upper = 0.0035e0_double, dSizeGrid
    real(double),allocatable::SizeGrid(:)
    real(double),allocatable::MagGrid(:)

    real(double),allocatable:: Joint_Size_Magnitude_Distribution(:,:), Magnitude_Distribution(:), TwoD_MagPrior_Distribution(:,:)

    real(double),allocatable::SizePrior_byMag(:,:) !-MagBin, GridValue-!
    real(double),allocatable::Aperture_Posterior_byMag(:,:,:) !-MagBin, Grid/Posterior, Value-! 

    !--Coarse Grid declarations
    integer::nCoarseGrid
    real(double):: Alpha_min, Alpha_Max
    real(double),allocatable:: Coarse_LikelihoodGrid(:)

    type(Binned_Catalogue)::BCat
    real(double),allocatable::Posterior_Single(:,:) !-Grid/Posterior, Value-!
    real(double),allocatable::Prior_Single(:,:) !!-Grid/Prior, Value-! 

    real(double)::Renorm

    character(2)::fmtString, apString

    !--Conversion to Mass--!
    real(double)::D_l, Area
    real(double),allocatable::Cluster_Mean(:), Cluster_Variance(:), Cluster_Mode(:), AntiSymm_Variance(:,:)
    real(double):: Scale_Mass, Scale_Mass_Error(2)
    character(10):: Mass_String, Error_Mass_String_Positive, Error_Mass_String_Negative

    !--MAgnitude Binning (by Absolute Magnitude)--!
    integer::nMag = 1
    real(double),allocatable::MagBins(:,:)
    integer::Magnitude_Binning_Type = 1 !-1:Absolute, 2:Apparent (MF606)-!

    character(500):: Output_File_Prefix

    type(catalogue)::TCat

    INTERFACE
         subroutine DM_Profile_Variable_Posteriors_CircularAperture(iCat, Ap_Pos, Ap_Radius, Lens_Redshift, Posteriors, Distribution_Directory, reproduce_Prior, Blank_Field_Catalogue)
           use Param_Types; use Catalogues
           !--Main routine that returns the Posteriors over all apertures.
           !--Ap_Radius in DEGREES
           type(Catalogue), intent(in)::iCat(:)
           real(double),intent(in)::Ap_Pos(:,:), Ap_Radius(:), Lens_Redshift(:)
           real(double),intent(out),allocatable::Posteriors(:,:,:) !-Aperture, Grid/Posterior, Value-! - Second dimension allows for a different grid for each aperture
           character(*), intent(in):: Distribution_Directory
           logical, intent(in):: Reproduce_Prior

           type(Catalogue), intent(in),optional::Blank_Field_Catalogue(:)
         end subroutine DM_Profile_Variable_Posteriors_CircularAperture
      end INTERFACE



    if(size(Ap_Radius)==1) then
       iAp_Radius = Ap_Radius(1)
    else
       iAp_Radius = Ap_Radius
    end if

    !---Recombine Source Sample
    call recombine_Source_Sample(iCat, Cat)

    !--Get Mean of size and magnitude distributions in masked apertures
    print *, '--__Applying mask to data as part of single cluster routine to get mean od size and magnitude (may be second time):'
    call Mask_Circular_Aperture(TCat, Core_Cut_Position, Core_Cut_Radius/60.e0_double)
    TCat = Cat

    print *, '**Global catalogue has mean [Core mask accounted for] (Size, Mag):', mean_discrete(TCat%Sizes), mean_discrete(TCat%MF606W)
    call Catalogue_Destruct(TCat)
!!$
    print *, ' '

    do i =1, size(Ap_Cats)
       Ap_Cats(i) = Select_Source_Sample(Cat, Ap_Pos, iAp_Radius, (/i/), Mag_Limits = Magnitude_Limits_byCluster, Size_Limits = Size_Limits_byCluster, Redshift_Limit = 1.2e0_double*Lens_Redshift)
    end do

!--Deprecated
!!$    !--Identify Reduced Catalogue for each aperture--!
!!$    do i =1, size(Ap_Cats)
!!$       print *, 'Searching for position of maxima in Size-Magnitude Shifts for Cluster:', i
!!$!       call Maximise_Convergenceb_yShifts_inAperture(Cat, Blank_Field_Catalogue, Ap_Pos(i,:), Ap_Radius(i))
!!$
!!$!       print *, 'Cutting on a core radius of:', Core_Cut_Radius(i), ' arcminutes for aperture:', i
!!$       call Identify_Galaxys_in_Circular_Aperture(Cat, Ap_Pos(i,:), iAp_Radius(i), Ap_Cats(i))!, Core_Radius = Core_Cut_Radius(i)/60.e0_double)
!!$       
!!$       print *, '* Ap ', i ,' has mean (Size,Mag):',  mean_discrete(Ap_Cats(i)%Sizes), mean_discrete(Ap_Cats(i)%MF606W)
!!$    end do
!!$    print *, ' '

    do i =1, size(Ap_Cats)
       print *, 'Aperture:', i, ' contains:', size(Ap_Cats(I)%RA), ' galaxies'
    end do

    print *, count(Cat%Redshift >= 0.e0_double), ' of ',size(Cat%Redshift), ' galaxies have redshift information'

    !--If a catalogue for the blank field is passed in, then use this catalogue to get the intrinsic distributions--!
    print *, 'RE-EVAL PRIOR?:', reproduce_Prior
    if(reproduce_Prior) then
       if(present(Blank_Field_Catalogue) == .false.) STOP 'DM_Profile_Variable_Posteriors_CircularAperture - Blankf Field Catalogue must be entered to allow for the production of the prior on a grid'
       write(*,'(A)') 'Producing Distribution from Catalogue'
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, TwoD_MagPrior_Distribution, Distribution_Directory, .true., Blank_Field_Catalogue)
    else
       write(*,'(A)') 'Reading in distribution from (here):', Distribution_Directory
       call return_Prior_Distributions(MagGrid, SizeGrid, Joint_Size_Magnitude_Distribution, Magnitude_Distribution, TwoD_MagPrior_Distribution, Distribution_Directory, .false., Blank_Field_Catalogue)
    end if

    if(size(magBins,1) > 1) STOP 'I have had to disable Magnitude Binning for now, youll have to edit the code to get this to work, stopping'

    !--Declarations for the use of a coarse likelihood grid
    !--Note that the Fitting Process used may include a search for the maximum, adding points to this grid
    select case(Surface_Mass_Profile)
    case(1) !-Flat-!
       nCoarseGrid = 100000
       Alpha_Min = -5.e-3_double; Alpha_Max = 1.e-2_double !-- SMD ~ Masses 10^12 -> 10^15 Msun/h, in Units of 10^18 Msun/h  
    case(2) !-SIS-!
       nCoarseGrid = 10000 !--This needs to be sigma_v^2--!
       Alpha_Min= -2.5e5_double; Alpha_Max = 9.e6_double    !--------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case(3) !-NFW-!
       nCoarseGrid = 100
       Alpha_Min= 0.05e0_double; Alpha_Max = 3.e0_double
    case default
       STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
    end select
    call equalscale(Alpha_min, Alpha_Max, nCoarseGrid, Coarse_LikelihoodGrid)

    Do Ap = 1, size(Ap_Cats)
       write(Output_File_Prefix,'(I2)') Ap
       Output_File_Prefix = trim(adjustl(Bayesian_Routines_Output_Directory))//'Aperture_'//trim(adjustl(Output_File_Prefix))//'_'

       !--Set input Posterior Grid
       allocate(Posterior_Single(2,size(Coarse_LikelihoodGrid))); Posterior_Single(1,:) = Coarse_LikelihoodGrid

       !--Evaluate posterior on grid
       call DM_Profile_Variable_Posterior_SingleFit(Ap_Cats(Ap), Surface_Mass_Profile, Lens_Redshift(Ap), Ap_Pos(Ap,:), Posterior_Single, Output_File_Prefix, .false., PriorMagGrid = MagGrid, PriorSizeGrid = SizeGrid, Prior = Joint_Size_Magnitude_Distribution, MagPrior = Magnitude_Distribution, TwoD_MagPrior_Distribution = TwoD_MagPrior_Distribution)  

       !--Allow the finalised posterior to be determined on a seperate grid to the input posterior, to account for the fact that the posterior routine may have added points as part of the search for the maximum of the posterior
       if(Ap == 1) then
          allocate(Posteriors(size(Ap_Cats), 2, maxval((/500, 2*size(Posterior_Single,2)/)))); Posteriors = 0.e0_double
          if(size(Posteriors,3) < size(Posterior_Single,2)) print *, 'WARNING - DM_Profile_Variable_Posteriors_CircularAperture - Posterior interpolated on COARSER grid than output'
          !--Set up Posterior Grid
          do i = 1, size(Posteriors,3)
             Posteriors(:,1,i) = minval(Posterior_Single(1,:)) + (i-1)*((maxval(Posterior_Single(1,:))-minval(Posterior_Single(1,:)))/(size(Posteriors,3)-1))
          end do
       end if
       Posteriors(Ap,2,:) = Linear_Interp(Posteriors(Ap,1,:), Posterior_Single(1,:), Posterior_Single(2,:))

       deallocate(Posterior_Single)
    end Do
    if(allocated(SizeGrid)) deallocate(SizeGrid)
    if(allocated(MagGrid)) deallocate(MagGrid)
    if(allocated(Joint_Size_Magnitude_Distribution)) deallocate(Joint_Size_Magnitude_Distribution)


    !--Output Posterior--!
    !--Moved to Single Run Routine
!!$    open(unit = 51, file = trim(Bayesian_Routines_Output_Directory)//'Posterior_per_Aperture.dat')
!!$    write(fmtstring,'(I2)') size(Posteriors,1)+1 !-Assumes all posteriors described by the same posterior grid-!
!!$    do j = 1, size(Posteriors,3)
!!$       write(51, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posteriors(1,1,j), Posteriors(:,2,j)
!!$    end do
!!$    close(51)
!!$    print *,'Output file to: ', trim(Bayesian_Routines_Output_Directory)//'Posterior_per_Aperture.dat'


  end subroutine DM_Profile_Variable_Posteriors_CircularAperture


  subroutine find_Maximum_by_Bisection(Likelihood, New_Point, tol, Keep_Going)
    use derivatives_lib, only:derivative_2pt; use Common_Functions, only: setNaN
    !--Assumes that the likelihood is well behaved and has a single maximum
    !--Returns the Likelihood with the addition of the bisected point, and the index of that point
    real(double), intent(inout),allocatable:: Likelihood(:,:)
    integer, intent(out):: New_Point
    real(double), intent(in)::tol
    logical, intent(out):: Keep_Going

    integer:: Maximum_Index
    integer:: i

    integer, save:: callcount, calltolerance = 100
    integer:: Bracket(3)
    logical::Bracket_Found

    real(double),dimension(size(Likelihood,1), size(Likelihood,2)+1):: tLikelihood
    real(double)::NewGridPoint

    !--Attemptd means using derivatives, but I am not convinced this works, therefroe obsolete
!!$    iMaximum_Boundary = (/1, size(Likelihood,2)/)
!!$    if(present(Maximum_Boundary)) then
!!$       iMaximum_Boundary = Maximum_Boundary
!!$    end if
!!$
!!$    !--Find maximum between boundary points - possibly do this using maxloc?
!!$    Maximum_Index = -1
!!$    do i = Maximum_Boundary(1), Maximum_Boundary(2)
!!$       if((derivative_2pt(Likelihood(2,i-1), Likelihood(2,i), Likelihood(1,i)-Likelihood(1,i-1)) > 0.e0_double) .and. (derivative_2pt(Likelihood(2,i+1), Likelihood(2,i+2), Likelihood(1,i+2)-Likelihood(1,i+1))< 0.e0_double)) Maximum_Index = i
!!$       !-Maximum_Index marks the lower boundary of the interval in which the maximum is exptect to occur
!!$    end do
!!$    if(Maximum_Index < 0) STOP 'find_Maximum_by_Bisection - FAILED TO FIND MAXIMUM'
!!$ 
!!$    print *, 'Maximum found to be between:', Likelihood(1,i), Likelihood(1,i+1)

    !--If call count too high, output warning, assign flag and output (flag can be used toeither stop the code or use the result as it stands)
    if(callcount > calltolerance) then
    print *, 'WARNING: CALL TOLERANCE FOR MAXIMUM FINDING IS TOO HIGH, EXITING'
       CALLCOUNT = 0
       Keep_Going = .false.
       return
    end if

    !--Find Bracket of maximum on original course grid
    !----This could potentially be sped up by keeping the previous braceting interval and only updating within some region of it
    Bracket_Found = .false.; Bracket = 0.
    do i = 2, size(Likelihood,2)-1
       if(Likelihood(2,i-1)<Likelihood(2,i) .and. Likelihood(2,i)>=Likelihood(2,i+1)) then
          if(Bracket_Found) then
             !--Compare the brackets to find the smallest (overcoming local minima
             print *, 'Multiple Mimina found!', Likelihood(:,i), ':', Likelihood(:,Bracket(2))
             if(Likelihood(2,i) > Likelihood(2, Bracket(2))) then
                cycle
             end if
             print *, 'NEW MAXIMA SET!'
          end if

          !--Bracket Found-!
          Bracket = (/i-1, i, i+1/) !--Middle of bracket is best guess at the minimum
          Bracket_Found = .true.
       end if
    end do

    if(Bracket(1) == 0 .or. Bracket(3) == size(Likelihood,2)) then
       print *, 'Maximum found to be intolerably close to the limits of the original grid, exiting without result'
       callcount = 0
       Keep_Going = .false.
       return
    end if

    !--Found within tolerance
    if( maxval((/Likelihood(1,Bracket(2))-Likelihood(1,Bracket(1)),Likelihood(1,Bracket(3))-Likelihood(1,Bracket(2))/)) <= tol ) then
       print *, 'Found Maximum to within a tolerance, expected to be:', Likelihood(1,Bracket(2))
       callcount = 0
       Keep_Going = .false.
       return
    end if

    !--Bisect the Grid and add this point onto the posterior 
!    NewGridPoint = 0.5e0_double*(Likelihood(1,Bracket(2))+Likelihood(1,Bracket(2+((-1)**callcount)))) !--1 allows for alternate sampling to left and right of entered bracket
 
    !--Add new point using golden ratio
    !--(which sets the new point at 0.38197 into the larger of the intervals)
    if( (Likelihood(1,Bracket(2))-Likelihood(1,Bracket(1))) >= (Likelihood(1,Bracket(3))-Likelihood(1,Bracket(2))) ) then
       NewGridPoint = Likelihood(1,Bracket(2)) - 0.38197e0_double*(Likelihood(1,Bracket(2))-Likelihood(1,Bracket(1)))
       Maximum_Index = Bracket(1)
    else
       NewGridPoint = Likelihood(1,Bracket(2)) + 0.38197e0_double*(Likelihood(1,Bracket(3))-Likelihood(1,Bracket(2)))
       Maximum_Index = Bracket(2)
    end if

!-- Use this if not using golden ratio--!
!!$    if( (-1)**callcount > 0) then
!!$       Maximum_Index = Bracket(2)
!!$    else
!!$       Maximum_Index = Bracket(1)
!!$    END if
    !--Assign new point
    tLikelihood(:,:Maximum_Index) = Likelihood(:,:Maximum_Index)
    tLikelihood(1,Maximum_Index+1) = NewGridPoint; tLikelihood(2,Maximum_Index+1) = setNaN()
    tLikelihood(:,Maximum_Index+2:) = Likelihood(:,Maximum_Index+1:)

    New_Point = Maximum_Index+1


    deallocate(Likelihood)
    allocate(Likelihood(size(tLikelihood,1), size(tLikelihood,2))); Likelihood = tLikelihood

  end subroutine find_Maximum_by_Bisection

  subroutine DM_Profile_Variable_Posterior_SingleFit(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorMagGrid, PriorSizeGrid, Prior, MagPrior, TwoD_MagPrior_Distribution)
    use cosmology, only:angular_diameter_distance_fromRedshift; use MC_Redshift_Sampling, only: Monte_Carlo_Redshift_Sampling_SigmaCritical; use Mass_Profiles; use Distributions, only: ch08_redshift_distribution_Array, CH08_Redshift_Distribution_Scalar; use Interpolaters, only: Linear_Interp; use Bayesian_Posterior_Evaluation, only: lnLikelihood_Evaluation_atVirialRadius_perSourceSample, get_Likelihood_Evaluation_Precursors
    use Integration, only:TrapInt, Integrate; use Matrix_methods, only: Determinant, Matrix_Invert; use Smoothing, only: KDE_BiVariate_Gaussian_Scalar, KDE_UniVariate_Gaussian; use Statistics, only:Discrete_Covariance, mean_discrete; use Common_Functions, only: setNaN
    !--Returns the Bayesian posterior for the DM Profile Variables, as defined by Mass_Profile
    !-For Mass_Profile = 1 (Flat): Sigma_0
    !-                   2 (SIS) : Velocity_Dispersion**2
    !-                   3 (NFW) : Virial Radius (r200)
    !--Prior is the prior distribution of p(m,R), or possibly p(m,lnR), m is *apparent magnitude*, R is *apparent size*
    !--Cat entered already assumed to be only the galaxies that are to be considered (e.g. in aperture), and well-described by posterior (e.g. same mag bin) 
    !--Mag-Prior entered seperately as in use we would require that no size cuts are used in the construction of this prior (including with KDE smoothing).

    !---TO DO:

    type(Catalogue)::Cat
    integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS, 3: NFW-!
    real(double),intent(in)::Lens_Redshift
    real(double),intent(in)::Lens_Position(2)
    !--Posterior(1,:) MUST contain the grid on which the posterior is to be evaluated
    real(double),allocatable,intent(InOut)::Posterior(:,:) !-Grid/Posterior, Value-!
    character(*), intent(in)::Output_Prefix
    logical,intent(in)::lnSize_Prior
    real(double),intent(in),optional:: PriorSizeGrid(:), PriorMagGrid(:), Prior(:,:), TwoD_MagPrior_Distribution(:,:) !-Magnitude, Size-!
    real(double), intent(inout), optional, allocatable:: MagPrior(:)

    !---INTERNAL DECLARATIONS
    integer::i, c, j, z, m
    integer:: Galaxy_Posterior_Method

    real(double),dimension(:,:),allocatable:: Survey_Renormalised_Prior
    real(double),dimension(:),allocatable:: Survey_Renormalised_MagPrior
    real(double),allocatable::Posterior_perGalaxy(:) !-Galaxy, Posterior-! - Uses Same Grid as overall Posterior -

    !--Cluster Model Delcarations
    real(double),allocatable::Sigma_Crit(:,:)
    real(double)::D_l, D_s, D_ls
    real(double),allocatable::Distance_from_Mass_Center(:) !-Galaxy-!

    logical:: Output_Posterior_Per_Galaxy = .true.
    character(7)::fmtstring
    character(500)::Filename
    logical::here

    character(3), allocatable:: Posterior_Flag(:,:)

    !--Redshift Distribution Declarations--!
!!$    integer,parameter:: nRedshift_Sampling = 50
!!$    real(double),parameter::Redshift_Lower = Lower_Redshift_Cut, Redshift_Higher = 4.e0_double 
    real(double), allocatable:: RedshiftGrid(:)

    real(double),dimension(size(Cat%RA),2):: Source_Positions

    !--Kappa dependant renormalisation--!
!DELETE    real(double),dimension(2):: Renormalisation_Magnitude_Limits, Renormalisation_Size_Limits
    real(double),allocatable:: MagnificationGrid(:), Renormalisation_by_Magnification(:,:), Convergence_Renorm_PerGalaxy(:,:),  MagOnly_Renormalisation_by_Magnification(:,:)
    !vv Must be set here vv!
    integer:: nMagnificationGrid = 2000
    real(double):: MagFactorGridLower = 1.e0_double, MagFactorGridHigher = 65.e0_double
    integer:: IntegrationFlag = -1000

!    integer:: nMagPosterior, nSizePosterior, nSizeMagPosterior

    !--Maximum Search Options
    logical:: Continue_To_Evaluate
    integer:: Alpha_Loop
    logical:: Search_For_Maximum = .true.
    !~~ Initial Grid Size and Tolerance to which the maximum are found set the error on confidence limits, and mode respectively, as well as setting minimum time for run
    integer:: Search__Coarse_Grid_Size = 100
    real(double):: Find_Maximum_Tolerance = 1.e-2_double

    !--Testing Declarations--!
    integer:: n_Default_Source_Redshift_Used, nGal_Ignored_MagLimits, nGal_Ignored_SizeLimits, nGal_Ignored_NaN

    INTERFACE
       subroutine DM_Profile_Variable_Posterior_SingleFit(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorMagGrid, PriorSizeGrid, Prior, MagPrior, TwoD_MagPrior_Distribution)
         use Param_Types; use Catalogues
         type(Catalogue)::Cat
         integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS, 3:NFW-!
         real(double),intent(in)::Lens_Redshift
         real(double),intent(in)::Lens_Position(2)
         real(double),allocatable,intent(Out)::Posterior(:,:) !-Grid/Posterior, Value-!
         character(*), intent(in)::Output_Prefix
         logical,intent(in)::lnSize_Prior

         real(double),intent(in),optional:: PriorSizeGrid(:), PriorMagGrid(:), Prior(:,:), TwoD_MagPrior_Distribution(:,:) !-Magnitude, Size-!
         real(double), intent(inout), optional,allocatable:: MagPrior(:)
       END subroutine DM_Profile_Variable_Posterior_SingleFit
    END INTERFACE
    
    print *, 'Called DM_Profile_Variable_Posterior'
    
    if(Analyse_with_Physical_Sizes) STOP 'DM_Profile_Variable_Posterior - I HAVE DISABLED THE ABILITY TO USE PHYISCAL SIZES AS UNNECESSARY, code still to be edited'
    
!!$    if((present(PriorMagGrid) == .false.) .or. (present(PriorSizeGrid)== .false.)) STOP 'DM_Profile_Variable_Posterior - Prior must be accompanied by grids'
!!$    allocate(Survey_Renormalised_Prior(size(Prior,1),size(Prior,2))); Survey_Renormalised_Prior = 0.e0_double
!!$    if(present(MagPrior)) then
!!$       allocate(Survey_Renormalised_MagPrior(size(PriorMagGrid))); Survey_Renormalised_MagPrior = 0.e0_double
!!$    end if
    print *, ' '
    print *, 'Attempting Prior Interpolation without Extrapolation'
    print *, ' '
    
    !--Check posterior grid set up correctly
    if(allocated(Posterior) == .false.) STOP 'DM_Profile_Variable_Posterior - Posterior Grid MUST be entered.'

    !--Start of Posterior Routines--!
    allocate(Posterior_perGalaxy(size(Cat%RA))); Posterior_perGalaxy = 1.e-100_double

    write(*,'(A)') '--------------------------------------------------------------------------------------------------'
    write(*,'(A)',advance = 'no') 'Getting Posterior for Cluster '
    if(Posterior_Method == 1) write(*,'(A)') 'using Sizes Only'
    if(Posterior_Method == 2) write(*,'(A)') 'using Sizes and Magnitudes'
    if(Posterior_Method == 3) write(*,'(A)') 'using Magnitudes Only'
    if(Posterior_Method == 4) write(*,'(A)') 'using Size and Magnitudes, and Magnitudes Only below the size limit'
    if(Enforce_Weak_Lensing) print *, '*** Weak lensing assumptions have been enforced'

    !--Set defaults for supplementary delarations
    n_Default_Source_Redshift_Used = 0; nGal_Ignored_MagLimits = 0; nGal_Ignored_SizeLimits = 0; nGal_Ignored_NaN = 0
!!DELETE    nMagPosterior = 0; nSizePosterior = 0; nSizeMagPosterior = 0

    !--Set up precursors
    call get_Likelihood_Evaluation_Precursors(Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, RedshiftGrid, Sigma_Crit, MagnificationGrid,Renormalisation_by_Magnification, MagOnly_Renormalisation_by_Magnification, PriorSizeGrid, PriorMagGrid, MagPrior, TwoD_MagPrior_Distribution, Prior, S1_Magnitude_Limits_byCluster, S1_Size_Limits_byCluster, S2_Magnitude_limits_byCluster, S2_Size_Limits_byCluster, (/Lens_Redshift/), Lower_Redshift_Cut, Output_Prefix, use_lnSize)

    Source_Positions(:,1) = Cat%RA; Source_Positions(:,2) = Cat%Dec

    print *, 'Allocating Posterior Flag'
    allocate(Posterior_Flag(size(Posterior,2), size(Cat%RA))); Posterior_Flag = '000'
    print *, 'Done'

    i = 0; Alpha_Loop = 0
    Continue_To_Evaluate = .true.

    do while (Continue_To_Evaluate)
       Posterior_perGalaxy = setNaN()
       !--Alpha_Loop counts the number of times the Alpha value has been looped over
       Alpha_Loop = Alpha_Loop+1
       i = i+1
       
       !--Set exit strategy in normal case
       if(Search_For_Maximum == .false. .and. Alpha_Loop > size(Posterior,2)) Continue_To_Evaluate = .false.
       if(Search_For_Maximum .and. Alpha_Loop>size(Posterior,2)) then
          !--Find the Maximum. In this case, i will be set to the new grid point (added to grid), and exit case will be set when either tolerance or call limits are met 
          call  find_Maximum_by_Bisection(Posterior, i, Find_Maximum_Tolerance, Continue_To_Evaluate)
       end if
       
       !--Exit here to remove need to evaluate an extra point when maximum has been found
       if(Continue_to_Evaluate == .false.) exit

       Posterior(2,i) = 1.e0_double
       if(Alpha_Loop == Size(Posterior,2)/2) print *, 'Approximately halfway done for this Aperture..'          
       
       Posterior(2,i) = lnLikelihood_Evaluation_atVirialRadius_perSourceSample(Posterior(1,i), Cat%Posterior_Method, Mass_Profile, Lens_Position, Lens_Redshift, Cat%Sizes, Cat%MF606W, Cat%Redshift, Source_Positions, PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, Survey_Renormalised_MagPrior, Size_Limits_byCluster, Magnitude_Limits_byCluster, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, Renormalisation_by_Magnification, RedshiftGrid, Sigma_Crit(1,:), use_lnSize, Posterior_Flag(i,:), Cat%RenormalisationGroup, output_Prefix)

    end do !--End of Posterior Loop

    if(any(dabs(Posterior) > huge(1.e0_double)) .or. any(isNaN(Posterior))) then
       print *, 'NaNs or infinities found in posterior:'
       do i =1, size(Posterior,2)
          print *, Posterior(:,i)
       end do
       STOP
    end if

    !--Convert from lnP to P for posterior, renormalise so that max(P) = 1, to avoid rounding errors
    Posterior(2,:) = dexp(Posterior(2,:) - maxval(Posterior(2,:)))

    if(n_Default_Source_Redshift_Used > 0) write(*,'(A,I3,A)') '****** Used the default redshift for:', n_Default_Source_Redshift_Used, ' galaxies ********'
    if(nGal_Ignored_SizeLimits > 0) write(*,'(A,I3)') '****** Number of galxies ignored as they fell outside the survey size limits:', nGal_Ignored_SizeLimits
    if(nGal_Ignored_NaN > 0) write(*,'(A,I3)') '****** Number of galaxies ignored as they were NaNs:', nGal_Ignored_NaN
    
    
    Filename= trim(adjustl(Output_Prefix))//'Posterior_Combined_Renormalised.dat'
    open(unit = 82, file = Filename)
    write(fmtstring,'(I1)') 2
    do i =1, size(Posterior,2)
       write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior(2,i)
    end do
    close(82)
    write(*,'(A)') 'Output Combined Posterior to: ', trim(adjustl(Filename))
    

    !----Note that posterior may not evaluate exactly to flags, if maximum by bisection was used
    Filename= trim(adjustl(Output_Prefix))//'Likelihood_Evaluation_Flags.dat'
    open(unit = 83, file = Filename)
    write(83, '((A,L1))') '# Note that posterior may not evaluate exactly to flags, if maximum by bisection was used: ', Search_For_Maximum
    write(fmtstring,'(I6)') size(Posterior_Flag,2)
    do i =1, size(Posterior_Flag,1)
       write(83, '(e14.7,x,'//trim(adjustl(fmtstring))//'(A,x))') Posterior(1,i), Posterior_Flag(i,:)
    end do
    close(83)
    write(*,'(A)') 'Output Posterior Flag to: ', trim(adjustl(Filename))

    !----Output a list of galaxy identifiers for this aperture
    Filename= trim(adjustl(Output_Prefix))//'Galaxy_Identifiers.dat'
    open(unit = 83, file = Filename)
    write(83, '((A))') '# The identifier is the number of the galaxy on read in'
    write(fmtstring,'(I1)') 1
    do i =1, size(Cat%RA,1)
       write(83, '('//trim(adjustl(fmtstring))//'(I8))') Cat%Galaxy_Number(i)
    end do
    close(83)
    write(*,'(2A)') 'Output Aperture Galaxy Idendtifiers to: ', trim(adjustl(Filename))


    if(any(isNAN(Posterior(2,:)))) then
       print *, 'Any NaNs in aperture posterior?', any(isNAN(Posterior(2,:))), count(isNAN(Posterior(2,:)) == .true.)
       print *, 'Stopping'
       STOP
    END if
    
    !--Remove any information less than a tolerance - Numerical Error--!
    where(Posterior(2,:) < 1.e-12_double)
       Posterior(2,:) = 0.e0_double
    end where
    
    if(allocated(Survey_Renormalised_Prior)) deallocate(Survey_Renormalised_Prior)
    if(allocated(Survey_Renormalised_MagPrior)) deallocate(Survey_Renormalised_MagPrior)
    
    !--On Successful Completion delete Poster per galaxy as large file--!
!!$    inquire(file = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat', exist = here)
!!$    if(here == .false.) STOP "Posterior per galaxy doesn't exist, stopping before accidental deletion"
!!$    call system('rm '//trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat')
    

  end subroutine DM_Profile_Variable_Posterior_SingleFit
  

  subroutine return_Prior_Distributions(MagGrid, SizeGrid, Dist, MagDist, TwoD_MagPrior_Dist, Dir, reEvaluate, BFCat, do_KDE)
    use IO, only: readin; use Distributions, only: produce_Joint_Size_Magnitude_Distribution, produce_Magnitude_Distribution; use Integration; use Interpolaters, only: linear_Interp
    !--Returns the joint size magnitude distribution (Dist) using sample with size and magnitude information (assumed to be first catalogue element), and amgnitude-only distribution (MagDist) produced for sample with only magnitude information (assumed to be the second array element)--!
    !-- If BFCat entered, then the Distribution is calculated from the catalogue, otherwise a read in is attempted --!
    !-- BFCat is an array of catalogues. This means that the secondary [magnitude] distribution can be constructed using a seperate sample to the size-magnitude distribution, e.g. in the case where size information is not present. One may alos want to consider the case where the magnitude-only distribution is produced using the full sample (may introduce bias if this includes galaxies not in the source sample)
    !--Dist is Size-Magnitude Distribution. MagDist is magnitude-only distribution, which may differ from the projected Dist due to size cuts on the prior. Magnitude_Distribution is constructed from a catalogue which only cuts between entered prior magnitude limits--!
    !--The Main outputs for this routine are Dist and TwoD_MagPrior_Dist, which correspond to the Size-Magnitude Distributions for Sample 1 and Sample 2 respectively. Previously, only the 1D magntiude distribution was output for sample 2, however it may be the case that the renormalisation for the magnitude case requires that size shifts are taken into account when sample 2 is constructed using a small size cut, even if only the magnitudes are used. The MagDist is therefore defunct, as it is recalcualted in the evualtion of the precursors on run time, but it has been left in here for now

    !!!!--- KNOWN BUG: Produces a SEGMENTATION FAULT if BFCat is not passed in. I have not got to the bottom of why, however as long as a catalogue is passed in, the re-evaluation of the prior can be avoided by settign reEvaluation to false

    !!!!--- POSSIBLE BUG: Distribution read in assumes that Mag Grid for the Size-Magnitude and Mag-Only Distributions are the same. This should be the case when this routine is used to get the distributions, however this may not always be the case.

    real(double),intent(out),allocatable:: MagGrid(:), SizeGrid(:), Dist(:,:), MagDist(:), TwoD_MagPrior_Dist(:,:)
    character(*),intent(in):: Dir
    logical:: reEvaluate
    type(Catalogue),intent(in),optional, dimension(:):: BFCat(:)
    logical,optional:: do_KDE

    character(500):: Filename, Output_Filename, Input_Filename, MagFilename, TDMagFilename
    real(double),allocatable:: Input_Array(:,:)

    character(5)::fmtstring
    integer::i,j
    logical:: iDo_KDE

    real(double),allocatable:: tMagGrid(:), tMagDist(:)
   
    real(double):: MLimit(2)
 
    !--Convergence Buffer used to extend prior slightly beyond maximum/minimum values of the sample
    real(double):: Mag_Limit_Convergence_Buffer = 0.3e0_double

    !--Cuts on the prior distribution--!
    real(double):: Redshift_Cuts(2)  = (/Lower_Redshift_Cut, 100.e0_double/) !-- If no cuts, then lower should still be < -1 to ensure only galaxies without redshift information are cut
    type(Catalogue)::Cut_Catalogue(size(BFCat))
    
   
    !--TESTING--!
    real(double),allocatable::Size_Only(:), Mag_Only(:)
    logical:: here
    

    INTERFACE
       subroutine return_Prior_Distributions(MagGrid, SizeGrid, Dist, MagDist, TwoD_MagPrior_Dist, Dir, reEvaluate, BFCat, do_KDE)
         use Catalogues, only: Catalogue; use Param_Types
         real(double),intent(out),allocatable:: MagGrid(:), SizeGrid(:), Dist(:,:), MagDist(:), TwoD_MagPrior_Dist(:,:)
         character(*),intent(in):: Dir
         logical, intent(in):: reEvaluate

         type(Catalogue),intent(in),optional,dimension(:):: BFCat
         logical,optional:: do_KDE
       end subroutine return_Prior_Distributions
    END INTERFACE

    write(*,'(A)') '_______________________________________________________________ Getting Priors _____________________________________________________________________________'
    print *, ' '

    print *, 'Started Distribution production subroutien'

    if(present(do_KDE)) then
       iDo_KDE = do_KDE
    else
       iDO_KDE = use_KDE_Smoothed_Distributions
    end if


 
    if(iDO_KDE) then
       Filename = 'Sample1_MagSize_Prior_Distribution_KDE.dat'
       MagFilename = 'Sample2_Magnitude_Prior_Distribution_KDE.dat'
       TDMagFilename = 'Sample2_MagSize_Prior_Distribution_KDE.dat'
    else
       Filename= 'Sample1_MagSize_Prior_Distribution_Histogram.dat'
       MagFilename = 'Sample2_Magnitude_Prior_Distribution_Histogram.dat'
       TDMagFilename = 'Sample2_MagSize_Prior_Distribution_Histogram.dat'
    end if


    if(reEvaluate) then

       if(present(BFCat) == .false.) then
          STOP 'return_Prior_Distributions - FATAL ERROR - Field Catalogue must be passed in to produce a-priori size and amgnitude distributions.'
       end if

       if(size(BFCat) == 0 .or. size(BFCat) > 2) then
          print *, 'return_Prior_Distributions: Blank-Field-Catalogue (BFCat) does not have the correct dimensions, FATAL. Size: ', size(BFCat)
          STOP
       END if



       print *, ' '
       print *, '--Producing Joint Size Magnitude Distribution from Catalogue.....'

       print *, 'Check on Field Catalogue:'
       print *, 'Size Catalogue:', Catalogue_Constructed(BFCat(1))
       print *, 'Size Range:', minval(BFCat(1)%Sizes), maxval(BFCat(1)%Sizes)
       print *, 'MF606W Range:', minval(BFCat(1)%MF606W), maxval(BFCat(1)%MF606W)
       print *, 'Magnitude Catalogue:', Catalogue_Constructed(BFCat(2))
       print *, 'Size Range:', minval(BFCat(2)%Sizes), maxval(BFCat(2)%Sizes)
       print *, 'MF606W Range:', minval(BFCat(2)%MF606W), maxval(BFCat(2)%MF606W)
       print *, ' '

       Cut_Catalogue = BFCat

       !--Removed: Should be done earlier
!!$       print *, '------Applying general cuts on the Catalogue from which the prior is constructed (Magnitude, Prior Size Cuts, Redshift)---------'
!!$       do i = 1, size(Cut_Catalogue)
!!$          !call Cut_By_Magnitude(Cut_Catalogue(i), Prior_Magnitude_Limits(1), Prior_Magnitude_Limits(2))
!!$          !call Cut_By_PixelSize(Cut_Catalogue(i), Prior_Size_Limits(1), Prior_Size_Limits(2))
!!$          call Cut_by_PhotometricRedshift(Cut_Catalogue(i), Redshift_Cuts(1), Redshift_Cuts(2))
!!$       end do
!!$       print *, '----- Finished Cuts on BF catalogue----------------------------------------------'


       !--Produce Joint Size Magnitude Distribution from the first catalogue
       MLimit = (/ minval((/(minval(Cut_Catalogue(i)%MF606W-2.17e0_double*Mag_Limit_Convergence_Buffer), i = 1, size(Cut_Catalogue))/)) , maxval((/(maxval(Cut_Catalogue(i)%MF606W+2.17e0_double*Mag_Limit_Convergence_Buffer), i = 1, size(Cut_Catalogue))/))/)
       !--Create Directory for intermediate storage--!
       inquire(directory = trim(Dir)//'Sample1_Prior/', exist = here)
       if(here == .false.) call system('mkdir '//trim(Dir)//'Sample1_Prior/')
       !---------------------------------------------!
       print *, 'Calling produce_Joint_Size_Magnitude_Distribution'
       call produce_Joint_Size_Magnitude_Distribution(SizeGrid, MagGrid, Dist, Cut_Catalogue(1), use_Physical_Sizes = Analyse_with_Physical_Sizes, Magnitude_Type = 2, Output_Dir = trim(Dir)//'Sample1_Prior/', ln_size_Distribution = .false., KDE_Smooth = ido_KDE, MagLimits = MLimit )
       print *, 'Finished Call'

       !--Output (matches the method of output used above)--!
       Output_Filename = trim(Dir)//trim(Filename)
       open(unit = 45, file = Output_Filename)
       
       write(fmtstring, '(I5)') size(MagGrid) + 1
       write(45, '('//trim(fmtstring)//'(e14.7,x))') 0.0e0_double, MagGrid
       do i= 1, size(SizeGrid)
          write(45, '('//trim(fmtstring)//'(e14.7,x))') SizeGrid(i), Dist(:,i)
       end do
       close(45)
       write(*,'(2(A))') 'Output to: ', trim(Output_Filename)
    
       !----Construct the prior magnitude distribution by integrating it between the entered entered size limits for the prior, NOT the survey limits
       print *, ' '
       write(*,'(A)') 'Beginning Magnitude-Only Distribution Construction:'
       if((size(Cut_Catalogue) >=2) .and. Catalogue_Constructed(Cut_Catalogue(2)) .and. (size(Cut_Catalogue(2)%RA) > 0)) then
          write(*,'(A)') 'Producing magnitude distribution using Mag-Only sub-Catalogue (Note: automatically uses KDE Smoothing'
          !--Two Dimsensional Case--!
          inquire(directory = trim(Dir)//'MagPrior/', exist = here)
          if(here == .false.) call system('mkdir '//trim(Dir)//'MagPrior/')

          !-- ** 2D case, usefult where the mag-only sample comes from size cuts. Disabeld as takes too much time
          !--call produce_Joint_Size_Magnitude_Distribution(SizeGrid, MagGrid, TwoD_MagPrior_Dist, Cut_Catalogue(2), use_Physical_Sizes = Analyse_with_Physical_Sizes, Magnitude_Type = 2, Output_Dir = trim(Dir)//'MagPrior/', ln_size_Distribution = .false., KDE_Smooth = ido_KDE, MagLimits = MLimit )
          allocate(TwoD_MagPrior_Dist(size(MagGrid), size(SizeGrid))); TwoD_MagPrior_Dist = 0.e0_double

          !--One-Dimensional Case--!
          !--What if MagGrid of Size-Magnitude Distribution does not cover this range - Should be ok with MagLimits set correctly above-!
          call produce_Magnitude_Distribution(MagGrid, MagDist, Cut_Catalogue(2), KDE_Smooth = .true.) !--Always KDE Smooth this distribution, since it is quick
       else
          print *, 'Producing magnitude distribution by integrating Size-Magnitude Distribution'
          allocate(MagDist(size(MagGrid))); MagDist = 0.e0_double
          do i= 1, size(SizeGrid)
             MagDist(i) = Integrate(SizeGrid, Dist(i,:), 2, lim = (/minval(SizeGrid), maxval(SizeGrid)/))
          end do
       end if

       !--2D Output (matches the method of output used above)--!
       Output_Filename = trim(Dir)//trim(TDMagFilename)
       open(unit = 45, file = Output_Filename)
       
       write(fmtstring, '(I5)') size(MagGrid) + 1
       write(45, '('//trim(fmtstring)//'(e14.7,x))') 0.0e0_double, MagGrid
       do i= 1, size(SizeGrid)
          write(45, '('//trim(fmtstring)//'(e14.7,x))') SizeGrid(i), TwoD_MagPrior_Dist(:,i)
       end do
       close(45)
       write(*,'(2(A))') 'Output to: ', trim(Output_Filename)

       !--1D Output (matches the method of input used below)--!
       Output_Filename = trim(Dir)//trim(MagFilename)
       open(unit = 46, file = Output_Filename)
       do i= 1, size(magGrid)
          write(46, '(2(e14.7,x))') magGrid(i), MagDist(i)
       end do
       close(46)
       write(*,'(2(A))') 'Output to: ', trim(Output_Filename)
       print *, ' '

       do i = 1, size(Cut_Catalogue)
          call Catalogue_Destruct(Cut_Catalogue(i))
       end do       

    else !---Prior Read in from file---!
       !--Read in 2D sample 1 distribution
       Input_Filename = trim(Dir)//trim(Filename)

       write(*,'(2A)') 'Reading In Size-Magnitude Distribution from:', trim(Input_Filename)

       inquire(file = Input_Filename, exist = here)
       if(here == .false.) STOP 'return_Prior_Distributions - Size-Magnitude Distribution File does not exist, stopping...'

       call ReadIn(Input_Array, filename  = trim(adjustl(Input_Filename)), tabbed = .false., header_label = '#')

       !--Output of ReadIn is (Col, Row)
       allocate(MagGrid(size(Input_Array,1)-1)); MagGrid = 0.e0_double
       allocate(SizeGrid(size(Input_Array,2)-1)); SizeGrid = 0.e0_double
       allocate(Dist(size(MagGrid), size(SizeGrid))); Dist = 0.e0_double

       MagGrid = Input_Array(2:,1)
       SizeGrid = Input_Array(1,2:)
       Dist = Input_Array(2:,2:)

       deallocate(Input_Array)

       !--Read in 2D sample 2 distribution
       Input_Filename = trim(Dir)//trim(TDMagFilename)

       write(*,'(2A)') 'Reading In Sample 2 Size-Magnitude Distribution from:', trim(Input_Filename)

       inquire(file = Input_Filename, exist = here)
       if(here == .false.) STOP 'return_Prior_Distributions - Size-Magnitude Distribution File does not exist, stopping...'

       call ReadIn(Input_Array, filename  = trim(adjustl(Input_Filename)), tabbed = .false., header_label = '#')
       !--Output of ReadIn is (Col, Row)

       !---ASSUMES THAT BOTH SAMPLE 1 AND SAMPLE 2 ARE ON THE SAME GRID
       allocate(TwoD_MagPrior_Dist(size(MagGrid), size(SizeGrid))); TwoD_MagPrior_Dist = 0.e0_double

       TwoD_MagPrior_Dist = Input_Array(2:,2:)

       deallocate(Input_Array)

       !---Read in 1D mag distribution for Sample 2
       Input_Filename =  trim(Dir)//trim(MagFilename)
       call ReadIn(Input_Array, filename  = trim(adjustl(Input_Filename)), tabbed = .false., header_label = '#')

       !--Assume Mag grid is the same as Size-Magnitude
       print *, ' '
       write(*,'(A)') 'NOTE: It is assumed that Mag-Only and Size-Mag distributions are on the same magnitude grid. I DO NOT CHECK FOR THIS'
       print *, ' '
       allocate(MagDist(size(Input_Array,2))); MagDist = Input_Array(2,:)

       if( (Input_Array(1,1) /= MagGrid(1)) .or. (Input_Array(1,size(Input_Array,2)) /= MagGrid(size(MagGrid))) .or. (size(Input_Array,2) /= size(MAgGrid))) then
          PRINT *, ' '
          write(*,'(A)') 'WARNING: Mag-Only Distribution may not be on the same grid as the Size_Mag Distribution on read in. This may introuce errors. Check this.'
          write(*,'(2A)') 'Mag_Only Distribution:', trim(Dir)//trim(Filename)
          write(*,'(2A)') 'Size-Mag Distribution:', trim(Dir)//trim(MagFilename)
          print *, ' '
       end if

       deallocate(Input_Array)

!!$       !--Magnitude Distribution-!
!!$       print *, 'Producing magnitude Distribution using Integration Method'
!!$       allocate(MagDist(size(MagGrid))); MagDist = 0.e0_double
!!$       do i= 1, size(SizeGrid)
!!$          MagDist(i) = Integrate(SizeGrid, Dist(i,:), 2, lim = (/minval(SizeGrid), maxval(SizeGrid)/))
!!$       end do
       
    end if

!    print *, 'Testing of Distribution Input:'
    allocate(Size_Only(size(SizeGrid))); Size_Only = 0.e0_double
    allocate(Mag_Only(size(MagGrid))); Mag_Only = 0.e0_double
    open(unit = 23, file = trim(Dir)//'Integrated_SizeOnlyDist.dat')
    open(unit =24, file = trim(Dir)//'Integrated_MagOnlyDist.dat')
    do i= 1, size(Size_Only)
       Size_Only(i) = Integrate(MagGrid, Dist(:,i),2, lim = (/minval(MagGrid), maxval(MagGrid)/))
       write(23,*) SizeGrid(i), Size_Only(i)
    end do
    do i= 1, size(Mag_Only)
       Mag_Only(i) = Integrate(SizeGrid, Dist(i,:),2, lim = (/minval(SizeGrid), maxval(SizeGrid)/))
       write(24,*) MagGrid(i),Mag_Only(i)
    end do
    close(23)
    close(24)

    write(*,'(A)') '_______________________________________________________________ Successfully got Priors _____________________________________________________________________________'
    print *, ' '

  end subroutine return_Prior_Distributions


  integer function Nearest_Neighbour_Index(Grid, Val)
    !--Returns the index of the point of the grid which is closest to the value - this requires the grid to be well sampled WHICH CAN'T BE TESTED--!
    !--Assumes ORDERED GRID--!
    real(double), intent(in):: Grid(:), Val

    integer:: i, index
    real(double)::Dist(2)

    !--Find indexes between which the value lies
    do i =1, size(Grid)
       if( (Grid(i) <= Val) .and. (Grid(i+1) > Val) ) then
          index = i
          exit
       end if
    end do

    !--Calculate Dist
    Dist = (/Val-Grid(index),Grid(index+1)-Val/)

    if(Dist(1) < Dist(2)) then
       Nearest_Neighbour_Index = index
    else
       Nearest_Neighbour_Index = index + 1
    end if

  end function Nearest_Neighbour_Index
    

  real(double) function Linear_Interpolation(Ref_Grid, Ref_Value, Grid_Wanted)
    real(double), intent(in)::Ref_Grid(:), Ref_Value(:)
    real(double),intent(in)::Grid_Wanted

    integer::i

    if(Grid_Wanted > maxval(Ref_Grid)) then
!       print *, 'Linear_Interpolation - Fatal Error - Grid value wanted is larger than refernce grid:', Grid_Wanted, maxval(Ref_Grid)
       Linear_Interpolation = 0.e0_double
!       STOP
    end if
    if(Grid_Wanted < minval(Ref_Grid)) then
!       print *,'Linear_Interpolation - Fatal Error - Grid value wanted is smaller than refernce grid, extrapolation needed:',Grid_Wanted, minval(Ref_Grid)
       Linear_Interpolation = 0.e0_double
!       STOP
    end if

    do i =1 , size(Ref_Grid)-1
       if(Ref_Grid(i) == Grid_Wanted) then
          Linear_Interpolation = Ref_Value(i)
          exit
       elseif(Ref_Grid(i+1) == Grid_Wanted) then
          Linear_Interpolation = Ref_Value(i+1)
          exit
       elseif( (Ref_Grid(i) < Grid_Wanted) .and. (Ref_Grid(i+1) > Grid_Wanted) ) then 
          Linear_Interpolation = Ref_Value(i) + ((Ref_Value(i+1)-Ref_Value(i))/(Ref_Grid(i+1)-Ref_Grid(i)))*(Grid_Wanted-Ref_Grid(i))
          exit
       end if
    end do

  end function Linear_Interpolation

end module Bayesian_Routines


!---------------POTENTIALLY OBSOLETE CODE---------------------------------------------------------------------------------------------------------------------------------------------------------------------!

!--10 Sept 2014 
!----DM_Profile_Variable_Posterior - Copy kept with sets up a grid on which the posterior is to be evaluated. Next step will edit to allow for a grid to be passed in, and will delete unecessary code, INCLUDING ANY KDE COVARAINCE

!!$  subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorCatalogue, PriorMagGrid, PriorSizeGrid, Prior, MagPrior)
!!$    use cosmology, only:angular_diameter_distance_fromRedshift; use MC_Redshift_Sampling, only: Monte_Carlo_Redshift_Sampling_SigmaCritical; use Mass_Profiles; use Distributions, only: ch08_redshift_distribution_Array, CH08_Redshift_Distribution_Scalar; use Interpolaters, only: Linear_Interp
!!$    use Integration, only:TrapInt, Integrate; use Matrix_methods, only: Determinant, Matrix_Invert; use Smoothing, only: KDE_BiVariate_Gaussian_Scalar, KDE_UniVariate_Gaussian; use Statistics, only:Discrete_Covariance, mean_discrete; use Common_Functions, only: setNaN
!!$    !--Returns the Bayesian posterior for the DM Profile Variables, as defined by Mass_Profile
!!$    !-For Mass_Profile = 1 (Flat): Sigma_0
!!$    !-                   2 (SIS) : Velocity_Dispersion**2
!!$    !-                   3 (NFW) : Virial Radius (r200)
!!$    !--Prior is the prior distribution of p(m,R), or possibly p(m,lnR), m is *apparent magnitude*, R is *apparent size*
!!$    !--Cat entered already assumed to be only the galaxies that are to be considered (e.g. in aperture), and well-described by posterior (e.g. same mag bin) 
!!$    !--Mag-Prior entered seperately as in use we would require that no size cuts are used in the construction of this prior (including with KDE smoothing).
!!$
!!$    !---TO DO:
!!$
!!$    type(Catalogue)::Cat
!!$    integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS, 3: NFW-!
!!$    real(double),intent(in)::Lens_Redshift
!!$    real(double),intent(in)::Lens_Position(2)
!!$    !--Posterior(1,:) MUST contain the grid on which the posterior is to be evaluated
!!$    real(double),allocatable,intent(InOut)::Posterior(:,:) !-Grid/Posterior, Value-!
!!$    character(*), intent(in)::Output_Prefix
!!$    logical,intent(in)::lnSize_Prior
!!$    type(Catalogue), intent(in), optional:: PriorCatalogue
!!$    real(double),intent(in),optional:: PriorSizeGrid(:), PriorMagGrid(:), MagPrior(:), Prior(:,:) !-Magnitude, Size-!
!!$
!!$    !---INTERNAL DECLARATIONS
!!$    integer::i, c, j, z, m
!!$    integer:: Galaxy_Posterior_Method
!!$
!!$    !-Size_Only_Prior contains the priors which depends only on size, which is evalutaed for each redshift and integrates over all magnitudes
!!$    real(double),dimension(:),allocatable:: Size_Only_Prior, Mag_Prior, Mag_Only_Prior
!!$    !--Size_Only_Mag_Prior contains the prior for size which evaluated over the mag grid, which will be integrated over. Contains p_[theta_0, m_0|z]*p[z|m_0] for all m in grid
!!$    real(double),dimension(:,:),allocatable::  Kappa_Renormalised_Prior, Survey_Renormalised_Prior, Size_Only_Mag_Prior
!!$    real(double),dimension(:),allocatable::Kappa_Renormalised_MagPrior, Survey_Renormalised_MagPrior
!!$    real(double),allocatable::Posterior_perGalaxy(:) !-Galaxy, Posterior-! - Uses Same Grid as overall Posterior -
!!$    real(double):: Effective_Magnification !- Only Model Dependant
!!$
!!$    logical:: Marginalise_Redshift_Distribution = .true.
!!$    logical::Known_Redshift
!!$
!!$    !--KDE_Smoothing Declarations--!
!!$    real(double),allocatable:: KDE_Gaussian_Covariance(:,:), Data_Vectors(:,:), KDE_Covariance_Inverse(:,:)
!!$    real(double):: KDE_Covariance_Determinant, KDE_Gaussian_Covariance_Reduction = 0.01e0_double !-How much is sig^2 which give KDE width reduced from measured covariance?--!
!!$    logical:: do_KDE_Extrapolation, do_KDE_OnTheFly
!!$    logical:: need_Extrapolate
!!$
!!$    real(double),allocatable::Sigma_Crit(:), Sigma_Crit_MC(:)
!!$    real(double):: Galaxy_Sigma_Critical
!!$    real(double)::D_l, D_s, D_ls
!!$    real(double),allocatable::Distance_from_Mass_Center(:) !-Galaxy-!
!!$
!!$    real(double)::Renorm
!!$
!!$    logical:: Output_Posterior_Per_Galaxy = .true.
!!$    character(7)::fmtstring
!!$    character(500)::Filename
!!$    logical::here
!!$
!!$    !--Redshift Distribution Declarations--!
!!$    integer,parameter:: nRedshift_Sampling = 50
!!$    real(double),parameter::Redshift_Lower = Lower_Redshift_Cut, Redshift_Higher = 4.e0_double !!!Edit to Lens_Redshift
!!$    real(double), dimension(nRedshift_Sampling):: RedshiftGrid
!!$    real(double)::RedshiftPDF
!!$    real(double),allocatable:: Posterior_perGalaxy_Redshift(:) !-Galaxy, Posterior, Redshift-!
!!$
!!$    !--Kappa dependant renormalisation--!
!!$!!    real(double),dimension(2):: Survey_Magnitude_Limits, Survey_Size_Limits !--User defined needs edit to below, these are set by default from distribution-!
!!$    real(double),dimension(2):: Renormalisation_Magnitude_Limits, Renormalisation_Size_Limits
!!$    real(double),allocatable:: MagnificationGrid(:), Renormalisation_by_Magnification(:), Convergence_Renorm_PerGalaxy(:,:),  MagOnly_Renormalisation_by_Magnification(:)
!!$    !vv Must be set here vv!
!!$    integer:: nMagnificationGrid = 2000
!!$    real(double):: MagFactorGridLower = 1.e0_double, MagFactorGridHigher = 65.e0_double
!!$    integer:: IntegrationFlag = -1000
!!$
!!$    integer:: nMagPosterior, nSizePosterior, nSizeMagPosterior
!!$
!!$    !--Maximum Search Options
!!$    logical:: Continue_To_Evaluate
!!$    integer:: Alpha_Loop
!!$    logical:: Search_For_Maximum = .true.
!!$    !~~ Initial Grid Size and Tolerance to which the maximum are found set the error on confidence limits, and mode respectively, as well as setting minimum time for run
!!$    integer:: Search__Coarse_Grid_Size = 100
!!$    real(double):: Find_Maximum_Tolerance = 1.e-2_double
!!$
!!$    !--Testing Declarations--!
!!$    real(double),dimension(3,size(Cat%RA)):: Convergence_per_Cluster
!!$    integer:: n_Default_Source_Redshift_Used, nGal_Ignored_MagLimits, nGal_Ignored_SizeLimits, nGal_Ignored_NaN
!!$    real(double), allocatable:: Aperture_Smoothed_Size_PDF(:), Aperture_Smoothed_Mag_PDF(:)
!!$    logical:: Produce_Shift_inAperture = .false.
!!$    real::Time1, Time2, Time3, Time4
!!$    integer::test
!!$
!!$    INTERFACE
!!$       subroutine DM_Profile_Variable_Posterior(Cat, Mass_Profile, Lens_Redshift, Lens_Position, Posterior, Output_Prefix, lnSize_Prior, PriorCatalogue, PriorMagGrid, PriorSizeGrid, Prior, MagPrior)
!!$         use Param_Types; use Catalogues
!!$         type(Catalogue)::Cat
!!$         integer,intent(in)::Mass_Profile !-1:Flat, 2:SIS, 3:NFW-!
!!$         real(double),intent(in)::Lens_Redshift
!!$         real(double),intent(in)::Lens_Position(2)
!!$         real(double),allocatable,intent(Out)::Posterior(:,:) !-Grid/Posterior, Value-!
!!$         character(*), intent(in)::Output_Prefix
!!$         logical,intent(in)::lnSize_Prior
!!$
!!$         type(Catalogue), intent(in), optional:: PriorCatalogue
!!$         real(double),intent(in),optional:: PriorSizeGrid(:), PriorMagGrid(:), MagPrior(:), Prior(:,:) !-Magnitude, Size-!
!!$       END subroutine DM_Profile_Variable_Posterior
!!$    END INTERFACE
!!$    
!!$    print *, 'Called DM_Profile_Variable_Posterior'
!!$    
!!$    if(Analyse_with_Physical_Sizes) STOP 'DM_Profile_Variable_Posterior - I HAVE DISABLED THE ABILITY TO USE PHYISCAL SIZES AS UNNECESSARY, code still to be edited'
!!$    
!!$    do_KDE_Extrapolation = .false.; do_KDE_OnTheFly = .false.
!!$    if(present(Prior) .or. present(MagPrior)) then
!!$       if((present(PriorMagGrid) == .false.) .or. (present(PriorSizeGrid)== .false.)) STOP 'DM_Profile_Variable_Posterior - Prior must be accompanied by grids'
!!$       allocate(Kappa_Renormalised_Prior(size(Prior,1), size(Prior,2))); Kappa_Renormalised_Prior = 0.e0_double
!!$       allocate(Survey_Renormalised_Prior(size(Prior,1),size(Prior,2))); Survey_Renormalised_Prior = 0.e0_double
!!$       if(present(MagPrior)) then
!!$          allocate(Survey_Renormalised_MagPrior(size(PriorMagGrid))); Survey_Renormalised_MagPrior = 0.e0_double
!!$          allocate(Kappa_Renormalised_MagPrior(size(PriorMagGrid))); Kappa_Renormalised_MagPrior =0.e0_double
!!$       end if
!!$       allocate(Size_Only_Mag_Prior(size(Prior,1),size(Prior,2))); Size_Only_Mag_Prior = 0.e0_double
!!$       if(allow_KDE_Extrapolation .and. present(PriorCatalogue)) then
!!$          do_KDE_Extrapolation = .true.
!!$          print *, ' '
!!$          print *, 'Attempting Prior Interpolation with KDE Extrapolation'
!!$          print *, ' '
!!$       else
!!$          print *, ' '
!!$          print *, 'Attempting Prior Interpolation without Extrapolation'
!!$          print *, ' '
!!$       end if
!!$    elseif(present(PriorCatalogue)) then
!!$       do_KDE_OnTheFly = .true.
!!$       print *, ' '
!!$       print *, 'Attempting KDE on the Fly - NOTE This does not include Kappa-Renormalisation (Seg fault straight after this..)'
!!$       print *, ' '
!!$    end if
!!$
!!$    !-Set Up Posterior Grid-!
!!$!!!$    select case(Mass_Profile)
!!$!!!$    case(1) !-Flat-!
!!$!!!$       nGrid = 100000
!!$!!!$       VGrid_Lower = -5.e-3_double; VGrid_Higher = 1.e-2_double !-- SMD ~ Masses 10^12 -> 10^15 Msun/h, in Units of 10^18 Msun/h  
!!$!!!$    case(2) !-SIS-!
!!$!!!$       nGrid = 10000 !--This needs to be sigma_v^2--!
!!$!!!$       VGrid_Lower= -2.5e5_double; VGrid_Higher = 9.e6_double    !--------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!$    case(3)
!!$!!!$       if(Search_for_Maximum) then
!!$!!!$          !--Set initial grid to a smaller value since the search for the maximum will reduce the number of points needed
!!$!!!$          nGrid = Search__Coarse_Grid_Size
!!$!!!$       else
!!$!!!$          nGrid = 500
!!$!!!$       end if
!!$!!!$       VGrid_Lower= 0.05e0_double; VGrid_Higher = 3.e0_double
!!$!!!$    case default
!!$!!!$       STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
!!$!!!$    end select
!!$!!!$    allocate(Posterior(2, nGrid)); Posterior = setNaN() !*!
!!$!!!$    do i =1, nGrid
!!$!!!$       Posterior(1,i) = VGrid_Lower + (i-1)*((VGrid_Higher-VGrid_Lower)/(nGrid-1))
!!$!!!$    end do
!!$    !----------------------!
!!$
!!$    if(allow_KDE_Extrapolation .or. KDE_OnTheFly) then
!!$       if(present(PriorCatalogue) == .false.) STOP 'DM_Profile_Variable_Posterior - Prior Catalogue needs to be entered to allow KDE Extrapolation'
!!$       !--Construct the Covariance that will be used for the KDE Smoothing--!
!!$       allocate(Data_Vectors(2,size(PriorCatalogue%Sizes))); Data_Vectors(1,:) = PriorCatalogue%MF606W; Data_Vectors(2,:) = PriorCatalogue%Sizes
!!$       call Discrete_Covariance(Data_Vectors, KDE_Gaussian_Covariance)
!!$       KDE_Gaussian_Covariance = KDE_Gaussian_Covariance_Reduction*KDE_Gaussian_Covariance
!!$       call Matrix_Invert(KDE_Gaussian_Covariance, KDE_Covariance_Inverse, 'S')
!!$       KDE_Covariance_Determinant = Determinant(KDE_Gaussian_Covariance)
!!$       deallocate(KDE_Gaussian_Covariance)
!!$    end if
!!$
!!$    if(produce_Shift_inAperture) then
!!$       !--Testing: output smoothed size and mag distributions within the aperture to see if there is any noticable shift
!!$       allocate(Data_Vectors(2,size(PriorCatalogue%Sizes))); Data_Vectors(1,:) = PriorCatalogue%MF606W; Data_Vectors(2,:) = PriorCatalogue%Sizes
!!$       call Discrete_Covariance(Data_Vectors, KDE_Gaussian_Covariance)
!!$       KDE_Gaussian_Covariance = KDE_Gaussian_Covariance_Reduction*KDE_Gaussian_Covariance
!!$       allocate(Aperture_Smoothed_Size_PDF(size(PriorSizeGrid))); Aperture_Smoothed_Size_PDF = 0.e0_double
!!$       allocate(Aperture_Smoothed_Mag_PDF(size(PriorMagGrid))); Aperture_Smoothed_Mag_PDF = 0.e0_double
!!$       call KDE_Univariate_Gaussian(Cat%Sizes, dsqrt(KDE_Gaussian_Covariance(2,2)), PriorSizeGrid, Aperture_Smoothed_Size_PDF)
!!$       call KDE_Univariate_Gaussian(Cat%MF606W, dsqrt(KDE_Gaussian_Covariance(1,1)), PriorMagGrid, Aperture_Smoothed_Mag_PDF)
!!$       open(unit = 31, file = trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Size.dat')
!!$       do i =1, size(PriorSizeGrid)
!!$          write(31, *) PriorSizeGrid(i), Aperture_Smoothed_Size_PDF(i)
!!$       end do
!!$       close(31)
!!$       open(unit = 31, file = trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Mag.dat')
!!$       do i =1, size(PriorSizeGrid)
!!$          write(31, *) PriorMagGrid(i), Aperture_Smoothed_Mag_PDF(i)
!!$       end do
!!$       close(31)
!!$       write(*, '(4A)') '**Output distributions in aperture to: ', trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Size.dat', ' : ', trim(adjustl(Output_Prefix))//'KDE_Distributions_in_Aperture_Mag.dat'
!!$       print *, '**with mean (Size, Mag):', mean_discrete(Cat%Sizes), mean_discrete(Cat%MF606W) 
!!$       deallocate(Aperture_Smoothed_Size_PDF, Aperture_Smoothed_Mag_PDF)
!!$    !---End of distributions in Aperture
!!$    end if
!!$
!!$    D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)
!!$    !--Set up the redshift grid--!
!!$    if(Marginalise_Redshift_Distribution) then
!!$       do z = 1, nRedshift_Sampling
!!$          RedshiftGrid(z) = Redshift_Lower + (z-1)*((Redshift_Higher-Redshift_Lower)/(nRedshift_Sampling-1))
!!$       end do
!!$    end if
!!$
!!$    !--Get Sigma_Critical for each point on the Redshift PDF grid--!
!!$    allocate(Sigma_Crit(size(RedshiftGrid))); Sigma_Crit = 0.e0_double
!!$    do z = 1, size(RedshiftGrid)
!!$       D_s = angular_diameter_distance_fromRedshift(0.e0_double, RedshiftGrid(z))
!!$       D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, RedshiftGrid(z))
!!$       Sigma_Crit(z) = 1.66492e0_double*(D_s/(D_l*D_ls)) !-(10^18 M_Sun/h)-!
!!$    end do
!!$    
!!$    !--Renormalise the prior within these size and magnitude limits--!
!!$    print *, 'Renormalisation of the intrinsic distribution:', Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)
!!$    Survey_Renormalised_Prior = Prior/Integrate(PriorMagGrid, PriorSizeGrid, Prior, 2, lim1 = Survey_Magnitude_Limits, lim2 = Survey_Size_Limits)
!!$
!!$    if(present(MagPrior)) then
!!$       !--Allow for seperate renormalisation of the magnitude prior. This is required if the magnitude prior is constructed from galaxies which are excluded from the joint size-magnitude analysis, e.g. due to size cuts--!
!!$       print *, 'Renormalisation of the intrinsic magnitude distribution:', Integrate(PriorMagGrid, MagPrior, 2, lim = Survey_Magnitude_Limits)
!!$       Survey_Renormalised_MagPrior = MagPrior/Integrate(PriorMagGrid, MagPrior, 2, lim = Survey_Magnitude_Limits)
!!$    end if
!!$    
!!$    !------This section renormalises the likelihood, taking into account kappa-dependent mag and size cuts, however this need to be implemented in the prior-----!
!!$    allocate(MagnificationGrid(nMagnificationGrid)); MagnificationGrid = 0.e0_double
!!$    allocate(Renormalisation_by_Magnification(size(MagnificationGrid))); Renormalisation_by_Magnification = 0.e0_double
!!$    if(present(MagPrior)) then
!!$       allocate(MagOnly_Renormalisation_by_Magnification(size(MagnificationGrid))); MagOnly_Renormalisation_by_Magnification = 0.e0_double
!!$    end if
!!$    allocate(Size_Only_Prior(size(Prior,2))); Size_Only_Prior = 0.e0_double
!!$    !Choose MagGrid lower to be the point where the full magnitude range is swept out:
!!$    MagFactorGridHigher = 1.05e0_double*(10.e0_double**((Survey_Magnitude_Limits(2)-Survey_Magnitude_Limits(1))/2.5e0_double)) !1.05 gives lee-way!
!!$    do i = 1, size(MagnificationGrid)
!!$       MagnificationGrid(i) = MagFactorGridLower + (i-1)*((MagFactorGridHigher- MagFactorGridLower)/(size(MagnificationGrid)-1))
!!$       Renormalisation_Size_Limits = Survey_Size_Limits/dsqrt(MagnificationGrid(i))
!!$       Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.5e0_double*dlog10(MagnificationGrid(i))
!!$
!!$       Renormalisation_by_Magnification(i) = Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
!!$       if(present(MagPrior)) then
!!$
!!$          !--Edit this code to calculate the magnitude renormalisation over the whole size grid if Posterior_Method == 3, and over (0, Survey_Size_Limit(1)) if Posterior_Method == 4
!!$          !-- vv This is only true if there are no cuts on the size-mag prior vv
!!$!MagRenorm          MagOnly_Renormalisation_by_Magnification(i) =  Integrate(PriorMagGrid, PriorSizeGrid, Survey_Renormalised_Prior, 2, lim1 = Renormalisation_Magnitude_Limits, lim2 = Renormalisation_Size_Limits)
!!$          MagOnly_Renormalisation_by_Magnification(i) = Integrate(PriorMagGrid, Survey_Renormalised_MagPrior, 2, lim = Renormalisation_Magnitude_Limits)
!!$       end if
!!$    end do
!!$    deallocate(Size_Only_Prior)
!!$       
!!$    open(unit = 53, file = trim(adjustl(Output_Prefix))//'Renormalisation_by_Magnification.dat')
!!$    do i = 1, size(Renormalisation_by_Magnification)
!!$       write(53, *) MagnificationGrid(i), Renormalisation_by_Magnification(i)
!!$    end do
!!$    write(*,'(2(A))') 'File output: ', trim(adjustl(Output_Prefix))//'Renormalisation_by_Magnification.dat'
!!$
!!$    !--Start of Posterior Routines--!
!!$    allocate(Posterior_perGalaxy(size(Cat%RA))); Posterior_perGalaxy = 1.e-100_double
!!$    Effective_Magnification = 0.e0_double
!!$    Convergence_Per_Cluster = 0.e0_double
!!$
!!$    write(*,'(A)') '--------------------------------------------------------------------------------------------------'
!!$    write(*,'(A)',advance = 'no') 'Getting Posterior for Cluster '
!!$    if(Posterior_Method == 1) write(*,'(A)') 'using Sizes Only'
!!$    if(Posterior_Method == 2) write(*,'(A)') 'using Sizes and Magnitudes'
!!$    if(Posterior_Method == 3) write(*,'(A)') 'using Magnitudes Only'
!!$    if(Posterior_Method == 4) write(*,'(A)') 'using Size and Magnitudes, and Magnitudes Only below the size limit'
!!$    if(Enforce_Weak_Lensing) print *, '*** Weak lensing assumptions have been enforced'
!!$
!!$
!!$    !--Get the distance from the mass centre for each galaxy
!!$    allocate(Distance_from_Mass_Center(size(Cat%RA))); Distance_from_Mass_Center = 0.e0_double
!!$    Distance_from_Mass_Center = dsqrt( (Cat%RA(:)-Lens_Position(1))**2.e0_double + (Cat%Dec(:)-Lens_Position(2))**2.e0_double ) !-in Degrees-!
!!$    Distance_from_Mass_Center = (D_l*Distance_from_Mass_Center*(3.142e0_double/(180.e0_double))) !-in Mpc/h-! 
!!$
!!$
!!$    n_Default_Source_Redshift_Used = 0; nGal_Ignored_MagLimits = 0; nGal_Ignored_SizeLimits = 0; nGal_Ignored_NaN = 0
!!$    nMagPosterior = 0; nSizePosterior = 0; nSizeMagPosterior = 0
!!$    Time1 = 0.; Time2 = 0.
!!$
!!$    i = 0; Alpha_Loop = 0
!!$    Continue_To_Evaluate = .true.
!!$    do while (Continue_To_Evaluate)
!!$       Posterior_perGalaxy = setNaN()
!!$       !--Alpha_Loop counts the number of times the Alpha value has been looped over
!!$       Alpha_Loop = Alpha_Loop+1
!!$       i = i+1
!!$
!!$       !--Set exit strategy in normal case
!!$       if(Search_For_Maximum == .false. .and. Alpha_Loop > size(Posterior,2)) Continue_To_Evaluate = .false.
!!$       if(Search_For_Maximum .and. Alpha_Loop>size(Posterior,2)) then
!!$          !--Find the Maximum. In this case, i will be set to the new grid point (added to grid), and exit case will be set when either tolerance or call limits are met 
!!$          call  find_Maximum_by_Bisection(Posterior, i, Find_Maximum_Tolerance, Continue_To_Evaluate)
!!$       end if
!!$
!!$       !--Exit here to remove need to evaluate an extra point when maximum has been found
!!$       if(Continue_to_Evaluate == .false.) exit
!!$
!!$       Posterior(2,i) = 1.e0_double
!!$       if(Alpha_Loop == Size(Posterior,2)/2) print *, 'Approximately halfway done for this Aperture..'          
!!$       do c = 1, size(Posterior_perGalaxy,1) !-Loop over galaxies-!
!!$          !             print *, 'Doing Galaxy Loop:', c, size(Posterior_pergalaxy,1)
!!$          
!!$          !~~Select the method of posterior reconstruction for that point based in input method   
!!$          select case(Posterior_Method)
!!$          case(1) !--SizeOnly--!
!!$             if( (Cat%Sizes(c) > Survey_Size_Limits(2)) .or. (Cat%Sizes(c) < Survey_Size_Limits(1))) then
!!$                nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
!!$                cycle
!!$             end if
!!$             !MagRenorm          iSurvey_size_Limits = Survey_Size_Limits
!!$             
!!$             nSizePosterior = nSizePosterior + 1
!!$             Galaxy_Posterior_Method = 1
!!$          case(2)!--SizeMag--!
!!$             if( (Cat%Sizes(c) > Survey_Size_Limits(2)) .or. (Cat%Sizes(c) < Survey_Size_Limits(1))) then
!!$                nGal_Ignored_SizeLimits = nGal_Ignored_SizeLimits + 1
!!$                cycle
!!$             end if
!!$             
!!$             !MagRenorm          iSurvey_size_Limits = Survey_Size_Limits
!!$             
!!$             nSizeMagPosterior = nSizeMagPosterior + 1
!!$             Galaxy_Posterior_Method = 2
!!$          case(3) !-Magnitude Only--!
!!$             nMagPosterior = nMagPosterior + 1
!!$             
!!$             !MagRenorm          iSurvey_size_Limits = (/0.e0_double, 1.e30_double/) !--Should encompase the whole data set
!!$             
!!$             if(present(MagPrior) == .false.) STOP 'Attempting to produce posterior using mag only, but no magnitde prior entered, stopping'
!!$             Galaxy_Posterior_Method = 3
!!$          case(4) !-Size mag above size data limit, Mag-Only below size data limit-!
!!$             if(present(MagPrior) == .false.) STOP 'Attempting to produce posterior using mag only under size limit, but no magnitde prior entered, stopping'
!!$             STOP 'Posterior Method 4 has been disabled as it needs further thought into its construction. In particular, with respect to the construction of the prior from the size-mag distriution, the effect of prior cuts, and correct renormalisation. The skeleton code has been added for this, but disabled for now. See commented code labeled MagRenorm'
!!$             
!!$             !--Note, in this case the data-renormalisation should take inot account that the mag only is constructed from a sample which takes small, faint galaxies
!!$             if(Cat%Sizes(c) < Survey_Size_Limits(1)) then
!!$                !MagRenorm             iSurvey_Size_Limits = (/0.e0_double, Survey_Size_Limits(1)/)
!!$                
!!$                Galaxy_Posterior_Method = 3
!!$                nMagPosterior = nMagPosterior + 1
!!$             else
!!$                !MagRenorm             iSurvey_Size_Limits = Survey_Size_Limits
!!$                
!!$                Galaxy_Posterior_Method = 2
!!$                nSizeMagPosterior = nSizeMagPosterior + 1
!!$             end if
!!$          end select
!!$          
!!$          !~~~ Cycle if the observed magnitude falls outside the survey limits
!!$          if( (Cat%MF606W(c) > Survey_Magnitude_Limits(2)) .or. (Cat%MF606W(c) < Survey_Magnitude_Limits(1))) then
!!$             nGal_Ignored_MagLimits = nGal_Ignored_MagLimits + 1
!!$             cycle
!!$          end if
!!$          if(isNaN(Cat%Sizes(c)) .or. isNaN(Cat%MF606W(c))) then
!!$             nGal_Ignored_NaN = nGal_Ignored_NaN + 1
!!$             cycle
!!$          end if
!!$          
!!$          !~~~Determine method of setting source redshift: Either Marginalise over distribution; Use Redshift of source if known; Use a default redshift
!!$          
!!$          !~~~Set Redshift PDF. If redshift is known, this is not needed. Alternatively, this could be determined externally and interpolated
!!$          Known_Redshift = .false.
!!$          !-Set to NaN as default (indicates not properly set
!!$          Galaxy_Sigma_Critical = dsqrt(-1.e0_double) 
!!$          if(Cat%Redshift(c) >= 0.e0_double) then
!!$             !--Use Galaxy Redshift if available--!
!!$             !--Ignore Galaxies with redshift less than the foreground-!
!!$             if(Cat%Redshift(c) < Lens_Redshift) cycle
!!$             
!!$             Galaxy_Sigma_Critical = Linear_Interp(Cat%Redshift(c), RedshiftGrid, Sigma_Crit)!1.66492e0_double*(D_s/(D_l*D_ls))
!!$             Known_Redshift = .true.
!!$             
!!$          elseif(Marginalise_Redshift_Distribution == .false.) then
!!$             !--Use Default Redshift if everything else fails (usually this should not be considered)--!
!!$             n_Default_Source_Redshift_Used = n_Default_Source_Redshift_Used + 1
!!$             
!!$             Galaxy_Sigma_Critical = Linear_Interp(Default_Source_Redshift, RedshiftGrid, Sigma_Crit)
!!$             Known_Redshift = .true.
!!$          elseif(Marginalise_Redshift_Distribution) then
!!$             !--Marginalising over the redshift distribution for that galaxy--!
!!$          else
!!$             STOP 'DM_Profile_Variable_Posterior - Both MC and Redshift Distribution methods set - this cannot be'
!!$          end if
!!$          if(Galaxy_Sigma_Critical < 0.e0_double) STOP 'DM_Profile_Variable_Posterior - Invalid Sigma Critical Entered, negative'
!!$          
!!$          if(Known_Redshift) then
!!$             allocate(Posterior_perGalaxy_Redshift(1)); Posterior_perGalaxy_Redshift = 0.e0_double
!!$          else
!!$             allocate(Posterior_perGalaxy_Redshift(nRedshift_Sampling)); Posterior_perGalaxy_Redshift = 0.e0_double
!!$          end if
!!$
!!$          do z = 1, size(Posterior_perGalaxy_Redshift,1)             
!!$
!!$             !--There must be a better way to do this, perhaps using interpolation? (i.e. it is very seperated from other methods
!!$             if(Marginalise_Redshift_Distribution .and. (Known_Redshift == .false.)) Galaxy_Sigma_Critical = Sigma_Crit(z)
!!$             
!!$             
!!$             if(isNaN(Galaxy_Sigma_Critical) .or. (Galaxy_Sigma_Critical > huge(1.e0_double))) STOP 'DM_Profile_Variable_Posterior - FATAL - Galaxy Sigma Critical not set correctly'
!!$             select case(Mass_Profile)
!!$             case(1) !-Flat-!
!!$                STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
!!$                Effective_Magnification = 1.e0_double+2.e0_double*(Posterior(1,i)/Galaxy_Sigma_Critical)
!!$             case(2) !-SIS-!
!!$                STOP 'The Strong Lensing approximation has not been implemented for this profile yet'
!!$                Effective_Magnification = 1.e0_double+2.e0_double*(SMD_SIS(Posterior(1,i), Distance_From_Mass_Center(c))/(Galaxy_Sigma_Critical*1.e18_double))
!!$             case(3) !-NFW-!
!!$                if(Enforce_Weak_Lensing) then
!!$                   Effective_Magnification = 1.e0_double + 2.e0_double*(SMD_NFW(Distance_From_Mass_Center(c), Lens_Redshift, Posterior(1,i))/(Galaxy_Sigma_Critical*1.e18_double))
!!$                else
!!$                   Effective_Magnification = Magnification_Factor(3, Distance_From_Mass_Center(c), Posterior(1,i),  Lens_Redshift, Galaxy_Sigma_Critical*1.e18_double)
!!$                end if
!!$             case default
!!$                STOP 'DM_Profile_Variable_Posterior - fatal error - Invalid Profile value entered'
!!$             end select
!!$
!!$             if(isNaN(Effective_Magnification)) then
!!$                print *, 'Magnification Factor is a NaN:', Distance_From_Mass_Center(c), Posterior(1,i), Lens_Redshift, Galaxy_Sigma_Critical*1.e18_double
!!$                print *, 'Loops, posterior, galaxy, redshift:', i, c, z
!!$                STOP
!!$             end if
!!$
!!$             if(Effective_Magnification < minval(MagnificationGrid) .or. Effective_Magnification > maxval(MagnificationGrid)) then
!!$                !--Skipping as outside limits on which magnification was evaluated--!
!!$                Posterior_perGalaxy_Redshift(z) = 1.e-100_double
!!$                cycle
!!$             end if
!!$
!!$             if(Effective_Magnification <= 0.e0_double) then
!!$                !--This won't be picked up unless the above cycle is relaxed--!
!!$                print *, Effective_Magnification, SMD_NFW_Scalar(Distance_From_Mass_Center(c), Lens_Redshift, Posterior(1,i)), Differential_SMD_Scalar(Distance_From_Mass_Center(c), Lens_Redshift, Posterior(1,i)), Galaxy_Sigma_Critical*1.e18_double, Posterior(1,i), Distance_From_Mass_Center(c)
!!$                STOP 'Effective Magnification invalid - negative. Stopping'
!!$             end if
!!$             
!!$             !--Get RedshiftPDF which depends on the unlensed magnitude - this could possibly be interpolated, which would hopefully be faster
!!$             if(Known_Redshift) then
!!$                RedshiftPDF = 1.e0_double
!!$             else
!!$                RedshiftPDF = CH08_Redshift_Distribution_Scalar(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification), RedshiftGrid(z))
!!$             end if
!!$
!!$             !--Construct Joint Size-Magnitde Prior which is renormalised according to the convergence value--!
!!$
!!$             Renormalisation_Size_Limits = Survey_Size_Limits/dsqrt(Effective_Magnification)
!!$             Renormalisation_Magnitude_Limits = Survey_Magnitude_Limits + 2.5e0_double*dlog10(Effective_Magnification)
!!$
!!$             if(Cuts_Renormalise_Likelihood) then
!!$                if(i == 1 .and. z ==1 .and. c == 1) print *, 'KAPPA RENORMALISATION TURNED ON'
!!$                Renorm = Linear_Interp(Effective_Magnification, MagnificationGrid, Renormalisation_by_Magnification, ExValue = 1.e30_double)
!!$                if(Renorm == 0) then
!!$                   Kappa_Renormalised_Prior = 0.e0_double
!!$                else
!!$                   Kappa_Renormalised_Prior = Survey_Renormalised_Prior/Renorm
!!$                end if
!!$                Renorm = 0.e0_double
!!$             else
!!$                if(i == 1 .and. z ==1 .and. c == 1) print *, 'KAPPA RENORMALISATION TURNED OFF'
!!$                Kappa_Renormalised_Prior =  Survey_Renormalised_Prior
!!$                Kappa_Renormalised_MagPrior = Survey_Renormalised_MagPrior
!!$             end if
!!$
!!$             select case (Galaxy_Posterior_Method)
!!$             case(1) !--Size Only--!
!!$                
!!$                !--Evaluate p_{theta_0, m_0}*p_{z|m_0} for the whole magnitude grid
!!$                do m =  1, size(PriorMagGrid)
!!$                   !--m_0, theta_0--!
!!$                   if((PriorMagGrid(m) < Renormalisation_Magnitude_Limits(1)) .or. (PriorMagGrid(m) > Renormalisation_Magnitude_Limits(2))) then
!!$                      !-Since Integrand does not extend over this region anyway-!
!!$                      Size_Only_Mag_Prior(m,:) = 1.e-100_double
!!$                   else
!!$                      if(Known_Redshift) then
!!$                         Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)
!!$                      else
!!$                         Size_Only_Mag_Prior(m,:) =  Kappa_Renormalised_Prior(m,:)*CH08_Redshift_Distribution_Scalar(PriorMagGrid(m), RedshiftGrid(z))
!!$                      end if
!!$                   end if
!!$                end do
!!$
!!$
!!$                if(lnSize_Prior) then
!!$                   STOP 'ln Size prior has been switched off since the conversion to strong lensing - check the theory carefully before use'
!!$                   if(dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification)) > maxval(PriorSizeGrid) .or. (dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification)) < minval(PriorSizeGrid))) then
!!$                      !--Extrapolation--!
!!$                      Posterior_perGalaxy_Redshift(z) = 0.e0_double
!!$                      cycle
!!$                   end if
!!$
!!$                   !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-tim by minimising the number of integrations required
!!$                   allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
!!$                   do j = 1, size(PriorSizeGrid)-1
!!$                      if( ( PriorSizeGrid(j)<= dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification))) .and. (  PriorSizeGrid(j+1) > dlog(Cat%Sizes(c)/dsqrt(Effective_Magnification)))) then
!!$!                      if( (PriorSizeGrid(j)<= dlog(Cat%Sizes(c))-Effective_Convergence(c,i)) .and. ( PriorSizeGrid(j+1) > dlog(Cat%Sizes(c))-Effective_Convergence(c,i)) ) then
!!$                         Size_Only_Prior(1) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
!!$                         Size_Only_Prior(2) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
!!$                         exit
!!$                      end if
!!$                   end do
!!$
!!$!!!$                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(dlog(Cat%Sizes(c))-Effective_Convergence(c,i), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue = 0.e0_double) !-Reinstate with SLensing--!
!!$                else
!!$                   if((Cat%Sizes(c)/dsqrt(Effective_Magnification) > maxval(PriorSizeGrid)) .or. (Cat%Sizes(c)/dsqrt(Effective_Magnification) < minval(PriorSizeGrid))) then
!!$                      !--Extrapolation--!
!!$                      Posterior_perGalaxy_Redshift(z) = 1.e-100_double
!!$                      cycle
!!$                   else
!!$                      
!!$                      !--Find Boundary Indexs for the Size at which the Prior will be evaluated - Aims to improve run-time by minimising the number of integrations required
!!$                      allocate(Size_Only_Prior(2)); Size_Only_Prior = 0.e0_double
!!$                      do j = 1, size(PriorSizeGrid)-1
!!$                         if( (PriorSizeGrid(j)<= Cat%Sizes(c)/dsqrt(Effective_Magnification)) .and. ( PriorSizeGrid(j+1) > Cat%Sizes(c)/dsqrt(Effective_Magnification)) ) then
!!$                            Size_Only_Prior(1) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j), 2, lim = Renormalisation_Magnitude_Limits)
!!$                            Size_Only_Prior(2) = Integrate(PriorMagGrid, Size_Only_Mag_Prior(:,j+1), 2, lim = Renormalisation_Magnitude_Limits)
!!$                            exit
!!$                         end if
!!$                      end do
!!$
!!$                      
!!$                      Posterior_perGalaxy_Redshift(z) = Linear_Interp(Cat%Sizes(c)/dsqrt(Effective_Magnification), PriorSizeGrid(j:j+1), Size_Only_Prior(:), ExValue =  1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification))
!!$
!!$                   end if
!!$                end if
!!$                deallocate(Size_Only_Prior)
!!$             case(2)!-Size and Magnitude-!
!!$                
!!$                if(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification) < 21.e0_double) then
!!$                   print *, 'Possible problem with de-lensed magnitude - falls outwith bright limit of Scrabback Fit (21)'
!!$                   print *, Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification)
!!$                end if
!!$                if(lnSize_Prior) then
!!$                   STOP 'DM_Profile_Variable_Posterior - lnSize with Size_Magnitude Method - I cannae do that captain!'
!!$                else
!!$                   !--Uses distributions of the apparent size--!
!!$                   !--Test for need for extrapolation---!
!!$                   need_Extrapolate = .false.
!!$                   if(( (Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification) > maxval(PriorMagGrid)) .or. (Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification) < minval(PriorMagGrid)) ) .or. ((Cat%Sizes(c)/dsqrt(Effective_Magnification) > maxval(PriorSizeGrid)) .or. (Cat%Sizes(c)/dsqrt(Effective_Magnification) > maxval(PriorSizeGrid)) )) then
!!$                      need_Extrapolate = .true.
!!$                   else
!!$                      need_Extrapolate = .false.
!!$                   end if
!!$      
!!$                   if((need_Extrapolate .and. do_KDE_Extrapolation) .or. do_KDE_OnTheFly) then !-.or. KDE_OnTheFly
!!$                      !--KDE_Extrapolation / KDE_OnTheFly(? - What about entry of prior?)
!!$                      Posterior_perGalaxy_Redshift(z) = KDE_BiVariate_Gaussian_Scalar(PriorCatalogue%MF606W, PriorCatalogue%Sizes, Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification), Cat%Sizes(c)/dsqrt(Effective_Magnification), Inverse_Covar = KDE_Covariance_Inverse, Det_Covar = KDE_Covariance_Determinant)*(1.e0_double/dsqrt(Effective_Magnification))*RedshiftPDF
!!$                   elseif(need_Extrapolate .and. (do_KDE_Extrapolation == .false.)) then
!!$                      !--Extrapolation, set to default, (zero)
!!$                      Posterior_perGalaxy_Redshift(z) = 1.e-100_double
!!$                   else
!!$                      !--Interpolation--!
!!$                      Posterior_perGalaxy_Redshift(z) = Linear_Interp(Cat%MF606W(c)+2.5e0_double*dlog10(Effective_Magnification), Cat%Sizes(c)/dsqrt(Effective_Magnification), PriorMagGrid, PriorSizeGrid, Kappa_Renormalised_Prior, ExValue = 1.e-100_double)*(1.e0_double/dsqrt(Effective_Magnification))*RedshiftPDF
!!$                   end if
!!$                   
!!$                   !--If no KDE Extrapolation, then this will set to a default value (effectively zero) outside the prior grid range
!!$                end if
!!$             case(3) !--Magnitude Only--!
!!$                !--Renormalise Magnitude Prior Distribution--!
!!$                
!!$                if(Cuts_Renormalise_Likelihood) then
!!$                   Renorm = Linear_Interp(Effective_Magnification, MagnificationGrid, MagOnly_Renormalisation_by_Magnification, ExValue = 1.e30_double)                      
!!$                   if(Renorm == 0) then
!!$                      Kappa_Renormalised_MagPrior = 0.e0_double
!!$                   else
!!$                      Kappa_Renormalised_MagPrior = Survey_Renormalised_MagPrior/Renorm
!!$                   end if
!!$                end if
!!$                
!!$                !--Could be edited for KDE Extrapolation, not done yet--!
!!$                
!!$                !--Construct p_[m_0] as the magnitude distribution for a sample of galaxies between cuts-corrected survey size limits 
!!$                !-MagRenorm
!!$!!!!$                   allocate(Mag_Only_Prior(2)); Mag_Only_Prior = 0.e0_double
!!$!!!$                   !--Find point where de-lensed magnitude lie on intrinsic distribution - could also use locate.
!!$!!!$                   do j = 1, size(PriorMagGrid)-1
!!$!!!$                      if( ( PriorMagGrid(j)<=  Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification) ) .and. (  PriorMagGrid(j+1) > Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification)) ) then
!!$!!!$                         Mag_Only_Prior(1) = Integrate(PriorSizeGrid, Kappa_Renormalised_Prior(j,:), 2, lim = Renormalisation_Size_Limits)
!!$!!!$                         Mag_Only_Prior(2) = Integrate(PriorSizeGrid, Kappa_Renormalised_Prior(j+1,:), 2, lim = Renormalisation_Size_Limits)
!!$!!!$                         exit
!!$!!!$                      end if
!!$!!!$                   end do
!!$!!!$                   Posterior_perGalaxy_Redshift(c,i,z) = Linear_Interp(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification), PriorMagGrid(j:j+1), Mag_Only_Prior, ExValue =  1.e-100_double)*RedshiftPDF(z)
!!$!!!$                   deallocate(Mag_Only_Prior)
!!$                !---------------------------------------------------------------------------------------------------------------------
!!$                
!!$                !--The following gives unbiased results if no size cuts are used, however is known to be biased in the presence of size cuts, I believe that this is due to the fact that the prior itself should depend on the size cuts in the presence of a size-magnitude correlation.
!!$                
!!$                Posterior_perGalaxy_Redshift(z) = Linear_Interp(Cat%MF606W(c) + 2.5e0_double*dlog10(Effective_Magnification), PriorMagGrid, Kappa_Renormalised_MagPrior, ExValue =  1.e-100_double)*RedshiftPDF
!!$
!!$             case default
!!$                STOP 'DM_Profile_Variable_Posterior - Error in choosing method of finding posterior'
!!$             end select
!!$          end do !--End of Redshift Loop
!!$          
!!$          if(size(Posterior_perGalaxy_Redshift,1) == 1) then
!!$             !--Redshift was taken to be exact
!!$             Posterior_perGalaxy(c) = Posterior_perGalaxy_Redshift(1)
!!$          else
!!$             !--Integrate over the redshift Information--!
!!$             Posterior_perGalaxy(c) = TrapInt(RedshiftGrid, Posterior_perGalaxy_Redshift(:))
!!$          end if
!!$          deallocate(Posterior_perGalaxy_Redshift)
!!$          
!!$       end do !--End of galaxy loop       
!!$
!!$       !~~Return lnP to ensure that we can renormalise in an alpha-independent way in the combined posterior, and to ensure that the PDF is correctly recovered even with ronding error (large negative lnP across all alphas). The renormalisation is done after the outer loop is finished
!!$       call Combine_Posteriors(Posterior(1,i), Posterior_perGalaxy(:), Combine_log_Posteriors, Renormalise = .false., Return_lnP = Combine_log_Posteriors, Combined_Posterior = Posterior(2,i))
!!$
!!$    end do !--End of Posterior Loop
!!$    deallocate(Posterior_perGalaxy)
!!$
!!$    if(any(dabs(Posterior) > huge(1.e0_double)) .or. any(isNaN(Posterior))) then
!!$       print *, 'NaNs or infinities found in posterior:'
!!$       do i =1, size(Posterior,2)
!!$          print *, Posterior(:,i)
!!$       end do
!!$       STOP
!!$    end if
!!$
!!$    !--Convert from lnP to P for posterior, renormalise so that max(P) = 1, to avoid rounding errors
!!$    if(Combine_log_Posteriors) Posterior(2,:) = dexp(Posterior(2,:) - maxval(Posterior(2,:)))
!!$
!!$    if(n_Default_Source_Redshift_Used > 0) write(*,'(A,I3,A)') '****** Used the default redshift for:', n_Default_Source_Redshift_Used, ' galaxies ********'
!!$    if(nGal_Ignored_SizeLimits > 0) write(*,'(A,I3)') '****** Number of galxies ignored as they fell outside the survey size limits:', nGal_Ignored_SizeLimits
!!$    if(nGal_Ignored_NaN > 0) write(*,'(A,I3)') '****** Number of galaxies ignored as they were NaNs:', nGal_Ignored_NaN
!!$    
!!$    
!!$    print *, '-------------------------------------------------------------'
!!$    print *, 'Constructed ', nSizePosterior, ' size-only posteriors'
!!$    print *, 'Constructed ', nSizeMagPosterior, ' size-magnitude posteriors'
!!$    print *, 'Constructed ', nMagPosterior, ' magnitude-only posteriors'
!!$    print *, '-------------------------------------------------------------'
!!$
!!$    !---DO NOT DELETE THIS. Instead, come up with a way to store the posterior per galaxy on a general grid for output
!!$!!!$    if(Debug_Mode .and. Output_Posterior_Per_Galaxy) then
!!$!!!$       
!!$!!!$       Filename = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat'
!!$!!!$       open(unit = 17, file = Filename)
!!$!!!$       write(fmtstring,'(I7)') size(Posterior_perGalaxy,1)+1
!!$!!!$       do i =1, size(Posterior,2)
!!$!!!$          write(17, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior_perGalaxy(:,i)
!!$!!!$       end do
!!$!!!$       close(17)
!!$!!!$       print *, 'Output Posterior per galaxy to: ', trim(adjustl(Filename))
!!$!!!$    end if
!!$    
!!$    
!!$    Filename= trim(adjustl(Output_Prefix))//'Posterior_Combined_Renormalised.dat'
!!$    open(unit = 82, file = Filename)
!!$    write(fmtstring,'(I1)') 2
!!$    do i =1, size(Posterior,2)
!!$       write(82, '('//trim(adjustl(fmtstring))//'(e14.7,x))') Posterior(1,i), Posterior(2,i)
!!$    end do
!!$    close(82)
!!$    print *, 'Output Combined Posterior to: ', trim(adjustl(Filename))
!!$    
!!$
!!$    if(any(isNAN(Posterior(2,:)))) then
!!$       print *, 'Any NaNs in aperture posterior?', any(isNAN(Posterior(2,:))), count(isNAN(Posterior(2,:)) == .true.)
!!$       print *, 'Stopping'
!!$       STOP
!!$    END if
!!$    
!!$    !--Remove any information less than a tolerance - Numerical Error--!
!!$    where(Posterior(2,:) < 1.e-12_double)
!!$       Posterior(2,:) = 0.e0_double
!!$    end where
!!$    
!!$    deallocate(Distance_From_Mass_Center)
!!$    if(allocated(Kappa_Renormalised_Prior)) deallocate(Kappa_Renormalised_Prior)
!!$    if(allocated(Survey_Renormalised_Prior)) deallocate(Survey_Renormalised_Prior)
!!$    if(allocated(Size_Only_Mag_Prior)) deallocate(Size_Only_Mag_Prior)
!!$    
!!$    !--On Successful Completion delete Poster per galaxy as large file--!
!!$!!!$    inquire(file = trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat', exist = here)
!!$!!!$    if(here == .false.) STOP "Posterior per galaxy doesn't exist, stopping before accidental deletion"
!!$!!!$    call system('rm '//trim(adjustl(Output_Prefix))//'Posterior_Per_Galaxy.dat')
!!$    
!!$  end subroutine DM_Profile_Variable_Posterior

