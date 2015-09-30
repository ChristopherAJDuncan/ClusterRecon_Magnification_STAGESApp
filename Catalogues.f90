module Catalogues_Declarations
  !--Modules encompasses declarations needed for referencesing in Interfaces--!
  use Param_Types; use RunTime_Input, only:verbose

end module Catalogues_Declarations

module Catalogues
  use Param_Types; use Catalogues_Declarations
  implicit none

  !--Directory and ReadIn Declarations--!
  character(50),private,parameter::Cat_Dir = 'Catalogues/'
!  character(60),private::Catalogue_Name = 'STAGES_shear.cat'

  !--Catalogue Derived Type--!
  type Catalogue
     character(60)::Name
     character(len = 10)::Mag_Label
     integer, dimension(:), allocatable:: Posterior_Method, RenormalisationGroup
     real(double),dimension(:),allocatable::flux, fluxerr, MF606W, magerr, Absolute_Magnitude, Surface_Brightness
     integer,dimension(:),allocatable::ntile
     integer,dimension(:), allocatable::Galaxy_Number
     real(double),dimension(:),allocatable::xpos, ypos
     real(double),dimension(:),allocatable::RA,Dec
     character(len = 4)::Sizes_Label
     logical::log_sizes
     real(double),dimension(:),allocatable::Sizes, Physical_Sizes!_FWHM, Size_FR, Size_KSB
     real(double),dimension(:), allocatable::Redshift
     real(double),dimension(:),allocatable:: g1,g2 !-These will contain the shear components.!
     !-Missing are observed ellipticity, anisotropy corrected ellipticity and shear responsivity correction-!
  end type Catalogue

  !-- Binned Catlogue contains a seperate Catalogue derived type for each reshift bin, which contains the subset of galaxies in that Bin--!
  type Binned_Catalogue_Mapping
     integer,allocatable::Map(:)
  end type Binned_Catalogue_Mapping

  type Binned_Catalogue
     character(5)::Label
     real(double),dimension(:,:),allocatable::Bin_Limits
     integer,dimension(:),allocatable::Occupation
     type(Catalogue),dimension(:),allocatable::Cat
     type(Binned_Catalogue_Mapping),dimension(:),allocatable::Index_Mapping
  end type Binned_Catalogue



  !--Global Parameters--!
  real(double)::ACSPixel_Size = 0.049 !-in arcseconds-!


  interface Catalogue_Assign_byGalaxy
     module procedure Catalogue_Assign_byGalaxy_byIndex, Catalogue_Assign_byGalaxy_byCatalogue 
  end interface Catalogue_Assign_byGalaxy

  INTERFACE Mask_Circular_Aperture
     MODULE PROCEDURE Mask_Multiple_Circular_Aperture, Mask_Single_Circular_Aperture
  end INTERFACE Mask_Circular_Aperture

  contains

    !---Catalogue Inforamtion Storage----!
    subroutine common_Catalogue_directories(Index, Directory, Filename, Columns)
      !-Index < 0 indeicates a blank field catalogue-!
      integer,intent(in)::Index
      character(*)::Filename, Directory
      integer,intent(out),allocatable::Columns(:) !-Galaxy_Number, ntile, RA, Dec, xpos, ypos, Absolute_Magnitude, MF606W, MagErr, Flux, FluxErr, Size, g1, g2, redshift-!    

      logical::here

      character(200)::Base_Directory != '/disk1/ps1/cajd/STAGES_SimultaneousFit/'!'/disk1/cajd/Size_Magnification/'

      call getcwd(Base_Directory)
      Base_Directory = trim(Base_Directory)//'/'

      if(allocated(Columns)) deallocate(Columns)
      allocate(Columns(15)); Columns = 0

      select case(Index)
      case(1) !-Full RRG+COMBO AND STAGES matched catalogue (source must occur in STAGES master and RRG measured)--!
         Directory = 'Catalogues/'
         Filename = 'STAGES_GALFIT.cat'
         Columns = (/-1,-1, 12, 13, 14, 15, -1, 3, -1, 1, 2, 9, -1, -1, 17/) !-MF606W: 7 (Galfit), 3: SEx
!!$
!!$         !!-- HST Extended - includes non-GALFIT option
!!$         Directory = 'Catalogues/'
!!$         Filename = 'STAGES_HSTExtended.cat'
!!$         Columns = (/-1,-1, 12, 13, 14, 15, -1, 3, -1, 1, 2, 9, -1, -1, 18/) !-MF606W: 7 (Galfit), 3: SEx

!!$         Directory = 'Catalogues/STAGES/'
!!$         Filename = 'RRG_AND_STAGES_GALFIT.cat'
!!$         Columns = (/-1,-1, 12, 13, 10, 11, -1,40, -1, 3, 4, 42, -1, -1, 44/) !--THis uses RRG Sex run for Magnitudes and Fluxs, and GALFIT RE from master where available - Mag : 5 SEX RRG, 38 SEX Master


         !Filename = 'STAGES_KSBf90_RRG+COMBO_AND_STAGES.cat.WCalib.SSNRCalib'
         !Columns = (/-1,1, 12, 13, 10, 11, -1, 41, -1, 3, 4, 31, 28, 29, 35/) !-Size: 31/32 RRG(TrQ/DetQ), 33/34 KSB(TrQ/DetQ), 21 SExtractor FR-! !--STAGES Master: 39: ST_FLUX_BEST, 40:ST_FLUXERR_BEST, 41:ST_MAG_BEST, 42: ST_MAGERR_BEST, 43: ST_MAG_GALFIT, 44: ST_MAGERR_GALFIT (DEFAULT MF606w: 5)

      case(-1) !--GEMS GALFIT
         Directory = 'Catalogues/'
         Filename = 'GEMS_GALFIT.cat'
         Columns = (/-1,-1, 12, 13, 14, 15, -1, 7, -1, 1, 2, 9, -1, -1, 17/) !-MF606W: 7 (Galfit), 3: SEx

      case(2) !--RRG+STAGES - All RRG sources, with STAGES matched where avaialble--!
!!$
!!$         Directory = 'Catalogues/STAGES/'
!!$         Filename = 'RRG_AND_STAGES_GALFIT.cat'
!!$         Columns = (/-1,-1, 12, 13, 10, 11, -1, 5, -1, 3, 4, 32, 28, 29, -1/)
!!$
         Directory = 'Catalogues/'
         Filename = 'RRG+STAGES_HSTExtended.cat'
         Columns = (/-1,-1, 12, 13, 10, 11, -1, 5, -1, 3, 4, 47, 28, 29, 35/) !--THis uses RRG Sex run for Magnitudes and Fluxs, and GALFIT RE from master where available
         !^^-- Mags:: 5: RRG SEx, 41: Master Mag

!!$         Filename = 'STAGES_KSBf90_RRG+COMBO.cat.WCalib.SSNRCalib.+SB'
!!$         Columns = (/-1,1, 12, 13, 10, 11, -1, 5, -1, 3, 4, 31, 28, 29, 35/) !-Size: 31/32 RRG(TrQ/DetQ), 33/34 KSB(TrQ/DetQ), 21 SExtractor FR-!
      case(3) !-- KSBf90 with RRG output--!
         Directory = 'Catalogues/'
         Filename =  'STAGES_KSBf90_RRG+COMBO.cat.WCalib.SSNRCalib.+SB'
         Columns = (/-1,-1, 12, 13, 10, 11, -1, 5, -1, 3, 4, 31, 28, 29, -1/) !-Size: 31/32 RRG(TrQ/DetQ), 33/34 KSB(TrQ/DetQ), 21 SExtractor FR-!
      case(4) !--STAGES-like Mock Catalogue--!
         Directory = 'Simulations/Output/'
         Filename = 'Mock_STAGES.cat'
         Columns = (/1,-1, 2, 3, -1, -1, 6, 7, -1, 11, 12, 4, -1, -1, -1/)
      case(5) !--COMBO17-like Mock Catalogue--!
         Directory = 'Simulations/Output/'
         Filename = 'Mock_COMBO.cat'
         Columns = (/1,-1, 2, 3, -1, -1, 6, 7, -1, 11, 12, 4, -1, -1, 5/)
      case(45) !--Combined STAGES and COMBO mock--!
         Directory = 'Simulations/Output/'
         Filename = 'Mock_STAGES+COMBO.cat'
         Columns = (/1,-1, 2, 3, -1, -1, 6, 7, -1, 11, 12, 4, -1, -1, 5/)
      case(6) !--COMBO-17 RRG Sizes--!
         Directory = 'Catalogues/'
         Filename = 'RRG_COMBO-17_gsMag.pzcat'
         Columns = (/-1,10, 1, 2, 13, 14, 8, 12, -1, 11, -1, 24, -1, -1, 5/)!-Sizes 24 (Tr), 25(det)-! 
      case(7) !--COMBO17 Catalogue: Unmatched, taken directly from COMBO-STAGES master :stages_20080529_v1pub.fits:--!
         Directory = 'Catalogues/'
         Filename = 'COMBO-17-Redshift_Cuts_21022014.pzcat'
         Columns = -1
      case(-4) !-STAGES-like Blank Field (Unlensed) Mock Catalogue-!
         Directory = 'Simulations/Output/'
         Filename = 'Mock_STAGES_Unlensed.cat'
         Columns = (/1,-1, 2, 3, -1, -1, 6, 7, -1, 11, 12, 4, -1, -1, -1/)         
      case(-5) !-COMBO-like Blank Field (Unlensed) Mock Catalogue-!
         Directory = 'Simulations/Output/'
         Filename = 'Mock_COMBO_Unlensed.cat'
         Columns = (/1,-1, 2, 3, -1, -1, 6, 7, -1, 11, 12, 4, -1, -1, 5/)
      case(-45) !--Combined STAGES and COMBO mock--!
         Directory = 'Simulations/Output/'
         Filename = 'Mock_STAGES+COMBO_Unlensed.cat'
         Columns = (/1,-1, 2, 3, -1, -1, 6, 7, -1, 11, 12, 4, -1, -1, 5/)
      case default
         print *, 'Catalogue Identifier:', Index
         STOP 'common_Catalogue_directories - Invalid Index entered: I do not have any information on this catalogue, retry with entry by hand'
      end select

      Directory = trim(Base_Directory)//trim(Directory)

!!$      inquire(file = trim(adjustl(Directory))//trim(adjustl(Filename)), exist = here)
!!$      if(here == .false.) then
!!$         print *, 'common_Catalogue_directories - Cannot find Catalogue, does it still exist?'
!!$         print *, trim(adjustl(Directory))//trim(adjustl(Filename))
!!$         STOP
!!$      END if

    end subroutine common_Catalogue_directories

    logical function check_Catalogue_Existence(Filename)
      character(*)::Filename

      logical::here

      inquire(File = Filename, exist = here)
      check_Catalogue_Existence = here
    end function check_Catalogue_Existence

    subroutine Concatonate_Catalogues(Master, Child)
      type(Catalogue):: Master, Child

      type(Catalogue):: tMaster
      integer::i

      tMaster = Master

      call Catalogue_Destruct(Master)

      call Catalogue_Construct(Master, size(tMaster%RA)+size(Child%RA))
      do i = 1, size(tMaster%RA)
         call Catalogue_Assign_byGalaxy_byCatalogue(Master, i, tMaster, i)
      end do

      do i = 1, size(Child%RA)
         call Catalogue_Assign_byGalaxy_byCatalogue(Master, size(tMaster%RA)+i, Child,i)
      end do

      call Catalogue_Destruct(tMaster)

    end subroutine Concatonate_Catalogues


    subroutine Mask_Single_Circular_Aperture(Cat, Ap_Pos, Ap_Radius)
      !--Cuts out galaxies which fall inside an aperture radius, where the radius is in DEGREES
      type(Catalogue), intent(inout)::Cat
      real(double), intent(in):: Ap_Pos(:)
      real(double), intent(in):: Ap_Radius

      real(double):: tAp_Pos(1,2)

      if(size(Ap_Pos) /= 2) STOP 'Mask_Single_Circular_Aperture - Aperture Position entered is not of the correct size'

      tAp_Pos(1,:) = Ap_Pos

      call Mask_Multiple_Circular_Aperture(Cat, tAp_Pos, (/Ap_Radius/))

    end subroutine Mask_Single_Circular_Aperture
      
    subroutine Mask_Multiple_Circular_Aperture(Cat, Ap_Pos, Ap_Radius)
      !--Cuts out galaxies which fall inside an aperture radius, where the radius is in DEGREES
      type(Catalogue), intent(inout)::Cat
      real(double), intent(in):: Ap_Pos(:,:)
      real(double), intent(in):: Ap_Radius(:)

      type(Catalogue)::Temp_Cat
      integer:: Ap, ac, c, Expected_Cut

      print *, '****Masking out Apertures:'
      do Ap = 1, size(Ap_Pos,1)
         print *, Ap, Ap_Pos(Ap,:), Ap_Radius(Ap)
      end do

      do Ap = 1, size(Ap_Pos,1)
         Expected_Cut = 0
         Expected_Cut = count( dsqrt( (Cat%RA-Ap_Pos(Ap,1))**2.e0_double +(Cat%Dec-Ap_Pos(Ap,2))**2.e0_double ) <= Ap_Radius(Ap) )
         print *, 'Cutting:', Expected_Cut, ' from Aperture:', Ap

         call Catalogue_Construct(Temp_Cat, size(Cat%RA)-Expected_Cut)

         ac = 0
         do c = 1, size(Cat%RA)
            if( (dsqrt( (Cat%RA(c)-Ap_Pos(Ap,1))**2.e0_double +(Cat%Dec(c)-Ap_Pos(Ap,2))**2.e0_double ) > Ap_Radius(Ap)) ) then
               if(ac > Size(Temp_Cat%RA)) STOP 'Identify_Galaxys_in_Circular_Aperture - Error in assigning aperture galaxy - GALAXY ASSINGATION IS LARGER THAN EXPECTED, stopping..'
               ac = ac + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(Temp_Cat, ac, Cat, c)
            end if
         end do
         Cat = Temp_Cat
         call Catalogue_Destruct(Temp_Cat)
      end do

    end subroutine Mask_Multiple_Circular_Aperture
      
    subroutine Identify_Galaxys_in_Circular_Aperture(Cat, Ap_Pos, Ap_Radius, Ap_Cat, Core_Radius)!, verbose)
      !--Returns a reduced catalogue containing only the galaxies in the aperture
      !--Ap_Pos, Ap_Radius and Core_Radius should be in DEGREES
      type(Catalogue), intent(in)::Cat
      real(double),intent(in)::Ap_Pos(:), Ap_Radius
      type(Catalogue), intent(out)::Ap_Cat
      real(double), intent(in), Optional:: Core_Radius
!      logical, intent(in), optional:: verbose

      real(double):: iCore
      integer::c, ac
      integer:: Expected_Number_In_Aperture
      

      if(present(Core_Radius)) then
         iCore = Core_Radius
      else
         iCore = 0.e0_double
      end if
      
      Expected_Number_in_Aperture = count( dsqrt( (Cat%RA-Ap_Pos(1))**2.e0_double +(Cat%Dec-Ap_Pos(2))**2.e0_double ) <= Ap_Radius) - count( dsqrt( (Cat%RA-Ap_Pos(1))**2.e0_double +(Cat%Dec-Ap_Pos(2))**2.e0_double ) <= iCore )
      
      write(*,'(A,2e9.5)') '----- Identifying sample in aperture palced at (RA,Dec):', Ap_Pos

      if(Expected_Number_in_Aperture == 0) STOP 'Identify_Galaxys_in_Circular_Aperture - There are no galaxies within the aperture, suggest increasing aperture size'
      
      call Catalogue_Construct(Ap_Cat, Expected_Number_in_Aperture)

!!$!      if(present(verbose) .and. verbose) then
         print *, ' '
         print *, 'Of ', count( dsqrt( (Cat%RA-Ap_Pos(1))**2.e0_double +(Cat%Dec-Ap_Pos(2))**2.e0_double ) <= Ap_Radius), ' galaxies within the aperture, ', count( dsqrt( (Cat%RA-Ap_Pos(1))**2.e0_double +(Cat%Dec-Ap_Pos(2))**2.e0_double ) <= iCore ), ' will be subtracted as they exisit in the core, leaving:', Expected_Number_in_Aperture
         print *, ' '
!!$!      end if

      ac = 0
      do c = 1, size(Cat%RA)
         if( (dsqrt( (Cat%RA(c)-Ap_Pos(1))**2.e0_double +(Cat%Dec(c)-Ap_Pos(2))**2.e0_double ) <= Ap_Radius) .and. (dsqrt( (Cat%RA(c)-Ap_Pos(1))**2.e0_double +(Cat%Dec(c)-Ap_Pos(2))**2.e0_double ) > iCore) ) then
            if(ac > Expected_Number_in_Aperture) STOP 'Identify_Galaxys_in_Circular_Aperture - Error in assigning aperture galaxy - GALAXY ASSINGATION IS LARGER THAN EXPECTED, stopping..'
            ac = ac + 1
            call Catalogue_Assign_byGalaxy_byCatalogue(Ap_Cat, ac, Cat, c)
         end if
      end do
      
      if(ac < Expected_Number_in_Aperture) print *, '**** WARNING - Identify_Galaxys_in_Circular_Aperture - Number assigned to catalogue does not equal the number expected'

  end subroutine Identify_Galaxys_in_Circular_Aperture

    !---Catalogue Binning Routines----!
    subroutine Calculate_Bin_Limits_by_equalNumber(Array, nBin, Limits)
      !--Returns the Bin Limits requires for nBin bins on Array, where the Bin is set so that there are an equal number of elements of Array in each bin (therefore suited to an array that contains information for a set of galaxies, and each element is a different galaxy--!
      use nr, only:sort
      real(double),intent(in)::Array(:)
      integer,intent(in)::nBin
      real(double), allocatable,intent(out):: Limits(:,:)

      integer:: nPerBin, i
      real(double),dimension(size(Array))::Sorted_Array

      if(allocated(Limits)) deallocate(Limits)
      allocate(Limits(nBin,2)); Limits = 0.e0_double

      if(nBin ==1) then
         Limits(1,:) = (/minval(Array), maxval(Array)/)
         return
      end if

      nPerBin = int(size(Array)/nBin + 1)

      Sorted_Array = Array
      !--Sort array--!
      call sort(Sorted_Array)
      
      Limits(1,:) = (/Sorted_Array(1), Sorted_Array(minval( (/nPerBin,size(Sorted_Array)/) ))/)
      do i = 2, nBin
         Limits(i,1) = Limits(i-1,2)
         Limits(i,2) = Sorted_Array(minval((/i*nPerBin,size(Sorted_Array)/)))
      end do
      Limits(nBin,2) = Sorted_Array(size(Sorted_Array)) 

    end subroutine Calculate_Bin_Limits_by_equalNumber

    subroutine bin_catalogue_by_magnitude(Cat,Limits,BCat, Mag_Type)
      USE NR, ONLY:SORT
      !--Bins the catalogue by magnitude, returning BCat wof type Binned_Catalogue, which contains a seperate catalogue for each magnitude bin entered in Limits
      type(Catalogue),intent(in)::Cat
      type(Binned_Catalogue)::BCat
      real(double), dimension(:,:),intent(in)::Limits
      integer,intent(in),optional::Mag_Type

      integer::Magnitude_Type !-1:Absolute; 2: MF606w-! !--Neds to be passed in eventually-!
      
      integer::c,i, j
      
      type(Catalogue)::tCat !-temporary allocation-!                                       
      integer::counter
      integer,allocatable::Expected_Occupation(:), Occupation(:)
      real(double),allocatable::Temporary_Magnitude(:)
      
      real(double),allocatable::Sorted_Array(:)
      logical::Stop_Flag
      logical::Index_Map_Found

      if(Verbose) print *, 'Binning Catalogue by Magnitude...'
      
      if(present(Mag_Type)) then
         Magnitude_Type = Mag_Type
      else
         Magnitude_Type = 1 !-Default-!
      end if

      allocate(Expected_Occupation(size(Limits,1))); Expected_Occupation = 0
      allocate(Occupation(size(Limits,1))); Occupation = 0

      !--Construct Array which stores the Magnitude with respect to which we are binning--!
      allocate(Temporary_Magnitude(size(Cat%RA)));
      select case(Magnitude_Type)
      case(1)
         print *, 'Binning Catalogue by Absolute_Magnitude'
         Temporary_Magnitude = Cat%Absolute_Magnitude
      case(2)
         print *, 'Binning Catalogue by MF606W'
         Temporary_Magnitude = Cat%MF606W
      case default
         STOP 'bin_catalogue_by_magnitude - Invalid choice of magnitude type , stopping...'
      end select
      
      if(all(Temporary_Magnitude == 0.e0_double) .or. (size(Limits,1) == 1 .and. (minval(Temporary_Magnitude)>= limits(1,1) .and. minval(Temporary_Magnitude)<= limits(1,2))) ) then
         print *, 'Bin_by_Magnitude, magnitude information empty, returning Binned Catalogue which contains all galaxies'
         call Binned_Catalogue_Construct(BCat, 'No Binning', 1)
         call Catalogue_Construct(BCat%Cat(1), size(Cat%RA))
         BCat%Cat(1) = Cat
         allocate(BCat%Index_Mapping(1)%Map(size(Cat%RA))); BCat%Index_Mapping(1)%Map = -1
         do c = 1, size(Cat%RA)
            BCat%Index_Mapping(1)%Map(c) = c
         end do
         BCat%Occupation = size(Cat%RA)
         return
      end if

      call Binned_Catalogue_Construct(BCat, 'By Magnitude', size(Limits,1))
      BCat%Bin_Limits = Limits
      
      !-Construct Expected Occupation-!                                                                                                                                                                           
      !--Set up the Binned Catalogue--!
      do i = 1, size(Limits,1)
         Expected_Occupation(i) = count(Temporary_Magnitude <= Limits(i,2)) - count(Temporary_Magnitude <= Limits(i,1))

         !--Add in galaxies at lower limit to the lowest bin
         if(i==1) Expected_Occupation(1) = Expected_Occupation(1) + count(Temporary_Magnitude==Limits(1,1))
         call Catalogue_Construct(BCat%Cat(i), Expected_Occupation(i))
         allocate(BCat%Index_Mapping(i)%Map(Expected_Occupation(i))); BCat%Index_Mapping(i)%Map = -1
      end do
      
      counter = 0
      do c = 1, size(Cat%Sizes)
         do i =1, size(Limits,1)
            if( (Temporary_Magnitude(c) > Limits(i,1)) .and. (Temporary_Magnitude(c) <= Limits(i,2)) ) then
               Occupation(i) = Occupation(i) + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(BCat%Cat(i), Occupation(i), Cat, c)
               BCat%Index_Mapping(i)%Map(Occupation(i)) = c
               exit !-Exit bin loop on success-!
            end if
            if(i == 1 .and. Temporary_Magnitude(c) == Limits(i,1)) then
               counter = counter + 1
               Occupation(i) = Occupation(i) + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(BCat%Cat(i), Occupation(i), Cat, c)
               BCat%Index_Mapping(i)%Map(Occupation(i))= c
               exit !-Exit bin loop on success-!
            end if
         end do
      end do



      
!--This error can be ignored if some galaxies are expected to fall outside the bin limits-!
!!$      do i =1, size(Limits,1)
!!$         if(any(BCat%Index_Mapping(i)%Map > size(Cat%Sizes)) .or. any(BCat%Index_Mapping(i)%Map < 0)) then
!!$            print *, 'Fatal Error - bin_catalogue_by_magnitude - Mapping to small or large:: Bin, Max/Min Mapping, nGal:', i, maxval(BCat%Index_Mapping(i)%Map), minval(BCat%Index_Mapping(i)%Map), size(Cat%Sizes)
!!$            STOP
!!$         end if
!!$      end do

      !--Check for Duplication of Mapping - i.e. there should be a monotonic mapping between catalogues--!
      allocate(Sorted_Array(sum(Occupation))); Sorted_Array = -100.e0_double
      counter = 0
      do i =1, size(Limits,1)
         Sorted_Array(counter + 1:counter+size(BCat%Index_Mapping(i)%Map)) = BCat%Index_Mapping(i)%Map
         counter = counter + size(BCat%Index_Mapping(i)%Map)
      end do
      call sort(Sorted_Array)
      do j = 2, size(Sorted_Array)-1
         if( (Sorted_Array(j-1) == Sorted_Array(j)) .or. (Sorted_Array(j+1) == Sorted_Array(j)) ) then
            print *, 'Duplication of Mapping, value:', Sorted_Array(j), j
            STOP_FLAG = .true.
         end if
      end do
      deallocate(Sorted_Array)
      if(STOP_FLAG) STOP
      !---------------------end duplication check-------------------------------------------------------!

      !--Search for galxys in the main catalogue which are between the binning limits, but have not been allocated to any bin--!
      !-Works, but problems is some galxies go missing, commented out for now-!
!(minval(Limits) <= minval(Temporary_Magnitude)) .and. (maxval(Limits) >= maxval(Temporary_Magnitude)) .and.
!!$      if((sum(Occupation) /= Size(Cat%RA))) THEN !!THIS IS ONLY TRUE IF LIMITS ARE SUCH THAT ALL SHOULD BE ACCOUNTED FOR!!
!!$         pRINT *, 'bin_catalogue_by_magnitude- FATAL ERROR - Sum of occupation is not equal to the size of the original array; not all galaxies accounted for:, Occ, Array', sum(Occupation), Size(Cat%RA)
!!$         do j = 1, size(Cat%RA)
!!$            Index_map_Found = .false.
!!$                  !-Ignore those that fall outside bin limits-!
!!$            if(Temporary_Magnitude(j) < minval(Limits)) cycle
!!$            if(Temporary_Magnitude(j) > maxval(Limits)) cycle
!!$            do i = 1, size(Limits,1)
!!$               if(j==1) print *, 'Occ, Expected:', Occupation(i), Expected_Occupation(i)
!!$               do c = 1, size(BCat%Cat(i)%RA)
!!$                  !---!
!!$                  if(BCat%Index_Mapping(i)%Map(c) == j) then
!!$                     !--Success--!
!!$                     Index_Map_Found = .true.
!!$                     exit
!!$                  end if
!!$                  if(Index_Map_Found== .false. .and. i == size(Limits,1) .and. c == size(BCat%Cat(i)%RA)) then
!!$                     print *, 'A Missing Galaxy was found:, galaxy:', j, 'with Magnitude:', Temporary_Magnitude(j), ' limits:', Limits(1,1), Limits(size(Limits,1),2)
!!$                  end if
!!$               end do
!!$               if(Index_map_found == .true.) exit
!!$            end do
!!$         end do
!!$
!!$         STOP
!!$      end if
      do i = 1, size(Limits,1)
         if(Occupation(i) /= Expected_Occupation(i)) then
            print *, 'Expected Error (FATAL): bin_catalogue_by_magnitude: Occupation of a bin is not equal to the expected occupation', i, Occupation(i), Expected_Occupation(i)
            STOP
         END if
      end do      

      if(Verbose) then
         print *, 'Sample of ', size(Cat%Sizes), ' galaxies split as:'
         do i =1, size(Limits,1)
            print *, 'Bin:', i, ' has ', size(BCat%Cat(i)%Sizes), ' galaxies. Bin is:', Limits(i,:)
         end do
      end if

      BCat%Occupation = Occupation
      
      deallocate(Temporary_Magnitude, Occupation, Expected_Occupation)
      if(Verbose) print *, 'Done.'
      
    end subroutine bin_catalogue_by_magnitude

    subroutine bin_catalogue_by_redshift(Cat,Limits,BCat)
      !--Bins the catalogue by redshift, returning BCat wof type Binned_Catalogue, which contains a seperate catalogue for each redshift bin entered in Limits
      !--On 13Jan2014, bin by magnitude was edited to account for the one galaxy that will fall into teh first bin, this will probably need done here - EDITED BUT NOT TESTED--!

      type(Catalogue),intent(in)::Cat
      type(Binned_Catalogue)::BCat
      real(double), dimension(:,:),intent(in)::Limits
      
      integer::c,i
      
      type(Catalogue)::tCat !-temporary allocation-!                                                                                                                                                                                                                          
      integer,allocatable::Expected_Occupation(:), Occupation(:)
      
      if(Verbose) print *, 'Binning Catalogue by Redshift...'
      
      allocate(Expected_Occupation(size(Limits,1))); Expected_Occupation = 0
      allocate(Occupation(size(Limits,1))); Occupation = 0
      
      call Binned_Catalogue_Construct(BCat, 'By Redshift', size(Limits,1))
      BCat%Bin_Limits = Limits
      
      !-Construct Expected Occupation-!                                                                                                                                                                                                                                       
      do i = 1, size(Limits,1)
!         Expected_Occupation(i) = count(Cat%Redshift <= Limits(i,2)) - count(Cat%Redshift <= Limits(1,1))  - sum(Expected_Occupation(:i))
         Expected_Occupation(i) = count(Cat%Redshift <= Limits(i,2)) - count(Cat%Redshift <= Limits(i,1))
         if(i==1) Expected_Occupation(1) = Expected_Occupation(1) + count(Cat%Redshift==Limits(1,1))
         !--Edit--!
!         Expected_Occupation(size(Expected_Occupation)) = Expected_Occupation(size(Expected_Occupation)) + count(Cat%Redshift==Limits(size(Limits,1),2))
         !--------!
         call Catalogue_Construct(BCat%Cat(i), Expected_Occupation(i))
         allocate(BCat%Index_Mapping(i)%Map(Expected_Occupation(i))); BCat%Index_Mapping(i)%Map = -1
      end do
      
      do c = 1, size(Cat%Sizes)
         do i =1, size(Limits,1)
            if( (Cat%Redshift(c) > Limits(i,1)) .and. (Cat%Redshift(c) <= Limits(i,2)) ) then
               Occupation(i) = Occupation(i) + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(BCat%Cat(i), Occupation(i), Cat, c)
               BCat%Index_Mapping(i)%Map(Occupation(i)) = c
            end if
            if(i == 1 .and. Cat%Redshift(c) == Limits(i,1)) then
               Occupation(i) = Occupation(i) + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(BCat%Cat(i), Occupation(i), Cat, c)
               BCat%Index_Mapping(i)%Map(Occupation(i)) = c
            end if
         end do
      end do
      !--Assign all the galaxies at the upper boundary to the lower redshift bin--!
      if(sum(Occupation) /= Size(Cat%RA)) THEN
         pRINT *, 'bin_catalogue_by_redshift- FATAL ERROR - Sum of occupation is not equal to the size of the original array; not all galaxies accounted for:, Occ, Array', sum(Occupation), Size(Cat%RA)
         do i = 1, size(Limits,1)
            print *, 'Occ, Expected:', Occupation(i), Expected_Occupation(i)
         end do
         STOP
      end if
         
      do i = 1, size(Limits,1)
         if(Occupation(i) /= Expected_Occupation(i)) print *, 'Expected Error (FATAL): bin_catalogue_by_redshift: Occupation of a bin is not equal to the expected occupation', Occupation(i), Expected_Occupation(i)
         
      end do
      
      if(Verbose) then
         print *, 'Sample of ', size(Cat%Sizes), ' galaxies split as:'
         do i =1, size(Limits,1)
            print *, 'Bin:', i, ' has ', size(BCat%Cat(i)%Sizes), ' galaxies. Bin is:', Limits(i,:)
         end do
      end if

      BCat%Occupation = Occupation
      
      if(Verbose) print *, 'Done.'
      
    end subroutine bin_catalogue_by_redshift

    subroutine bin_catalogue_by_Redshift_and_Magnitude(Cat, RedLimits, MagLimits, BBCat, Bin_Occupation)
      type(Binned_Catalogue),allocatable,intent(out)::BBCat(:) !-Size nZ, i.e. one Binned_Catalogue for each redshift bin-!

      real(double),intent(in)::MagLimits(:,:), RedLimits(:,:)
      type(Catalogue),intent(in)::Cat
      integer,allocatable,intent(out),optional::Bin_Occupation(:,:)

      type(Binned_Catalogue)::BCat
      integer::BinZ, BinM

      !--Bin first by redshift using pre-defined Redshift--!                                                                                                                                                                                 
      call bin_catalogue_by_redshift(Cat, RedLimits, BCat)

      allocate(BBCat(size(RedLimits,1)))
      do BinZ = 1, size(BBCat)
         call Bin_Catalogue_by_magnitude(BCat%Cat(BinZ), MagLimits, BBCat(BinZ))
         if(present(Bin_Occupation)) then
            allocate(Bin_Occupation(size(RedLimits,1), size(MagLimits,1))); Bin_Occupation = 0.e0_double
            do BinM = 1, size(MagLimits,1)
               Bin_Occupation(BInZ, BinM) = BBCat(BinZ)%Occupation(BinM)
            end do
         end if
      end do


    end subroutine bin_catalogue_by_Redshift_and_Magnitude

    subroutine unBin_Binned_Catalogue(BCat, Cat)
      !--Unpacks the binned catalogue into the full Catalogue, in the order in which it was binned--!
      type(Binned_Catalogue)::BCat
      type(Catalogue)::Cat

      integer::b, c, bc

      integer::nGal, counter

      !--Tested 13Jan2014--!
      integer,allocatable:: Removed_Galaxy_Mapping(:) !-As some galaxies may have been removed from the original-!
      logical:: Galaxy_Found
      integer::Galaxy_Index_Adjustment, Max_Mapping

      !--Testing--!
      integer::Index

      nGal = 0
      do b = 1, size(BCat%Cat)
         nGal = nGal + size(BCat%Cat(b)%RA)
      end do

      call Catalogue_Destruct(Cat)
      call Catalogue_Construct(Cat, nGal)

!--Unbinning by mapping index only works if NO galaxies where thrown away from the original, binned catalogue, thus there must be a mapping between the full original catalogue (when binned) and the full catalogue now, accounting for discarded galaxies--!
      Max_Mapping = 0
      do b = 1, size(BCat%Cat)
         Max_Mapping = maxval( (/Max_Mapping, maxval(BCat%Index_Mapping(b)%Map)/) )
      end do
      allocate(Removed_Galaxy_Mapping(Max_Mapping)); Removed_Galaxy_Mapping = -1

      do c = 1, Max_Mapping
         Galaxy_Found = .false.
         !Search through the mappings in BCat to determine whether galaxy is present or cut!
         do b = 1, size(BCat%Cat)
            do bc = 1, size(BCat%Index_Mapping(b)%Map)
               if(BCat%Index_Mapping(b)%Map(bc) == c) then
                  !--Success--!
                  Galaxy_Found = .true.
                  exit
               end if
               if(b == size(BCat%Cat) .and. bc == size(BCat%Index_Mapping(b)%Map) .and. Galaxy_Found==.false.) then
                  !-Final Binned Galaxy checked and galxy not found-!
                  Galaxy_Index_Adjustment = Galaxy_Index_Adjustment + 1
               end if
            end do
            !-Don't contine to search bins for galaxy-!
            if(Galaxy_Found) exit
         end do
         Removed_Galaxy_Mapping(c) = c - Galaxy_Index_Adjustment
         if(Galaxy_Found == .false.) then
            Removed_Galaxy_Mapping(c) = -1000!sqrt(-1.e0_double)!1/0 !-Convert to NaN for missing galaxy so if referenced then error, and adjust all following-!
         end if
      end do

      counter = 0
      do b = 1, size(BCat%Cat)
         if(allocated(BCat%Index_Mapping(b)) == .false.) STOP 'unBin_Binned_Catalogue - Mapping not allocated'
         do c = 1, size(BCat%Cat(b)%RA)
            counter = counter + 1
            if(BCat%Index_Mapping(b)%Map(c) < 0) then
               print *, 'For bin:', b, ' galaxy:', c, ' of:', size(BCat%Index_Mapping(b)%Map),':'
               STOP 'unBin_Binned_Catalogue - Mapping is Negative, stopping'
            end if
            !call Catalogue_Assign_byGalaxy_byCatalogue(Cat, counter, BCat%Cat(b), c) if order is not important
            if(Removed_Galaxy_Mapping(BCat%Index_Mapping(b)%Map(c)) <= 0) STOP 'unBin_Binned_Catalogue - NaN references in removed mapping, stopping'
            call Catalogue_Assign_byGalaxy_byCatalogue(Cat, Removed_Galaxy_Mapping(BCat%Index_Mapping(b)%Map(c)), BCat%Cat(b), c)
            Index = BCat%Index_Mapping(b)%Map(c)
         end do
      end do

      call Binned_Catalogue_Destruct(BCat)

      deallocate(Removed_Galaxy_Mapping)

    end subroutine unBin_Binned_Catalogue

    !------------------END BINNING ROUTINES-------------------------!
    real(double) function Physical_Size_from_Pixel_Size(Redshift, Pixel_Size)
      use Cosmology, only:angular_diameter_distance_fromRedshift
      real(double),intent(in)::Redshift, Pixel_Size
      real(double)::Pixel_Size_Radians

      Pixel_Size_Radians = (ACSPixel_Size*3.141592654e0_double)/(180.e0_double*3600.e0_double)

      if(Redshift > 0.e0_double) Physical_Size_from_Pixel_Size = (Pixel_Size_Radians*Pixel_Size*angular_diameter_distance_fromRedshift(0.e0_double,Redshift))

    end function Physical_Size_from_Pixel_Size

    subroutine convert_Size_from_Physical_to_Pixel(Cat)
      use Cosmology
      !-Converts from a phyiscal galaxy size to pixel size. Can only be done if there is redshift information for each galaxy!!!-!
      !--Physical Size in Mpc/h, and is in proper size co-ords (ie extra factor of a)
      type(Catalogue)::Cat

      integer::g

      real(double)::Pixel_Size_Radians

      if(Verbose) print *, 'Calculating Physical Size using Pixel Size = ', ACSPixel_Size, ' arcseconds'

      if(all(Cat%Redshift < 0.e0_double)) STOP 'convert_Size_from_Pixel_to_Physical - All Galaxies in catalogue have not been assigned a redshift, exiting...'

    if(allocated(Cat%Sizes)) deallocate(Cat%Sizes)
    allocate(Cat%Sizes(size(Cat%Physical_Sizes))); Cat%Sizes = -1.e0_double

    Pixel_Size_Radians = (ACSPixel_Size*3.141592654e0_double)/(180.e0_double*3600.e0_double)
!!$    do g = 1, size(Cat%Sizes)
!!$       if(Cat%Redshift(g) > 0.e0_double) Cat%Physical_Sizes(g) = Pixel_Size_Radians*Cat%Sizes(g)*angular_diameter_distance_fromRedshift(Cat%Redshift(g))
!!$    end do

    where(Cat%Redshift > 0.e0_double)
!       Cat%Sizes = (Cat%Physical_Sizes*(1.e0_double+Cat%Redshift))/(Pixel_Size_Radians*angular_diameter_distance_fromRedshift(0.e0_double,Cat%Redshift)) !-Use if angular diameter distance is comoving-!
       Cat%Sizes = (Cat%Physical_Sizes)/(Pixel_Size_Radians*angular_diameter_distance_fromRedshift(0.e0_double,Cat%Redshift)) !-Use if angular diameter distance is proper-!
      
    end where

    if(Verbose) print *, 'Done'

  end subroutine convert_Size_from_Physical_to_Pixel


    subroutine convert_Size_from_Pixel_to_Physical(Cat)
      use Cosmology
      !-Converts from a pixel size to physical size. Can only be done if there is redshift information for each galaxy!!!-!
      !--PHysical Size in Mpc/h, and is in proper size coo-ords (ie extra factor of a)
      type(Catalogue)::Cat

      integer::g

      real(double)::Pixel_Size_Radians

      if(Verbose) print *, 'Calculating Physical Size using Pixel Size = ', ACSPixel_Size, ' arcseconds'

      if(all(Cat%Redshift < 0.e0_double)) STOP 'convert_Size_from_Pixel_to_Physical - All Galaxies in catalogue have not been assigned a redshift, exiting...'

    if(allocated(Cat%Physical_Sizes)) deallocate(Cat%Physical_Sizes)
    allocate(Cat%Physical_Sizes(size(Cat%Sizes))); Cat%Physical_Sizes = -1.e0_double

    Pixel_Size_Radians = (ACSPixel_Size*3.141592654e0_double)/(180.e0_double*3600.e0_double)
!!$    do g = 1, size(Cat%Sizes)
!!$       if(Cat%Redshift(g) > 0.e0_double) Cat%Physical_Sizes(g) = Pixel_Size_Radians*Cat%Sizes(g)*angular_diameter_distance_fromRedshift(Cat%Redshift(g))
!!$    end do

    where(Cat%Redshift > 0.e0_double)
!       Cat%Physical_Sizes = (Pixel_Size_Radians*Cat%Sizes*angular_diameter_distance_fromRedshift(0.e0_double,Cat%Redshift))/(1.e0_double+Cat%Redshift) !-Use if angular diameter distance is comoving-!
       Cat%Physical_Sizes = (Pixel_Size_Radians*Cat%Sizes*angular_diameter_distance_fromRedshift(0.e0_double,Cat%Redshift))  !-Use if angular diameter distance is proper-!
    end where

    if(Verbose) print *, 'Done'

    end subroutine convert_Size_from_Pixel_to_Physical

    subroutine match_Redshift_byCatalogue_byPosition(Cat, Ref_Cat)
      !-For each galaxy in Cat, finds the *exact* match (within small tolerance for rounding errors) in RA and Dec in teh Ref_Cat, and sets the Redshift of that galaxy to the redshift stored in ref Cat
      type(Catalogue)::Cat, Ref_Cat

      integer::nFail, nPass
      real(double)::tol = 1.e-6_double

      integer::i,j

      nPass = 0.e0_double; nFail = 0.e0_double
      do i = 1, size(Cat%Redshift)
         do j = 1, size(Ref_Cat%Redshift)
            if( (dabs(Cat%RA(i) - Ref_Cat%RA(j)) <= tol) .and. (dabs(Cat%Dec(i) - Ref_Cat%Dec(j)) <= tol) ) then
               Cat%Redshift(i) = Ref_Cat%Redshift(j)
               nPass = nPass + 1
               exit
            end if
            if(j==size(Ref_Cat%Redshift)) nFail = nFail + 1
         end do
      end do
      if(nPass + nFail /= size(Cat%Redshift)) print *, 'match_Redshift_byCatalogue_byPosition - Error in assigning redshifts: nPass and nFail do not add up to original length'
      
      if(Verbose) print *, 'Matching Redshifts:', nPass, ' galaxies out of', size(Cat%Redshift),' were matched'

    end subroutine match_Redshift_byCatalogue_byPosition

    subroutine Cut_By_PhotometricRedshift(Cat, LowerCut, UpperCut, Discard_Gals_wout_Redshift)
      !-Galaxies with Redshifts can be isolated by setting LowerCut <= 0.e0_double and Discard_Gals_wout_Redshift = .true.-!
      !-if Discard_Gals_wout_Redshift is present and true, then only galaxies with redhift information (defined as having redshift > 0) are returned-!
      type(Catalogue)::Cat
      real(double), intent(in),optional::LowerCut,UpperCut
      logical, intent(in),optional::Discard_Gals_wout_Redshift

      real(double)::ilower, iupper
      type(catalogue)::Temp_Cat
      logical::Discard_Gals
      integer::nPass,i, nGal

!!$      INTERFACE
!!$         subroutine Cut_By_PhotometricRedshift(Cat, LowerCut, UpperCut, Discard_Gals_wout_Redshift)
!!$           use Param_Types; use Catalogues_Declarations
!!$           type(Catalogue)::Cat
!!$           
!!$           real(double), intent(in),optional::LowerCut,UpperCut
!!$           logical, intent(in),optional::Discard_Gals_wout_Redshift
!!$         end subroutine Cut_By_PhotometricRedshift
!!$      END INTERFACE
      
      if(all(Cat%Redshift < 0) .or. allocated(Cat%Redshift) == .false.) then
         print *, 'Cut_By_PhotometricRedshift - Redshift data is empty or unassigned, RETURNING WITHOUT ALTERATION'
         return
      end if

      if((present(lowerCut)==.false.) .and. (present(upperCut)==.false.)) then
         print *, 'Error - Cut_By_PhotometricRedshift - Either a lower or upper cut needs to be entered, returning without cutting'
         return
      end if

      Discard_Gals = .false.
      if(present(Discard_Gals_wout_Redshift)) Discard_Gals = Discard_Gals_wout_Redshift

      if(present(lowerCut)) then
         ilower =lowerCut
      else
         ilower = minval(Cat%Redshift)
      end if
      ilower = maxval((/0.e0_double, ilower/)) !-Ensures lower never goes below zero-!
      if(present(upperCut)) then
         iupper =upperCut
      else
         iupper = maxval(Cat%Redshift)
      end if
  
      print *, 'Number of assigned galaxy redshifts:', count(Cat%Redshift >= 0.e0_double)
      print *, 'Cutting Catalogue to Phot-Z between limits:', ilower, iupper

      if(ilower >= iupper) then
         print *, 'lower/upper:', ilower, iupper
         STOP 'FATAL ERROR - Cut_By_PhotometricRedshift - lower cut larger than upper cut'
      end if
      
      nGal = size(Cat%RA)

      Temp_Cat = Cat

      call Catalogue_Destruct(Cat)

      if(Discard_Gals) then
         nPass = count(Temp_Cat%Redshift <= iupper) - count(Temp_Cat%Redshift < iLower) 
      else
         nPass = count(Temp_Cat%Redshift <= iupper) - count(Temp_Cat%Redshift < iLower) + count(Temp_Cat%Redshift < 0.e0_double)
      end if

      call Catalogue_Construct(Cat, nPass)

      nPass = 0
      if(Discard_Gals) then
         do i = 1, size(Temp_Cat%Redshift)
            if((Temp_Cat%Redshift(i) >= ilower) .and.(Temp_Cat%Redshift(i) <= iupper)) then
               nPass = nPass + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(Cat, nPass, Temp_Cat, i)
            end if
         end do
      else
         do i = 1, size(Temp_Cat%Redshift)
            if( (Temp_Cat%Redshift(i)< 0.e0_double) .or. ((Temp_Cat%Redshift(i) >= ilower) .and.(Temp_Cat%Redshift(i) <= iupper)) ) then
               nPass = nPass + 1
               call Catalogue_Assign_byGalaxy_byCatalogue(Cat, nPass, Temp_Cat, i)
            end if
         end do
      end if
      if(nPass /= size(Cat%RA)) then
         print *, 'Error - Cut_By_PhotometricRedshift - Error in Assigning galxies that pass cuts', nPass, size(Cat%RA)
         STOP
      end if

      print *, 'Out of:', size(Temp_Cat%Redshift),' orginal galaxies, ', nGal-nPass, ' were cut by redshift' 
      print *, '       leaving:', size(Cat%RA), ' galaxies in catalogue'
!      print *, 'Out of:', size(Temp_Cat%Redshift),' orginal galaxies, ', size(Cat%Redshift), ' galaxies passed the redshift cut' 

      call Catalogue_Destruct(Temp_Cat)
            
    end subroutine Cut_By_PhotometricRedshift


    subroutine Cut_By_Magnitude(Cat, lowerCut, upperCut)
      !--Only Set up to use MF606W, i.e. apparent magnitude--!
      type(Catalogue)::Cat
      real(double),intent(in),optional::lowerCut, upperCut

      real(double)::ilower, iupper
      type(Catalogue)::Temp_Cat
      integer::i
      integer::nPass

      if((present(lowerCut)==.false.) .and. (present(upperCut)==.false.)) then
         print *, 'Error - Cut_By_PixelSize - Either a lower or upper cut needs to be entered, returning without cutting'
         return
      end if

      if(all(Cat%MF606W <= 0.e0_double) .or. allocated(Cat%MF606W)==.false.) then
         print *, 'Cut_by_Magnitude: MF606W is not allocated or is empty, returning WITHOUT ALTERATION'
         return
      end if
      
      if(present(lowerCut)) then
         ilower =lowerCut
      else
         ilower = minval(Cat%MF606W)
      end if
      if(present(upperCut)) then
         iupper =upperCut
      else
         iupper = maxval(Cat%MF606W)
      end if

      if(ilower >= iupper) then
         print*,  'FATAL ERROR - Cut_By_Magnitude - lower cut larger than upper cut'
         print *, ilower, iupper
         STOP
      end if

      
      Temp_Cat = Cat

      call Catalogue_Destruct(Cat)

      nPass = 0
      do i = 1, size(Temp_Cat%MF606W)
         if(isNaN(Temp_Cat%MF606W(i))) cycle
         if((Temp_Cat%MF606W(i) >= ilower) .and. (Temp_Cat%MF606W(i) <= iupper)) nPass = nPass + 1
      end do

      !--This does not account for NaNs yet, which may result from negative magnification
!      nPass = count(Temp_cat%MF606W <= iupper) - count(Temp_Cat%MF606W < ilower)!count(Temp_cat%MF606W <= iupper) + count(Temp_Cat%MF606W >= ilower) - size(Temp_Cat%MF606W)


      call Catalogue_Construct(Cat, nPass)

      nPass = 0
      do i = 1, size(Temp_Cat%MF606W)
         if(isNaN(Temp_Cat%MF606W(i))) cycle
         if((Temp_Cat%MF606W(i) >= ilower) .and.(Temp_Cat%MF606W(i) <= iupper)) then
            nPass = nPass + 1
            call Catalogue_Assign_byGalaxy_byCatalogue(Cat, nPass, Temp_Cat, i)
         end if
      end do
      if(nPass /= size(Cat%MF606W)) then
         print *, 'Error - Cut_By_magnitude - Error in Assigning galxies that pass cuts', nPass, size(Cat%MF606W)
         STOP
      end if

      print *, 'Cut Catalogue by magnitude between limits:', ilower, iupper, '; nCut = ', size(Temp_Cat%RA) - nPass

      call Catalogue_Destruct(Temp_Cat)

    end subroutine Cut_By_Magnitude

    subroutine Cut_By_SNR(Cat, lowerCut, upperCut, Discarded)
      type(Catalogue)::Cat
      real(double),intent(in),optional::lowerCut, upperCut
      type(Catalogue), intent(out), optional:: Discarded

      real(double)::ilower, iupper
      type(Catalogue)::Temp_Cat
      integer::i
      integer::nPass, nFail
      real(double), dimension(size(Cat%RA)):: SNR

      if((present(lowerCut)==.false.) .and. (present(upperCut)==.false.)) then
         print *, 'Error - Cut_By_SNR - Either a lower or upper cut needs to be entered, returning without cutting'
         return
      end if

      if(all(Cat%flux <= 0.e0_double) .or. all(Cat%fluxerr <= 0.e0_double)) then
         print *, '--------- Cut_By_SNR - Flux or Flux_Err not set, returning without cut'
         print *, ' '
         return
      end if

      !---Construct mSNR
      SNR = Cat%flux/Cat%fluxerr

      if(present(lowerCut)) then
         ilower =lowerCut
      else
         ilower = 0.0
      end if
      if(present(upperCut)) then
         iupper =upperCut
      else
         iupper = 10000.
      end if

      if(ilower >= iupper) then
         STOP 'FATAL ERROR - Cut_By_SNR - lower cut larger than upper cut'
      end if
      
      Temp_Cat = Cat

      call Catalogue_Destruct(Cat)

      nPass = 0
      do i = 1, size(Temp_Cat%RA)
         if(isNaN(SNR(i))) cycle
         if((SNR(i) >= ilower) .and. (SNR(i) <= iupper)) nPass = nPass + 1
      end do

      call Catalogue_Construct(Cat, nPass)
      if(present(Discarded)) then
         call Catalogue_Destruct(Discarded)
         call Catalogue_Construct(Discarded, size(SNR)-nPass)
      end if

      nPass = 0; nFail = 0
      do i = 1, size(SNR)
         if(isNaN(SNR(i))) cycle
         if((SNR(i) >= ilower) .and.(SNR(i) <= iupper)) then
            nPass = nPass + 1
            call Catalogue_Assign_byGalaxy_byCatalogue(Cat, nPass, Temp_Cat, i)
         elseif(present(Discarded)) then
            nFail = nFail + 1
            call Catalogue_Assign_byGalaxy_byCatalogue(Discarded, nFail, Temp_Cat, i)
         end if
      end do
      if(nPass /= size(Cat%Sizes)) then
         print *, 'Error - Cut_By_SNR - Error in Assigning galxies that pass cuts', nPass, size(Cat%Sizes)
         STOP
      end if

      print *, 'Cut Catalogue by SNR between limits:', ilower, iupper, '; nCut = ', size(Temp_Cat%RA) - nPass

      call Catalogue_Destruct(Temp_Cat)

    end subroutine Cut_By_SNR
    
    subroutine Cut_By_PixelSize(Cat, lowerCut, upperCut, Discarded)
      type(Catalogue)::Cat
      real(double),intent(in),optional::lowerCut, upperCut
      type(Catalogue), intent(out), optional:: Discarded

      real(double)::ilower, iupper
      type(Catalogue)::Temp_Cat
      integer::i
      integer::nPass, nFail

      if((present(lowerCut)==.false.) .and. (present(upperCut)==.false.)) then
         print *, 'Error - Cut_By_PixelSize - Either a lower or upper cut needs to be entered, returning without cutting'
         return
      end if

      if(present(lowerCut)) then
         ilower =lowerCut
      else
         ilower = minval(Cat%Sizes)
      end if
      if(present(upperCut)) then
         iupper =upperCut
      else
         iupper = maxval(Cat%Sizes)
      end if

      if(ilower >= iupper) then
         STOP 'FATAL ERROR - Cut_By_PixelSize - lower cut larger than upper cut'
      end if
      
      Temp_Cat = Cat

      call Catalogue_Destruct(Cat)

      nPass = 0
      do i = 1, size(Temp_Cat%Sizes)
         if(isNaN(Temp_Cat%Sizes(i))) cycle
         if((Temp_Cat%Sizes(i) >= ilower) .and. (Temp_Cat%Sizes(i) <= iupper)) nPass = nPass + 1
      end do

!      nPass = count(Temp_cat%Sizes <= iupper) + count(Temp_Cat%Sizes >= ilower) - size(Temp_Cat%Sizes) - count(isNaN(Temp_cat%Sizes))


      call Catalogue_Construct(Cat, nPass)
      if(present(Discarded)) then
         call Catalogue_Destruct(Discarded)
         call Catalogue_Construct(Discarded, size(Temp_Cat%RA)-nPass)
      end if
      

      nPass = 0; nFail = 0
      do i = 1, size(Temp_Cat%Sizes)
         !--TESTING< THIS NEEDS REMOVED
         !if(i<=10) print *, 'cut_by_pixel size: removing D = -1 from both samples: THIS NEEDS DEALT WITH'
         !if(isNaN(Temp_Cat%Sizes(i))) cycle

         if((Temp_Cat%Sizes(i) >= ilower) .and.(Temp_Cat%Sizes(i) <= iupper)) then
            nPass = nPass + 1
            call Catalogue_Assign_byGalaxy_byCatalogue(Cat, nPass, Temp_Cat, i)
         elseif(present(Discarded)) then
            nFail = nFail + 1
            call Catalogue_Assign_byGalaxy_byCatalogue(Discarded, nFail, Temp_Cat, i)
         end if
      end do
      if(nPass /= size(Cat%Sizes)) then
         print *, 'Error - Cut_By_PixelSize - Error in Assigning galxies that pass cuts', nPass, size(Cat%Sizes)
         STOP
      end if
!!$      if(present(Discarded) .and. nFail /= size(Discarded%Sizes)) then
!!$         !--Previously this has been due to the presence of NaNs
!!$         print *, 'FATAL Error - Cut_By_PixelSize - Error in Assigning galxies that fail cuts', nFail, size(Discarded%RA)
!!$         print *, 'nFail, nFail+nPass, Catalogue Size:', nFail, nPass+nFail, size(Temp_Cat%RA)
!!$         STOP
!!$      end if


      print *, 'Cut Catalogue by size between limits:', ilower, iupper, '; nCut = ', size(Temp_Cat%RA) - nPass

      call Catalogue_Destruct(Temp_Cat)

    end subroutine Cut_By_PixelSize
    
    subroutine PSF_Correction(Cat, Correction_Type)
      !--Applies size correction by correction type--!
      type(Catalogue)::Cat
      integer::Correction_Type !-0:None, 1:Pixel Size Cut, 2:Subtract PSF size in quadrature-!

      real(double)::Pixel_Size_Cut(2)
      type(Catalogue)::Temporary_Cat
      integer::ii, i

      if(Correction_Type == 0) return
      if(any((/0,1,2/) == Correction_Type)==.false.) print *, 'WARNING - PSF_Correction - Correction Type not supported' 
      
      !--Main body, applies correction--!
      select case(Correction_Type)
      case(1) !-Cut by Pixel_Size_Cut-!
         !-Set pixel size cut, lower to reduce PSF Corrections, Upper to remove outliers-!
         Pixel_Size_Cut = (/4.e0_double, 20.e0_double/)
         call Clip_Sizes(Cat, Pixel_Size_Cut)
      end select

    end subroutine PSF_Correction

    subroutine Clip_Sizes(Cat, SizeClip)
      type(catalogue),intent(in)::Cat
      real(double),intent(in),dimension(2)::SizeClip !-Lower, Upper-!
      
      print *, 'Clipping Size by:', SizeClip

      call Cut_By_PixelSize(Cat, SizeClip(1), SizeClip(2))
    end subroutine Clip_Sizes

    function global_mean_size(Cat, by_Size_Type)
      !-Returns the mean size of galaxies in the catalogue, as 1/N*sum(size)-!
      type(Catalogue)::Cat
      real(double)::global_mean_size
      character(*), intent(in)::by_Size_Type

      if(by_Size_Type == 'Physical' .or. by_Size_Type == 'physical' .or. by_Size_Type == 'Phys'.or. by_Size_Type == 'phys') then
         global_mean_size = sum(Cat%Physical_Sizes)/size(Cat%Physical_Sizes)
      elseif(by_Size_Type == 'Pixel' .or. by_Size_Type == 'pixel' .or. by_Size_Type == 'Pix'.or. by_Size_Type == 'pix') then
         global_mean_size = sum(Cat%Sizes)/size(Cat%Sizes)
      else
         STOP 'global_mean_size - invalid size type entered'
      end if

    end function global_mean_size

    function size_variance(Cat, by_Size_Type, global_mean)
      !-Returns the mean size of galaxies in the catalogue, as 1/N*sum(size)-!                                                                                                                                                               
      use Statistics, only:variance_discrete
      type(Catalogue)::Cat
      real(double)::size_variance
      character(*), intent(in)::by_Size_Type
      real(double), intent(in),optional::global_mean

      real(double):: mean_size_global

      if(present(global_mean)) then
         mean_size_global = global_mean
      else
!         print *, 'Size Varaince: Calculating mean from Catalogue entered'
         mean_size_global =  global_mean_size(Cat, by_Size_Type)
      end if

      if(by_Size_Type == 'Physical' .or. by_Size_Type == 'physical' .or. by_Size_Type == 'Phys'.or. by_Size_Type == 'phys') then
         size_variance =  variance_discrete(Cat%Physical_Sizes, Cat%Physical_Sizes, mean_size_global, mean_size_global)
      elseif(by_Size_Type == 'Pixel' .or. by_Size_Type == 'pixel' .or. by_Size_Type == 'Pix'.or. by_Size_Type == 'pix') then
         size_variance =  variance_discrete(Cat%Sizes, Cat%Sizes, mean_size_global, mean_size_global)
      else
         STOP 'global_mean_size - invalid size type entered'
      end if

    end function size_variance

    !--Catalogue IO----!

    subroutine Catalogue_ReadIn(Cat, Cat_Filename, Size_label, Cols, get_lnSize)
      use IO, only: ReadIn, IO_Verbose, skip_header; use Strings_lib, only:remove_FileExtension_fromString
      !-Direct ReadIn of Catalogue File-!
      !--if get_lnSize is present and true, then return a catalogue where Sizes = lnT
      type(Catalogue)::Cat
      character(*),intent(in)::Cat_Filename
      logical, intent(in),optional:: get_lnSize

      real(double),allocatable::Cat_2D(:,:)

      integer,intent(in)::Cols(:) !- Number, ntile, RA, Dec, xpos, ypos, Absolute_Magnitude, MF606W, MagErr, Flux, FluxErr, Size, g1, g2, redshift-!
      integer::Cols_Length = 15
      integer, allocatable::iCols(:)
      integer::Column_Index, i
      character(500)::header_Filename
      logical::here

      character(*)::Size_Label

      !--Cols describes where the code is to find the relelvant information--!
      if( all(Cols == 0) ) then
         !--Attempt read in from header--!
         header_Filename = remove_FileExtension_fromString(Cat_Filename)
         if(len(header_filename) == 1) STOP 'Catalogue_ReadIn - Error with constructing header filename, length of 1'
         header_Filename = trim(adjustl(header_Filename(1:len(trim(adjustl(header_filename)))-1)))//'_header.dat'
         inquire(file = header_filename, exist = here)
         if(here) then
            !-Read in cols-!
            open(37, file = header_filename)
            call skip_header(37, '#')
            allocate(iCols(Cols_Length)); iCols = -100
            read(37, *) iCols
            if(any(iCols == -100)) STOP 'Catalogue_ReadIn - Error  in iCols header readin, not all columns are accounted for'
         else
            print *, 'Header_Filename:', Header_Filename
            STOP 'Catalogue_ReadIn - Header_Filename does not exist'
         end if
      else
         if(size(Cols) /= Cols_Length) then 
            print *, 'Catalogue Read In - Cols Input is not of the correct size, please input (any < 0 will be ignored):'
            print *, 'ntile, RA, Dec, xpos, ypos, Mag, MagErr, Flux, FluxErr, Size, g1, g2'
            allocate(iCols(Cols_Length)); iCols = -1
            read(*,*) iCols
         else
            allocate(iCols(size(Cols))); iCols = Cols
         end if
      end if

      if(check_Catalogue_Existence(trim(adjustl(Cat_Filename))) == .false.) then
         print *, 'Catalogue:', trim(adjustl(Cat_Filename))
         STOP 'Catalogue readin: Catalogue does not exist'
      end if

      print *, 'Creating Catalogues with Cols:', iCols

      !--Read into 2D format--!
      call ReadIn(Cat_2D, filename  = trim(adjustl(Cat_Filename)), tabbed = .false., header_label = '#')
      !--This outputs as Cat2D(cols, rows), so 2nd dimension lists all galaxies in catalogue--!
      print *, 'Catalogue read in from:', trim(adjustl(Cat_Filename))

      print *, 'Constructing Catalogue with:', size(Cat_2d,2), ' entries'

      !-Split 2D format into permanent Catalogue Storage--!
      call Catalogue_Construct(Cat, size(Cat_2D,2))
      Cat%Name = trim(adjustl(Cat_Filename))

      PRINT *, 'Using Cols:', icols

      do Column_Index = 1, size(iCols)
         if(iCols(Column_Index) < 0) cycle
         select case (Column_Index)
         case(1)
            Cat%Galaxy_Number = Cat_2D(iCols(Column_Index),:)
         case(2) !-ntile-!
            Cat%ntile = Cat_2D(iCols(Column_Index),:)
         case(3) !-RA-!
            Cat%RA = Cat_2D(iCols(Column_Index),:)
         case(4) !-Dec-!
            Cat%Dec = Cat_2D(iCols(Column_Index),:)
         case(5) !-xpos-!
            Cat%xpos = Cat_2D(iCols(Column_Index),:)
         case(6) !-ypos-!
            Cat%ypos = Cat_2D(iCols(Column_Index),:)
         case(7) !-Absolute_Mag-!
            Cat%Absolute_Magnitude = Cat_2D(iCols(Column_Index),:)
         case(8)
            Cat%MF606W = Cat_2D(iCols(Column_Index),:)
         case(9) !-MagErr-!
            Cat%MagErr = Cat_2D(iCols(Column_Index),:)
         case(10) !-Flux-!
            Cat%Flux = Cat_2D(iCols(Column_Index),:)
         case(11) !-FluxErr-!
            Cat%FluxErr = Cat_2D(iCols(Column_Index),:)
         case(12) !-Size-!
            Cat%Sizes_Label= Size_Label
            Cat%Sizes = Cat_2D(iCols(Column_Index),:)
         case(13) !-g1-!
            Cat%g1 = Cat_2D(iCols(Column_Index),:)
         case(14) !-g2-!
            Cat%g2 = Cat_2D(iCols(Column_Index),:)
         case(15)
            Cat%Redshift = Cat_2D(iCols(Column_Index),:)
         case default
            STOP 'iCols constructed to a size which I cannot deal with yet'
         end select
      end do

      !--Calculate Physical_sizes for the Catalogue--!
!      call convert_Size_from_Pixel_to_Physical(Cat)

      Cat_2D = 0.e0_double; deallocate(Cat_2D)

      !--- Set Galaxy Number [id] (for now this overwrites any input)
      Cat%Galaxy_Number = (/(i, i = 1, size(Cat%RA))/)

      !print *, 'Alpplying log-Size?:', (present(get_lnSize) .and. get_lnSize), minval(Cat%Sizes), maxval(Cat%Sizes)

      if(present(get_lnSize)) then

         if(get_lnSize) then
            print *, ' '
            print *, '****************************'
            print *, '~~~ Catalogue uses log-Sizes'
            print *, '****************************'
            print *, ' '
            where (Cat%Sizes <= 0.e0_double)
               Cat%Sizes = 0.e0_double/0.e0_double
            elsewhere
               Cat%Sizes = dlog(Cat%Sizes)
            end where
            Cat%log_Sizes = .true.
         end if

      end if

      if(count(isNaN(Cat%MF606W)) > 0) then
         print *, 'Error - Catalogue_ReadIn - ', count(isNaN(Cat%MF606W)), ' NaNs found in magnitude'
         STOP
      END if

      write(*,'(A)') '_______ Catalogue Read in Successful ___________'
      print *, ' '

    end subroutine Catalogue_ReadIn

    subroutine Catalogue_Output(Cat, Output_Filename, Pad)
      type(Catalogue),intent(in)::Cat
      character(len = *),intent(in)::Output_filename
      logical,intent(in),optional::Pad

      logical::iPad !--Pad with zeros?--!
      integer::i

      if(present(Pad)) then
         iPad = Pad
      else
         iPad  = .false.
      end if

      open(unit = 30, file =  trim(adjustl(Output_Filename)))
      write(30, '(A)') '## 1:nTile; 2:RA; 3:Dec; 4:x; 5:y; 6:Magnitude; 7:MagError; 8: Flux; 9: FluxError; 10: Size; 11:g1; 12:g2'
      do i = 1, size(Cat%Sizes)
         if(iPad) then
            if(trim(adjustl(Cat%Sizes_Label))=='FWHM') then
               write(30, '(I2, x, 20(e14.7, x))') Cat%ntile(i), Cat%RA(i), Cat%Dec(i), Cat%xpos(i), Cat%ypos(i), Cat%MF606W(i), Cat%MagErr(i), Cat%Flux(i), Cat%FluxErr(i), Cat%Sizes(i), 0., 0., 0., 0., 0., 0., 0., 0., Cat%g1(i), Cat%g2(i)
            elseif(trim(adjustl(Cat%Sizes_Label))=='FR') then
               write(30, '(I2, x, 20(e14.7,x))') Cat%ntile(i), Cat%RA(i), Cat%Dec(i), Cat%xpos(i), Cat%ypos(i), Cat%MF606W(i), Cat%MagErr(i), Cat%Flux(i), Cat%FluxErr(i),0., Cat%Sizes(i),  0., 0., 0., 0., 0., 0., 0., Cat%g1(i), Cat%g2(i)
            elseif(trim(adjustl(Cat%Sizes_Label))=='KSB') then
               write(30, '(I2, x, 20(e14.7,x))') Cat%ntile(i), Cat%RA(i), Cat%Dec(i), Cat%xpos(i), Cat%ypos(i), Cat%MF606W(i), Cat%MagErr(i), Cat%Flux(i), Cat%FluxErr(i), 0., 0., Cat%Sizes(i), 0., 0., 0., 0., 0., 0., Cat%g1(i), Cat%g2(i)
            end if
         else
            write(30, '(I2, x, 11(e14.7,x))') Cat%ntile(i), Cat%RA(i), Cat%Dec(i), Cat%xpos(i), Cat%ypos(i), Cat%MF606W(i), Cat%MagErr(i), Cat%Flux(i), Cat%FluxErr(i), Cat%Sizes(i), Cat%g1(i), Cat%g2(i)
         end if
      end do
      close(30)

    end subroutine Catalogue_Output

    !----------------CATALOGUE DERIVED TYPE SUB/FUNC----------------------------!
    subroutine Binned_Catalogue_Construct(BCat, Label, nBin)
      type(Binned_Catalogue),intent(out)::BCat
      character(*),intent(in)::Label
      integer, intent(in)::nBin!#Bins

      integer::i

      BCat%Label = Label
      allocate(BCat%Bin_limits(nBin,2)); BCat%Bin_Limits = 0.e0_double
      allocate(BCat%Occupation(nBin)); BCat%Occupation = 0
      allocate(BCat%Cat(nBin))
      allocate(BCat%Index_Mapping(nBin))

    end subroutine Binned_Catalogue_Construct

    subroutine Binned_Catalogue_Destruct(BCat)
      type(Binned_Catalogue)::BCat

      integer::i

      BCat%Label = ' '
      if(allocated(BCat%Bin_Limits)) deallocate(BCat%Bin_Limits)
      if(allocated(BCat%Occupation)) deallocate(BCat%Occupation)
      if(allocated(BCat%Index_Mapping)) deallocate(BCat%Index_Mapping)
      if(allocated(BCat%Cat)) then
         do i = 1, size(BCat%Cat)
            call Catalogue_Destruct(BCat%Cat(i))
         end do
         deallocate(BCat%Cat)
      end if

    end subroutine Binned_Catalogue_Destruct


    subroutine Catalogue_Destruct(Cat)
      type(Catalogue)::Cat

      Cat%Name = ' '

      if(allocated(Cat%Posterior_Method)) deallocate(Cat%Posterior_Method)
      if(allocated(Cat%RenormalisationGroup)) deallocate(Cat%RenormalisationGroup)
      if(allocated(Cat%Galaxy_Number)) deallocate(Cat%Galaxy_Number)
      if(allocated(Cat%ntile)) deallocate(Cat%ntile)
      if(allocated(Cat%flux)) deallocate(Cat%flux)
      if(allocated(Cat%fluxerr)) deallocate(Cat%fluxerr)
      if(allocated(Cat%Surface_Brightness)) deallocate(Cat%Surface_Brightness)
      if(allocated(Cat%MF606W)) deallocate(Cat%MF606W)
      if(allocated(Cat%magerr)) deallocate(Cat%magerr)
      if(allocated(Cat%Absolute_Magnitude)) deallocate(Cat%Absolute_Magnitude)
      if(allocated(Cat%xpos)) deallocate(Cat%xpos)
      if(allocated(Cat%ypos)) deallocate(Cat%ypos)
      if(allocated(Cat%RA)) deallocate(Cat%RA)
      if(allocated(Cat%Dec)) deallocate(Cat%Dec)
      if(allocated(Cat%Sizes)) deallocate(Cat%Sizes)
      if(allocated(Cat%Physical_Sizes)) deallocate(Cat%Physical_Sizes)
      if(allocated(Cat%Redshift)) deallocate(Cat%Redshift)
      if(allocated(Cat%g1)) deallocate(Cat%g1)
      if(allocated(Cat%g2)) deallocate(Cat%g2)
!      if(allocated(Cat%Size_FR)) deallocate(Cat%Size_KSB)
!      if(allocated(Cat%Size_KSB)) deallocate(Cat%Size_KSB)

    end subroutine Catalogue_Destruct

    subroutine Catalogue_Construct(Cat, nObj)
      type(Catalogue)::Cat
      integer,intent(in)::nObj
       
      !-If Catalogue contains information then destroy-!
      if(Catalogue_Constructed(Cat)) call Catalogue_Destruct(Cat)

      Cat%Name = ' '

      allocate(Cat%Posterior_Method(nObj)); Cat%Posterior_Method = 0
      allocate(Cat%RenormalisationGroup(nObj)); Cat%RenormalisationGroup = 0
      allocate(Cat%Galaxy_Number(nObj)); Cat%Galaxy_Number = 0
      allocate(Cat%flux(nObj)); Cat%flux = -100.e0_double
      Cat%Mag_Label = ' '
      !allocate(Cat%Surface_Brightness(nObj)); Cat%Surface_Brightness = -100.e0_double
      allocate(Cat%MF606W(nObj)); Cat%MF606W = 0.e0_double
      allocate(Cat%Absolute_Magnitude(nObj)); Cat%Absolute_Magnitude = 0.e0_double
      allocate(Cat%xpos(nObj)); Cat%xpos = 0.e0_double
      allocate(Cat%ypos(nObj)); Cat%ypos = 0.e0_double 
      Cat%Sizes_Label = ' '
      allocate(Cat%Sizes(nObj)); Cat%Sizes = -100.e0_double
      allocate(Cat%Physical_Sizes(nObj)); Cat%Physical_Sizes = 0.e0_double
 
      allocate(Cat%Fluxerr(nObj)); Cat%Fluxerr = -100.e0_double
      allocate(Cat%Magerr(nObj)); Cat%Magerr = 0.e0_double
      allocate(Cat%ntile(nObj)); Cat%ntile = 0
      allocate(Cat%RA(nObj)); Cat%RA = 0.e0_double
      allocate(Cat%Dec(nObj)); Cat%Dec = 0.e0_double
      allocate(Cat%g1(nObj)); Cat%g1 = 0.e0_double
      allocate(Cat%g2(nObj)); Cat%g2 = 0.e0_double
      allocate(Cat%Redshift(nObj)); Cat%Redshift = -100.e0_double
!      allocate(Cat%Size_FR(nObj)); Cat%Size_FR = 0.e0_double
!      allocate(Cat%Size_KSB(nObj));  Cat%Size_KSB = 0.e0_double

    end subroutine Catalogue_Construct

    logical function Catalogue_Constructed(Cat)
      type(Catalogue),intent(in)::Cat


      Catalogue_Constructed = .false.

      if(Cat%Name /= ' ') Catalogue_Constructed = .true.

      Catalogue_Constructed = .false.
      if(allocated(Cat%Posterior_Method)) Catalogue_Constructed = .true.
      if(allocated(Cat%RenormalisationGroup)) Catalogue_Constructed = .true.
      if(allocated(Cat%Galaxy_Number)) Catalogue_Constructed = .true.
      if(allocated(Cat%flux)) Catalogue_Constructed = .true.
      if(allocated(Cat%Surface_Brightness)) Catalogue_Constructed = .true.
      if(allocated(Cat%MF606W)) Catalogue_Constructed = .true.
      if(allocated(Cat%Absolute_Magnitude)) Catalogue_Constructed = .true. 
      if(allocated(Cat%xpos)) Catalogue_Constructed = .true.
      if(allocated(Cat%ypos)) Catalogue_Constructed = .true.
      if(allocated(Cat%Sizes)) Catalogue_Constructed = .true.
      if(allocated(Cat%Physical_Sizes)) Catalogue_Constructed = .true.
!      if(allocated(Cat%Size_FR)) Catalogue_Constructed = .true.
!      if(allocated(Cat%Size_KSB)) Catalogue_Constructed = .true.
      if(allocated(Cat%ntile)) Catalogue_Constructed = .true.
      if(allocated(Cat%fluxerr)) Catalogue_Constructed = .true.
      if(allocated(Cat%magerr)) Catalogue_Constructed = .true.
      if(allocated(Cat%RA)) Catalogue_Constructed = .true.
      if(allocated(Cat%Dec)) Catalogue_Constructed = .true.
      if(allocated(Cat%g1)) Catalogue_Constructed = .true.
      if(allocated(Cat%g2)) Catalogue_Constructed = .true.
      if(allocated(Cat%Redshift)) Catalogue_Constructed = .true.

    end function Catalogue_Constructed

    subroutine Catalogue_Assign_byGalaxy_byIndex(Cat, Galaxy_Index, flux, mag, xpos, ypos, Sizes,Status)
      !--DEPRECATED CODE - OUT OF DATE--!
      !-Assigns a single value by entered index. Function is overloaded with equivalent which assigns by Catalogue ()-!
      type(Catalogue)::Cat
      integer,intent(in)::Galaxy_Index
      real(double),intent(in),optional::flux, mag, xpos,ypos, Sizes
      integer,optional::Status ! "-1 Failed", "0 Unresolved", "1 Success" !

      integer::iStatus

!!$  INTERFACE
!!$     subroutine Catalogue_Assign_byGalaxy_byIndex(Cat, Galaxy_Index, flux, mag, xpos,ypos, Size_FWHM, Size_FR, Size_KSB, Status)
!!$       use Param_Types
!!$       type(Catalogue)::Cat
!!$       integer,intent(in)::Galaxy_Index
!!$       
!!$       real(double),intent(in),optional::flux, mag, xpos,ypos, Size_FWHM, Size_FR, Size_KSB
!!$       integer,optional::Status
!!$     end subroutine Catalogue_Assign_byGalaxy_byIndex
!!$  END INTERFACE


      if(Galaxy_Index > size(Cat%MF606W)) then; print *, 'Catalogue_Assign_byGalaxy - Index Entered not valid - Too Large'; iStatus = -1; return;
      elseif(Galaxy_Index <= 0) then; print *, 'Catalogue_Assign_byGalaxy - Index Entered not valid - Too Small'; iStatus = -1; return;
      end if

      if(present(Status)) Status = iStatus
      if(iStatus < 0) return

      if(present(flux)) Cat%flux(Galaxy_Index) = flux
      if(present(mag)) Cat%MF606W(Galaxy_Index) = mag
      if(present(xpos)) Cat%xpos(Galaxy_Index) =xpos
      if(present(ypos)) Cat%ypos(Galaxy_Index) =ypos
      if(present(Sizes)) Cat%Sizes(Galaxy_Index) =Sizes
      !if(present(flux)) Cat%Size_FR(Galaxy_Index) =Size_FR
      !if(present(flux)) Cat%Size_KSB(Galaxy_Index) =Size_KSB

    end subroutine Catalogue_Assign_byGalaxy_byIndex

    subroutine Catalogue_Assign_byGalaxy_byCatalogue(Cat, Index, Cat_Ref, Index_Ref)
      !--Equates a galaxy at position: Index_Ref in Catalogue: Cat_Ref to: Index in Cat
      !--NEEDS TO BE KEPT UP TO DATE WITH ADDITIONS/SUBTRACTIONS TO THE CATALOGUE DERIVED TYPE DEFINITION--!
      type(Catalogue), intent(inout)::Cat
      integer,intent(in)::Index
      type(Catalogue), intent(in)::Cat_Ref
      integer,intent(in)::Index_Ref
      
      integer::iStatus

      if(Index > size(Cat%RA)) then; print *, 'Catalogue_Assign_byGalaxy_byCatalogue - Index Entered not valid - Too Large', Index, size(Cat%RA); iStatus = -1; return;
      elseif(Index <= 0) then; print *, 'Catalogue_Assign_byGalaxy - Index Entered not valid - Too Small'; iStatus = -1; return;
      end if

      if(Index_Ref > size(Cat_Ref%RA)) then; print *, 'Catalogue_Assign_byGalaxy - Index_Ref Entered not valid - Too Large'; iStatus = -1; return;
      elseif(Index_Ref <= 0) then; print *, 'Catalogue_Assign_byGalaxy - Index_Ref Entered not valid - Too Small'; iStatus = -1; return;
      end if

      Cat%Galaxy_Number(Index) = Cat_Ref%Galaxy_Number(Index_Ref)
      Cat%Posterior_Method(Index) = Cat_Ref%Posterior_Method(Index_Ref)
      Cat%RenormalisationGroup(Index) = Cat_Ref%RenormalisationGroup(Index_Ref)
      Cat%flux(Index) = Cat_Ref%flux(Index_Ref)
      Cat%fluxerr(Index) = Cat_Ref%fluxerr(Index_Ref)
      !Cat%Surface_Brightness(Index) = Cat_Ref%Surface_Brightness(Index_Ref)
      Cat%Absolute_Magnitude(Index) = Cat_Ref%Absolute_Magnitude(Index_Ref)
      Cat%MF606W(Index) = Cat_Ref%MF606W(Index_Ref)
      Cat%magerr(Index) = Cat_Ref%magerr(Index_Ref)
      Cat%xpos(Index) = Cat_Ref%xpos(Index_Ref)
      Cat%ypos(Index) = Cat_Ref%ypos(Index_Ref)
      Cat%RA(Index) = Cat_Ref%RA(Index_Ref)
      Cat%Dec(Index) = Cat_Ref%Dec(Index_Ref)
      Cat%Sizes(Index) = Cat_Ref%Sizes(Index_Ref)
      Cat%Physical_Sizes(Index) = Cat_Ref%Physical_Sizes(Index_Ref)
      !Cat%Size_FR(Index) = Cat_Ref%Size_FR(Index_Ref)
      !Cat%Size_KSB(Index) = Cat_Ref%Size_KSB(Index_Ref)
      Cat%g1(Index) = Cat_Ref%g1(Index_Ref)
      Cat%g2(Index) = Cat_Ref%g2(Index_Ref)
      Cat%Redshift(Index) = Cat_Ref%Redshift(Index_Ref)

    end subroutine Catalogue_Assign_byGalaxy_byCatalogue

    function get_Total_Corrected_Ellipticity_g(Cat)
      type(Catalogue)::Cat
      
      real(double),dimension(:),allocatable::get_Total_Corrected_Ellipticity_g

      integer::i

      if(Catalogue_Constructed(Cat)==.false.) STOP 'FATAL ERROR - et_Total_Corrected_Ellipticity_g - Catalogue entered not constructed'

      allocate(get_Total_Corrected_Ellipticity_g(size(Cat%g1))); get_Total_Corrected_Ellipticity_g = 0.e0_double

      do i = 1, size(Cat%g1)
         get_Total_Corrected_Ellipticity_g(i) =dsqrt( Cat%g1(i)*Cat%g1(i) + Cat%g2(i)*Cat%g2(i) )
      end do
    
    end function get_Total_Corrected_Ellipticity_g

    subroutine get_redshift_Information_combine_Catalogues(Cat, RCat, Default_Redshift)
      !--Appends the redshift information in RCat onto Cat using RA/Dec Information--!
      type(Catalogue),intent(in)::RCat
      type(Catalogue),intent(inout)::Cat
      real(double),optional::Default_Redshift

      real(double)::iDz!-Internal default Redshift
      real(double)::Position_Tol = 1.0 !--In arcseconds--!
      real(double)::Pos_Tol_Degree != Position_Tol/(3600.e0_double)

      integer::c, rc
      integer::nMatch

      integer,dimension(Size(Cat%Redshift)):: Index_Mapping

      print *, 'Matching Redshift Information to Catalogue'

      Pos_Tol_Degree = Position_Tol/(3600.e0_double)

      !--Set all redshifts to the default--!
      if(present(Default_Redshift)) then
         print *, 'Setting Unassigned Redshifts to the Default Redshift of:', Default_Redshift
         Cat%Redshift = Default_Redshift
      end if
     
      Index_Mapping = -1
      nMatch = 1
      do rc = 1, size(RCat%Redshift)
         do c = 1, size(Cat%Redshift)
            if( (dsqrt( (Cat%RA(c)-Rcat%RA(rc))*(Cat%RA(c)-Rcat%RA(rc)) + (Cat%Dec(c)-Rcat%Dec(rc))*(Cat%Dec(c)-Rcat%Dec(rc)) ) <= Pos_Tol_Degree) .and. (Index_Mapping(c) < 0)) then
               Cat%Redshift(c) = RCat%Redshift(rc)
               Cat%Absolute_Magnitude(c) = RCat%Absolute_Magnitude(rc)
               nMatch = nMatch + 1
               Index_Mapping(c) = rc
               exit
            end if
         end do
      end do

      print *, 'Of ', size(RCat%Redshift), 'galaxies with redshifts, ', nMatch, ' have been matched'

    end subroutine get_redshift_Information_combine_Catalogues


    !"---------------------------------------DEPRECATED CODE------------------------------------------------------------------------!



    !--COMBO17 Redshift Routines---!

    subroutine match_COMBO17_Redshifts(Cat, return_only_Matches)
      !-Matches redshift information from COMBO17 file to those galaxies contained in the STAGES catalogue-!
      !-30 Aug 13 -- this version uses Catherine's pre-matched file. Deprecated version which uses the COMBO17 data in Deprecated_Code File
      !-if "return_only_Matches" is present and true, the catalogue is reduced to only those glaxies with redshift information-!
      type(Catalogue)::Cat
      logical,optional::return_only_Matches

      logical::here
      character(120)::COMBO17_Filename = trim(adjustl(Cat_Dir))//'STAGES_shear_pz_matched.dat'
      real(double)::COMBO17_Input(10788, 15)

      integer::nMatch_Single, nMatch_Total
      integer::n_Multiple_Matches, nFailed_Matches
      integer::i,j, COMBO_Index, STAGES_Index, counter
      integer::COMBO_z_Index = 6
      character(200)::COMBO17_Line

      inquire(file = COMBO17_Filename, exist = here)
      if(here == .false.) then
         print *, 'Filename:', COMBO17_Filename
         STOP 'match_COMBO17_Redshifts - COMBO17 redshift file not found - check existence, stopping'
      end if

      !--Read in the COMBO17 file--!
      print *, 'Reading in COMBO17 file...'
      open(unit = 35, file = trim(adjustl(COMBO17_Filename)))
      !-Read Line by Line-!
      i = 0; j = 0; counter = 0
      do while (j < size(COMBO17_Input,1))
         counter = counter + 1
         if(counter <= 16) then
            read(35, *) COMBO17_Line
            cycle
         end if
         j = j+1
         read(35, *, IOSTAT=i) COMBO17_Input(j,:)
         if(i < 0) exit
      end do
      if(j /= size(COMBO17_Input,1)) STOP 'ERROR - COMBO17 input - Input File Length not as expected'
      i = 0; j = 0
      close(35)
      print *, 'Done.'

      if(allocated(Cat%Redshift)) deallocate(Cat%Redshift)
      allocate(Cat%Redshift(size(Cat%Sizes))); Cat%Redshift = -1.e0_double !-Default < 0, meaning unassigned-!
      nMatch_Total = 0; n_Multiple_Matches = 0; nFailed_Matches =  0
      do COMBO_Index = 1, size(COMBO17_Input,1) !-Loop through all COMBO17 Galaxies-!
         nMatch_Single = 0
         do STAGES_Index = 1, size(Cat%RA)
            !-Match on both RA and Dec to within a tolerance-!
            if( (dabs(Cat%RA(STAGES_Index)-COMBO17_Input(COMBO_Index,11)) <= 1.e-8_double) .and. (dabs(Cat%Dec(STAGES_Index)-COMBO17_Input(COMBO_Index,12)) <= 1.e-8_double)) then
               nMatch_Single = nMatch_Single + 1
               if(Cat%Redshift(STAGES_Index) < 0.e0_double) nMatch_Total = nMatch_Total + 1
               Cat%Redshift(STAGES_Index) = COMBO17_Input(COMBO_Index,COMBO_z_Index)
               
               !--Test shear measurements to ensure they match those held in my code--!
               if(dabs(Cat%g1(STAGES_Index) - COMBO17_Input(COMBO_Index,14)) > 1.e-0_double) PRINT *, 'WARNING - SHEARS (g1) DO NOT MATCH FOR MATHCED GALAXY (combo_index, stages_index):', COMBO_Index, STAGES_Index, Cat%g1(STAGES_Index), COMBO17_Input(COMBO_Index,14)
!               if(dabs(Cat%g2(STAGES_Index) - COMBO17_Input(COMBO_Index,15)) > 1.e-0_double) PRINT *,  'WARNING - SHEARS (g2) DO NOT MATCH FOR MATHCED GALAXY (combo_index, stages_index):', COMBO_Index, STAGES_Index
               exit
            end if
         end do
         if(nMatch_Single == 0) then
!            print *, 'WARNING - A COMBO17 Galaxy has not been matched to a STAGES galaxy'
            nFailed_Matches = nFailed_Matches + 1
            print *, 'Catherines mathced galaxy, line:', COMBO_Index+16, ' failed to find a match', COMBO17_Input(COMBO_Index,:)
         end if
         if(nMatch_Single > 1) then
!            print *, 'WARNING - A COMBO17 Galaxy has been matched to MULTIPLE STAGES galaxies'
            n_Multiple_Matches = n_Multiple_Matches + 1
         end if
      end do
      if(nFailed_Matches > 0) print *, 'WARNING:', nFailed_Matches,' COMBO17 galaxies failed to find a match'
      if(n_Multiple_Matches > 0) print *, 'WARNING:', n_Multiple_Matches,' COMBO17 galaxies found multiple matches'
      print *, nMatch_Total, ' of ', size(COMBO17_Input,1), ' galaxies where matched to STAGES galaxies'
      STOP 'STOPPING FOR TESTING'

      if((present(return_only_Matches)) .and. (return_only_Matches == .true.)) call reduce_Catalogue_toGalaxieswithRedshift(Cat)

    end subroutine match_COMBO17_Redshifts

    subroutine reduce_Catalogue_toGalaxieswithRedshift(Cat)
      !--Reduces the input catalogue to only those with redshift information. Matching is already assumed, and those with Redshift > 0 are considered as having redshift information--!
      type(Catalogue)::Cat
    

    end subroutine reduce_Catalogue_toGalaxieswithRedshift




  end module Catalogues
