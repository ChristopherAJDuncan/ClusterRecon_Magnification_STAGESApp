module Data_Tests
  !---Contains routines that test the sample and data, which will be used in a given analysis---!
  use Catalogues; use Param_Types
  implicit none
  
  
contains

!!$  subroutine Foreground_Contamination_NumberDensity_MultipleGrid(Cat, Ap_Pos, Output_Directory, MaskedSurveyArea)
!!$    use Data_Types
!!$    !--Improvements:
!!$    !~~~~ Edit code to output a global number density, calculated from the entered masked survey area, and source number of catalogue using 3' apertures around the clusters. Be sure to edit the survey area to subtract the area of the masks. --- This will reduce contamination from the clusters, but this may be expected to be small
!!$    type(Catalogue), intent(in):: Cat
!!$    real(double), intent(in):: Ap_Pos(:,:) !--in DEGREES; Apeture, RA/Dec-!
!!$    character(*), intent(in):: Output_Directory
!!$    real(double), intent(in):: MaskedSurveyArea
!!$
!!$    integer:: Ap, i
!!$
!!$    real(double):: Core_Radius = 0.e0_double
!!$    type(Array_of_Arrays1D):: NumberDensity
!!$    type(Array_of_Arrays1D)
!!$!    real(double), allocatable:: NumberDensity(:,:) !-Ap, Seperation-!
!!$    real(double), allocatable:: SeperationGrid(:), SeperationBins(:,:)
!!$    integer:: nSep = 30
!!$    real(double),parameter:: Seperation_Limits(2) = (/0.e0_double, 3.e0_double/) !--In ArcMinutes--!
!!$    real(double),allocatable:: Annulus_Area(:)
!!$
!!$    real(double):: global_Number_Density
!!$
!!$    type(Catalogue):: Annulus_Cat
!!$
!!$    character(3)::fmtstring
!!$
!!$    print *, '***Foreground Contamination tests taken on:', size(Cat%RA), 'galaxies'
!!$
!!$    !--Set up number density grid--!
!!$    allocate(SeperationGrid(nSep)); SeperationGrid = 0.e0_double
!!$    allocate(SeperationBins(nSep, 2)); SeperationBins = 0.e0_double
!!$    allocate(Annulus_Area(size(SeperationGrid))); Annulus_Area = 0.e0_double
!!$    do i =1 , nSep
!!$       if(i == 1) then
!!$          SeperationBins(i,1) = Seperation_Limits(1)
!!$       else
!!$          SeperationBins(i,1) = SeperationBins(i-1,2)
!!$       end if
!!$       SeperationBins(i,2) = SeperationBins(i,1) + ((Seperation_Limits(2)-Seperation_Limits(1))/nSep)
!!$
!!$       SeperationGrid(i) = 0.5e0_double*( SeperationBins(i,2) + SeperationBins(i,1) )
!!$       Annulus_Area(i) = 3.14159e0_double*( (SeperationBins(i,2)**2.e0_double) - (SeperationBins(i,1)**2.e0_double) ) !--Pi*r^2--!
!!$    end do
!!$
!!$    allocate(NumberDensity(size(Ap_Pos,1), size(SeperationGrid))); NumberDensity = 0.e0_double
!!$    do Ap = 1, size(Ap_Pos,1)
!!$       do i = 1, size(SeperationGrid)
!!$          
!!$          NumberDensity(Ap,i) = count( dsqrt( (Cat%RA-Ap_Pos(Ap,1))**2.e0_double +(Cat%Dec-Ap_Pos(Ap,2))**2.e0_double ) <= SeperationBins(i,2)/60.e0_double) - count( dsqrt( (Cat%RA-Ap_Pos(Ap,1))**2.e0_double +(Cat%Dec-Ap_Pos(Ap,2))**2.e0_double ) <= SeperationBins(i,1)/60.e0_double )
!!$!          call Identify_Galaxys_in_Circular_Aperture(Cat, Ap_Pos(Ap,:), SeperationBins(i,2)/60.e0_double, Annulus_Cat, Core_Radius = SeperationBins(i,1)/60.e0_double)
!!$
!!$          NumberDensity(Ap,i) = (1.e0_double*NumberDensity(Ap,i))/Annulus_Area(i)
!!$
!!$          call Catalogue_Destruct(Annulus_Cat)
!!$       end do
!!$    end do
!!$
!!$    !---Output---!
!!$    open(unit = 31, file = trim(Output_Directory)//'Foreground_Contamination_NumberDensity.dat')
!!$    !--Header--!
!!$    write(31, '(A)') '#1,2:Seperation Bin [Annulus Limits], 3+: Number Density of galaxies in annulus'
!!$    write(31, '(A, x, I6, x, A, x, e8.1)') '# nSource: ', size(Cat%Dec)
!!$    write(fmtstring, '(I2)') size(Ap_Pos)+2
!!$    do i = 1, size(SeperationGrid)
!!$       write(31, '('//trim(fmtstring)//'(e9.3,x))') SeperationBins(i,:), NumberDensity(:,i)
!!$    end do
!!$    close(31)
!!$
!!$    write(*,'(2A)') '*** Cluster Contamination Test output to: ', trim(Output_Directory)//'Foreground_Contamination_NumberDensity.dat'
!!$
!!$
!!$    global_Number_Density = size(Cat%Dec)/MaskedSurveyArea
!!$    open(unit = 31, file = trim(Output_Directory)//'Foreground_Contamination_Fraction.dat')
!!$    !--Header--!
!!$    write(31, '(A)') '#1,2:Seperation Bin [Annulus Limits], 3+: Cantamination Fraction (ngal_annulus/ngal_field) of galaxies in annulus'
!!$    write(31, '(A, x, I6, x, A, x, e8.1)') '# nSource: ', size(Cat%Dec), ' global_Number_Density:', global_Number_Density
!!$    write(fmtstring, '(I2)') size(Ap_Pos)+2
!!$    do i = 1, size(SeperationGrid)
!!$       write(31, '('//trim(fmtstring)//'(e9.3,x))') SeperationBins(i,:), NumberDensity(:,i)/global_Number_Density
!!$    end do
!!$    close(31)
!!$
!!$    write(*,'(2A)') '*** Cluster Contamination Test output to: ', trim(Output_Directory)//'Foreground_Contamination_Fraction.dat'
!!$
!!$    
!!$    deallocate(Annulus_Area,SeperationGrid,SeperationBins,NumberDensity)
!!$
!!$  end subroutine Foreground_Contamination_NumberDensity_MultipleGrid


  subroutine Foreground_Contamination_NumberDensity(Cat, Ap_Pos, Output_Directory, MaskedSurveyArea)
    !--Improvements:
    !~~~~ Edit code to output a global number density, calculated from the entered masked survey area, and source number of catalogue using 3' apertures around the clusters. Be sure to edit the survey area to subtract the area of the masks. --- This will reduce contamination from the clusters, but this may be expected to be small
    type(Catalogue), intent(in):: Cat
    real(double), intent(in):: Ap_Pos(:,:) !--in DEGREES; Apeture, RA/Dec-!
    character(*), intent(in):: Output_Directory
    real(double), intent(in):: MaskedSurveyArea

    integer:: Ap, i

    real(double):: Core_Radius = 0.e0_double
    real(double), allocatable:: NumberDensity(:,:) !-Ap, Seperation-!
    real(double), allocatable:: SeperationGrid(:), SeperationBins(:,:)
    integer:: nSep = 40
    real(double),parameter:: Seperation_Limits(2) = (/0.e0_double, 4.e0_double/) !--In ArcMinutes--!
    real(double),allocatable:: Annulus_Area(:)

    real(double):: global_Number_Density

    type(Catalogue):: Annulus_Cat

    character(3)::fmtstring

    print *, '***Foreground Contamination tests taken on:', size(Cat%RA), 'galaxies'

    !--Set up number density grid--!
    allocate(SeperationGrid(nSep)); SeperationGrid = 0.e0_double
    allocate(SeperationBins(nSep, 2)); SeperationBins = 0.e0_double
    allocate(Annulus_Area(size(SeperationGrid))); Annulus_Area = 0.e0_double
    print *, 'Allocated', nSep
    do i =1 , nSep
       if(i == 1) then
          SeperationBins(i,1) = Seperation_Limits(1)
       else
          SeperationBins(i,1) = SeperationBins(i-1,2)
       end if
       SeperationBins(i,2) = SeperationBins(i,1) + ((Seperation_Limits(2)-Seperation_Limits(1))/nSep)

       SeperationGrid(i) = 0.5e0_double*( SeperationBins(i,2) + SeperationBins(i,1) )
       Annulus_Area(i) = 3.14159e0_double*( (SeperationBins(i,2)**2.e0_double) - (SeperationBins(i,1)**2.e0_double) ) !--Pi*r^2--!
    end do

    print *, 'Set 1'

    allocate(NumberDensity(size(Ap_Pos,1), size(SeperationGrid))); NumberDensity = 0.e0_double
    print *, 'Looping over Ap_Pos and SepGrid:', allocated(SeperationGrid), size(Ap_Pos), size(Cat%RA)
    do Ap = 1, size(Ap_Pos,1)
       do i = 1, size(SeperationGrid)
          
          NumberDensity(Ap,i) = count( dsqrt( (Cat%RA-Ap_Pos(Ap,1))**2.e0_double +(Cat%Dec-Ap_Pos(Ap,2))**2.e0_double ) <= SeperationBins(i,2)/60.e0_double) - count( dsqrt( (Cat%RA-Ap_Pos(Ap,1))**2.e0_double +(Cat%Dec-Ap_Pos(Ap,2))**2.e0_double ) <= SeperationBins(i,1)/60.e0_double )
!          call Identify_Galaxys_in_Circular_Aperture(Cat, Ap_Pos(Ap,:), SeperationBins(i,2)/60.e0_double, Annulus_Cat, Core_Radius = SeperationBins(i,1)/60.e0_double)

          NumberDensity(Ap,i) = (1.e0_double*NumberDensity(Ap,i))/Annulus_Area(i)

          call Catalogue_Destruct(Annulus_Cat)
       end do
    end do

    print *, 'Set 2'

    !---Output---!
    open(unit = 31, file = trim(Output_Directory)//'Foreground_Contamination_NumberDensity.dat')
    !--Header--!
    write(31, '(A)') '#1,2:Seperation Bin [Annulus Limits], 3+: Number Density of galaxies in annulus'
    write(31, '(A, x, I6, x, A, x, e8.1)') '# nSource: ', size(Cat%Dec)
    write(fmtstring, '(I2)') size(Ap_Pos)+2
    do i = 1, size(SeperationGrid)
       write(31, '('//trim(fmtstring)//'(e9.3,x))') SeperationBins(i,:), NumberDensity(:,i)
    end do
    close(31)

    write(*,'(2A)') '*** Cluster Contamination Test output to: ', trim(Output_Directory)//'Foreground_Contamination_NumberDensity.dat'


    global_Number_Density = size(Cat%Dec)/MaskedSurveyArea
    open(unit = 31, file = trim(Output_Directory)//'Foreground_Contamination_Fraction.dat')
    !--Header--!
    write(31, '(A)') '#1,2:Seperation Bin [Annulus Limits], 3+: Cantamination Fraction (ngal_annulus/ngal_field) of galaxies in annulus'
    write(31, '(A, x, I6, x, A, x, e8.1)') '# nSource: ', size(Cat%Dec), ' global_Number_Density:', global_Number_Density
    write(fmtstring, '(I2)') size(Ap_Pos)+2
    do i = 1, size(SeperationGrid)
       write(31, '('//trim(fmtstring)//'(e9.3,x))') SeperationBins(i,:), NumberDensity(:,i)/global_Number_Density
    end do
    close(31)

    write(*,'(2A)') '*** Cluster Contamination Test output to: ', trim(Output_Directory)//'Foreground_Contamination_Fraction.dat'

    
    deallocate(Annulus_Area,SeperationGrid,SeperationBins,NumberDensity)

  end subroutine Foreground_Contamination_NumberDensity


end module Data_Tests
