module Size_Histograms
  use catalogues; use param_types 
  implicit none

contains

  !--Plotter routine--!
  subroutine Size_Histogram_Plotter(filename, filename2)
    character(*)::filename
    character(*),optional::filename2
    character(120)::Plotter_Name = 'Size_Histograms_Plotter.py'

    logical::here

    inquire(file = Plotter_Name, exist = here)
    if(here == .false.) STOP 'Size_Histogram_Plotter - FATAL ERROR - Plotter does not exist'

    inquire(file = filename, exist = here)
    if(here == .false.) STOP 'Size_Histogram_Plotter - FATAL ERROR - File to be plotted does not exist'

    if(present(filename2)) then
       call system('python '//trim(adjustl(Plotter_Name))//' '//trim(adjustl(filename))//' '//trim(adjustl(filename2)))
    else
       call system('python '//trim(adjustl(Plotter_Name))//' '//trim(adjustl(filename)))
    end if

  end subroutine Size_Histogram_Plotter



  !--Size Histogram across all catalogue--! 
  subroutine Size_Histogram_Catalogue(Sizes, filename)
    real(double),intent(in)::Sizes(:)
    character(*),intent(in):: filename

    integer,parameter::nBin = 30
    real(double)::Bin_Limits(nBin,2) 
    integer,allocatable::Number_in_Bin(:)

    real(double)::minSize, maxSize

    integer::i,j,counter
    
    minSize = minval(Sizes); maxSize = maxval(Sizes)
    print *, 'Plotting Size histrogram across the whole catalogue, between range:', minsize, maxsize
    !--Set up limits of histogram--!
    do i = 1, nBin
       Bin_Limits(i,:) = (/ minSize+(i-1)*((maxSize-minSize)/(nBin)) , minSize+(i)*((maxSize-minSize)/(nBin)) /)
    end do

    !--Set histogram bin values--!
    allocate(Number_in_Bin(nBin)); Number_in_Bin = 0.e0_double
    do i =1, size(Sizes)
       if(Sizes(i) < 0.e0_double) cycle
       do j =1, nBin
          if( (Sizes(i)>= Bin_Limits(j,1)) .and. (Sizes(i)<Bin_Limits(j,2)) ) then
             Number_in_Bin(j) = Number_in_Bin(j) + 1
             exit
          end if
       end do
    end do

    !--Output--!
    open(file = filename, unit  = 40)
    do j = 1, nBin
       write(40, *) Bin_Limits(j,:), Number_In_Bin(j)
    end do
    close(40)
    print *, 'Output histogram of sizes within circular aperture to: ', trim(adjustl(filename))



  end subroutine Size_Histogram_Catalogue


  subroutine Size_Histogram_Circular_Aperture(Cat, Ap_Position, Ap_Radius, filename)
    type(Catalogue)::Cat
    real(double), intent(in)::Ap_Position(:) !-In RA/Dec-!
    real(double), intent(in)::Ap_Radius !-In the same units as RA and Dec entry-!
    character(*), intent(in)::filename

    real(double),allocatable::Sizes_in_Ap(:) !-Lists all the sizes that fall within the aperture-!
    
    integer,parameter::nBin = 10
    real(double)::Bin_Limits(nBin,2) 
    integer,allocatable::Number_in_Bin(:)

    real(double)::minSize, maxSize

    integer::i,j,counter

    integer::Expected_Number_in_Aperture
    real(double)::Mean_of_Catalogue, Mean_of_Aperture

    print *, 'Attempting to get a size histogram in filename:', trim(filename)

    Mean_of_Catalogue = sum(Cat%Physical_sizes)/size(Cat%Physical_sizes)

    minSize = 1.e10_double; maxSize = 0.e0_double

    Expected_Number_in_Aperture = 0
    do i = 1,  size(Cat%Physical_Sizes)
       if( dsqrt( (Cat%RA(i)-Ap_Position(1))**2.e0_double +(Cat%Dec(i)-Ap_Position(2))**2.e0_double ) <= Ap_Radius ) then
          Expected_Number_in_Aperture = Expected_Number_in_Aperture + 1
       end if
    end do

!    Expected_Number_in_Aperture = count( dsqrt( (Cat%RA-Ap_Position(1))**2.e0_double +(Cat%Dec-Ap_Position(2))**2.e0_double ) <= Ap_Radius)
    print *, 'Expected Number in Aperture = ', Expected_Number_in_Aperture

    allocate(Sizes_in_Ap(Expected_Number_in_Aperture)); Sizes_In_Ap = -100.e0_double

    counter = 1
    do i =1, size(Cat%Physical_Sizes)
       if( Cat%Physical_Sizes(i) <= 0.e0_double ) cycle
       if( dsqrt( (Cat%RA(i)-Ap_Position(1))**2.e0_double +(Cat%Dec(i)-Ap_Position(2))**2.e0_double ) <= Ap_Radius ) then !--Assumes that both are in the same units-!
          Sizes_in_Ap(counter) = Cat%Physical_Sizes(i)
          counter = counter + 1
          if( Cat%Physical_Sizes(i) < minSize ) minSize = Cat%Physical_Sizes(i)
          if( Cat%Physical_Sizes(i) > maxSize ) maxSize = Cat%Physical_Sizes(i)
       end if
    end do

    if((counter-1) < size(Sizes_in_Ap)-1) then !-1 to allow for rounding errors-!
       print *, counter-1, size(Sizes_in_Ap)
       STOP 'Size_Histogram_Circular_Aperture - FATAL ERROR - error in assigning sizes within aperture'
    end if
!    if(any(Sizes_in_Ap < 0.e0_double)) STOP 'Size_Histogram_Circular_Aperture - FATAL ERROR - negatives in the list of sizes in aperture'

    Mean_of_Aperture = sum(Sizes_in_Ap)/size(Sizes_in_Ap)

    print *, '!---------------------------------------------------------------------!'
    print *, 'Mean (of Catalogue/ in Aperture):', Mean_of_Catalogue, Mean_of_Aperture
    print *, 'Expected Convergence:', (Mean_of_Aperture/Mean_of_Catalogue) - 1.e0_double
    print *, '!---------------------------------------------------------------------!'

    !--Set up limits of histogram--!
    do i = 1, nBin
       Bin_Limits(i,:) = (/ minSize+(i-1)*((maxSize-minSize)/(nBin)) , minSize+(i)*((maxSize-minSize)/(nBin)) /)
    end do
    !--Adjust top limit to account for the fact that some galaxies will lie exaclty on this limit-!
    Bin_Limits(nBin,2) = Bin_Limits(nBin,2) + 0.1*(Bin_Limits(nBin,2)-Bin_Limits(nBin,1))

    allocate(Number_in_Bin(nBin)); Number_in_Bin = 0.e0_double
    do i =1, size(Sizes_in_Ap)
       if(Sizes_in_Ap(i) < 0.e0_double) cycle
       do j =1, nBin
          if( (Sizes_In_Ap(i)>= Bin_Limits(j,1)) .and. (Sizes_In_Ap(i)<Bin_Limits(j,2)) ) then
             Number_in_Bin(j) = Number_in_Bin(j) + 1
             exit
          end if
       end do
    end do
    Number_in_Bin(size(Number_in_Bin)) = Number_in_Bin(size(Number_in_Bin)) + count(Sizes_In_Ap == Bin_Limits(size(Bin_Limits,1),2))


!!$    if(Sum(Number_in_Bin) /= Expected_Number_in_Aperture) then 
!!$       print *, 'Expected, Actual:', Expected_Number_in_Aperture, Sum(Number_in_Bin)
!!$       STOP 'Size_Histogram_Circular_Aperture - FATAL ERROR - differences in the expected number of galaxies in the aperture'
!!$    end if

    !--Output--!
    open(file = filename, unit  = 40)
    do j = 1, nBin
       write(40, *) Bin_Limits(j,:), Number_In_Bin(j)
    end do
    close(40)
    print *, 'Output histogram of sizes within circular aperture to: ', trim(adjustl(filename))

    deallocate(Sizes_in_Ap)

  end subroutine Size_Histogram_Circular_Aperture

end module Size_Histograms
