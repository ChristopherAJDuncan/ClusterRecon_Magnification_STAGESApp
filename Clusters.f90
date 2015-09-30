module Foreground_Clusters
  use Param_Types; use Cosmology, only: angular_diameter_distance_fromRedshift
  implicit none

!"#################################################################################!
 !- Contains the code used to position the clusters, as well as assign them an aperture radius etc -!

  !____________________________________________________________________________________________________!
  !-Derived type cluster constains all the information on the clusters/apertures.
  !-Used both to define cluster (simulations) and define position of aperture on field (during fitting)
  !-Fitting group groups clusters together. Each in a group will be fit simultaneously
  !____________________________________________________________________________________________________!
  type Foreground
     real(double),allocatable::Position(:,:) !-2nd Dimension is RA(1), Dec(2)--!
     real(double),allocatable::DM_Profile_Parameter(:)
     real(double),allocatable::Redshift(:)
     integer,allocatable:: Fitting_Group(:)
     integer:: Surface_Mass_Profile
  end type Foreground


  contains

    !--Getter/Setters--!
    subroutine Get_Clusters(Cl, Input_File)
      !--Constructor - TYPICALLY NOT A GETTER--!
      !--Reads in Cluster/Aperture information from a specified input file
      !--Stores Cluster information--!                                                                                                                                     
      type(Foreground),intent(out)::Cl
      character(*), intent(in)::Input_File
      
      integer::c, i
      logical::Here
      
      integer::nClusters, SMD_Profile_Type
      real(double),allocatable:: Positions(:)
      real(double),allocatable:: DM_Parameters(:)
      integer,allocatable:: Cluster_Fitting_Group(:)
      real(double), dimension(100)::Redshift
      
      namelist/Number/nClusters, SMD_Profile_Type
      namelist/Cluster_Variables/Positions, DM_Parameters, Redshift, Cluster_Fitting_Group

      call Foreground_Destruct(Cl)
      
      inquire(file = Input_File, exist = Here)
      if(here == .false.) then 
         print *, 'Input_File:', Input_File
         STOP 'Cluster Input File not present'
      end if
      open(unit = 87, file = Input_File)
      
      read(87, nml = Number)
      
      allocate(Positions(2*nClusters));
      allocate(DM_Parameters(nClusters))
      allocate(Cluster_Fitting_Group(nClusters)); Cluster_Fitting_Group = (/ (i,i=1,nClusters) /)
      Redshift = -1.

      read(87, nml = Cluster_Variables)

      close(87)

      !--Unpack--!                                                                      
                                                                                                                             
      allocate(Cl%Position(nClusters,2))
      allocate(Cl%DM_Profile_Parameter(nClusters))
      allocate(Cl%Fitting_Group(nClusters)); 
      do i = 1, nClusters
         Cl%Position(i,:) = (/ Positions(2*i-1),Positions(2*i)/)
         Cl%DM_Profile_Parameter(i) = DM_Parameters(i)
         Cl%Fitting_Group(i) = Cluster_Fitting_Group(i)
      end do

      allocate(Cl%Redshift(nClusters)); Cl%Redshift = -1.      
      if(count(Redshift>=0) == 1) then
         Cl%Redshift = Redshift(1)
      else
         Cl%Redshift = Redshift
      end if

      Cl%Surface_Mass_Profile = SMD_Profile_Type

      deallocate(Positions, DM_Parameters, Cluster_Fitting_Group)
    end subroutine Get_Clusters

    subroutine Foreground_Destruct(Cl)
      type(Foreground)::Cl

      Cl%Surface_Mass_Profile = 0
      if(allocated(Cl%Position)) deallocate(Cl%Position)
      if(allocated(Cl%DM_Profile_Parameter)) deallocate(Cl%DM_Profile_Parameter)
      if(allocated(Cl%Redshift)) deallocate(Cl%Redshift)
    
    end subroutine Foreground_Destruct

    !------------------------End Of Getter/Setters-----------------------------------------------------------------------------!
    function get_distance_between_Clusters(Positions, Redshift)
      real(double), intent(in):: Positions(:,:)
      real(double),intent(in),optional::Redshift
      real(double), dimension(size(Positions,1), size(Positions,1)):: get_distance_between_Clusters
      
      integer:: i,j

      get_distance_between_Clusters = 0.e0_double
      do i = 1, size(Positions,1)
         do j = i+1, size(Positions,1)
            if(i == j) cycle
            get_distance_between_Clusters(i,j) = dsqrt( (Positions(i,1)-Positions(j,1))**2.e0_double + (Positions(i,2)-Positions(j,2))**2.e0_double )
            get_distance_between_Clusters(j,i) = get_distance_between_Clusters(i,j)
         end do
      end do

      print *, 'Returned distance between clusters:', get_distance_between_Clusters

    end function get_distance_between_Clusters

    subroutine print_distance_between_Clusters(Positions, Redshift)
      real(double), intent(in):: Positions(:,:)
      real(double),intent(in)::Redshift

      integer::i, j 
      real(double)::D_l, angular_distance

      D_L = angular_diameter_distance_fromRedshift(0.e0_double, Redshift)
      do i= 1, size(Positions,1)
         do j = i+1, size(Positions,1)
            angular_distance = dsqrt( (Positions(i,1)-Positions(j,1))**2.e0_double + (Positions(i,2)-Positions(j,2))**2.e0_double )
            write(*,'(A,I1,A,I1,A,e12.5)') 'Distance Between Clusters ',i,' and ',j,' (on Sky) is: Angular:', angular_distance
            write(*,'(A,e12.5)') '                                                 and Physical:', D_l*(angular_distance)*(3.142e0_double/180.e0_double)
         end do
      end do

    end subroutine print_distance_between_Clusters

  end module Foreground_Clusters
