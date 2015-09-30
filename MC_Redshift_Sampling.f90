module MC_Redshift_Sampling
  use Catalogues; use Distributions
  implicit none

  contains

    subroutine Monte_Carlo_Redshift_Sampling_SigmaCritical(Cat, Lens_Redshift, Sigma_Crit)
      !--Returns Sigma Critical for each galaxy in the Catalogue, where Sigma_Critical is returned using MC technique of CH08 is there exisits not redshift information for the galaxy--!
      !¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬!
      !--SIGMA CRITICAL RETURNED IN UNITS OF 10^18 MSun/h--!
      !¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬!
      use Cosmology, only:angular_diameter_distance_fromRedshift
      type(Catalogue),intent(in)::Cat
      real(double), intent(in)::Lens_Redshift
      real(double), intent(out),allocatable:: Sigma_Crit(:)

      real(double)::D_l, D_s, D_ls

      integer,parameter::N_mc = 100

      integer::N_mcd, N_Fail_Magnitude

      integer::c, i, mc
      real(double), allocatable::PDF(:), ZGrid(:)
      real(double)::Z_Lower = 0.e0_double, Z_Upper = 4.e0_double, dZ
      integer::nz = 1000
      real(double)::z_med

      !-Temporary storage of MC'd redshifts and varaible values-!
      real(double),dimension(N_mc)::Redshift_Sample, Variable_To_Calculate 
      real(double):: Calculated_Variable

      if(allocated(Sigma_Crit)) deallocate(Sigma_Crit)
      allocate(Sigma_Crit(size(Cat%RA))); Sigma_Crit = -1.e0_double

      D_l = angular_diameter_distance_fromRedshift(0.e0_double, Lens_Redshift)

      !--Set up redshift grid which will be used to get PDF for each galaxy--!                                                                                                                                                             
      allocate(ZGrid(nZ)); ZGrid = 0.e0_double
      dZ = (Z_Upper-Z_Lower)/(1.e0_double*(nZ-1))
      do i = 1, nZ
         ZGrid(i) = Z_Lower + (i-1)*dZ
      end do

      print *, 'Assigning Sigma_Critical information by redshift information using MC sampling to :', count(Cat%Redshift < 0.e0_double), ' galaxies. This is expected to take:', 2.85e-4_double*count(Cat%Redshift < 0.e0_double)*(N_mc/100.e0_double), ' minutes......'

      N_mcd = 0; n_Fail_Magnitude = 0
      do c = 1, size(Cat%RA)
         if(Cat%Redshift(c) >= 0.e0_double) then
            !--Pre-assigned redshift for that galaxy should be used--!
            !*!
            D_s = angular_diameter_distance_fromRedshift(0.e0_double, Cat%Redshift(c))
            D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Cat%Redshift(c))
            
            Calculated_Variable = 1.66492e0_double*(D_s/(D_l*D_ls))

         else !-MC Sampling-!
            N_mcd = N_mcd + 1
            !--Get a smail et al PDF                                                            
            z_med = 0.29e0_double*(Cat%MF606W(c) - 22.e0_double) + 0.31e0_double
            
            
            if(z_med <= 0.e0_double) then
               z_med = 0.02e0_double !--These have m<~21, which *SHOULD* be cut from the sample anyway--!!!                                                                                                                                      
               n_Fail_Magnitude = n_Fail_Magnitude + 1
            end if
            call Analytic_Source_Redshift_Distribution(2.e0_double, 1.5e0_double, z_med, ZGrid, PDF)
            
            !--Assign Redshift by Randomly sampling from that distribution--!
            call randomly_Sample_From_Distribution(ZGrid, PDF, Redshift_Sample)
            
            do mc = 1, N_mc
               !--Calculate Variable for each redshift sample--!
               if(Redshift_Sample(mc) < 0.e0_double) STOP 'Assign Redshifts - CH08 - Sampled Redshift Negative, fatal'
               !*!
               D_s = angular_diameter_distance_fromRedshift(0.e0_double, Redshift_Sample(mc))
               D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Redshift_Sample(mc))

               Variable_To_Calculate(mc) = 1.66492e0_double*(D_s/(D_l*D_ls))
            end do
            
            
            Calculated_Variable = sum(Variable_To_Calculate)/N_mc
            Variable_To_Calculate = 0.e0_double; Redshift_Sample = -1.e0_double            
         end if
         
         Sigma_Crit(c) = Calculated_Variable
         
         Calculated_Variable = 0.e0_double
      end do

      print *, 'Of ', size(Cat%RA), 'galaxies, ', N_mcd, " where Monte Carlo'd"
      if(n_Fail_Magnitude > 0) print *, 'Of these, ', n_Fail_Magnitude, ' where below an acceptable magnitude limit for the catalogue'

    end subroutine Monte_Carlo_Redshift_Sampling_SigmaCritical


    subroutine Monte_Carlo_Redshift_Sampling_Catalogue(Cat, Variable)
      use Distributions
      !-Populates *Physical Size* for a sample of galaxies without redshift information (tested by looking for z<0) in the Catalogue
      !--Uses the MC method described in CH08 using:
      !--------------------------------------------: Smail et al nz for magnitude, with a z_m(mag) given by z_m(m) = 0.29(m-22)+0.31
      !--TO DO:
      !-------:
      !-- Lines which would need to be changed dependingon the variable to be calculated are marked with !*!

      type(Catalogue)::Cat
      integer, intent(in),optional:: Variable !-1: Physical Size-!


      integer::iVariable
      integer,parameter::N_mc = 100

      integer::N_mcd, N_Fail_Magnitude

      integer::c, i, mc
      real(double), allocatable::PDF(:), ZGrid(:)
      real(double)::Z_Lower = 0.e0_double, Z_Upper = 4.e0_double, dZ
      integer::nz = 1000
      real(double)::z_med

      !-Temporary storage of MC'd redshifts and varaible values-!
      real(double),dimension(N_mc)::Redshift_Sample, Variable_To_Calculate 
      real(double):: Calculated_Variable

      INTERFACE
         subroutine Monte_Carlo_Redshift_Sampling_Catalogue(Cat, Variable)
           use Param_Types; use Catalogues
           type(Catalogue)::Cat
           
           integer, intent(in),optional:: Variable !-1: Physical Size-!
         end subroutine Monte_Carlo_Redshift_Sampling_Catalogue
      END INTERFACE

      !--Set up redshift grid which will be used to get PDF for each galaxy--!                                                                                                                                                             
      allocate(ZGrid(nZ)); ZGrid = 0.e0_double
      dZ = (Z_Upper-Z_Lower)/(1.e0_double*(nZ-1))
      do i = 1, nZ
         ZGrid(i) = Z_Lower + (i-1)*dZ
      end do

      iVariable = 1 !-Physical Sizes as default-!
      if(present(Variable)) iVariable = Variable

      print *, 'Assigning variable information by redshift information using MC sampling to :', count(Cat%Redshift < 0.e0_double), ' galaxies. This is expected to take:', 2.85e-4_double*count(Cat%Redshift < 0.e0_double), ' minutes......'

      N_mcd = 0; n_Fail_Magnitude = 0
      do c = 1, size(Cat%RA)
         if(Cat%Redshift(c) >= 0.e0_double) then
            !--Pre-assigned redshift for that galaxy should be used--!
            !*!
            Calculated_Variable =  Physical_Size_from_Pixel_Size(Cat%Redshift(c), Cat%Sizes(c))
         else !-MC Sampling-!
            N_mcd = N_mcd + 1
            !--Get a smail et al PDF                                                            
            z_med = 0.29e0_double*(Cat%MF606W(c) - 22.e0_double) + 0.31e0_double
            
            
            if(z_med <= 0.e0_double) then
               z_med = 0.02e0_double !--These have m<~21, which *SHOULD* be cut from the sample anyway--!!!                                                                                                                                      
               n_Fail_Magnitude = n_Fail_Magnitude + 1
            end if
            call Analytic_Source_Redshift_Distribution(2.e0_double, 1.5e0_double, z_med, ZGrid, PDF)
            
            !--Assign Redshift by Randomly sampling from that distribution--!
            call randomly_Sample_From_Distribution(ZGrid, PDF, Redshift_Sample)
            
            do mc = 1, N_mc
               if(Redshift_Sample(mc) < 0.e0_double) STOP 'Assign Redshifts - CH08 - Sampled Redshift Negative, fatal'
               !*!
               Variable_To_Calculate(mc) = Physical_Size_from_Pixel_Size(Redshift_Sample(mc), Cat%Sizes(c))
            end do
            
            
            Calculated_Variable = sum(Variable_To_Calculate)/N_mc
            Variable_To_Calculate = 0.e0_double; Redshift_Sample = -1.e0_double            
         end if
         
         select case(iVariable)
         case(1)!--Physical Sizes-!
            Cat%Physical_Sizes(c) = Calculated_Variable
         case default 
            STOP 'Monte_Carlo_Redshift_Sampling_Catalogue - Variable to be returned not supported, stopping...'
         END select
         
         Calculated_Variable = 0.e0_double
      end do

      print *, 'Of ', size(Cat%RA), 'galaxies, ', N_mcd, " where Monte Carlo'd"
      if(n_Fail_Magnitude > 0) print *, 'Of these, ', n_Fail_Magnitude, ' where below an acceptable magnitude limit for the catalogue'

    end subroutine Monte_Carlo_Redshift_Sampling_Catalogue



end module MC_Redshift_Sampling
