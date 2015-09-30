!--#Version that uses r200 defined as the radius on which density contained in halo is 200rho_crit - See Wright and Brainerd
module Mass_Profiles
  use Param_Types
  implicit none

  INTERFACE get_Lensing_Quantities
     module procedure get_Lensing_Quantities_SingleSource
  END INTERFACE get_Lensing_Quantities

  INTERFACE SMD_SIS
     module procedure SMD_SIS_Array, SMD_SIS_Scalar
  END INTERFACE SMD_SIS
  INTERFACE SMD_NFW
     module procedure SMD_NFW_Array, SMD_NFW_Scalar
  END INTERFACE SMD_NFW
  
  INTERFACE get_NFW_VirialRadius_from_VirialMass
     module procedure get_NFW_VirialRadius_from_VirialMass_Scalar
  END INTERFACE

  contains

    !--------------------------------GENERAL ROUTINES---------------------------------------------------------!
    real(double) function Magnification_Factor(Profile, Radius, Param, Redshift, Sigma_Critical)
      integer, intent(in):: Profile
      real(double), intent(in):: Param
      real(double), intent(in):: Radius, Redshift, Sigma_Critical

      real(double):: Gamma, Kappa

      select case(Profile)
      case(1)
         STOP 'I cannot calculate the Magnification Factor for this profile yet'
      case(2) !-SIS-!
         !--As SIS has Kappa \propto 1/theta, |gamma| = kappa, so only Kappa needs calculated
         Kappa = SMD_SIS_Scalar(Param, Radius)/Sigma_Critical
         Magnification_Factor = 1.e0_double/(1.e0_double-(2.e0_double*Kappa))
      case(3) !-NFW-!
         Kappa = SMD_NFW_Scalar(Radius, Redshift, Param)/Sigma_Critical
         Gamma = Differential_SMD_Scalar(Radius, Redshift, Param)/Sigma_Critical
         Magnification_Factor = 1.e0_double/( ((1.e0_double-Kappa)**2.e0_double) - Gamma*Gamma)
      end select

    end function Magnification_Factor
      
    real(double) function Total_MagnificationFactor_MultipleClusters(Profile, Position, Cluster_Pos, Param, Redshift, Source_Redshift, Sigma_Crit)
      !--Tests and works on single cluster case, and indications that it works on multiple cluster cases
      !-This could be edited to work with SIS by using the fact that for Kappa \propto 1/\theta, |Gamma| = Kappa, mu^-1 = 1-2Kappa
      use Cosmology, only: angular_diameter_distance_fromRedshift
      integer, intent(in):: Profile
      real(double), intent(in)::Param(:), Redshift(:)
      real(double),intent(in):: Cluster_Pos(:,:) !-Cluster, Position (RA,Dec)-!
      real(double), intent(in):: Position(2) !-RA, Dec-!
      real(double), intent(in),optional:: Source_Redshift, Sigma_Crit(:)

      real(double):: dRA, dDec, Theta, Radius,  D_l, D_ls, D_s
      real(double), dimension(size(Cluster_Pos,1)):: Sigma_Critical
      integer:: C, nCl
      real(double):: Gamma1, Gamma2, Kappa, Gamma
      real(double):: KappaT, Gamma1T, Gamma2T, GammaT
      
      real(double):: mag_Factor

      integer, save:: callcount = 0

      INTERFACE 
         real(double) function Total_MagnificationFactor_MultipleClusters(Profile, Position, Cluster_Pos, Param, Redshift, Source_Redshift, Sigma_Crit)
           use Param_Types
           integer, intent(in):: Profile
           real(double), intent(in)::Param(:), Redshift(:)
           real(double),intent(in):: Cluster_Pos(:,:) !-Cluster, Position (RA,Dec)-!                                                                                                                      
           real(double), intent(in):: Position(2) !-RA, Dec-!                                                                                                                                                
           
           real(double), intent(in),optional:: Source_Redshift, Sigma_Crit(:)
         end function Total_MagnificationFactor_MultipleClusters
      END INTERFACE

      if(present(Sigma_Crit)) then
         call get_Lensing_Quantities_SingleSource(Profile, Position, Cluster_Pos, Param, Redshift, Sigma_Crit = Sigma_Crit, magnification_Factor = mag_Factor)
      elseif(present(Source_Redshift)) then
         call get_Lensing_Quantities_SingleSource(Profile, Position, Cluster_Pos, Param, Redshift, Source_Redshift = Source_Redshift, magnification_Factor = mag_Factor)
      else
         STOP 'Total_MagnificationFactor_MultipleClusters - Sigma_Critical or Source_Redshift must be entered'
      end if

      Total_MagnificationFactor_MultipleClusters = mag_Factor

    end function Total_MagnificationFactor_MultipleClusters

    subroutine get_Lensing_Quantities_SingleSource(Profile, Position, Cluster_Pos, Param, Redshift, Source_Redshift, Sigma_Crit, reduced_Shear, shear, convergence, magnification_Factor, tangential_Shear)
      !----THIS SHOULD BE THE GO-TO ROUTINE FOR LENSING QUANTITIES.
      !--Returns all lensing quantities for a single source and multiple lens system.
      !--Author: cajd
      !--Date: 23Feb2015
      !--Tested 23Feb2015 and good to numerical error (amplitude of shear signal untested. Mag Factor (including shear amplitude) agrees with previous routines [which this has replaced], and they were well tested)

      use Cosmology, only: angular_diameter_distance_fromRedshift
      integer, intent(in):: Profile
      real(double), intent(in)::Param(:), Redshift(:)
      real(double),intent(in):: Cluster_Pos(:,:) !-Cluster, Position (RA,Dec)-!                                                                                                                                   
      real(double), intent(in):: Position(2) !-RA, Dec-!
                                                                                                                                                          
      real(double), intent(in),optional:: Source_Redshift, Sigma_Crit(:)
      real(double), intent(out), optional:: reduced_Shear(2), shear(2), convergence, magnification_Factor, tangential_Shear

      real(double):: dRA, dDec, Theta, Radius,  D_l, D_ls, D_s, DR
      real(double):: sin2Theta, cos2Theta
      real(double), dimension(size(Cluster_Pos,1)):: Sigma_Critical
      integer:: C, nCl
      real(double):: Gamma1, Gamma2, Kappa, Gamma
      real(double):: KappaT, Gamma1T, Gamma2T, GammaT, MGammaT !-T = Total Quantities

      integer, save:: callcount = 0

      if(Profile /=3) STOP 'get_Lensing_Quantities_SingleSource - I can only do this for NFW for now'

      if(size(Redshift) /= size(Cluster_Pos,1)) STOP 'Total_MagnificationFactor_MultipleClusters - Error in redshift - not correct size'
      if(size(Param) /= size(Cluster_Pos,1)) STOP 'Total_MagnificationFactor_MultipleClusters - Error in DM Profile Param - not correct size'

      KappaT = 0.e0_double; Gamma1T = 0.e0_double; Gamma2T = 0.e0_double; MGammaT = 0.e0_double; GammaT = 0.e0_double;
      nCl = size(Cluster_Pos,1)

      !--Get Sigma_Critical
      if(present(Sigma_Crit)) then
         if(size(Sigma_Crit) /= nCl) STOP 'Total_MagnificationFactor_MultipleClusters - Sigma_Crit entered not of the correct size'
         Sigma_Critical = Sigma_Crit
      else
         if(present(Source_Redshift) == .false.) STOP 'Total_MagnificationFactor_MultipleClusters - Sigma_Crit or Source Redshift Must be entered'
         Sigma_Critical = dsqrt(-1.e0_double)
         D_s = angular_diameter_distance_fromRedshift(0.e0_double, Source_Redshift)
         do C = 1, nCl
            if(Source_Redshift <= Redshift(C)) cycle
            if(Source_Redshift < 0.e0_double) STOP 'Total_MagnificationFactor_MultipleClusters - Source Redshift Entered invalid (negative)'
            D_l = angular_diameter_distance_fromRedshift(0.e0_double, Redshift(C))
            D_ls = angular_diameter_distance_fromRedshift(Redshift(C), Source_Redshift)
            Sigma_Critical(C) = 1.66492e18_double*(D_s/(D_l*D_ls))
         end do
      end if

      callcount = callcount + 1

      do C = 1, nCl
         !--Skip if Source Redshift is less than the lens redshift
         if(isNaN(Sigma_Critical(C)) .or. dabs(Sigma_Critical(C)) > huge(1.e0_double) .or. (Sigma_Critical(C) < 0.e0_double)) cycle

         dRA = (Position(1)-Cluster_Pos(C, 1))
         dDec = (Position(2) - Cluster_Pos(C,2))
         dR = dsqrt( dRA*dRA + dDec*dDec ) !-in Degrees-!

         D_l = angular_diameter_distance_fromRedshift(0.e0_double, Redshift(C))
         Radius = D_l*dR*(3.14159265359e0_double/180.e0_double) !-In Mpc/h-!

         !--Get angle wrt centre of cluster, on cartesian co-ordinate frame
         if(present(reduced_Shear) .or. present(shear) .or. present(magnification_Factor) .or. present(tangential_Shear)) then

            !-Get Gamma Tangential
            Gamma = Differential_SMD_Scalar(Radius, Redshift(C), Param(C))
            Gamma = Gamma/Sigma_Critical(C)

            !--Get angular transformation
            sin2Theta = 2.e0_double*((dDec*dRA)/(dR*dR))
            cos2Theta = 2.e0_double*((dRA*dRA)/(dR*dR)) - 1.e0_double
                        
            !--Convert to Shear Components on Cartesian Frame, which is universal to all lensing clusters
            Gamma1 = -1.e0_double*Gamma*cos2Theta
            Gamma2 = -1.e0_double*Gamma*sin2Theta

            if(isNaN(Gamma1)) print *, 'Gamma1 is a NaN:', Gamma, Redshift(C), Sigma_Critical(C), Param(C), Radius, D_l
            if(isNaN(Gamma2)) print *, 'Gamma2 is a NaN:', Gamma, Redshift(C), Sigma_Critical(C), Param(C), Radius, D_l
            
            !--Get Summed Quantities (sum over all clusters)
            Gamma1T = Gamma1T + Gamma1
            Gamma2T = Gamma2T + Gamma2

            GammaT = GammaT + Gamma
         end if

         if(present(reduced_Shear) .or. present(convergence) .or. present(magnification_Factor)) then
            !--Get Convergence Part--!
            Kappa = SMD_NFW_Scalar(Radius, Redshift(C), Param(C))
            Kappa = Kappa/Sigma_Critical(C)
            KappaT = KappaT + Kappa
         end if
      end do

      !---Set Outputs
      if(present(shear)) then
         shear = (/Gamma1T, Gamma2T/)
      end if

      if(present(reduced_Shear)) then
         reduced_Shear = (/Gamma1T, Gamma2T/)/(1.e0_double- KappaT)
      end if

      if(present(convergence)) then
         convergence = KappaT
      end if

      if(present(magnification_Factor)) then
         MGammaT = dsqrt(Gamma1T*Gamma1T + Gamma2T*Gamma2T)
         magnification_Factor = 1.e0_double/( ((1.e0_double-KappaT)**2.e0_double) - MGammaT*MGammaT)
      end if

      if(present(tangential_Shear)) then
         tangential_Shear = GammaT
      end if

    end subroutine get_Lensing_Quantities_SingleSource

!!$    function get_Cartesian_Shear_NFW(Profile, Position, Cluster_Pos, Param, Redshift, Source_Redshift, Sigma_Crit, get_Reduced_Shear)
!!$      !--Gets the tangential shear analytically, rotates to cartesian co-ordiantes.
!!$      !--Returns 2-element real(double) array containing gamma_1, gamma_2, (or g_1, g_2 if get_Reduced_Shear == True)
!!$      integer, intent(in):: Profile
!!$      real(double), intent(in)::Param(:), Redshift(:)
!!$      real(double),intent(in):: Cluster_Pos(:,:) !-Cluster, Position (RA,Dec)-!                                                                                                                                   
!!$      real(double), intent(in):: Position(2) !-RA, Dec-!                                                                                                                                                          
!!$      real(double), intent(in),optional:: Source_Redshift, Sigma_Crit(:)
!!$
!!$      logical:: optional:: get_Reduced_Shear
!!$
!!$
!!$      if(Profile /=3) STOP 'Total_MagnificationFactor_MultipleClusters - I can only do this for NFW for now'
!!$
!!$      if(size(Redshift) /= size(Cluster_Pos,1)) STOP 'Total_MagnificationFactor_MultipleClusters - Error in redshift - not correct size'
!!$      if(size(Param) /= size(Cluster_Pos,1)) STOP 'Total_MagnificationFactor_MultipleClusters - Error in DM Profile Param - not correct size'
!!$
!!$      KappaT = 0.e0_double; Gamma1T = 0.e0_double; Gamma2T = 0.e0_double
!!$      nCl = size(Cluster_Pos,1)
!!$
!!$
!!$      !--Get Sigma_Critical
!!$      if(present(Sigma_Crit)) then
!!$         if(size(Sigma_Crit) /= nCl) STOP 'Total_MagnificationFactor_MultipleClusters - Sigma_Crit entered not of the correct size'
!!$         Sigma_Critical = Sigma_Crit
!!$      else
!!$         if(present(Source_Redshift) == .false.) STOP 'Total_MagnificationFactor_MultipleClusters - Sigma_Crit or Source Redshift Must be entered'
!!$         Sigma_Critical = dsqrt(-1.e0_double)
!!$         D_s = angular_diameter_distance_fromRedshift(0.e0_double, Source_Redshift)
!!$         do C = 1, nCl
!!$            if(Source_Redshift <= Redshift(C)) cycle
!!$            if(Source_Redshift < 0.e0_double) STOP 'Total_MagnificationFactor_MultipleClusters - Source Redshift Entered invalid (negative)'
!!$            D_l = angular_diameter_distance_fromRedshift(0.e0_double, Redshift(C))
!!$            D_ls = angular_diameter_distance_fromRedshift(Redshift(C), Source_Redshift)
!!$            Sigma_Critical(C) = 1.66492e18_double*(D_s/(D_l*D_ls))
!!$         end do
!!$      end if
!!$
!!$      callcount = callcount + 1
!!$
!!$      do C = 1, nCl
!!$         !--Skip if Source Redshift is less than the lens redshift
!!$         if(isNaN(Sigma_Critical(C)) .or. dabs(Sigma_Critical(C)) > huge(1.e0_double) .or. (Sigma_Critical(C) < 0.e0_double)) cycle
!!$
!!$         !--Get angle wrt centre of cluster, on cartesian co-ordinate frame
!!$         dRA = dabs(Position(1)-Cluster_Pos(C, 1))
!!$         dDec = dabs(Position(2) - Cluster_Pos(C,2))
!!$         !--Theta needs edited to account for four quadrants (ok for just the magnification part)
!!$         Theta = atan(dDec/dRA)
!!$         
!!$         D_l = angular_diameter_distance_fromRedshift(0.e0_double, Redshift(C))
!!$         Radius = dsqrt( dRA*dRA + dDec*dDec ) !-in Degrees-!
!!$         Radius = D_l*Radius*(3.142e0_double/180.e0_double) !-In Mpc/h-!
!!$
!!$         !-Get Gamma Tangential
!!$         Gamma = Differential_SMD_Scalar(Radius, Redshift(C), Param(C))
!!$         Gamma = Gamma/Sigma_Critical(C)
!!$         
!!$         !--Convert to Shear Components on Cartesian Frame, which is universal to all lensing clusters
!!$         Gamma1 = -1.e0_double*Gamma*dcos(2.e0_double*Theta)
!!$         Gamma2 = -1.e0_double*Gamma*dsin(2.e0_double*Theta)
!!$
!!$         if(isNaN(Gamma1)) print *, 'Gamma1 is a NaN:', Gamma, dcos(2.e0_double*Theta), Theta, Redshift(C), Sigma_Critical(C), Kappa, Param(C), Radius, D_l
!!$         if(isNaN(Gamma2)) print *, 'Gamma2 is a NaN:', Gamma, dsin(2.e0_double*Theta), Theta, Redshift(C), Sigma_Critical(C), Kappa, Param(C), Radius, D_l
!!$
!!$         !--Get Summed Quantities
!!$         Gamma1T = Gamma1T + Gamma1
!!$         Gamma2T = Gamma2T + Gamma2
!!$
!!$         if(get_Reduced_Shear) then
!!$            Kappa = SMD_NFW_Scalar(Radius, Redshift(C), Param(C))
!!$            Kappa = Kappa/Sigma_Critical(C)
!!$            KappaT = KappaT + Kappa
!!$         end if
!!$      end do
!!$
!!$      if(get_Reduced_Shear) then
!!$         Gamma1T = Gamma1T/(1.e0_double- KappaT)
!!$         Gamma2T = Gamma2T/(1.e0_double- KappaT)
!!$      end if
!!$
!!$      get_Cartesian_Shear_NFW = (/Gamma1T, Gamma2T/)
!!$
!!$    end function get_Cartesian_Shear_NFW

    subroutine Halo_Mass(Profile, Param, Param_Error, Integrated_Mass, Mass_Error, Redshift)
      !--Gets Virial Radius for SMD type and returns Halo Mass--!
      real(double), intent(in):: Param, Param_Error(:)
      real(double), intent(out):: Integrated_Mass
      real(double),intent(out):: Mass_Error(:)
      real(double), intent(in),optional::Redshift
      integer::Profile
      
      real(double):: Virial_Radius

      Virial_Radius = virial_Radius_from_ProfileFreeParameter(Profile, Param)
      
      if(present(Redshift)) then
         call Integrated_Mass_Within_Radius(Profile, Virial_Radius, Param, Param_Error, Integrated_Mass, Mass_Error, Redshift)
      else
         call Integrated_Mass_Within_Radius(Profile, Virial_Radius, Param, Param_Error, Integrated_Mass, Mass_Error)
      end if

    end subroutine Halo_Mass

    real(double) function virial_Radius_from_ProfileFreeParameter(Profile, Param)
      integer, intent(in):: Profile
      real(double), intent(in):: Param

      select case(Profile)
      case(1) !-Flat-!
         STOP 'virial_Radius_from_ProfileFreeParameter - I cannot do FLAT yet'
      case(2) !--SIS--!
         STOP 'virial_Radius_from_ProfileFreeParameter - I cannot do SIS yet'
      case(3) !--NFW, r200 only--!
         virial_Radius_from_ProfileFreeParameter = Param
      case default
         STOP 'virial_Radius_from_ProfileFreeParameter - Incorrect profile entered'
      end select

    end function virial_Radius_from_ProfileFreeParameter

    subroutine Integrated_Mass_Within_Radius(Profile, Scale, Param, Param_Error, Integrated_Mass, Mass_Error, Redshift)
      !-Param is the free parameter for the profile: Flat->Sigma_0; SIS->Sigma_v^2; 3: NFW -> r_200, (c fixed)
      !--Scale must be in Mpc/h --!                                                                        
      !-Assumes Spherically Symmetric-!
      real(double), intent(in)::Scale, Param, Param_Error(:)
      real(double), intent(out):: Integrated_Mass
      real(double),intent(out):: Mass_Error(:)
      real(double), intent(in),optional::Redshift
      integer::Profile

      integer::i
      real(double)::iScale

      if(size(Param_Error) /= size(Mass_Error)) STOP 'Integrated_Mass_Within_Radius - Param Error and Mass Error MUST be the same size'

      if(Scale <= 0.0e0_double) then
         !--Use Virial Radius--!
         iScale = virial_Radius_from_ProfileFreeParameter(Profile, Param)
      else
         iScale = Scale
      end if

      select case(Profile)
      case(1) !-Flat-!
         Integrated_Mass = Flat_Mass_withinRadius(iScale, Param)
         do i = 1, size(Param_Error)
            Mass_Error(i) = Error_Flat_Mass_withinRadius(iScale, Param_Error(i))
         end do
      case(2) !-SIS-!
         Integrated_Mass = SIS_Mass_withinRadius(iScale, Param)
         do i = 1, size(Param_Error)
            Mass_Error(i) = SIS_Mass_withinRadius(iScale, Param_Error(i))
         end do
      case(3) !-NFW-!
         if(present(Redshift) == .false.) STOP 'Integrated_Mass_Within_Radius - Cluster Redshift MUST be present'
         Integrated_Mass = NFW_Mass_withinRadius(iScale, Redshift,r_200 = Param)
         do i = 1, size(Param_Error)
            Mass_Error(i) = Error_NFW_Mass_withinRadius(iScale, Param, Param_Error(i), Redshift)
         end do
      case default 
         STOP 'Integrated_Mass_Within_Radius - Incorrect profile entered'
      end select

    end subroutine Integrated_Mass_Within_Radius



    !------------------------------------NFW------------------------------------------------------------------!
    real(double) function Error_NFW_Mass_withinRadius(R, r200, Er200, z, r_s)
      use Cosmology, only: Normalised_Hubble_Parameter
      !--Returns the error on the integrated NFW mass within a radius--!
      !--Uses an analytic for for dDelcdr200 given by Wolfram Alpha, but could alos be done analytically--!
      !--Er200 is the measured error on r200-!
      !--TESTED with values 18Mar2014--!
      real(double), intent(in):: R, r200, Er200, z
      real(double), intent(in), optional::r_s

      real(double)::rho_c, M200, delta_c, c, rs

      real(double):: dDelcdr200, dMdr200, MR

      real(double):: Omega_Matter, rho_c0

      if(present(r_s)) then
         STOP 'Error_NFW_Mass_withinRadius - I cannot deal with rs yet'
      end if

      call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c)
      rs = r200/c

      !--Use normalise Hubble Parameter as 4.3017e-13 labels G/H_0^2
      rho_c = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(z)
      rho_c0 = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(0.0e0_double)
!      Omega_matter = 1.e0_double !-Use for halo defined wrt critical density
      Omega_matter = (rho_c0*0.3e0_double*((1+z)**3.e0_double)/rho_c) !-Use for halo defined wrt mass density

      !--Following is analytic form for the derivative, taken from Wolfram Alpha--!
      dDelcdr200 = Omega_Matter*(200.e0_double/3.e0_double)*((c*c)/rs)*( (3.e0_double*((1.e0_double+c)**2.e0_double)*dlog(1.e0_double+c) - c*(3.e0_double+4.e0_double*c))/((c-(1.e0_double+c)*dlog(1.e0_double+c))**2.e0_double) )

      MR = NFW_Mass_withinRadius(R, z, r_200 = r200)
      dMdr200 = (MR/delta_c)*dDelcdr200 !--Isolates the parts of the integrated mass which do not depend on delta_c, and therefore r200--!

      Error_NFW_Mass_withinRadius = dMdr200*Er200
           
  end function Error_NFW_Mass_withinRadius

    real(double) function NFW_Mass_withinRadius(R, z,r_200, M_200, r_s)
      !--Returns the mass contained within radius R--!
      !--TESTED with values 18Mar2014--!
      !--Can be used to obtain M200 and eM200 directly by entering R = r200--!
      real(double), intent(in):: R, z
      real(double), intent(in),optional:: M_200, r_200
      real(double), intent(in), optional::r_s

      real(double)::rho_c, M200, delta_c, c, rs, r200

      INTERFACE
         real(double) function NFW_Mass_withinRadius(R, z,r_200, M_200, r_s)
           use Param_Types
           real(double), intent(in):: R, z

           real(double), intent(in),optional:: M_200, r_200
           real(double), intent(in), optional::r_s
         end function NFW_Mass_withinRadius
      END INTERFACE

      if(present(r_200) .and. present(M_200)) then
         STOP 'NFW_Mass_withinRadius - ONLY ONE OF virial radius or virial mass must be entered, stopping...'
      elseif(present(r_200)) then
         r200 = r_200
      elseif(present(m_200)) then
         r200 = get_NFW_VirialRadius_from_VirialMass(z, m_200)
      else
         STOP 'NFW_Mass_withinRadius - Either virial radius or virial mass must be entered, stopping...'
      END if

      if(r200 == 0.e0_double) then
         NFW_Mass_withinRadius = 0.e0_double
         return
      end if

      if(present(r_s)) then
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c, r_s)
         rs = r_s
      else
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c)
         rs = r200/c
      end if
      !Check entered M200 against returned?!

      NFW_Mass_withinRadius = 4.e0_double*3.142e0_double*delta_c*rho_c*rs*rs*rs*(dlog(1.e0_double +(r/rs)) - ((r/rs)/(1.e0_double+(r/rs))) )


    end function NFW_Mass_withinRadius

    subroutine rescale_NFW_byMass(M_scale, r_scale, z, r200)
      !-M_Scale is the mass wanted within the scale r_scale
      !--Returns r200 so that the mass contained within a scale r_scale is equal to M_scale
      !--Assumes r200 is the only free parameter, and concentration is set by mass (fit) so that rs is set--!
      use Interpolaters, only:Linear_Interp
      real(double), intent(in)::M_scale, r_scale, z
      real(double),intent(out):: r200
      
      integer::i
      real(double), dimension(1000)::Mass_Storage, r200Grid
      real(double)::r200_l = 0.1e0_double, r200_u = 20.e0_double, dr200

      real(double)::rho_c, M200, delta_c, c, rs

      !--Loop over reasonable r200 values until the mass within the scale is as expected--!
      dr200 = (r200_u-r200_l)/size(Mass_Storage)
      do i = 1, size(Mass_Storage)
         r200Grid(i) = r200_l+(i-1)*dr200
         call get_NFW_Parameters(z, r200Grid(i), rho_c, M200, delta_c, c)
         rs = r200Grid(i)/c
         Mass_Storage(i) = NFW_Mass_withinRadius(R_Scale, Z, r_200 = r200Grid(i), r_s = rs)
      end do      

      if(M_Scale > maxval(Mass_Storage)) then
         print *, 'Mass_Scale, Max Tabulated:', M_Scale, maxval(Mass_Storage)
         STOP 'rescale_NFW_byMass - Mass Scale too large for tabulated values'
      end if
      if(M_Scale < minval(Mass_Storage)) STOP 'rescale_NFW_byMass - Mass Scale too small for tabulated values'

      r200 = Linear_Interp(M_Scale, Mass_Storage,r200Grid)

      !--Output--!
!!$      open(unit = 62, file = 'NFW_rescale_byr200.dat')
!!$      do i = 1, size(Mass_Storage)
!!$         write(62, *) r200Grid(i), Mass_Storage(i)
!!$      end do
!!$      close(62)

    end subroutine rescale_NFW_byMass

    function SMD_NFW_Array(r, z, r200, r_s)
      real(double),intent(in)::r(:), r200, z
      real(double), intent(in),optional::r_s
      
      real(double),dimension(size(r))::SMD_NFW_Array
      
      integer::i

      do i = 1, size(r)
         if(present(r_s)) then
            SMD_NFW_Array(i) = SMD_NFW_Scalar(r(i), z, r200, r_s)
         else
            SMD_NFW_Array(i) = SMD_NFW_Scalar(r(i), z, r200)
         end if
      end do

    end function SMD_NFW_Array

    real(double) function SMD_NFW_Scalar(r, z, r200, r_s)
      real(double),intent(in)::r, r200, z
      real(double), intent(in),optional::r_s

      real(double)::delta_c, c, rs, M200
      real(double)::omega_matter !-Evaluated at the redshift of the Halo-!
      real(double)::rho_c
      real(double)::x

      if(present(r_s)) then
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c,  r_s)
      else
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c)
      end if

      if(r200 == 0.e0_double) then
         SMD_NFW_Scalar = 0.e0_double
      end if

      rs = r200/c

      x = r/rs
      !--Following from Wright, Brainerd--!
      if(x < 1) then
         SMD_NFW_Scalar = ( (2.e0_double*rs*delta_c*rho_c)/(x*x - 1.e0_double) )*(1.e0_double - (2.e0_double/dsqrt(1.e0_double-x*x))*atanh(dsqrt((1.e0_double-x)/(1.e0_double+x))) )
      elseif(x == 1) then
         SMD_NFW_Scalar = (2.e0_double*rs*delta_c*rho_c)/3.e0_double
      else
         SMD_NFW_Scalar = ( (2.e0_double*rs*delta_c*rho_c)/(x*x - 1.e0_double) )*(1.e0_double - (2.e0_double/dsqrt(x*x-1.e0_double))*atan(dsqrt((x-1.e0_double)/(1.e0_double+x))) )
      end if



    end function SMD_NFW_Scalar

    function Differential_SMD_Array(r, z, r200, r_s)
      real(double),intent(in)::r(:), r200, z
      real(double), intent(in),optional::r_s

      real(double), dimension(size(r)):: Sigma, SigmaBar, Differential_SMD_Array

      if(present(r_s)) then
         Sigma = SMD_NFW_Array(r,z,r200,r_s)
         SigmaBar = Mean_SMD_NFW_Array(r, z, r200, r_s)
      else
         Sigma = SMD_NFW_Array(r,z,r200)
         SigmaBar = Mean_SMD_NFW_Array(r, z, r200)
      end if

      Differential_SMD_Array = SigmaBar-Sigma

    end function Differential_SMD_Array

    real(double) function Differential_SMD_Scalar(r, z, r200, r_s)
      !--Can be converted into a tangential shear by dividing by 1/Sigma_{Crit}
      real(double),intent(in)::r, r200, z
      real(double), intent(in),optional::r_s

      real(double):: Sigma, SigmaBar

      if(present(r_s)) then
         Sigma = SMD_NFW_Scalar(r,z,r200,r_s)
         SigmaBar = Mean_SMD_NFW_Scalar(r, z, r200, r_s)
      else
         Sigma = SMD_NFW_Scalar(r,z,r200)
         SigmaBar = Mean_SMD_NFW_Scalar(r, z, r200)
      end if

      Differential_SMD_Scalar = SigmaBar-Sigma

    end function Differential_SMD_Scalar

    real(double) function Shear_NFW(r,z,r200, Sigma_Critical, r_s)
      !--Eqn 14 Wright and Brainerd. Tested against Differential_SMD/Sigma_Critical and works
      real(double),intent(in)::r, r200, z, Sigma_Critical
      real(double), intent(in),optional::r_s

      real(double)::delta_c, c, rs, M200
      real(double)::omega_matter !-Evaluated at the redshift of the Halo-!
                                                                                                                                                                                                                  
      real(double)::rho_c, g
      real(double)::x

      if(present(r_s)) then
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c,  r_s)
      else
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c)
      end if

      rs = r200/c

      x = r/rs
      !--Following from Wright, Brainerd--!                                                                                                                                                                       
      if(x < 1) then
         g = ( (8.e0_double*atanh(dsqrt((1.e0_double-x)/(1.e0_double+x))))/((x*x)*dsqrt(1.e0_double-(x*x)))  ) + ((4.e0_double/(x*x))*dlog(0.5e0_double*x)) - (2.e0_double/((x*x)-1.e0_double)) + ( (4.e0_double*atanh(dsqrt((1.e0_double-x)/(1.e0_double+x))))/( ((x*x)-1.e0_double)*((1.e0_double-(x*x))**0.5e0_double)  )  )
         Shear_NFW = (rs*delta_c*rho_c*g)/Sigma_Critical
      elseif(x==1) then
         Shear_NFW = ((rs*delta_c*rho_c)/Sigma_Critical)*( (10.e0_double/3.e0_double) + 4.e0_double*dlog(0.5e0_double) )
      else
         g = ( (8.e0_double*atan(dsqrt((x-1.e0_double)/(1.e0_double+x))))/((x*x)*dsqrt((x*x)-1.e0_double))  ) + ((4.e0_double/(x*x))*dlog(0.5e0_double*x)) - (2.e0_double/((x*x)-1.e0_double)) + ((4.e0_double*atan(dsqrt((x-1.e0_double)/(1.e0_double+x))))/(((x*x)-1.e0_double)**1.5e0_double)   )
         Shear_NFW = (rs*delta_c*rho_c*g)/Sigma_Critical
      end if

    end function Shear_NFW

    function Mean_SMD_NFW_Array(r, z, r200, r_s)
      real(double),intent(in)::r(:), r200, z
      real(double), intent(in),optional::r_s
      
      real(double),dimension(size(r))::Mean_SMD_NFW_Array

      integer::i

      do i = 1, size(r)
         if(present(r_s)) then
            Mean_SMD_NFW_Array(i) = Mean_SMD_NFW_Scalar(r(i), z, r200, r_s)
         else
            Mean_SMD_NFW_Array(i) = Mean_SMD_NFW_Scalar(r(i), z, r200)
         end if
      end do

    end function Mean_SMD_NFW_Array

    real(double) function Mean_SMD_NFW_Scalar(r, z, r200, r_s)
      real(double),intent(in)::r, r200, z
      real(double), intent(in),optional::r_s
      
      real(double)::delta_c, c, rs, M200
      real(double)::omega_matter !-Evaluated at the redshift of the Halo-!                                                                                                                                         
      real(double)::rho_c
      real(double)::x

      if(present(r_s)) then
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c,  r_s)
      else
         call get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c)
      end if

      rs = r200/c

      x = r/rs
      !--Following from Wright, Brainerd--!
      if(x < 1) then
         Mean_SMD_NFW_Scalar = (4.e0_double/(x*x))*(rs*delta_c*rho_c)*( (2.e0_double/dsqrt(1.e0_double-(x*x)))*atanh(dsqrt((1.e0_double-x)/(1.e0_double+x))) + dlog(0.5e0_double*x) )
!         Mean_SMD_NFW_Scalar = (4.e0_double/(x*x))*rs*delta_c*rho_c*( (2.e0_double/dsqrt(1.e0_double-x*x))*atanh(dsqrt((1.e0_double-x)/(1.e0_double+x))) +dlog(0.5e0_double*x) )
      elseif(x == 1) then
         Mean_SMD_NFW_Scalar = 4.e0_double*rs*delta_c*rho_c*(1.e0_double + dlog(0.5e0_double))
!         Mean_SMD_NFW_Scalar =  4.e0_double*rs*delta_c*rho_c*(1.e0_double+dlog(0.5e0_double))
      else
         Mean_SMD_NFW_Scalar = (4.e0_double/(x*x))*(rs*delta_c*rho_c)*( (2.e0_double/dsqrt((x*x)-1.e0_double))*atan(dsqrt((x-1.e0_double)/(1.e0_double+x))) + dlog(0.5e0_double*x) )
!         Mean_SMD_NFW_Scalar = (4.e0_double/(x*x))*rs*delta_c*rho_c*( (2.e0_double/dsqrt(x*x-1.e0_double))*atan(dsqrt((x-1.e0_double)/(1.e0_double+x))) + dlog(0.5e0_double*x) )
      end if

    end function Mean_SMD_NFW_Scalar

    function get_NFW_VirialMass_from_VirialRadius(z, VirialRadius)
      use Cosmology, only: Normalised_Hubble_Parameter
      real(double), intent(in)::VirialRadius(:), z
      real(double),dimension(size(VirialRadius))::get_NFW_VirialMass_from_VirialRadius

      real(double)::rho_c, Omega_matter, rho_c0

      rho_c = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(z)
      rho_c0 = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(0.0e0_double)
      Omega_matter = (rho_c0*0.3e0_double*((1+z)**3.e0_double)/rho_c)
      !--Use Omega_Matter = 1 to use definition of r200 as where density of matter contained in halo is ** 200 times the critical density ** (e.g. Wright, Brainerd 1999) In this case, the fit of Dolag et al may not be correct
      !--Use Omega_Matter = 1 to use definition of r200 as where density of matter contained in halo is ** 200 times the mean matter density ** (e.g. Dolag 2004, Heymans 2008)
      !Omega_matter = 1.e0_double


      get_NFW_VirialMass_from_VirialRadius = ( (800e0_double*3.142e0_double)/(3.e0_double) )*(VirialRadius**3.e0_double) * rho_c * Omega_Matter

    end function get_NFW_VirialMass_from_VirialRadius

    function get_NFW_VirialRadius_from_VirialMass_Scalar(z,VirialMass)
      use Cosmology, only: Normalised_Hubble_Parameter
      real(double), intent(in)::VirialMass, z
      real(double)::get_NFW_VirialRadius_from_VirialMass_Scalar

      real(double)::rho_c, Omega_matter, rho_c0

      rho_c = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(z)
      rho_c0 = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(0.0e0_double)

      !--Use Omega_Matter = 1 to use definition of r200 as where density of matter contained in halo is ** 200 times the critical density ** (e.g. Wright, Brainerd 1999)
      !--Use Omega_Matter \= 1 to use definition of r200 as where density of matter contained in halo is ** 200 times the mean matter density ** (e.g. Dolag 2004, Heymans 2008)
      !Omega_matter = 1.e0_double
      Omega_matter = (rho_c0*0.3e0_double*((1+z)**3.e0_double)/rho_c)

      get_NFW_VirialRadius_from_VirialMass_Scalar = ( (VirialMass*3.e0_double)/(800e0_double*3.142e0_double*rho_c*Omega_Matter) )**(1.e0_double/3.e0_double)

    end function get_NFW_VirialRadius_from_VirialMass_Scalar


    subroutine get_NFW_Parameters(z, r200, rho_c, M200, delta_c, c, r_s)
      !--TESTED with values 18Mar2014--!
      use Cosmology, only: Normalised_Hubble_Parameter
      
      real(double),intent(in)::r200, z
      real(double), intent(in),optional::r_s
      real(double),intent(out)::delta_c, c, M200, rho_c

      real(double)::omega_matter,rs, rho_c0 !-Evaluated at the redshift of the Halo-!                                                                                                                                        

      !--Get M200 for entered r200--!
      rho_c = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(z) !-in Units [M_Sun/h][h/Mpc]^3. Uses G/H_0 = 4.3017e-13 (Mpc/h)^3(h/M_sun)-!
      rho_c0 = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(0.0e0_double)

      !--Use Omega_Matter = 1 to use definition of r200 as where density of matter contained in halo is ** 200 times the critical density ** (e.g. Wright, Brainerd 1999)
      !--Use Omega_Matter = 1 to use definition of r200 as where density of matter contained in halo is ** 200 times the mean matter density ** (e.g. Dolag 2004, Heymans 2008)
      !Omega_matter = 1.e0_double
      Omega_matter = (rho_c0*0.3e0_double*((1+z)**3.e0_double)/rho_c)

      M200 = ( (800e0_double*3.142e0_double)/(3.e0_double) )*(r200**3.e0_double) * rho_c * Omega_Matter

      !--Get Concentration/rs
      if(present(r_s)) then
         rs = r_s      
         c = r200/rs
      else
         c = NFW_get_concentration_fit(z, M200)
         rs = r200/c
      end if

      delta_c =  ((200.e0_double*Omega_Matter*c*c*c)/3.e0_double)/(dlog(1.e0_double+c) - c/(1.e0_double+c))
      
    end subroutine get_NFW_Parameters

    real(double) function NFW_get_Concentration_fit(z, M)
      !--Uses the fit of Dolag et al 2004 (see CH08) to get the concentration from the halo virial mass (M, usually M200) and the redshift of the halo--! 
      !--M must be in M_Sun/h--!
      real(double), intent(in)::z, M

      !--Fit Parameters--!
      real(double), parameter:: alpha = -0.102e0_double, c_0 = 9.59e0_double

      NFW_get_Concentration_fit = (c_0/(1.e0_double+z))*((M/10.e14_double)**alpha)

    end function NFW_get_Concentration_fit

    !--------------------------_End NFW Routines_-------------------------------------------------------------!

    !--------------------------------- Singular Isothermal Sphere -----------------------------------------------!
    real(double) function SIS_Mass_withinRadius(R, Sigv2)
      !--Free Parameter here is velocity disperion squared--! 
      real(double), intent(in)::R, Sigv2

      SIS_Mass_withinRadius = (6.2832e0_double*Sigv2*R)/(2.e0_double*4.3017e-9_double) !-(2*Pi*Sigv^2*R)/(2G)-!

    end function SIS_Mass_withinRadius

    real(double) function Error_SIS_Mass_withinRadius(R, eSigv2)
      !--Free Parameter here is velocity dispersion squared--!
      real(double), intent(in)::R, eSigv2

      Error_SIS_Mass_withinRadius = (6.2832e0_double*eSigv2*R)/(2.e0_double*4.3017e-9_double)

    end function Error_SIS_Mass_withinRadius

    real(double) function VelocityDispersionSquared_SIS_fromHaloMass(M200, Lens_Redshift)
      !@--Sets sigma_v^2, which sets the scale of the profile, from the halo mass M200
      use Cosmology, only: Normalised_Hubble_Parameter
      real(double), intent(in):: M200, Lens_Redshift

      real(double):: rho_c, rho_c0, Omega_Matter
      real(double):: G = 4.302e-9_double !-in units [Mpc/M_{sun}][km/s]^2
      
      !--Use normalise Hubble Parameter as 4.3017e-13 labels G/H_0^2
      !--rho_c then in units [M_{Sun}/h][h/Mpc]^3
      rho_c = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(Lens_Redshift)
      rho_c0 = (3.e0_double/(8.e0_double*3.142e0_double*4.3017e-13_double))*Normalised_Hubble_Parameter(0.0e0_double)
      Omega_matter = 1.e0_double !-Use for halo defined wrt critical density
      !Omega_matter = (rho_c0*0.3e0_double*((1+z)**3.e0_double)/rho_c) !-Use for halo defined wrt mass density
      
      VelocityDispersionSquared_SIS_fromHaloMass = (( ((4.e0_double*200.e0_double)/3.e0_double)*rho_c*Omega_Matter*((M200/3.14159265359e0_double)**2.e0_double))**0.3333e0_double)*G

    end function VelocityDispersionSquared_SIS_fromHaloMass

    real(double) function get_Einstein_Radius_SIS(velocity_dispersion2, Lens_Redshift, Source_Redshift)
      use cosmology, only: angular_diameter_distance_fromRedshift
      !--Utilises two methods: 1 searches for point where kappa = 0.5, the other uses an analytic form for the Einstein Radius
      real(double), intent(in):: velocity_dispersion2, Lens_Redshift, Source_Redshift

      integer:: Method = 2

      !--Method 1 declarations--!
      real(double):: D_ls, D_s

      !--Method 2 Declarations--!

      select case(Method)
      case(1)
         !--Search for kappa = 0.5
         STOP 'get_Einstein_Radius_SIS - Search Method not coded yet'
      case(2)
         !-- \theta_E = \frac{4\pi D_{ds}}{D_s}\left(\frac{\sigma_v}{c}\right)^2
         D_s = angular_diameter_distance_fromRedshift(0.e0_double, Source_Redshift)
         D_ls = angular_diameter_distance_fromRedshift(Lens_Redshift, Source_Redshift)

         !-- Velocity Dispersion should be km/s
         get_Einstein_Radius_SIS = (4.e0_double*3.141592e0_double)*(D_ls/D_s)*(velocity_dispersion2/((2.99792e5_double)**2.e0_double))
      case default
         STOP 'get_Einstein_Radius_SIS - Invalid method entered'
      end select

    end function get_Einstein_Radius_SIS

    function SMD_SIS_Array(velocity_dispersion2, radius, core, truncation)
      real(double), intent(in):: velocity_dispersion2, radius(:)
      real(double), intent(in),optional:: core, truncation
      
      real(double), dimension(size(radius)):: SMD_SIS_Array

      real(double)::Core_Radius, Truncation_Radius

      integer::i

      Core_Radius = 0.e0_double; Truncation_Radius = 1.e50_double
      if(present(Core)) Core_Radius = Core
      if(present(Truncation)) Truncation_Radius = Truncation

      do i = 1, size(radius)
         SMD_SIS_Array(i) = SMD_SIS_Scalar(velocity_dispersion2, radius(i), Core_Radius, Truncation_Radius)
      end do

    end function SMD_SIS_Array


    real(double) function SMD_SIS_Scalar(velocity_dispersion2, radius, core, truncation)
      !--Velocity Dispersion^2 should be (km/s)^2
      real(double), intent(in):: velocity_dispersion2, radius
      real(double), intent(in),optional:: core, truncation

      real(double)::Core_Radius, Truncation_Radius

      Core_Radius = 0.e0_double; Truncation_Radius = 1.e50_double
      if(present(Core)) Core_Radius = Core
      if(present(Truncation)) Truncation_Radius = Truncation

      !--Newtons Constant here is taken to be 4.3017e-9 [M_Sun/Mpc][km/s]^2
      SMD_SIS_Scalar =  ( (Velocity_Dispersion2)/(2.e0_double*4.3017e-9_double) )*( 1.e0_double/(dsqrt(Radius*Radius + Core_Radius*Core_Radius)) - 1.e0_double/(dsqrt(Radius*Radius + Truncation_Radius*Truncation_Radius)) )


    end function SMD_SIS_Scalar

    real(double) function SIS_velocity_dispersion_byMass(Mass, Mass_Scale)
      !--Does not work if there is a core radius-!
      !--Mass_Scale needs to be in Mpc/h--!
      !--Returns Km/s--!
      real(double), intent(in):: Mass, Mass_Scale

      SIS_velocity_dispersion_byMass = dsqrt((4.3017e-9_double*Mass)/(3.142e0_double*Mass_Scale))

    end function SIS_velocity_dispersion_byMass

    !--------Flat Profile---------------------------!

    real(double) function Flat_Mass_withinRadius(R, Sigma_0)
      !--Scale must be in Mpc/h --!
      !-Assumes Spherically Symmetric-!
      real(double), intent(in)::R, Sigma_0

      Flat_Mass_withinRadius = 3.142e0_double*R*R*Sigma_0

    end function Flat_Mass_withinRadius


    real(double) function Error_Flat_Mass_withinRadius(R, eSigma_0)
      !--Scale must be in Mpc/h --!
      !-Assumes Spherically Symmetric-!
      real(double), intent(in)::R, eSigma_0

      Error_Flat_Mass_withinRadius = 3.142e0_double*R*R*eSigma_0

    end function Error_Flat_Mass_withinRadius


end module Mass_Profiles
