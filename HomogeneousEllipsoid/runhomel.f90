program runhomel
    use intreal_types
    use params
    use cosmoparams
    use input_parameters
    use HomogeneousEllipsoid
    implicit none
    character(len=:), allocatable :: filename
    integer len_filename, idynax
    real Omcurv, zdynax, Ddynax, fcdynax

    ! Open logfile
    read(*,*) len_filename
    allocate( character(len=len_filename) :: filename )
    read(*,*) filename
    open(unit=10, file=filename)

    ! Reads in homel parameters from command line
    call readhomel(Omcurv) ! subroutine below

    ! Make tables for linear-order approximation to Zeldovich D(t)
    !call Dlinear_cosmology(Omx,Omb,Omvac,h,iamcurved,dcurv) ! subroutine in src/hpkvd/psubs)_Dlinear.f90
    !call Dlinear_setup ! subroutine in src/hpkvd/psubs)_Dlinear.f90

    idynax=1
    !call evolve_ellipse_full(idynax,zdynax,Ddynax,fcdynax)

    !write(*,*) Frho, e_v, p_v
    close(10)
end program runhomel

subroutine readhomel(Omcurv)
    use intreal_types
    use params
    use cosmoparams
    use input_parameters
    use HomogeneousEllipsoid
    implicit none
    real Omcurv

    ! Declared in params.f90
    read(*,*) iwant_evmap! integer, 1) turn z_vir vs e_v, 2) z vs p_v for fixed e_v, 3) table
    read(*,*) nstepmax   ! integer, stops neverending looping, using same value as hpkvd
    ! nyp:integer assigned 6 in evolve_ellipse_full
    read(*,*) iwant_rd   ! integer, set to 1 to use Carlson's elliptic integrals
    read(*,*) tfac   ! real, fraction of local 1-axis Hubble time for dt
    read(*,*) zinit  ! real, initial redshift of homel simulation
    read(*,*) dcrit  ! real \
    read(*,*) e_vmax ! real  |
    read(*,*) de_v   ! real  |
    read(*,*) p_vbar ! real  } all 0.0 unless iwant_evmap != 3
    read(*,*) dp_v   ! real  |
    read(*,*) e_vbar ! real  |
    read(*,*) Fbar   ! real /
    read(*,*) fcoll_3! real \
    read(*,*) fcoll_2! real  } radial freeze out factors
    read(*,*) fcoll_1! real /
    ! a_3eq,a_2eq,a_1eq:real, assigned in evolve_ellipse_full
    read(*,*) ivir_strat  ! integer, 1) a_jeq=fcoll_3 a_b3, 2) a_jeq=fcoll_j a_bj
    read(*,*) iforce_strat! integer, 0) no bg, 1) sbg, 3) bg+NLstrain, 4) stbg+Lstrain, 5) Lstrain, 6) SW b_i

    ! Declared in cosmoparams.f90
    read(*,*) Omb  ! real, Omega_baryon, the baryonic matter density fraction
    read(*,*) Omx  ! real, Omega_x, the DM density fraction
    ! Omt real, Omega_total, the total density fraction (excluding curvature), assigned Omt=Omnr+Omvac in Dlinear_cosmology
    ! Omer real, assigned 0.0 in Dlinear_cosmology
    ! Omnr real, Omega_nr, the non-relativistic energy density fraction, assigned Omnr=Omb+Omx+Omhdm in Dlinear_cosmology
    ! Omhdm real, Omega_hotDM, the hot DM energy density fraction, assigned in Dlinear_cosmology
    read(*,*) Omvac! real, Omega_Lambda, DE density fraction
    read(*,*) Omcur! real, Omega_k, the curvature parameter
    Omcurv = Omcur
    read(*,*) h    ! real, little h, the dimensionless Hubble constant, h = H_0/(100 km s^-1 Mpc^-1)
    ! Rvac_nr:real, assigned Omvac/Omnr in homel
    ! Rcurv_nr:real, assigned Omcurv/Omnr in homel
    ! iamcurved:integer, assigned 1 or 0 in Dlinear_cosmology
    ! dcurv!:real, assigend in Dlinear_cosmology

    ! Declared in input_parameters.f90
    read(*,*) ihard ! integer, used in formatting output
    ! fhdmclus:real, assigned in psubs_Dlinear

    read(*,*) Frho,e_v,p_v

end subroutine readhomel
