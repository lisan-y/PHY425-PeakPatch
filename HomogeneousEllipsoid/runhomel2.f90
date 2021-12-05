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
    !call readhomel(Omcurv) ! subroutine below
    open(unit=6, file='home_params.txt')
    read(6,*) iwant_evmap, iwant_rd, ihard, ivir_strat, iforce_strat,&
    omx, omb, h, omvac, omnr, omcurv, omt, zinit, tfac, fcoll_3, fcoll_2,&
    fcoll_1, idynax, nstepmax, nout, dcrit, Frho, e_v, p_v, e_vmax
    close(6)
    ! Make tables for linear-order approximation to Zeldovich D(t)
    call Dlinear_cosmology(Omx,Omb,Omvac,h,iamcurved,dcurv) ! subroutine in src/hpkvd/psubs)_Dlinear.f90
    call Dlinear_setup ! subroutine in src/hpkvd/psubs)_Dlinear.f90
    idynax=1
    call evolve_ellipse_full(idynax,zdynax,Ddynax,fcdynax)

    write(*,*) Frho, e_v, p_v
    close(10)
end program runhomel

