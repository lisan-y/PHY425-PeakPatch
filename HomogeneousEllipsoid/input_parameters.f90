module input_parameters

  use cosmoparams ! src/modules/GlobalVariables/cosmoparams.f90

  character*128 &
       filein, fileout, densfilein, pkfile, filterfile, fielddir

  integer &
       ireadfield,ibatch,ioutshear,&
       num_redshifts,&
       iseedFFT,&
       nlx,nly,nlz,&
       nbuff,next,ievol,&
       ihard,&
       debug,wsmooth,ioutfield,NonGauss,&
       ilpt,iwrap,iwant_field_part,largerun

  real &
       global_redshift,maximum_redshift,&
       fhdmclus,&
       dcore_box,dL_box,cenx,ceny,cenz,&
       zinit_fac_ell,dFrho,&
       rmax2rs,fNL

  ! Other parameters

  real    fbuffer_box
  real    biasold,biasnew

  integer nboxes, idocore, idoboxf, ifmt
  integer min_core_box, max_core_box
  logical verbose
  integer*8 ncxm,ncxp,ncym,ncyp,nczm,nczp

end module input_parameters
