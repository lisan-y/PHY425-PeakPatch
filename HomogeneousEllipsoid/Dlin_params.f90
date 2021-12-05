module Dlin_params

  integer, parameter :: itabmax    =  2001
  
  real D_ainfty,D_a0,HD_Hainfty,HD_Ha0,chininfty,chin0

  integer nn,ilm,nmax

  integer imaketab_Dlin, ntab_DLin
  real amintab,amaxtab,chimintab,chimaxtab,atab(itabmax),D_atab(itabmax),dlnD_dlnatab(itabmax),chitab(itabmax),&
       ttab(itabmax),Dtab(itabmax),Dmintab,Dmaxtab,tmintab,tmaxtab


  !decay
  integer ldecay
  real taud,zzd,ad,Omnrdi,Omnrd,Omerd

  real aeq,t0,t1


  ! qtime
  real alogn,a,zz,b,fb


end module Dlin_params

