module homogeneous_ellipse

  use Solvers

  real bvec_2(3)
  real aLam_1,aLam_2,aLam_3,Frho,e_v,p_v

contains

  subroutine evolve_ellipse_full(idynax,zdynax,Ddynax,fcdynax)
    
    !------------------------------------------------------------
    ! INTEGRATION VARIABLE t = scale factor
    ! COSMIC TIME UNIT= (8\pi G/3 rho_b*)^{-1/2}
    
    ! a_idot=[p_i](Ha)^-1 
    ! p_idot=-(3/4)[delta1_e b_i + (2/3 - b_i)]a_b^-3 a_i (Ha)^-1
    ! delta1_e = (a_b/a_1)(a_b/a_2)(a_b/a_3)
    ! dot wrt a_b  .... 
    ! Ha=sqrt(1+OMvac/Omnra_b^3+Omcurv/Omnr a_b) /a_b^{1/2}
    ! b_idot=[-2p_i/a_i + SUM_{j.ne.i} [b_i-b_j]/[a_i^2-a_j^2] 
    !     (a_ip_i-a_jp_j) + b_i SUM_j p_j/a_j](Ha)^-1
    ! INSTEAD WE USE THE BINNEY-TREMAINE FORMULA FOR b_i  
    !------------------------------------------------------------
    
    USE intreal_types
    use input_parameters
    use params
    
    PARAMETER (ny=6)
    PARAMETER (one_third=1.0/3.0)
    parameter (nzdynmax=7)
    dimension y(ny),dy(ny)
    
    dimension zdynv(nzdynmax),delta1_ev(nzdynmax),aturnv(3),dcritv(3)
    integer(i4b) ldynv(nzdynmax),ldyn_pv(nzdynmax)
    external get_derivs
    no_vir=3
    no_turn=3
    no_dens=1
    nzdyn=no_vir+no_turn+no_dens
    ivir0=1
    ivir1=no_vir
    iturn0=no_vir+1
    iturn1=no_vir+no_turn
    idens0=no_vir+no_turn+1
    idens1=no_vir+no_turn+no_dens
    do i=1,3
       aturnv(i)=0.0
    enddo
    dcritv(1)=dcrit
    do ii=1,nzdyn
       ldynv(ii)=0
       delta1_ev(ii)=-1.0
       zdynv(ii)=-1.0
    enddo
    
    nyp=ny
    a_3eq=0.0
    a_2eq=0.0
    a_1eq=0.0
    lvirv(1)=0
    lvirv(2)=0
    lvirv(3)=0
    zdynax=-1
    Ddynax=-1
    fcdynax=-1
    if(iwant_evmap.eq.1) then
       Frho=Fbar
       e_v=e_v+de_v
       if(e_v.gt.e_vmax) then
          stop
       endif
       p_v=p_vbar
    elseif(iwant_evmap.eq.2) then
       Frho=Fbar
       e_v=e_vbar
       p_v=p_v+dp_v
       if(p_v.gt.e_v) then
          stop
       endif
    elseif(iwant_evmap.eq.0) then
       if(e_v.lt.0.0) then
          stop
       endif
    elseif(iwant_evmap.eq.-1) then
       iwant_evmap=0
    endif
    
    aLam_3 = Frho/3.0*(1.0+3.0*e_v+p_v)
    aLam_2 = Frho/3.0*(1.0-2.0*p_v)
    aLam_1 = Frho/3.0*(1.0-3.0*e_v+p_v)
    
    if(aLam_3.lt.aLam_2.or.aLam_3.lt.aLam_1.or.aLam_2.lt.aLam_1) then
       write(*,*) 'problem with aLams, exiting...',aLam_1,aLam_2,aLam_3
       call mpi_finalize(ierr)
       stop
    endif
    
    !C ROUGH ESTIMATE OF COLLAPSE zzc FROM aLam_1
    zzc_Fsph=(aLam_1+aLam_2+aLam_3)/1.686
    
    ! used to be in ic_set
    Rvac_nr=Omvac/Omnr
    Rcurv_nr=Omcurv/Omnr
    zz=zinit  
    a_b=1.0/zz
    t=a_b
    !C THIS ASSUMES THAT WE START IN THE nr DOMINATED REGIME
    
    dlin=Dlinear(a_b,chi,HD_Ha,D_a)
    !C RATIO IS TAKEN HERE wrt H_0 Omnr^{1/2}
    Ha_b_nr=sqrt((1.0+Rvac_nr*a_b**3+Rcurv_nr*a_b)/a_b)
    Ha_b_nrinv=1.0/Ha_b_nr
    y(3)=a_b*(1.0-dlin*aLam_3)
    y(2)=a_b*(1.0-dlin*aLam_2)
    y(1)=a_b*(1.0-dlin*aLam_1)
    y(6)=Ha_b_nr*(1.0-(1.0+HD_Ha)*dlin*aLam_3)
    y(5)=Ha_b_nr*(1.0-(1.0+HD_Ha)*dlin*aLam_2)
    y(4)=Ha_b_nr*(1.0-(1.0+HD_Ha)*dlin*aLam_1)
    !  call ic_set(zinit,aLam_1,aLam_2,aLam_3,t,nyp,y)
    istep=0
111 continue
    istep=istep+1
    a_3=y(3)
    dtstep=tfac*a_3*sqrt(a_3)*Ha_b_nr
    
    !C NOTE FUNNY RK4 OF PRESS REQUIRES FIRST CALL TO derivs
    call get_derivs(t,y,dy)
    call rk4(y,dy,nyp,t,dtstep,y,get_derivs)
    t=t+dtstep
    delta1_e=t**3/(y(1)*y(2)*y(3))
    
    if(ivir1.ge.ivir0) then
       do ii=ivir0,ivir1
          if(lvirv(ii).ne.ldynv(ii)) then
             zdynv(ii)=1.0/t
             ldynv(ii)=1
             delta1_ev(ii)=delta1_e
          endif
       enddo
    endif
    if(iturn1.ge.iturn0) then
       do ii=iturn0,iturn1
          if(ldynv(ii).ne.1) then
             iax=ii-ivir1
             if(y(iax).lt.aturnv(iax)) then
                zdynv(ii)=1.0/t
                ldynv(ii)=1
                delta1_ev(ii)=delta1_e
             else
                aturnv(iax)=y(iax)
             endif
          endif
       enddo
    endif
    
    if(idens1.ge.idens0) then
       do ii=idens0,idens1
          if(ldynv(ii).ne.1) then
             if(delta1_e.ge.dcritv(ii-iturn1)) then
                zdynv(ii)=1.0/t
                ldynv(ii)=1
                delta1_ev(ii)=delta1_e
             endif
          endif
       enddo
    endif
    if(iwant_evmap.eq.0) then
       if((istep.eq.(istep/nout)*nout).or.lvirv(1).eq.1.or.istep.ge.nstepmax) then
          zz=1.0/t
          a_b=1.0/zz
          delta1_e=t**3/(y(1)*y(2)*y(3))
          tau_1=(bvec_2(1)-one_third) 
          tau_2=(bvec_2(2)-one_third) 
          tau_3=(bvec_2(3)-one_third) 
11        format(i5,6(1pe11.3),2x,2i2)
115       format(' tau321=(b_i/2-1/3) ',3(1pe11.3))
       endif
    endif
    if(lvirv(1).eq.1.or.istep.ge.nstepmax.or.t.ge.2.0) then
       if(iwant_evmap.eq.0) write(*,*) 'last axis (1) is virialized '
67     format(' F,e_v,p_v,Lam123: ',6f9.3)
       !c         write(*,66) zturn3,zvir3,zturn2,zvir2,zturn1,zvir1,zdcrit
66     format('zvir123,zturn,dcr: ',7(f9.3))
68     format('zc180sph,zcF: ',2(F9.3),' dvir,dcr: ',4(F9.3))
       if(ihard.eq.1.and.iwant_evmap.ne.3) then
          write(10,67) Frho,e_v,p_v,aLam_1,aLam_2,aLam_3
          write(10,66) (zdynv(ii),ii=1,nzdyn)
          write(10,68) zzc_3sph,zzc_Fsph,delta1_ev(3),delta1_ev(2),delta1_ev(1),delta1_ev(7)
          write(10,169) (delta1_ev(ii),ii=1,nzdyn)
       elseif(ihard.eq.2.and.iwant_evmap.ne.3) then
          write(10,69) Frho,e_v,p_v,zdynv(3),zdynv(2),zdynv(1),zdynv(7),zzc_Fsph
69        format(8F9.3)
169       format('densvir,turn,cut: ',7(f9.3))
       endif
       zdynax=zdynv(idynax)
       if(zdynax.eq.-1) then
          Ddynax=-1
          fcdynax=-1
       else
          Ddynax=Dfnofa(1.0/zdynax)
          fcdynax=Frho*Ddynax
       endif
       return
    endif
    goto 111
  end subroutine evolve_ellipse_full
  
  !----------------------------------------------------------------------
  subroutine get_derivs(t,y,dy)
    USE intreal_types
    use cosmoparams
    use params
   
    PARAMETER (one_third=1.0/3.0)
    PARAMETER (ny=6)
    dimension y(ny),dy(ny)
    dimension avec(3)
    
    do j=1,3
       avec(j)=y(j)
    enddo
    !C THESE b_i ARE WHITE-SILK alpha_i THAT THEY SOLVE BY ODE
    a_b=t
    a_b3=a_b**3
    Ha_b_nr=sqrt((1.0+Rvac_nr*a_b3+Rcurv_nr*a_b)/a_b)
    Ha_b_nrinv=1.0/Ha_b_nr
    
    if(lvirv(3).eq.0.and.(avec(3).le.fcoll_3*a_b)) then
       lvirv(3)=1
       a_3eq=fcoll_3*a_b
       if(ivir_strat.eq.1) then
          a_3eq2=a_3eq*1.001
          a_3eq1=a_3eq*1.0003
       endif
    endif
    if(ivir_strat.eq.2.and.lvirv(2).eq.0) a_3eq2=fcoll_2*a_b
    if(ivir_strat.eq.2.and.lvirv(1).eq.0) a_3eq1=fcoll_1*a_b
    
    if(lvirv(3).eq.1.and.(avec(2).le.a_3eq2.and.lvirv(2).eq.0)) then
       lvirv(2)=1
       a_2eq=a_3eq2
    endif
    if(lvirv(3).eq.1.and.(avec(1).le.a_3eq1.and.lvirv(1).eq.0)) then
       lvirv(1)=1
       a_1eq=a_3eq1
    endif
    
    if(lvirv(3).eq.1) then
       avec(3)=a_3eq
       y(3)=avec(3)
       y(6)=0.0
       dy(3)=0.0
       dy(6)=0.0
    endif
    if(lvirv(2).eq.1) then
       avec(2)=a_2eq
       y(2)=avec(2)
       y(5)=0.0
       dy(2)=0.0
       dy(5)=0.0
    endif
    if(lvirv(1).eq.1) then
       avec(1)=a_1eq
       y(1)=avec(1)
       y(4)=0.0
       dy(1)=0.0
       dy(4)=0.0
    endif
    delta1_e=a_b3/avec(1)/avec(2)/avec(3)
    
    if(iforce_strat.eq.0) then
       d0=-Rvac_nr*Ha_b_nrinv
       d1_int=1.5*delta1_e/a_b3*Ha_b_nrinv
    else
       d0=(0.5/a_b3-Rvac_nr)*Ha_b_nrinv
       d1_int=1.5*(delta1_e-1.0)/a_b3*Ha_b_nrinv
       if(iforce_strat.eq.4.or.iforce_strat.eq.5.or.iforce_strat.eq.6) then
          dlin=Dlinear(a_b,chi,HD_Ha,D_a)
          d1_ext=1.5/a_b3*Ha_b_nrinv*dlin
          Frho_3=Frho/3.0
          !C I THINK 3 IS FUNDAMENTALLY FLAWED
       elseif(iforce_strat.eq.3) then
          d1_ext=1.5*2.5/a_b3*Ha_b_nrinv
       endif
    endif
    if(iforce_strat.ne.5.and.iforce_strat.ne.6) then
       call get_b_2(avec,bvec_2)
       !C THESE b_i ARE WHITE-SILK alpha_i THAT THEY SOLVE BY ODE
    elseif(iforce_strat.eq.6) then
       sum=0.0
       do j=1,1,3
          sum=sum+avec(j)
       enddo
       sum_3=sum/3.0
       do j=1,3
          bvec_2(j)=1.0/3.0
          !c +0.4*(sum_3-avec(j))/sum_3
       enddo
       bvec_2(1)=bvec_2(1)+0.4*dlin*(aLam_1-Frho_3)
       bvec_2(2)=bvec_2(2)+0.4*dlin*(aLam_2-Frho_3)
       bvec_2(3)=bvec_2(3)+0.4*dlin*(aLam_3-Frho_3)
       !C WS approximation (SUM a^-1) doesn't work, SUM a doesn't work,
       !C     AND linear CORRECTION ISN'T QUITE ASYMMETRIC ENOUGH
       !C     BUT ISN'T BAD
    endif
    if(lvirv(3).eq.0) then
       dy(3)=y(6)*Ha_b_nrinv
       if(iforce_strat.eq.0) then
          dy(6)=-y(3)*d1_int*bvec_2(3)
       elseif(iforce_strat.eq.1) then
          dy(6)=-y(3)*(d1_int*bvec_2(3)+d0)
       elseif(iforce_strat.eq.3) then
          dy(6)=-y(3)*(d1_int*bvec_2(3)+d0+d1_ext*(bvec_2(3)-one_third))
       elseif(iforce_strat.eq.4.or.iforce_strat.eq.6) then
          dy(6)=-y(3)*(d1_int*bvec_2(3)+d0+d1_ext*(aLam_3-Frho_3))
       elseif(iforce_strat.eq.5) then
          dy(6)=-y(3)*d0-d1_ext*a_b*aLam_3
       endif
    endif
    if(lvirv(2).eq.0) then
       dy(2)=y(5)*Ha_b_nrinv
       if(iforce_strat.eq.0) then
          dy(5)=-y(2)*d1_int*bvec_2(2)
       elseif(iforce_strat.eq.1) then
          dy(5)=-y(2)*(d1_int*bvec_2(2)+d0)
       elseif(iforce_strat.eq.3) then
          dy(5)=-y(2)*(d1_int*bvec_2(2)+d0+d1_ext*(bvec_2(2)-one_third))
       elseif(iforce_strat.eq.4.or.iforce_strat.eq.6) then
          dy(5)=-y(2)*(d1_int*bvec_2(2)+d0+d1_ext*(aLam_2-Frho_3))
       elseif(iforce_strat.eq.5) then
          dy(5)=-y(2)*d0-d1_ext*a_b*aLam_2
       endif
    endif
    if(lvirv(1).eq.0) then
       dy(1)=y(4)*Ha_b_nrinv
       if(iforce_strat.eq.0) then
          dy(4)=-y(1)*d1_int*bvec_2(1)
       elseif(iforce_strat.eq.1) then
          dy(4)=-y(1)*(d1_int*bvec_2(1)+d0)
       elseif(iforce_strat.eq.3) then
          dy(4)=-y(1)*(d1_int*bvec_2(1)+d0+d1_ext*(bvec_2(1)-one_third))
       elseif(iforce_strat.eq.4.or.iforce_strat.eq.6) then
          dy(4)=-y(1)*(d1_int*bvec_2(1)+d0+d1_ext*(aLam_1-Frho_3))
       elseif(iforce_strat.eq.5) then
          dy(4)=-y(1)*d0-d1_ext*a_b*aLam_1
       endif
    endif
    return
  end subroutine get_derivs
  
  !----------------------------------------------------------------------
  subroutine get_b_2(avec,bvec_2)
    USE intreal_types
    use params
    
    PARAMETER (one_third=1.0/3.0)
    dimension avec(3),bvec_2(3)
    !C THIS ROUTINE RETURNS 
    !C   b_i/2=0.5*a_1a_2a_3 \int du [a_i^2+u] Prod_j [a_j^2+u]^{1/2}
    !C   SO THE ELLIPSOID GRAV POT IS 4PI G RHO_E (b_i/2.0) X_I^2/2  
    !C   IN THE PRINCIPAL AXIS SYSTEM
    r_2=avec(2)/avec(1)
    r_3=avec(3)/avec(1)
    if(r_3.ge.r_2) then
       if(r_3.ge.0.999) then
          bvec_2(1)=one_third
          bvec_2(2)=one_third
          bvec_2(3)=one_third
          return
       elseif(r_3.le.0.001) then
          ecc2=1.0-r_3**2
          one_ecc2=r_3**2
          ecc=sqrt(ecc2)
          one_ecc=0.5*r_3**2
          bvec_2(3)=0.5*((one_ecc2)/ecc2*(1.0/(one_ecc2)-log((1.0+ecc)/(one_ecc))*0.5/ecc))
          bvec_2(2)=bvec_2(3)
          bvec_2(1)=1.0-(bvec_2(2)+bvec_2(3))
          return
       else
          ecc2=1.0-r_3**2
          ecc=sqrt(ecc2)
          bvec_2(3)=0.5*((1.0-ecc2)/ecc2*(1.0/(1.0-ecc2)-log((1.0+ecc)/(1.0-ecc))*0.5/ecc))
          bvec_2(2)=bvec_2(3)
          bvec_2(1)=1.0-(bvec_2(2)+bvec_2(3))
          return
       endif
    endif
    if(iwant_rd.eq.1) then
       r_22=r_2*r_2
       r_32=r_3*r_3
       qrat=one_third*r_2*r_3
       bvec_2(3)=qrat*elliptic_rd(r_22,1.0,r_32)
       bvec_2(2)=qrat*elliptic_rd(r_32,1.0,r_22)
       bvec_2(1)=1.0-(bvec_2(2)+bvec_2(3))
    else
       sinth=sqrt(1.0-r_3**2)
       x = sinth/r_3
       rk2=(1.0-r_2**2)/(1.0-r_3**2)
       rkappa2=1.0-rk2
       rkappa=sqrt(rkappa2)
       FLegell_1 = el2(x,rkappa,1.0,1.0)
       ELegell_2 = el2(x,rkappa,1.0,rkappa2)
       qrat=r_2*r_3/sinth**3
       bvec_2(1)=qrat*(FLegell_1-ELegell_2)/rk2
       bvec_2(3)=qrat*((r_2/r_3)*sinth-ELegell_2)/rkappa2
       bvec_2(2)=1.0-(bvec_2(1)+bvec_2(3))
    endif
    return
  end subroutine get_b_2
  
  !----------------------------------------------------------------------

end module homogeneous_ellipse
