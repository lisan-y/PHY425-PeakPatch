module Solvers
  
contains
  
  !----------------------------------------------------------------------
  
  SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
    USE intreal_types
    PARAMETER (NMAX=10)
    DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
    EXTERNAL DERIVS
    HH=H*0.5
    H6=H/6.
    XH=X+HH
    DO I=1,N
       YT(I)=Y(I)+HH*DYDX(I)
    ENDDO
    CALL DERIVS(XH,YT,DYT)
    DO I=1,N
       YT(I)=Y(I)+HH*DYT(I)
    ENDDO
    CALL DERIVS(XH,YT,DYM)
    DO I=1,N
       YT(I)=Y(I)+H*DYM(I)
       DYM(I)=DYT(I)+DYM(I)
    ENDDO
    CALL DERIVS(X+H,YT,DYT)
    DO I=1,N
       YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
    ENDDO
    RETURN
  END SUBROUTINE RK4
  
  !----------------------------------------------------------------------

  FUNCTION elliptic_rd(x,y,z)
    USE intreal_types
    real(sp)  elliptic_rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
    PARAMETER (ERRTOL=0.05,TINY=1.0e-25,BIG=4.5e21,C1=3.0/14.0,C2=1.0/6.0,&
         C3=9.0/22.0,C4=3.0/26.0,C5=0.25*C3,C6=1.5*C4)
    REAL(SP) alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,&
         sqrtz,sum,xt,yt,zt
    !c     COMPUTES CARLSON'S ELLIPTIC INTEGRAL 
    !C     R_D(x,y,z)=1.5*int_0^infty dt/[(t+z)^3/2 (t+x)^1/2 (t+y)^1/2]
    !C     x,y >= 0, AND AT MOST ONE ZERO
    if(min(x,y).lt.0.0.or.min(x+y,z).lt.TINY.or.max(x,y,z).gt.BIG) then
       write(*,*) 'invalid arguments in rd '
    endif
    xt=x
    yt=y
    zt=z
    sum=0.0
    fac=1.0
1   continue
    sqrtx=sqrt(xt)
    sqrty=sqrt(yt)
    sqrtz=sqrt(zt)
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
    sum=sum+fac/(sqrtz*(zt+alamb))
    fac=0.25*fac
    xt=0.25*(xt+alamb)
    yt=0.25*(yt+alamb)
    zt=0.25*(zt+alamb)
    ave=0.2*(xt+yt+3.0*zt)
    delx=(ave-xt)/ave
    dely=(ave-yt)/ave
    delz=(ave-zt)/ave
    if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL) go to 1
    ea=delx*dely
    eb=delz*delz
    ec=ea-eb
    ed=ea-6.0*eb
    ee=ed+ec+ec
    elliptic_rd=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
    return
  END FUNCTION elliptic_rd
  
  !----------------------------------------------------------------------
  
  FUNCTION EL2(X,QQC,AA,BB)
    USE intreal_types
    PARAMETER(PI=3.14159265, CA=.0003, CB=1.E-9)
    IF(X.EQ.0.)THEN
       EL2=0.
    ELSE IF(QQC.NE.0.)THEN
       QC=QQC
       A=AA
       B=BB
       C=X**2
       D=1.+C
       P=SQRT((1.+QC**2*C)/D)
       D=X/D
       C=D/(2.*P)
       Z=A-B
       EYE=A
       A=0.5*(B+A)
       Y=ABS(1./X)
       F=0.
       L=0
       EM=1.
       QC=ABS(QC)
1      B=EYE*QC+B
       E=EM*QC
       G=E/P
       D=F*G+D
       F=C
       EYE=A
       P=G+P
       C=0.5*(D/P+C)
       G=EM
       EM=QC+EM
       A=0.5*(B/EM+A)
       Y=-E/Y+Y
       IF(Y.EQ.0.)Y=SQRT(E)*CB
       IF(ABS(G-QC).GT.CA*G)THEN
          QC=SQRT(E)*2.
          L=L+L
          IF(Y.LT.0.)L=L+1
          GO TO 1
       ENDIF
       IF(Y.LT.0.)L=L+1
       E=(ATAN(EM/Y)+PI*L)*A/EM
       IF(X.LT.0.)E=-E
       EL2=E+C*Z
    ELSE
       write(*,*) 'failure in EL2'
    ENDIF
    RETURN
  END FUNCTION EL2
  
end module Solvers
