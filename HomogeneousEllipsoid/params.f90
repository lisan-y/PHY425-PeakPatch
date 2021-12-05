module params

  integer iwant_evmap,nstepmax,nyp,nout,iwant_rd
  real tfac,zinit,dcrit,e_vmax,de_v,p_vbar,dp_v,e_vbar,a_bcrit,Fbar

  real fcoll_3,fcoll_2,fcoll_1,a_3eq,a_2eq,a_1eq
  integer lvirv(3),ivir_strat,iforce_strat

end module params

!c  iforce_strat / 0 no bgnd, 1 std bgnd, 
!c        3 bgnd + nonlinear strain, 4 std bgnd + linear strain /'
!c   ivir_strat / 1 a_jeq=fcoll_3 a_b3,2 a_jeq=fcoll_j a_bj/'


