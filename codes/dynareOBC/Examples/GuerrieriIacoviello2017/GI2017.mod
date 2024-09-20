// This is derived from code kindly provided by Matteo Iacoviello here: https://www2.bc.edu/matteo-iacoviello/research.htm

var aj ak ap arr aw az b bnot c c1 data_ctot data_dp data_dwtot data_ik data_ntot data_q data_r data_rnot dp dp1 dp2 dp3 dw dw1 h1 ik k lm maxlev n n1 q qk r rk rnot uc uc1 uh uh1 w w1 xp xw xw1 y z_j;
varexo eps_j eps_k eps_p eps_r eps_w eps_z;

parameters BETA BETA1 EC EH ETA JEI M ALPHA PHIK DK LAGP LAGW PIBAR INDTR SIGMA TAYLOR_P TAYLOR_Q TAYLOR_R TAYLOR_Y TETAP TETAW XP_SS XW_SS RHO_J RHO_K RHO_P RHO_R RHO_W RHO_Z STD_J STD_K STD_P STD_R STD_W STD_Z RHO_J2 RHOD ITAYLOR_W;

load_params_and_steady_state( 'Steady.txt' );

model;

    # llr = 1 / BETA;
    # llrk = llr - (1-DK);
    # llxp = XP_SS;
    # llxw = XW_SS;
    # llxw1 = XW_SS;
    # lllm = (1 - BETA1/BETA) / (1 - BETA1*RHOD/PIBAR);
    # QHTOC = JEI/(1-BETA);
    # QH1TOC1 = JEI/(1-BETA1-lllm*M*(1-RHOD));
    # KTOY = ALPHA/(llxp*llrk);
    # BTOQH1 = M*(1-RHOD)/(1-RHOD/PIBAR);
    # C1TOY = (1-ALPHA)*SIGMA/(1+(1/BETA-1)*BTOQH1*QH1TOC1)*(1/llxp);
    # CTOY = (1-C1TOY-DK*KTOY);
    # lln = ((1-SIGMA)*(1-ALPHA)/(llxp*llxw*CTOY))^(1/(1+ETA));
    # lln1 = (SIGMA*(1-ALPHA)/(llxp*llxw1*C1TOY))^(1/(1+ETA));
    # lly = KTOY^(ALPHA/(1-ALPHA))*(lln^(1-SIGMA))*lln1^SIGMA;
    # llctot = lly-DK*KTOY*lly;
    # llik = DK*KTOY*lly; 
    # llk = KTOY*lly; 
    # llq = QHTOC*CTOY*lly + QH1TOC1*C1TOY*lly;

    c + c1 + ik = y;
    uc = BETA * r/dp(1)*uc(1);
    w*uc/xw = az*n^ETA;
    q*uc = uh + BETA*q(+1)*uc(1);
    c1 + q*(h1-h1(-1)) + r(-1)*b(-1)/dp = w1*n1 + b + INDTR*log(ap);
    uc1*(1-lm) = BETA1 * (r/dp(1)-RHOD*lm(+1)/dp(1))*uc1(1);
    w1*uc1/xw1 = az*n1^ETA;
    q*uc1 = uh1 + BETA1*q(+1)*uc1(1) + lm*(1-RHOD)*uc1*M*q;
    y = n^((1-ALPHA)*(1-SIGMA))*n1^((1-ALPHA)*SIGMA)*k(-1)^ALPHA;
    (1-ALPHA)*(1-SIGMA)*y = xp*w*n;
    (1-ALPHA)*(SIGMA)*y = xp*w1*n1;
    log(dp/PIBAR) - LAGP*log(dp(-1)/PIBAR) = BETA*(log(dp(1)/PIBAR)-LAGP*log(dp/PIBAR)) - ((1-TETAP)*(1-BETA*TETAP)/TETAP)*(log(xp/XP_SS)) + (1-INDTR)*log(ap);
    log(dw/PIBAR) - LAGW*log(dw(-1)/PIBAR) = BETA*(log(dw(+1)/PIBAR)-LAGW*log(dw/PIBAR)) - ((1-TETAW)*(1-BETA*TETAW)/TETAW)*log(xw/XW_SS) + log(aw);
    log(dw1/PIBAR) - LAGW*log(dw1(-1)/PIBAR) = BETA*(log(dw1(+1)/PIBAR)-LAGW*log(dw1/PIBAR)) - ((1-TETAW)*(1-BETA*TETAW)/TETAW)*log(xw1/XW_SS) + log(aw);
    log(rnot) = TAYLOR_R*log(r(-1)) + (1-TAYLOR_R)*(TAYLOR_P)*(0.25*log(dp/PIBAR)+0.25*log(dp1/PIBAR)+0.25*log(dp2/PIBAR)+0.25*log(dp3/PIBAR)) + (1-TAYLOR_R)*TAYLOR_Y*(log(y/lly)) + (1-TAYLOR_R)*TAYLOR_Q/4*(log(q/q(-1))) + (1-TAYLOR_R)*log(PIBAR/BETA) + log(arr);
    dp1=dp(-1);
    dp2=dp1(-1);
    dp3=dp2(-1);
    uc = (1-EC)/(1-BETA*EC)*(az/(c-EC*c(-1))-BETA*EC*az(+1)/(c(+1)-EC*c));
    uc1 = (1-EC)/(1-BETA1*EC)*(az/(c1-EC*c1(-1))-BETA1*EC*az(+1)/(c1(+1)-EC*c1));
    uh = (1-EH)/(1-BETA*EH)*JEI*(az*aj/(1-h1-EH*(1-h1(-1)))-BETA*EH*az(+1)*aj(+1)/((1-h1(+1))-EH*(1-h1)));
    uh1 = (1-EH)/(1-BETA1*EH)*JEI*(az*aj/(h1-EH*h1(-1))-BETA1*EH*az(+1)*aj(+1)/(h1(+1)-EH*h1));
    //%uh = az*JEI*aj/(1-h1);
    //%uh1 = az*JEI*aj/h1;
    uc*qk*(1-PHIK*(ik-ik(-1))/llik) = uc - BETA*uc(+1)*qk(+1)*PHIK*(ik(+1)-ik)/llik;
    uc*qk/ak = BETA*uc(+1)*(rk(+1)+qk(+1)*(1-DK)/ak(+1));
    k/ak = ik + (1-DK)*k(-1)/ak;
    ALPHA*y = xp*rk*k(-1);
    dw = w*dp/w(-1);
    dw1 = w1*dp/w1(-1);

    log(aj) = RHO_J*log(aj(-1)) + z_j;
    z_j = RHO_J2*z_j(-1) + eps_j;
    log(ak) = RHO_K*log(ak(-1)) + eps_k;
    log(ap) = RHO_P*log(ap(-1)) + eps_p;
    log(aw) = RHO_W*log(aw(-1)) + eps_w;
    log(arr) = RHO_R*log(arr(-1)) + eps_r;
    log(az) = RHO_Z*log(az(-1)) + eps_z;

    data_ctot = (c + c1)/llctot-1;
    data_ntot = SIGMA*n/lln+(1-SIGMA)*n1/lln1-1;
    data_ik = (ik/llik)-1; 
    data_q = q/llq-1;
    data_r = r-PIBAR/BETA; 
    data_rnot = rnot-PIBAR/BETA; 
    data_dp = log(dp/PIBAR);
    data_dwtot =log((SIGMA*dw1 + (1-SIGMA)*dw)/PIBAR);

    0 = min( bnot - b, lm ); // b = bnot OR lm=0, b <= bnot, lm >= 0
    bnot = (1-RHOD)*M*q*h1+RHOD*b(-1)/dp;
    maxlev = b-bnot;
    log(r) = max( 0, log(rnot) );
    
    #irf1 = ( q - STEADY_STATE( q ) ) / STEADY_STATE( q );
    #irf2 = ( c + c1 - STEADY_STATE( c + c1 ) ) / STEADY_STATE( c + c1 );
    #irf3 = data_ntot - STEADY_STATE( data_ntot );
    #irf4 = lm / 100;

end;

shocks;
    var eps_j; stderr STD_J;
    var eps_k; stderr STD_K;
    var eps_p; stderr STD_P;
    var eps_r; stderr STD_R;
    var eps_w; stderr STD_W;
    var eps_z; stderr STD_Z;
end;

steady_state_model;

    r = PIBAR / BETA ;
    rk = 1/BETA - (1-DK) ;
    xp = XP_SS ;
    xw = XW_SS ;
    xw1 = XW_SS ;
    lm = (1 - BETA1/BETA) / (1 - BETA1*RHOD/PIBAR) ;
    QHTOC_ = JEI/(1-BETA);
    QH1TOC1_ = JEI/(1-BETA1-lm*(1-RHOD)*M);
    KTOY_ = ALPHA/(xp*rk);
    BTOQH1_ = M*(1-RHOD)/(1-RHOD/PIBAR) ;
    C1TOY_ = (1-ALPHA)*SIGMA/(1+(1/BETA-1)*BTOQH1_*QH1TOC1_)*(1/xp) ;
    CTOY_ = (1-C1TOY_-DK*KTOY_) ;
    n = ((1-SIGMA)*(1-ALPHA)/(xp*xw*CTOY_))^(1/(1+ETA));
    n1 = (SIGMA*(1-ALPHA)/(xp*xw1*C1TOY_))^(1/(1+ETA));
    y = KTOY_^(ALPHA/(1-ALPHA))*(n^(1-SIGMA))*n1^SIGMA ;
    c = CTOY_*y;
    c1 = C1TOY_*y;
    k = KTOY_*y;
    ik = DK*k;
    w = xw*c*n^ETA;
    w1 = xw1*c1*n1^ETA;
    q = QHTOC_*c + QH1TOC1_*c1 ;
    h = QHTOC_*c/q ;
    h1 = QH1TOC1_*c1/q ;
    b = BTOQH1_*q*h1 ;
    uc = 1/c;
    uc1 = 1/c1;
    uh = JEI/h;
    uh1 = JEI/h1;
    dp = PIBAR ;
    dp1 = dp;
    dp2 = dp;
    dp3 = dp;
    dw = PIBAR ;
    dw1 = PIBAR ;
    aj = 1 ;
    arr =1;
    az = 1;
    ak=1;
    ap=1;
    aw=1;
    qk = 1;
    rnot = r;
    bnot = b;
    maxlev = 0;
    data_ctot = 0 ;
    data_ik = 0;   
    data_q = 0;
    data_r = 0; 
    data_rnot = 0;
    data_dp = 0 ;
    data_dwtot =0;
    data_ntot=0;
    z_j=0;

end;

steady;
check;

stoch_simul( order=1, irf=0, periods=38 ) irf1 irf2 irf3 irf4;
