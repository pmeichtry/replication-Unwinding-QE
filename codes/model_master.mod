%% Replication files for Cantore and Meichtry (2024)
% Unwinding Quantitative Easing: State Dependency and Household Heterogeneity
%
% This file contains the basic model structure of the different 
% specifications covered in the paper.
%
% Available shock specifications (defined in model-specific mod-files):
% - QE_shock_on  	With (=1) or without (=0) QE shock
% - QE_shock_neg  	Negative (=1) or positive (=0) shock to CB long-term bond holdings
% - pref_shock_on  	With (=1) or without (=0) preference shock
%
% C. Cantore and P. Meichtry, September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DECLARATION OF ENDOGENOUS VARIABLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var 
    UCS           $u^{\prime}(c^S)$            (long_name='Marginal Utility of Consumption, Savers')
    UCB           $u^{\prime}(c^B)$            (long_name='Marginal Utility of Consumption, Borrowers')
    UNS           $u^{\prime}(N^S)$            (long_name='Marginal Utility of Leisure, Savers')
    UNB           $u^{\prime}(N^B)$            (long_name='Marginal Utility of Leisure, Borrowers')
    C             $c$                          (long_name='Aggregate Consumption')
    CS            $c^S$                        (long_name='Consumption, Savers')
    CB            $c^B$                        (long_name='Consumption, Borrowers')
    N             $N$                          (long_name='Aggregate Hours')
    NS            $N^S$                        (long_name='Hours, Savers')
    NB            $N^B$                        (long_name='Hours, Borrowers')
    A             ${\theta}$                   (long_name='Preference shock process')
    W             $w$                          (long_name='Real Wage')
    Y             $y$                          (long_name='Real Output')
    R             $r$                          (long_name='Real Interest Rate')
    Rn            $R$                          (long_name='Nominal Interest Rate')
    Rn_unc        $R_{unc}$                    (long_name='Unconstrained Nominal Interest Rate, Taylor Rule')
    PIE           $\Pi$                        (long_name='Inflation')
    MC            $mc$                         (long_name='Real Marginal Cost')
    Z             $Z$                          (long_name='Labor Augmenting shock process')
    profits       $d$                          (long_name='Aggregate Profits')
    T             $t$                          (long_name='Lump Sum Tax')
    G             $g$                          (long_name='Government Spending')
    Btot          $b^{tot}$                    (long_name='Total Debt')
    B             $b$                          (long_name='Short Term Bonds Supply')
    BB            $b^{H,B}$                    (long_name='Short Term Bonds Demand, Borrowers')
    BS            $b^{H,S}$                    (long_name='Short Term Bonds Demand, Savers')
    BL            $b^{L}$                      (long_name='Long Term Bonds Supply')
    BHL           $b^{H,L}$                    (long_name='Long Term Bonds Demand, All Households')
    BHLB          $b^{H,L,B}$                  (long_name='Long Term Bonds Demand, Borrowes')
    BHLS          $b^{H,L,S}$                  (long_name='Long Term Bonds Demand, Savers')
    RL            $r^{L}$                      (long_name='Long Term Bonds Real Return')
    RLn           $R^{L}$                      (long_name='Long Term Bonds Nominal Return')
    V             $V$                          (long_name='Long Term Bonds Price')
    ZZ            ${\Omega}$                   (long_name='Net Purchases of Long Term Bonds by Central Bank')
    Q             $b^{CB,L}$                   (long_name='Value of Long Term Bond Purchases by Central Bank')
    QY            ${\frac{b^{CB,L}}{y}}$       (long_name='Central Bank Purchases to Annualized GDP Ratio')
    q             $q$                          (long_name='Fraction of Market Value of Long Term Bond Purchases by Central Bank')
    PSI           $\Psi$                       (long_name='Lagrange Multiplier on Borrowing Constraint')
 % Variables transformed for log deviations
    Y_obs C_obs N_obs W_obs PIE_obs R_obs Rn_obs CS_obs CB_obs NS_obs NB_obs P_obs
    G_obs T_obs Btot_obs B_obs BL_obs BHL_obs Q_obs QY_obs q_obs RL_obs RLn_obs V_obs
    CSl_obs CBl_obs NSl_obs NBl_obs
    ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DECLARATION OF EXOGENOUS VARIABLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo 
    epsZ           ${\varepsilon^{z}}$       (long_name='Technology shock')
    epsM           ${\varepsilon^{m}}$       (long_name='Monetary Policy shock')
    epsq           ${\varepsilon^{q}}$       (long_name='QE shock')
    epsG           ${\varepsilon^{g}}$       (long_name='Government Spending shock')
    epsA           ${\varepsilon^{\theta}}$  (long_name='Preference shock')
    epsBL          ${\varepsilon^{{b}^L}}$   (long_name='Long Term Bonds shock')
    eps_common                               (long_name='Auxiliary shock for multiple-shock simulation')
    ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DECLARATION OF PARAMETERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters 
   betta          ${\beta^S}$             (long_name='Discount factor, Savers')
   bettaH         ${\beta^H}$             (long_name='Discount factor, Borrowers')
   sigma_c        ${\sigma}$              (long_name='Intertemporal elasticity of substitution')
   varphi         ${\varphi}$             (long_name='Inverse of Frisch elasticity of labor supply')
   //psi_p          ${\psi_p}$              (long_name='Rotemberg price adjustment cost')
   eepsil         ${\epsilon}$            (long_name='Elasticity of substitution between intermediate goods varieties')
   rho_Z          ${\rho_{Z}}$            (long_name='Autoregressive parameter for technology shock')
   rho_G          ${\rho_{G}}$            (long_name='Autoregressive parameter for government spending shock')
   rho_A          ${\rho_{\theta}}$       (long_name='Autoregressive parameter for preference shock')
   rho_BL         ${\rho_{BL}}$           (long_name='Autoregressive parameter for long term bonds shock')
   phi_pie        ${\phi_{\pi}}$          (long_name='Taylor rule coefficient on inflation')
   rho_r          ${\rho_{r}}$            (long_name='Interest rate smoothing')
   zetaS          ${\zeta^S}$             (long_name='Weight on hours in utility, Savers')
   zetaH          ${\zeta^H}$             (long_name='Weight on hours in utility, Borrowers')
   lambda         ${\lambda}$             (long_name='Share of Borrowers')
   tauD           ${\tau^D}$              (long_name='Tax on profits')
   tauS           ${\tau^S}$              (long_name='Production subsidy')
   rho_tauT       ${\rho^{\tau t}}$       (long_name='Tax smoothing in fiscal rule')
   rho_tauT_B     ${\rho^{\tau B}}$       (long_name='Tax response to total debt')
   rho_tauT_G     ${\rho^{\tau G}}$       (long_name='Tax response to government spending')
   chi            $\chi$                  (long_name='Long term bond coupon decay rate')
   deltatilde     $\tilde{\delta}$        (long_name='Steady state ratio of long term to short term bonds')
   nu             $\nu$                   (long_name='Portfolio share adjustment cost')
   qbar           ${\bar{q}}$             (long_name='Steady state long term bond holdings by central bank')
   rho_q          ${\rho_{q}}$            (long_name='QE smoothing')
   Nss            ${\bar{N}}$             (long_name='Steady state Hours')
   PIEss          ${\bar{\Pi}}$           (long_name='Steady state Inflation')
   GYss           ${\frac{\bar{G}}{\bar{Y}}}$  (long_name='Steady state gov spending to GDP ratio')
   BYss           ${\frac{\bar{B}}{\bar{Y}}}$  (long_name='Steady state debt to GDP ratio')
   tr             $tr$                    (long_name='Steady state transfer')
   Dbar           $\bar{D}$               (long_name='Borrowing limit')
   price_duration                         (long_name='Frequency of price reset')
   weight_QEshock                         (long_name='Relative size of QE shock')
   weight_Ashock                          (long_name='Relative size of preference shock')
   ;

%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER VALUES %%
%%%%%%%%%%%%%%%%%%%%%%
load param_values;

set_param_value('betta',betta);
set_param_value('bettaH',bettaH);
set_param_value('rho_Z',rho_Z);
set_param_value('rho_G',rho_G);
set_param_value('rho_A',rho_A);
set_param_value('weight_QEshock',weight_QEshock);
set_param_value('weight_Ashock',weight_Ashock);
set_param_value('lambda',lambda);
set_param_value('sigma_c',sigma_c);
set_param_value('varphi',varphi);
set_param_value('rho_r',rho_r);
set_param_value('phi_pie',phi_pie);
set_param_value('rho_BL',rho_BL);
set_param_value('rho_tauT',rho_tauT);
set_param_value('rho_tauT_B',rho_tauT_B);
set_param_value('rho_tauT_G',rho_tauT_G);
set_param_value('chi',chi);
set_param_value('deltatilde',deltatilde);
set_param_value('nu',nu);
set_param_value('qbar',qbar);
set_param_value('rho_q',rho_q);
set_param_value('eepsil',eepsil);
set_param_value('price_duration',price_duration);
set_param_value('tauD',tauD);
set_param_value('tauS',tauS);
set_param_value('PIEss',PIEss);
set_param_value('Nss',Nss);
set_param_value('GYss',GYss);
set_param_value('BYss',BYss);
set_param_value('Dbar',Dbar);

%%%%%%%%%%%%%%%%%%%%%%
%% MODEL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
model;
#delta=deltatilde*(1-STEADY_STATE(q));
#calvo=1-1/price_duration; %implied Calvo sticky prices parameter
#psi_p=calvo*(eepsil-1)/((1-calvo)*(1-betta*calvo)); %implied Rotemberg sticky prices parameter (exploit first-order equivalence)

%% Households
[name='Marginal Utility of Consumption, Savers']
UCS=A*(CS)^(-1/sigma_c);
[name='Marginal Utility of Leisure, Savers']
UNS=-A*zetaS*NS^varphi;
[name='Labor Supply, Savers']
W=-UNS/UCS;

% [name='Budget constraint, Savers']
% CS+BS+BHLS=R(-1)*BS(-1)+RL(-1)*BHLS(-1)+W*NS+(1-tauD)/(1-lambda)*profits-T-nu/2*(delta*BS/BHLS-1)^2-tr/(1-lambda);

[name='Euler Equation short term bonds, Savers']
UCS*(1+delta*nu/BHLS*(delta*BS/BHLS-1))=betta*R*UCS(+1);
[name='Euler Equation long term bonds, Savers']
UCS*(1-(delta*nu*BS/BHLS^2*(delta*BS/BHLS-1)))=betta*RL*UCS(+1);

[name='Marginal Utility of Consumption, Borrowers']
UCB=A*(CB)^(-1/sigma_c);
[name='Marginal Utility of Leisure, Borrowers']
UNB=-A*zetaH*NB^varphi;
[name='Labor Supply, Borrowers']
W=-UNB/UCB;

[name='Budget constraint, Borrowers']
CB+BB+BHLB=R(-1)*BB(-1)+RL(-1)*BHLB(-1)+W*NB+tauD/lambda*profits-T+tr/lambda-nu/2*(delta*BB/BHLB-1)^2;
[name='Borrowing constraint']
BB+BHLB=-Dbar;  %0=min(PSI, BB+BHLB+Dbar);

[name='Euler Equation short term bonds, Borrowers']
UCB*(1+delta*nu/BHLB*(delta*BB/BHLB-1))=bettaH*R*UCB(+1)+PSI;
[name='Euler Equation long term bonds, Borrowers']
UCB*(1-(delta*nu*BB/BHLB^2*(delta*BB/BHLB-1)))=bettaH*RL*UCB(+1)+PSI;


%% Government & Monetary Policy
[name='Government budget constraint']
B+BL=R(-1)*B(-1)+RL(-1)*BL(-1)+ZZ-T+G;

[name='Total debt']
Btot=B+BL;

[name='Tax rule']
log(T/STEADY_STATE(T))=rho_tauT*log(T(-1)/STEADY_STATE(T))+rho_tauT_B*log(Btot/STEADY_STATE(Btot))+rho_tauT_G*log(G/STEADY_STATE(G));
[name='Supply long term bonds']
log(BL)-log(STEADY_STATE(BL))=rho_BL*(log(BL(-1))-log(STEADY_STATE(BL)))+epsBL;

[name='Long term bonds nominal return']
RLn=(1+chi*V(+1))/V;
[name='Long term bonds real return']
RL=RLn/PIE(+1);
[name='CB long term bonds net purchases']
ZZ=Q-RL(-1)*Q(-1);
[name='QE']
Q=q*BL;
[name='CB purchases to annualized GDP ratio']
QY=Q/(4*Y);

[name='Fisher equation'] 
R=Rn/PIE(+1);
[name='Taylor rule'] 
log(Rn_unc/STEADY_STATE(Rn_unc))=rho_r*(log(Rn_unc(-1)/STEADY_STATE(Rn_unc)))+(1-rho_r)*phi_pie*log(PIE/STEADY_STATE(PIE))+epsM;
Rn=max(1,Rn_unc); %log(Rn)=max(0,log(Rn_unc));

%% Firms
[name='Prodution function']
Y=Z*N;
[name='Real wage']
W=MC*Z;
[name='Profits']
profits=(1-MC-psi_p/2*(PIE-1)^2)*Y;

[name='Non-Linear Phillips Curve'] 
PIE*psi_p*(PIE-1)=(1+tauS)*(1-eepsil)+eepsil*MC+betta*(UCS(+1)/UCS*PIE(+1)*psi_p*(PIE(+1)-1)*Y(+1)/Y);


%% Aggregation & Shocks
[name='Aggregation Consumption']
C=lambda*CB+(1-lambda)*CS;
[name='Aggregation Labor']
N=lambda*NB+(1-lambda)*NS;

[name='Resource constraint'] 
Y=C+G+psi_p/2*(PIE-1)^2*Y;

[name='Bonds short, clearing']
B=lambda*BB+(1-lambda)*BS;

[name='Aggregation long term bonds']
BHL=lambda*BHLB+(1-lambda)*BHLS;
[name='Bonds long, clearing']
BL=BHL+Q; %BHL=BL*(1-q);

[name='QE shock']
@#if QE_shock_on == 1
    @#if QE_shock_neg == 1
        log(q/STEADY_STATE(q))=rho_q*(log(q(-1)/STEADY_STATE(q)))-epsq/STEADY_STATE(q) - weight_QEshock*eps_common/STEADY_STATE(q);
    @#else
        log(q/STEADY_STATE(q))=rho_q*(log(q(-1)/STEADY_STATE(q)))+epsq/STEADY_STATE(q) + weight_QEshock*eps_common/STEADY_STATE(q);
    @#endif
@#else
    log(q/STEADY_STATE(q))=rho_q*(log(q(-1)/STEADY_STATE(q)))+epsq/STEADY_STATE(q);
@#endif

[name='Labor augmenting shock'] 
log(Z)-log(STEADY_STATE(Z))=rho_Z*(log(Z(-1))-log(STEADY_STATE(Z)))+epsZ;

[name='Government spending shock'] 
log(G)-log(STEADY_STATE(G))=rho_G*(log(G(-1))-log(STEADY_STATE(G)))+epsG;

[name='Preference shock']
@#if pref_shock_on == 1
    log(A)-log(STEADY_STATE(A))=rho_A*(log(A(-1))-log(STEADY_STATE(A)))-epsA - weight_Ashock*eps_common;
@#else
    log(A)-log(STEADY_STATE(A))=rho_A*(log(A(-1))-log(STEADY_STATE(A)))-epsA;
@#endif


%% Variables in log (for log-linearization)
Y_obs=log(Y);
C_obs=log(C);
N_obs=log(N);
W_obs=log(W);
PIE_obs=log(PIE);
R_obs=log(R);
Rn_obs=log(Rn);
CS_obs=log(CS);
CB_obs=log(CB);
NS_obs=log(NS);
NB_obs=log(NB);
P_obs=profits/STEADY_STATE(Y);
G_obs=log(G);
T_obs=log(T);
Btot_obs=log(Btot);
B_obs=log(B);
BL_obs=log(BL);
BHL_obs=log(BHL);
Q_obs=log(Q);
QY_obs=log(QY);
q_obs=log(q);
RL_obs=log(RL);
RLn_obs=log(RLn);
V_obs=log(V);

%% Total-population variables (across each household type)
CBl_obs=lambda*CB_obs;
CSl_obs=(1-lambda)*CS_obs;
NBl_obs=lambda*NB_obs; 
NSl_obs=(1-lambda)*NS_obs; 
end;

shocks;
var epsZ; stderr 1;
var epsM; stderr 1;
var epsq; stderr 1;
var epsG; stderr 1;
var epsA; stderr 1;
var epsBL; stderr 1;
var eps_common; stderr 1e-8; %set small enough not to affect model covariance structure
end;

steady;
check;
resid(1);

/*
write_latex_parameter_table;
write_latex_dynamic_model;
write_latex_definitions;
collect_latex_files;
*/

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STOCHASTIC SIMULATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
stoch_simul(order=1,irf=20,periods=0, irf_shocks = ( eps_common ), nograph, noprint);
