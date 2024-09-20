function [ys,check] = model_master_steadystate_base(ys,exo)
% Base file for computing the steady state, called during the simulations 
% of the various model specifications.
% Code runs with Dynare 4.5. For Dynare 4.6, see commented out sections.
%
% C. Cantore and P. Meichtry, September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
global M_ 

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end

% initialize indicator
check = 0;

%%%% If Dynare 4.6, use this section instead %%%%
% %% Initialization
% function [ys,params,check] = model_master_steadystate_base(ys,exo,M_,options_)
% 
% % read out parameters to access them with their name
% NumberOfParameters = M_.param_nbr;
% for ii = 1:NumberOfParameters
%   paramname = M_.param_names{ii};
%   eval([ paramname ' = M_.params(' int2str(ii) ');']);
% end
% % initialize indicator
% check = 0;


%% Model equations
PIE=PIEss;
Z=1;
A=1;

R=1/betta;
RL=R;
Rn=PIE/betta;
RLn=Rn;
Rn_unc=Rn;

N=Nss;
NB=N;
NS=N;

Y=Z*N;
MC=(1+tauS)*(eepsil-1)/eepsil;
W=MC*Z;
profits=Y-W*N; 

V=1/(RLn - chi);
q=qbar;

B=BYss*4*Y*1/(1+deltatilde);
BL=B*deltatilde; %BL=B*delta/(1-q);
Btot=B+BL;

BB=-Dbar/(1+deltatilde*(1-q)); %BB=-Dbar/(1+delta);
BHLB=-Dbar-BB;
BS=(B-lambda*BB)/(1-lambda);

Q=q*BL;
QY=Q/(4*Y);
ZZ=Q*(1-RL);
BHL=BL-Q;
BHLS=(BHL-lambda*BHLB)/(1-lambda);

G=GYss*Y;
C=(1-GYss)*Y;
T=G+ZZ-B*(1-R)-BL*(1-RL);

CB=C;
CS=C;
UCB=A*(CB)^(-1/sigma_c);
UCS=A*(CS)^(-1/sigma_c);

zetaH=(UCB)*(W)/(A*NB^varphi);
UNB=-A*zetaH*NB^varphi;
zetaS=(UCS)*(W)/(A*NS^varphi);
UNS=-A*zetaS*NS^varphi;

PSI=UCB*(1-bettaH/betta);
tr=lambda*(CB+(1-R)*BB+(1-RL)*BHLB-W*NB-tauD/lambda*profits+T);

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
P_obs=profits/Y;
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

CBl_obs=lambda*CB_obs;
CSl_obs=(1-lambda)*CS_obs;
NBl_obs=lambda*NB_obs; 
NSl_obs=(1-lambda)*NS_obs; 


%% Update parameters and variables

for iter = 1:length(M_.params)  %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr;
for ii = 1:NumberOfEndogenousVariables  %auxiliary variables are set automatically
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end


%%%% If Dynare 4.6, use this section instead %%%%

% params=NaN(NumberOfParameters,1);
% for iter = 1:length(M_.params)  %update parameters set in the file
%   eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
% end
% 
% NumberOfEndogenousVariables = M_.orig_endo_nbr;
% for ii = 1:NumberOfEndogenousVariables  %auxiliary variables are set automatically
%   varname = M_.endo_names{ii};
%   eval(['ys(' int2str(ii) ') = ' varname ';']);
% end

end