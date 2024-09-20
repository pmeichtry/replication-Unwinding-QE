%% Replication files for Cantore and Meichtry (2024)
% Unwinding Quantitative Easing: State Dependency and Household Heterogeneity
%
% This file sets up the baseline model with a QT shock, switching off the 
% preference shock.
%
% C. Cantore and P. Meichtry, September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% With QE shock 
@#define QE_shock_on = 1

% Negative shock to CB long-term bond holdings (=QT)
@#define QE_shock_neg = 1

% No preference shock
@#define pref_shock_on = 0


% Include base structure of model
@#include "model_master.mod"
