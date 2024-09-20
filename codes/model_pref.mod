%% Replication files for Cantore and Meichtry (2024)
% Unwinding Quantitative Easing: State Dependency and Household Heterogeneity
%
% This file sets up the baseline model with a preference shock, switching
% off the QE/QT shock.
%
% C. Cantore and P. Meichtry, September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% No QE/QT shock 
@#define QE_shock_on = 0
@#define QE_shock_neg = 0

% With preference shock
@#define pref_shock_on = 1


% Include base structure of model
@#include "model_master.mod"
