%% Replication files for Cantore and Meichtry (2024)
% Unwinding Quantitative Easing: State Dependency and Household Heterogeneity
%
% This file creates Figures 1, 3, 4, and 5 in the paper.
%
% C. Cantore and P. Meichtry, September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

addpath('functions');

%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IRF settings
opt.periods_plot    = 20;   %Horizon for IRFs
opt.rescale         = 100;  %Rescaling data to plot

% Subplot dimensions
opt.nrows_subplots  = 3;
opt.ncols_subplots  = 3;
opt.plot_size       = [0.2 0.2 0.6 0.65];

% Variables
opt.var_names   = char('Y_obs', 'C_obs', 'PIE_obs', ...
    'CSl_obs', 'CBl_obs', ...
    'W_obs', 'Rn_obs', ...
    'R_obs', 'RL_obs');
opt.var_labels  = char('Output', 'Consumption', 'Inflation', ...
    'Consumption S (total)', 'Consumption B (total)', ...
    'Real Wage', 'Short Nominal Rate', ...
    'Short Real Rate', 'Long Real Rate');
opt.shock_names	= char('eps_common');

% Legend                     
opt.nolegend    = 0;
opt.legendLoc  = 'southeast';
opt.legendPos  = [0.515 0 0 0.04];
opt.legendOrient = 'horizontal';

% Axes
opt.xlabel      = {'Quarters'};
opt.ylabel      = {'% dev. from SS','','',...
    '% dev. from SS','','',...
    '% dev. from SS','',''};
opt.xtickstep   = 5;

% Font settings
opt.fontsizeTitle   = 11;
opt.fontsizeAxis    = 10;
opt.fontsizeLegend  = 10;
opt.fontType        = 'times';

% Save
opt.savePlot        = 1;    %Export figures as PDF files
opt.savePath        = '../outputs/';


%% Create figures

% Paths IRF data
path_TA             = '../simul_results/baseTA/';
path_TA_closeZLB    = '../simul_results/baseTA_closeZLB/';
path_RA             = '../simul_results/baseRA/';

% Figure 1
opt.shock_labels = char('QE/QT (non-binding ZLB, no pref shock) + QT close to ZLB');
opt.results_files = { [path_TA 'irfs_QE_nobound']; [path_TA 'irfs_QT_nobound']; [path_TA_closeZLB 'irfs_QT'] };
opt.model_names = { 'QE, no ZLB'; 'QT, no ZLB'; 'QT, close to ZLB' };
opt.linestyle = {'--', '-', ':'};
opt.linecolor = {[0 0.3570 0.6810], [0.8500, 0.2, 0.0980], [0.4060, 0.6040, 0.1880]};
opt.linewidth = {1.5, 1.5, 1.5};
opt.saveFilename = 'Figure1';
IRF_multi_plot( opt );

% Figure 3
opt.shock_labels = char('QE net of preference shock (binding ZLB) and QT (non-binding ZLB)');
opt.results_files = { [path_TA 'irfs_QE_net_of_pref']; [path_TA 'irfs_QT_nobound'] };
opt.model_names = { 'QE, ZLB'; 'QT, no ZLB' };
opt.linestyle = {'-.', '-'};
opt.linecolor = {[0.35 0.35 0.35], [0.8500, 0.2, 0.0980]};
opt.linewidth = {1.5, 1.5};
opt.saveFilename = 'Figure3';
IRF_multi_plot( opt );

% Figure 4
opt.shock_labels = char('QT (non-binding ZLB, no pref shock): RANK vs. TANK-BS');
opt.results_files = { [path_TA 'irfs_QT_nobound']; [path_RA 'irfs_QT_nobound'] };
opt.model_names = { 'TANK-BS'; 'RANK' };
opt.linestyle = {'-', '--'};
opt.linecolor = {[0.8500, 0.2, 0.0980], [1, 0.5, 0.35]};
opt.linewidth = {1.5, 1.5};
opt.saveFilename = 'Figure4';
IRF_multi_plot( opt );

% Figure 5
opt.shock_labels = char('QE net of preference shock (binding ZLB): RANK vs. TANK-BS');
opt.results_files = { [path_TA 'irfs_QE_net_of_pref']; [path_RA 'irfs_QE_net_of_pref'] };
opt.model_names = { 'TANK-BS'; 'RANK' };
opt.linestyle = {'-.', '--'};
opt.linecolor = {[0.35 0.35 0.35], [0.7 0.7 0.7]};
opt.linewidth = {1.5, 1.5};
opt.saveFilename = 'Figure5';
IRF_multi_plot( opt );
