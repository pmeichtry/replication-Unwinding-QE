%% Replication files for Cantore and Meichtry (2024)
% Unwinding Quantitative Easing: State Dependency and Household Heterogeneity
%
% This file creates Figure 2 in the paper.
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
opt.nrows_subplots  = 2;
opt.ncols_subplots  = 4;
opt.plot_size       = [0.2 0.2 0.65 0.55];

% Variables (will be created below)
opt.var_names   = char('CS_obs', 'BondIncS_obs', 'LaborIncS_obs', 'ProfitIncS_obs', ...
     'CB_obs', 'BondIncB_obs', 'LaborIncB_obs', 'ProfitIncB_obs');
opt.var_labels  = char('Consumption S', 'Bond demand/interest S', 'Net labor income S', 'Profit income S', ...
    'Consumption B', 'Bond demand/interest B', 'Net labor income B', 'Profit income B');
opt.shock_names	= char('eps_common');

% Legend                     
opt.nolegend    = 0;
opt.legendLoc  = 'southeast';
opt.legendPos  = [0.515 0 0 0.04];
opt.legendOrient = 'horizontal';

% Axes
opt.xlabel      = {'Quarters'};
opt.ylabel      = {'% dev. from SS','','','', '% dev. from SS','','',''};
opt.xtickstep   = 5;

% Font settings
opt.fontsizeTitle   = 11;
opt.fontsizeAxis    = 10;
opt.fontsizeLegend  = 10;
opt.fontType        = 'times';

% Save
opt.savePlot        = 1;  	%Export figures as PDF files
opt.savePath        = '../outputs/';


%% Data preparation

% Calculate new variables
shockName = strtrim(opt.shock_names);

path_bcc_results = { '../simul_results/baseTA/irfs_QE_nobound';...
    '../simul_results/baseTA/irfs_QT_nobound';...
    '../simul_results/baseTA_closeZLB/irfs_QT';...
    '../simul_results/baseTA/irfs_QE_net_of_pref';...
    '../simul_results/baseRA/irfs_QT_nobound';...
    '../simul_results/baseRA/irfs_QE_net_of_pref'
    };

tauD=0; %no redistribution in baseline calibration

for i = 1:length(path_bcc_results)
    select_data = load(path_bcc_results{i});

    if contains(path_bcc_results{i},'RANK')
        lambda=0;
    else
        lambda=0.35;
    end

    % Borrower budget constraint components (in % deviations from SS)
    data_bcc.CB_obs.(shockName) = select_data.irfs.CB_obs.(shockName);
    data_bcc.OtherIncomeB_obs.(shockName) = select_data.irfs.W_obs.(shockName) + select_data.irfs.NB_obs.(shockName) ...
        - select_data.irfs.T_obs.(shockName) + tauD/lambda*select_data.irfs.P_obs.(shockName);
    data_bcc.BondIncB_obs.(shockName) = data_bcc.CB_obs.(shockName) - data_bcc.OtherIncomeB_obs.(shockName); %computed as residual
    data_bcc.LaborIncB_obs.(shockName) = select_data.irfs.W_obs.(shockName) + select_data.irfs.NB_obs.(shockName) ...
        - select_data.irfs.T_obs.(shockName);
    data_bcc.ProfitIncB_obs.(shockName) = tauD/lambda*select_data.irfs.P_obs.(shockName);

    % Saver budget constraint components (in % deviations from SS)
    data_bcc.CS_obs.(shockName) = select_data.irfs.CS_obs.(shockName);
    data_bcc.OtherIncomeS_obs.(shockName) = select_data.irfs.W_obs.(shockName) + select_data.irfs.NS_obs.(shockName) ...
        - select_data.irfs.T_obs.(shockName) + (1-tauD)/(1-lambda)*select_data.irfs.P_obs.(shockName);
    data_bcc.BondIncS_obs.(shockName) = data_bcc.CS_obs.(shockName) - data_bcc.OtherIncomeS_obs.(shockName); %computed as residual
    data_bcc.LaborIncS_obs.(shockName) = select_data.irfs.W_obs.(shockName) + select_data.irfs.NS_obs.(shockName) ...
        - select_data.irfs.T_obs.(shockName);
    data_bcc.ProfitIncS_obs.(shockName) = (1-tauD)/(1-lambda)*select_data.irfs.P_obs.(shockName);

    % Save
    irfs_bcc.irfs = data_bcc;
    save([path_bcc_results{i} '_bdg_constr.mat'], '-struct', 'irfs_bcc')
    
end

% Set borrower-specific IRFs to zero in mat-files of RANK model
varsB = { 'CB_obs'; 'OtherIncomeB_obs'; 'BondIncB_obs'; 'LaborIncB_obs'; 'ProfitIncB_obs'};
for i = 1:length(path_bcc_results)
    if contains(path_bcc_results{i},'RANK')
        fileList = {[path_bcc_results{i} '_bdg_constr.mat']};

        select_data = load(fileList);
        for k = 1:length(varsB)  %borrower-specific variable
            select_data.irfs.(varsB{k}).(shockName) = zeros(1,length(select_data.irfs.(varsB{k}).(shockName))); 
        end
        
        save(fileList, '-struct', 'select_data')
    end
end


%% Create figure

% Paths IRF data
path_TA             = '../simul_results/baseTA/';
path_TA_closeZLB    = '../simul_results/baseTA_closeZLB/';

% Figure 2
opt.shock_labels = char('Budget constraint in % dev: QE/QT (non-binding ZLB, no pref shock) + QT close to ZLB');
opt.results_files = { [path_TA 'irfs_QE_nobound_bdg_constr']; [path_TA 'irfs_QT_nobound_bdg_constr']; [path_TA_closeZLB '/irfs_QT_bdg_constr'] };
opt.model_names = { 'QE, no ZLB'; 'QT, no ZLB'; 'QT, close to ZLB' };
opt.linestyle = {'--', '-', ':'};
opt.linecolor = {[0 0.3570 0.6810], [0.8500, 0.2, 0.0980], [0.4060, 0.6040, 0.1880]};
opt.linewidth = {1.5, 1.5, 1.5};
opt.saveFilename = 'Figure2';
IRF_multi_plot( opt );
