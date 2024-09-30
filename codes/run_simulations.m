%% Replication files for Cantore and Meichtry (2024)
% Unwinding Quantitative Easing: State Dependency and Household Heterogeneity
%
% This file runs the different model simulations in the paper.
% Simulation results are saved in ../simul_results
%
% Required: 
%   - Dynare, version 4.4 or later. Code written and tested with version 4.5.7.
%   - dynareOBC toolkit. Code tested with v3.30.54.1968
% Recommended: install a mixed integer linear programming solver to get
% exact results. See dynareOBC documentation for more details.
%
% C. Cantore and P. Meichtry, September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

addpath('functions')
addpath ('dynareOBC') %dynareOBC toolbox (add without subfolders)

addpath('C:\Program Files (x86)\Dynare\4.5.7\matlab') %Dynare
addpath('C:\Produits\gurobi1003\win64\matlab') %MILP solver


%% GENERAL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model specifications
spec.modelNames = {'TANK', 'TANK_closeZLB', 'RANK', 'RANK_closeZLB'};
    %TANK: TANK-BS model
    %TANK_closeZLB: TANK-BS model with economy close to ZLB
    %RANK: RANK model
    %RANK_closeZLB: RANK model with economy close to ZLB

% Shock specifications
spec.shockLabels = {'QE', 'QT', 'pref', 'QE_pref', 'QT_pref'};
    %QE: Model with QE shock and ZLB
    %QT: Model with QT shock and ZLB
    %pref: Model with preference shock and ZLB (no QE/QT shock)
    %QE_pref: Model with QE shock, preference shock and ZLB
    %QT_pref: Model with QT shock, preference shock and ZLB
    
% IRF extraction
spec.save_results = 1; %extract and save results (=1) or stop after simulation (=0)


save('specification.mat', 'spec');


%% Run model simulations

for m = 1:length(spec.modelNames)
    
    %% SET PARAMETER VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Households
    if contains(spec.modelNames{m}, 'RANK')
        lambda 	= 0.0001;   %Share of Borrowers
    else
        lambda 	= 0.35;
    end
    sigma_c  	= 1;    	%Intertemporal elasticity of substitution
    varphi      = 1;    	%Inverse Frisch elasticity of labor supply

    if contains(spec.modelNames{m}, 'ZLB')
        betta 	= 0.99955;   %Discount factor, Savers
    else
        betta  	= 0.99;
    end
    bettaB      = 0.95;    	%Discount factor, Borrowers
    Dbar        = 0.5;      %Borrowing limit (defined as D/Y)

    rho_A       = 0.8;     	%AR coeff preference shock

    % Firms
    eepsil      = 6;        %Elasticity of substitution between differentiated goods 
    tauD        = 0;      	%Tax on profits
    price_duration = 3.5;   %Number of quarters required to reset prices

    tauS        = (eepsil-1)^(-1); %Production subsidy ensuring zero profits in SS
    rho_Z       = 0.75;   	%AR coeff technology shock (labor augmenting)

    % Government
    rho_G       = 0.9;     	%AR coeff government spending shock
    rho_tauT    = 0.7;      %Tax smoothing in fiscal rule
    rho_tauT_B  = 0.33;     %Tax response to total debt
    rho_tauT_G  = 0.1;    	%Tax response to government spending

    % Central Bank
    phi_pie     = 1.5;   	%Taylor rule coefficient on inflation
    rho_r       = 0.8;      %Interest rate smoothing

    qbar       	= 0.25;     %Steady-state long-term bond holdings by central bank
    rho_q     	= 0.9;      %QE Smoothing
    rho_BL      = 0.9;  	%AR coeff long term bonds shock

    % Long-term bonds
    chi         = 0.975;  	%Long-term bond coupon decay rate
    deltatilde  = 0.3;      %Steady-state ratio of long-term to short-term bonds
    nu          = 0.05; 	%Portfolio share adjustment cost

    % Other steady-state paramaters
    Nss         = 1;       	%SS hours
    PIEss       = 1;       	%SS inflation
    GYss        = 0.2;    	%SS gov spending to GDP ratio
    BYss        = 0.8;      %SS debt to GDP ratio

    % Relative shock sizes within auxiliary shock 'eps_common'
    weight_QEshock = 0.11;  %relative size QE shock
    weight_Ashock  = 0.42;  %relative size pref shock


    save param_values lambda sigma_c varphi betta bettaB Dbar rho_A ...
        eepsil tauD price_duration tauS rho_Z ...
        rho_G rho_tauT rho_tauT_B rho_tauT_G ...
        phi_pie rho_r qbar rho_q rho_BL ...
        chi deltatilde nu ...
        Nss PIEss GYss BYss weight_QEshock weight_Ashock


    %% Run simulation and save results

    for s = 1:length(spec.shockLabels)
        save('modelCount.mat', 'm');
        save('shockCount.mat', 's');
        
        load('specification.mat');
        
        % Steady state base file
        ststFile = 'model_master_steadystate_base.m';
        
        % Run simulation for individual shocks using dynareOBC toolbox
        if strcmp(spec.shockLabels{s}, 'QE')
            copyfile(ststFile, 'model_QE_steadystate.m')
            dynareOBC model_QE.mod ShockScale=5e6 IRFsAroundZero %Bypass MILPsolver='gurobi'

        elseif strcmp(spec.shockLabels{s}, 'QT')
            copyfile(ststFile, 'model_QT_steadystate.m')
            dynareOBC model_QT.mod ShockScale=5e6 IRFsAroundZero %Bypass MILPsolver='gurobi'

        elseif strcmp(spec.shockLabels{s}, 'pref')
            copyfile(ststFile, 'model_pref_steadystate.m')
            dynareOBC model_pref.mod ShockScale=5e6 IRFsAroundZero %Bypass MILPsolver='gurobi'

        elseif strcmp(spec.shockLabels{s}, 'QE_pref')
            copyfile(ststFile, 'model_QE_pref_steadystate.m')
            dynareOBC model_QE_pref.mod ShockScale=5e6 IRFsAroundZero %Bypass MILPsolver='gurobi'

        elseif strcmp(spec.shockLabels{s}, 'QT_pref')
            copyfile(ststFile, 'model_QT_pref_steadystate.m')
            dynareOBC model_QT_pref.mod ShockScale=5e6 IRFsAroundZero %Bypass MILPsolver='gurobi'

        end

        delete *_steadystate.m

        %%%% OPTIONS (see dynareOBC documentation for full list) %%%%
        % ShockScale=[#]            Scale shock up/down (+/-)
        % IRFsAroundZero            Center IRFs around zero (default: IRFs centered around risky SS)
        % Bypass                    Ignore ZLB (equal to standard Dynare simulation). Useful for debugging.
        % MILPsolver=[STRING]       Solver for linear programming problem


        %% Extract and save IRF data
        load('specification.mat');
        
        if spec.save_results ~= 1
            return; %exit script
        end
        
        load('modelCount.mat');
        load('shockCount.mat');

        % Define storage path and variables
        if strcmp(spec.modelNames{m}, 'TANK')
            path_results = '../simul_results/baseTA/';
        elseif strcmp(spec.modelNames{m}, 'TANK_closeZLB')
            path_results = '../simul_results/baseTA_closeZLB/';
        elseif strcmp(spec.modelNames{m}, 'RANK')
            path_results = '../simul_results/baseRA/';
        elseif strcmp(spec.modelNames{m}, 'RANK_closeZLB')
            path_results = '../simul_results/baseRA_closeZLB/';
        end
        
        varNames = char(dynareOBC_.EndoVariables(1:end-1));
    %     varNames = char('Y', 'R', 'PIE', 'C', 'CS', 'CB',    'Y_obs', 'R_obs', 'PIE_obs', 'C_obs', 'CS_obs', 'CB_obs', ...
    %             'N', 'NS', 'NB', 'W', 'profits',    'N_obs', 'NS_obs', 'NB_obs', 'W_obs', 'P_obs', ...
    %             'G', 'T', 'Btot', 'B', 'BL', 'BHL',    'G_obs', 'T_obs', 'Btot_obs', 'B_obs', 'BL_obs', 'BHL_obs', ...
    %             'Q', 'QY', 'q', 'Rn', 'RL', 'RLn', 'V',    'Q_obs', 'QY_obs', 'q_obs', 'Rn_obs', 'RL_obs', 'RLn_obs', 'V_obs', ...
    %             'BB', 'BS', 'BHLB', 'BHLS', 'ZZ', 'PSI', ...
    %             'CSl_obs', 'CBl_obs', 'NBl_obs', 'NSl_obs');
        shockNames = char('eps_common');

        % Extract IRF data for plotting
        [irfsOBC,irfsNoBounds,IRFoffset] = extract_irfs_dynareOBC( dynareOBC_, oo_, varNames , shockNames );
        idx_epscommon = find(contains(string(M_.exo_names),'eps_common'));
        stderr_epscommon = sqrt(M_.Sigma_e(idx_epscommon,idx_epscommon));
        shockScale_dynareOBC = dynareOBC_.ShockScale;

        % Save
        if ~exist('irfs','var')
            % IRFs with binding OBC
            irfs = irfsOBC;
            save([path_results 'irfs_' spec.shockLabels{s} '.mat'], 'irfs','IRFoffset',...
                'weight_QEshock','weight_Ashock','stderr_epscommon','shockScale_dynareOBC')
            clear irfs

            % IRFs with non-binding OBC
            irfs = irfsNoBounds;
            save([path_results 'irfs_' spec.shockLabels{s} '_nobound.mat'], 'irfs','IRFoffset',...
                'weight_QEshock','weight_Ashock','stderr_epscommon','shockScale_dynareOBC')
        else
            error('Variable name irfs already in workspace.')
        end
    end
    
    
    %% Additional data preparation and manipulation
    shockNames = char('eps_common');

    % Calculate pure effects (i.e. net of negative preference shock)
        % File paths of input data
        path_QE_pref = [path_results 'irfs_QE_pref.mat'];
        path_QE_pref_nobound = [path_results 'irfs_QE_pref_nobound.mat'];
        path_QT_pref = [path_results 'irfs_QT_pref.mat'];
        path_QT_pref_nobound = [path_results 'irfs_QT_pref_nobound.mat'];
        path_pref = [path_results 'irfs_pref.mat'];
        path_pref_nobound = [path_results 'irfs_pref_nobound.mat'];
        
        % File paths for saving output data 
        path_QE_net = [path_results 'irfs_QE_net_of_pref.mat'];
        path_QE_nobound_net = [path_results 'irfs_QE_nobound_net_of_pref.mat'];
        path_QT_net = [path_results 'irfs_QT_net_of_pref.mat'];
        path_QT_nobound_net = [path_results 'irfs_QT_nobound_net_of_pref.mat'];
        
        % Compute net IRFs and save
        % QE (with binding ZLB)
        irfs_net_temp = calculate_irfs_net( path_QE_pref, path_pref, shockNames );
        save(path_QE_net, '-struct', 'irfs_net_temp')

        % QE (with non-binding ZLB)
        irfs_net_temp = calculate_irfs_net( path_QE_pref_nobound, path_pref_nobound, shockNames );
        save(path_QE_nobound_net, '-struct', 'irfs_net_temp')

        % QT (with binding ZLB)
        irfs_net_temp = calculate_irfs_net( path_QT_pref, path_pref, shockNames );
        save(path_QT_net, '-struct', 'irfs_net_temp')

        % QT (with non-binding ZLB)
        irfs_net_temp = calculate_irfs_net( path_QT_pref_nobound, path_pref_nobound, shockNames );
        save(path_QT_nobound_net, '-struct', 'irfs_net_temp')


    % Set borrower-specific IRFs to zero in mat-files of RANK model
    if contains(spec.modelNames{m}, 'RANK')
        varsB = { 'UCB'; 'UNB'; 'CB'; 'NB'; 'BB'; 'BHLB';...
            'CB_obs'; 'NB_obs'; 'CBl_obs'; 'NBl_obs';
            };

        % Get file names
        fileList = dir([path_results '*.mat']);
        fileList = {fileList.name};

        % Set zeros and save
        for i = 1:length(fileList)  %file
            temp_data = load([path_results fileList{i}]);

            for j = 1:length(varsB)  %borrower-specific variable
                temp_data.irfs.(varsB{j}).(shockNames) = zeros(1,length(temp_data.irfs.(varsB{j}).(shockNames))); 
            end

            save([path_results fileList{i}], '-struct', 'temp_data')
        end    
    end
end

delete('*Count.mat')
delete('param_values.mat')
delete('specification.mat')
