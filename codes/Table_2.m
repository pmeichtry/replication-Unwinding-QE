%% Replication files for Cantore and Meichtry (2024)
% Unwinding Quantitative Easing: State Dependency and Household Heterogeneity
%
% This file replicates the numbers of Table 2 in the paper.
%
% C. Cantore and P. Meichtry, September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables
var_names   = {'Y_obs'; 'PIE_obs'; 'C_obs';};
var_labels  = {'Output'; 'Inflation'; 'Consumption'};
shock_name  = {'eps_common'};

% Parameters
betta           = 0.99;     %Discount factor to compute PDVs
periods_cumul   = 4;        %PDV horizon

rescale         = 100;      %Rescale factor (set =100 for %)
round_digits    = 2;        %Round numbers to # digits

% Save
saveTable   = 1;            %Write output tables to xlsx files
savePath   	= '../outputs/';


%% Prepare and generate table

% Set path
path_results = {'../simul_results/baseTA/'; ...
    '../simul_results/baseRA/'};

for mm = 1:length(path_results)
    
    % Load data from mat files
    irfs_QE_pref    = load([path_results{mm} 'irfs_QE_pref' '.mat']);
    irfs_QT         = load([path_results{mm} 'irfs_QT' '.mat']); %'irfs_QT_nobound.mat'
    irfs_pref       = load([path_results{mm} 'irfs_pref' '.mat']);

    % Check for same underlying parameters in data
    param = char('shockScale_dynareOBC');
    if ~isequal(irfs_QE_pref.(param), irfs_QT.(param), irfs_pref.(param))
       error('Shock scales of simulated IRF data to match are unequal.')
    end

    param = char('stderr_epscommon');
    if ~isequal(irfs_QE_pref.(param), irfs_QT.(param), irfs_pref.(param))
       error('Std errors of simulated IRF data to match are unequal.')
    end

    % Extract IRF data
    for j = 1:length(var_names)

        % pure QE effect at ZLB (net of preference shock effect)
        data.QE.(var_names{j}) = irfs_QE_pref.irfs.(var_names{j}).(char(shock_name)) - irfs_pref.irfs.(var_names{j}).(char(shock_name));

        % QT effect off-ZLB (w/o preference shock)
        data.QT.(var_names{j}) = irfs_QT.irfs.(var_names{j}).(char(shock_name));

    end

    fnames = fieldnames(data);


    % Cumulate data
    bettaSeries = betta.^(0:(periods_cumul-1));
    for i = 1:length(fnames) %loop over cases
        fName = fnames{i};
        for j = 1:length(var_names) %loop over variables
            select_data = data.(fName).(var_names{j});
            weighted_data = bettaSeries .* select_data(1:periods_cumul);
            dataCumul.(fName).(var_names{j}) = sum(weighted_data);

            clear select_data weighted_data
        end
    end

    % Create table
    for i = 1:length(var_names) %loop over variables
        for j = 1:length(fnames) %loop over cases
            fName = fnames{j};

            colName = [var_names{i}(1:end-4) '_' fName];

            select_data = data.(fName).(var_names{i});
            select_dataCumul = dataCumul.(fName).(var_names{i});
            tbl.(colName) = [select_data(1,1); select_dataCumul]; %[on impact; cumulative]

            clear select_data
        end
    end
    
    TBL = struct2table(tbl);
    TBL{:,:} = round(TBL{:,:}*rescale,round_digits); %round and scale


    % Print and save
    if mm == 1
        disp(['Multipliers on impact and cumulated (PDV for ' num2str(periods_cumul) ' quarters)'])
        TBL.Model_multiplier = {'TANK-BS (impact)    '; 'TANK-BS (cumulative)'};
        TBL = movevars(TBL,"Model_multiplier",'Before',1);
        disp(TBL)
        if saveTable == 1
            writetable(TBL,[savePath 'Table2.txt'],'Delimiter','\t')
        end
    elseif mm == 2
        TBL.Model_multiplier = {'RANK (impact)    '; 'RANK (cumulative)'};
        TBL = movevars(TBL,"Model_multiplier",'Before',1);
        disp(TBL)
        if saveTable == 1
            writetable(TBL,[savePath 'Table2.txt'],'WriteMode','Append',...
                'Delimiter','\t','WriteVariableNames',false) 
        end
    end
end
