function [irfs_net] = calculate_irfs_net( pathM, pathS, shockName )
% This function computes net IRFs. It takes the difference between two IRF 
% data series that are outputs from a dynareOBC simulation.
%
% Inputs:
%   - pathM         File path of structure of IRFs to be used as minuend
%   - pathS         File path of structure of IRFs to be deducted as subtrahend
%  	- shockName     char with the shock name of interest (one shock only)
%
% Outputs:
% 	- irfs_net      [structure] Net IRFs
%
% C. Cantore and P. Meichtry, September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
irfM = load(pathM);
irfS = load(pathS);

% Check if data series have same underlying parameters
stderr = char(['stderr_' shockName(isstrprop(shockName,'alpha'))]);
if ~isequal(irfM.shockScale_dynareOBC, irfS.shockScale_dynareOBC) || ~isequal(irfM.(stderr), irfS.(stderr))
    error('Shock scales or std errors of simulated IRF data to match are unequal.')
end

% Calculate net effect
fnames = fieldnames(irfM.irfs);
for j = 1:length(fnames)
    irfs_diff.(fnames{j}).(shockName) = irfM.irfs.(fnames{j}).(shockName) - irfS.irfs.(fnames{j}).(shockName);
end

% Update structure
irfs_net = irfM;
irfs_net = rmfield(irfs_net,'irfs');
irfs_net.irfs = irfs_diff;

end
