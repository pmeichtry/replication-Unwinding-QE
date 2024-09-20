function [ ] = IRF_multi_plot( opt )
% This function plots IRF simulations for multiple models
% 
% Input: Structure 'opt' with fields:
%   -periods_plot       # periods to plot
%   -nrows_subplots     # rows for subplots
%   -ncols_subplots     # columns for subplots
%   -plot_size          Figure size ([x1 x2 y1 y2])
%   -results_files      Cell array with mat-file names containing IRFs
%   -model_names        Cell array with model names (used for legend)
%   -var_names          Char array with variable names
%   -var_labels         Char array with variable labels (used for titles)
%   -shock_names        Char array with shock names
%   -shock_labels       Char array with shock labels (used for figure name)
%   -irf_around_zero    Center IRFs around zero (=1) or not. Default: 0
%   -nolegend           Without legend (=1) or with legend (=0). Default: 0
%   -legendLoc          Char array with location of legend wrt. axes. Default: 'northeast'
%   -legendPos          Legend position ([x1 x2 y1 y2]). Default: []
%   -legendOrient       Orientation of legend. Default: 'vertical'
%   -xlabel             Cell array with string label x-axis. Default: 'Periods'
%   -ylabel             Cell array with one string label for y-axis or multiple
%                       labels for each variable. Default: 'Deviation from SS'.
%   -xtickstep          x-axis tick step value. Default: 5
%   -linestyle          Cell array with line styles
%   -linecolor          Cell array with line colors
%   -linewidth          Cell array with line widths
%   -fontsizeTitle      Font size of subplot titles
%   -fontsizeAxis       Font size of subplot axes
%   -fontsizeLegend     Font size of legend
%   -fontType           Char array with font type of titles and axes
%   -rescale            Rescaling factor for data to plot
%
% C. Cantore and P. Meichtry, based on code by J. Swarbrick, September 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default values
if ~isfield(opt,'xlabel')
    opt.xlabel = {'Periods'};
end
if ~isfield(opt,'ylabel')
    opt.ylabel = {'Deviation from SS'};
end
if ~isfield(opt,'xtickstep')
    opt.xtickstep = 5;
end
if ~isfield(opt,'irf_around_zero')
    opt.irf_around_zero = 0;
end
if ~isfield(opt,'nolegend')
    opt.nolegend = 0;
end
if ~isfield(opt,'legendLoc')
    opt.legendLoc = 'northeast';
end
if ~isfield(opt,'legendPos')
    opt.legendPos = [];
end
if ~isfield(opt,'legendOrient')
    opt.legendOrient = 'vertical';
end
if ~isfield(opt,'var_labels')
    opt.var_labels = opt.var_names;
end
if ~isfield(opt,'shock_labels')
    opt.shock_labels = opt.shock_names;
end


%% Data Preparation
opt.size_var        = size(opt.var_names);
opt.num_var         = opt.size_var(1);
opt.size_shock_names = size(opt.shock_names);
opt.number_shocks   = opt.size_shock_names(1);
opt.num_datasets    = length(opt.results_files);
opt.dim_fig         = opt.nrows_subplots*opt.ncols_subplots;
opt.num_figs        = ceil(opt.num_var/opt.dim_fig);

index = 0;
for ii=1:opt.num_datasets
    clearvars -except ii opt data index
    index = index+1;
    load([cell2mat(strcat((opt.results_files(index,:)),'.mat'))]);
    for jj = 1:opt.number_shocks
    curr_shock = strtrim(opt.shock_names(jj,:));
        for kk=1:opt.num_var
            try
                if max(abs(irfs.(strtrim(opt.var_names(kk,:))).(curr_shock)))>1e-9
                    if opt.irf_around_zero == 1
                        irfs.plot.(strtrim(opt.var_names(kk,:)))(jj,:) = irfs.(strtrim(opt.var_names(kk,:))).(curr_shock)(:,1:opt.periods_plot)+IRFoffset.(strtrim(opt.var_names(kk,:))).(curr_shock);
                    else
                        irfs.plot.(strtrim(opt.var_names(kk,:)))(jj,:) = irfs.(strtrim(opt.var_names(kk,:))).(curr_shock)(:,1:opt.periods_plot);
                    end
                else
                    eval( strcat('irfs.',strtrim(opt.var_names(kk,:)),'(',num2str(jj),',:) = zeros(1,opt.periods_plot);'));
                    irfs.plot.(strtrim(opt.var_names(kk,:)))(jj,:) = zeros(1,opt.periods_plot);
                end
            catch
                disp(['Problem storing IRF for ',strtrim(opt.var_names(kk,:))])
                irfs.plot.(strtrim(opt.var_names(kk,:)))(jj,:) = zeros(1,opt.periods_plot);
            end
        end
    end
    for jj = 1:opt.num_var
        curr_var = strtrim(opt.var_names(jj,:));
        data(:,:,index,jj) = irfs.plot.(curr_var)(:,1:opt.periods_plot);
    end
end

data(abs(data)<1e-12)=0; %set very small number =0

clearvars -except opt data


%% Create figure
for shocks = 1:opt.number_shocks  
    for figs = 1:opt.num_figs
        figure('Name',strtrim(opt.shock_labels(shocks,:)));
        for pos_subplot = 1:opt.dim_fig
            vars = pos_subplot + opt.dim_fig*(figs-1);
            if vars <= opt.num_var
                subplot(opt.nrows_subplots,opt.ncols_subplots,pos_subplot)
                for jj = 1:opt.num_datasets
                    plot((opt.rescale).*data(shocks,1:opt.periods_plot,jj,vars), opt.linestyle{jj}, 'color',opt.linecolor{jj}, 'LineWidth',opt.linewidth{jj}); hold on;
                    plot(zeros(1,opt.periods_plot),'-','HandleVisibility','off','Color','k','Linewidth',0.05); hold on;
                end
                set(gca,'XTick',0:opt.xtickstep:opt.periods_plot,'FontSize',opt.fontsizeAxis,'FontName',opt.fontType);
                if numel(opt.ylabel) == 1
                    ylabel(opt.ylabel,'FontSize',opt.fontsizeAxis,'FontName',opt.fontType);
                elseif numel(opt.ylabel) > 1
                    ylabel(opt.ylabel{vars},'FontSize',opt.fontsizeAxis,'FontName',opt.fontType);
                end
                grid off
                titlename=strtrim(opt.var_labels(vars,:));
                title(titlename,'FontSize',opt.fontsizeTitle,'FontName',opt.fontType,'FontWeight','normal');
                axis tight;
            end
            if pos_subplot >= opt.dim_fig-(opt.ncols_subplots-1)
                xlabel(opt.xlabel,'FontSize',opt.fontsizeAxis,'FontName',opt.fontType);
            end
        end
        
        if ~opt.nolegend
            Lgnd = legend(opt.model_names,'FontSize',opt.fontsizeLegend,'FontName',opt.fontType); %'interpreter','latex'
            Lgnd.Location = opt.legendLoc;
            if ~isempty(opt.legendPos)
                Lgnd.Position = opt.legendPos;
                Lgnd.Orientation = opt.legendOrient;
            end
        end
        
        if isfield(opt,'plot_size')
            set(gcf, 'units','normalized', 'position',opt.plot_size)
        end
        
        % Save
        if opt.savePlot == 1
            set(gcf,'PaperUnits','normalized','PaperPosition',[0.2 0.2 0.6 0.65])
            if opt.num_figs == 1
                plotName = [opt.savePath opt.saveFilename '.pdf'];
            elseif opt.num_figs > 1
                plotName = [opt.savePath opt.saveFilename '_' num2str(figs) '.pdf'];
            end
            exportgraphics(gcf, plotName, 'ContentType','vector','Resolution',300)
        end
        
    end
end

end
