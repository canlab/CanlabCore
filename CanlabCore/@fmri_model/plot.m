function plot(obj)
% Plot an fmri_model object
%
% :Usage:
% ::
%
%     plot(obj)

% ..
%    Setup
% ..

TR = obj.xY.RT;

bf = obj.xBF.bf;

% Define sessions and number of conditions
nsess = length(obj.Sess);



% ----------------------------------------------
% Image
% ----------------------------------------------
create_figure('design matrix image');
imagesc(obj.xX.X); 
title('Design matrix (X)');
set(gca, 'YDir', 'reverse');
axis tight;
xlabel('Regressors');
ylabel('Images (TRs)');
colormap gray
drawnow

% ----------------------------------------------
% Basis sets
% ----------------------------------------------
nconds = length(obj.xBF);
create_figure('basis sets', 1, nconds);
for i = 1:nconds
    
    subplot(1, nconds, i);
    xx = linspace(0, obj.xBF(i).length - 1, size(obj.xBF(i).bf, 1));
    
    plot(xx, obj.xBF(i).bf, 'LineWidth', 2);
    title(sprintf('Basis set, Cond. %3.0f, %s', i, obj.xX.name{i}));
    axis tight
    xlabel('Time (sec)');
    
end


% ----------------------------------------------
% Plot
% ----------------------------------------------


switch obj.build_method
% ----------------------------------------------
    case 'Separate sessions'
% ----------------------------------------------
        nr = 2; nc = ceil(nsess) ./ 2;
        
        create_figure('design matrix', nr, nc);
        
        for i = 1:nsess
            subplot(nr, nc, i)
            set(gca, 'FontSize', 10);
            
            [Xs, dummy, C, dummy, names] = get_session_X(obj, i);
            
            plot_subfcn(obj, Xs, C, names);
            title(sprintf('Session %3.0f', i));
            
        end
        
        
        
        
% ----------------------------------------------
    case 'Concatenated sessions'
% ----------------------------------------------
        create_figure('design matrix');
    
        plot_subfcn(obj, obj.xX.X(:, obj.xX.iH), obj.xX.X(:, obj.xX.iC), obj.xX.name);
       
    otherwise
        error('Unknown build method for fmri_model object. See help for valid strings.');
        
end


end % main function




% ----------------------------------------------
% ----------------------------------------------

% SUB-functions

% ----------------------------------------------
% ----------------------------------------------


function plot_subfcn(obj, Xs, C, names)

% ----------------------------------------------
% SETUP
% ----------------------------------------------

% conditions for unique color codes
nconds = size(obj.xX.cond_assignments, 2);

% bf = obj.xBF.bf;
% nbf = size(bf, 2);

colors = {'r' 'g' 'b' 'm' [1 .5 0], [0 1 .5], [0 .5 1], [1 0 .5]};
while length(colors) < nconds, colors = [colors colors]; end

disp('Plotting design matrix cols');

% ----------------------------------------------
% Plot
% ----------------------------------------------

han = plot_matrix_cols([Xs C], 'horiz');
axis tight
drawnow

% ***need to fix for Separate Session design

for i = 1:nconds
%     whstart = nbf * (i - 1) + 1;
%     whend = nbf * i;
    wh = find(obj.xX.cond_assignments(:, i));
    set(han(wh), 'Color', colors{i});
end

if ~isempty(C)
    names = [names {'Covariates'}];
end

yvals = linspace(1, size([Xs C], 2), nconds + 2);
yvals = yvals(2:end-1);

set(gca, 'YTick', yvals, 'YTickLabel', names)

xlabel('Time (TRs)')

title(obj.build_method);

drawnow

end




