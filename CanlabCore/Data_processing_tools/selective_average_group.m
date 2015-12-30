function [group_avgs, group_stes, subject_avgs, subject_stes, braindata] = selective_average_group(V, onsets, xyz_mm_pos, varargin)
% ..
%    DOCUMENTATION NOT COMPLETE!
% ..
% uses selective_average.m
% used in selective_average_interactive_view_init.m
%
% :Examples:
% ::
%
%    [group_avgs, group_stes, subject_avgs, subject_stes] = selective_average_group(V, onsets, vox, 'basepts', 1:2, 'plotstes', 0);
%
%    % Format onsets from onsets2 (NSF study format) into correct format for this function
%    N = length(imgs); n_conditions = size(eventdesign{1}, 2);
%    onsets = cell(1, N);
%    for i = 1:N
%       for j = 1:n_conditions
%           onsets{i}{j} = onsets2{i}(find(eventdesign{i}(:, j)));
%       end
%    end

    basepts = [];
    t = 20;
    doplot = 1; % not used
    plotstes = 1;
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}

                % functional commands

                case {'baseline' 'basepts' 'basepoints'}, basepts = varargin{i+1};
                case {'t', 'timepoints'}, t = varargin{i+1};
                case 'plot', doplot = 1;

                case 'plotstes', plotstes = varargin{i+1};


                case {'S'}, S = varargin{i+1};
                case {'scans'}, scans_per_sess = varargin{i+1};
                case {'I'}, I = varargin{i+1};

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    
    % Get voxel coordinates from mm coords for each subject
    % Subjects may be in different space!
    % -----------------------------------------------------
     N = length(V);
     vox = zeros(3, N);
     
     for i = 1:N

         %pos = spm_orthviews('Pos');
         vox(:, i) = (mm2voxel(xyz_mm_pos, V{i}(1)))';
     end

     vox(4, :) = 1;
     
     meanvox = mean(vox, 2);
     stdvox = std(vox')';
     inreg = all(stdvox == 0);
     if inreg, regstr = 'All subjects in same space.'; else regstr = 'Subjects are in different spaces.'; end

     fprintf('Position: %3.0f %3.0f %3.0f , Voxel mean: %3.2f %3.2f %3.2f, Std: %3.2f %3.2f %3.2f  %s\n', xyz_mm_pos(1), xyz_mm_pos(2), xyz_mm_pos(3), ...
         meanvox(1), meanvox(2), meanvox(3), stdvox(1), stdvox(2), stdvox(3), regstr);
     
     go = input('Plot selective averages for this voxel? (1 or 0) : ');
    
     if ~go
         group_avgs = [];
         group_stes = [];
         subject_avgs = [];
         subject_stes = [];
         braindata = [];
        
         return
     end

    % conditions per subject; These are assumed to be the same for all
    % subjects!
    n_conditions = length(onsets{1});

    % Load data for this voxel
    % -----------------------------------------------------
   

    fprintf('Loading voxel data for all subjects: ');
    tic
    braindata = cell(1, N);
    
    for i = 1:N
        fprintf('%3.0f ', i); 
        braindata{i} = spm_get_data(V{i}, vox(:, i)); 
    end
    fprintf(' %3.0f s\n', toc)

    
    % Should high-pass filter here!
    % -----------------------------------------------------
    fprintf('High-pass filtering\n');
    for i = 1:N
        braindata{i} = hpfilter(braindata{i}, [], S, scans_per_sess, I);
    end
    
    
    % Get selective averages for each subject
    % -----------------------------------------------------

    nx = ceil(sqrt(N));
    ny = floor(sqrt(N));
    create_figure('Selective average individual plot', nx, ny);

    subject_avgs = cell(1, n_conditions);
    subject_stes = cell(1, n_conditions);

    for i = 1:N

        subplot(nx, ny, i);
        [averages, stderrs] = selective_average(braindata{i}, onsets{i}, 't', t, 'plot');
        xlabel('Time'); ylabel('Response');

        for j = 1:n_conditions
            subject_avgs{j}(:, i) = averages{j};
            subject_stes{j}(:, i) = stderrs{j};
        end

    end

    % Adjust data, if asked for: subtract baseline points
    % -----------------------------------------------------
    if ~isempty(basepts)

        fprintf('Subtracting baseline timepoints:')
        fprintf(' %3.0f', basepts);
        fprintf('\n')
        
        for j = 1:n_conditions
            basemean = nanmean(subject_avgs{j}(basepts, :), 1);

            basemean = basemean(ones(t, 1), :);

            subject_avgs{j} = subject_avgs{j} - basemean;

        end
    end


    % Plot
    % -----------------------------------------------------
    create_figure('Selective Average: Group Plot');

    colors = {'ro-' 'go-' 'bo-' 'yo-' 'co-' 'mo-' 'r^:' 'g^:' 'b^:' 'y^:' 'c^:' 'm^:'};

    x = (1:t)' - 1; 
    
    group_avgs = cell(1, n_conditions);
    group_stes = cell(1, n_conditions);
    
    for j = 1:n_conditions

        group_avgs{j} = nanmean(subject_avgs{j}, 2);
        group_stes{j} = ste(subject_avgs{j}');

        if plotstes, tor_fill_steplot(subject_avgs{j}', colors{j}, 0, [], x); end

        plot(x, group_avgs{j}, colors{j}, 'LineWidth', 3);

    end


end
