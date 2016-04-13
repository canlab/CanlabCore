function selective_average_interactive_view_init(imgs, onsets, scans_per_sess, TR, hp_length, t, basepts, plotstes)
% :Usage:
% ::
%
%    selective_average_interactive_view_init(imgs, onsets, 20, 1:2, 1);
%
% Initialize point-and-click data extraction and plotting of selective
% averages
%
% :Inputs:
%
%   **imgs:**
%        cell array of image files, one cell per subject
%
%        imgs is cell array of time series image names for each subject
%
%   **onsets:**
%        cell array of onsets
%
%        Each cell contains its own cell array, with one cell per
%        condition.
%
%   **t:**
%        time points in average to estimate
%
%   **basepts:**
%        indices of baseline points to subtract from individual averages
%
%   **plotstes:**
%        plot standard errors: 1 or 0


    disp('Initializing Selective Average interactive viewer') % need to do the stuff below only once

    disp('-------------------------------------------------')

    %     disp('Loading mediation_SETUP.mat')
    %     [SETUP, imgs, wh_is_image, name] = load_mediation_setup;


    N = length(imgs);



    % Setting up filtering
    % ------------------------------------------------------
    fprintf('Setting up high-pass filtering (dummy params for first 2 time points per run): ');
    % %     I = cell(1, N);
    % %     S = cell(1, N);
    tic

    for i = 1:N
        [tmp, I, S] = hpfilter(rand(sum(scans_per_sess), 1), TR, hp_length, scans_per_sess, [], 1:2);
    end

    fprintf(' %3.0f s\n', toc)

    % ------------------------------------------------------
    V = cell(1, N);
    fprintf('Mapping image volumes to memory: ');
    tic
    for i = 1:N
        fprintf('%3.0f ', i);
        V{i} = spm_vol(imgs{i});
    end

    fprintf(' %3.0f s\n', toc)




    % ------------------------------------------------------

    disp('Registering graphics callback with SPM registry')

    callback_handle = @(str, pos, reg, hReg) selavg_interactive_callback_wrapper(str, pos, reg, hReg);

    hSpmFig = spm_figure('GetWin', 'Graphics');

    hReg = uicontrol(hSpmFig, 'Style', 'Text', 'String', 'InteractiveViewer hReg', ...
        'Position', [100 200 100 025], 'Visible', 'Off', ...
        'FontName', 'Times', 'FontSize', 14, 'FontWeight', 'Bold', ...
        'HorizontalAlignment', 'Center');

    hReg = spm_XYZreg('InitReg', hReg, V{1}(1).mat, V{1}(1).dim(1:3)');
    spm_XYZreg('Add2Reg', hReg, 0, callback_handle);
    spm_orthviews('Register', hReg);


    %fh = findobj('Tag', 'Graphics');
    %set(fh, 'WindowButtonUpFcn', callback_handle);

    disp('Ready!')


    % inline

    function selavg_interactive_callback_wrapper(str, pos, reg, hReg)

        switch str
            case 'SetCoords'
                selavg_interactive_callback(pos, onsets, V, t, basepts, plotstes, S, scans_per_sess, I);

            otherwise
                disp('Unknown callback command from spm_XYZreg');
        end

    end

end





function selavg_interactive_callback(pos, onsets, V, t, basepts, plotstes, S, scans_per_sess, I)
    %% Done each time you click

    % %     %pos = spm_orthviews('Pos');
    % %     vox = (mm2voxel(pos', V{1}(1)))';
    % %     vox(4) = 1;
    % %
    % %     fprintf('Position: %3.0f %3.0f %3.0f\n', pos(1), pos(2), pos(3));
    % %     go = input('Plot selective averages for this voxel? (1 or 0) : ');

    % %     if go
    % % %         % Load data for this voxel
    % % %         % -----------------------------------------------------
    % % %         N = length(V);
    % % %
    % % %         fprintf('Loading voxel data for all subjects: ');
    % % %         for i = 1:N, fprintf('%3.0f ', i); braindata{i} = spm_get_data(V{i}, vox); end

    % Load data and get selective averages for each subject
    % -----------------------------------------------------
    [group_avgs, group_stes, subject_avgs, subject_stes, braindata] = selective_average_group( ...
        V, onsets, pos, 't', t, 'basepts', basepts, 'plotstes', plotstes, 'S', S, 'scans', scans_per_sess, 'I', I);


    if ~isempty(group_avgs)

        disp('-------------------------------------------------------');
        disp('Assigning group_avgs and stes in base workspace.');
        disp('(and subject_avgs and subject_stes)');
        disp('Assigning braindata to base workspace: Timeseries for each subject');
        disp('-------------------------------------------------------');
        disp(' ')

        assignin('base', 'group_avgs', group_avgs);
        assignin('base', 'group_stes', group_stes);
        assignin('base', 'subject_avgs', subject_avgs);
        assignin('base', 'subject_stes', subject_stes);
        assignin('base', 'braindata', braindata);

    end
    % %     end
end



function [SETUP, imgs, wh_is_image, name] = load_mediation_setup

    SETUP = [];
    imgs = [];

    fname = [pwd filesep 'mediation_SETUP.mat'];
    if exist(fname,'file')
        load(fname);

        % try to find names (single level)
        if exist('SETUP','var') && isfield(SETUP, 'M') && ischar(SETUP.M)
            imgs = SETUP.M;
            name = 'From mediation_SETUP SETUP.M';
            wh_is_image = 'M';

        elseif exist('SETUP','var') && isfield(SETUP, 'X') && ischar(SETUP.X)
            imgs = SETUP.X;
            name = 'From mediation_SETUP SETUP.X';
            wh_is_image = 'X';

        elseif exist('SETUP','var') && isfield(SETUP, 'Y') && ischar(SETUP.Y)
            imgs = SETUP.Y;
            name = 'From mediation_SETUP SETUP.Y';
            wh_is_image = 'Y';

        end

        % try to find names: multi-level
        if isfield(SETUP, 'data')
            switch SETUP.cmdstring
                case 'Search for mediators'
                    imgs = SETUP.data.M;
                    name = 'Multilevel, from SETUP.data.M';
                    wh_is_image = 'M';

                case 'Search for indirect influences'
                    imgs = SETUP.data.X;
                    name = 'Multilevel, from SETUP.data.X';
                    wh_is_image = 'X';

                case 'Search for mediated outcomes'
                    imgs = SETUP.data.Y;
                    name = 'Multilevel, from SETUP.data.Y';
                    wh_is_image = 'Y';

                otherwise
                    error('Unknown cmdstring: "%s".', cmdstring);
            end

            if ~iscell(imgs) || ~ischar(imgs{1})
                imgs = []; % invalid data here
            end



        end



        if isempty(imgs)
            fprintf(1,'Could not find image list.\n');
        end

    else
        fprintf(1,'Go to valid mediation directory with SETUP.mat to use interactive plotting.\n');
    end

end

