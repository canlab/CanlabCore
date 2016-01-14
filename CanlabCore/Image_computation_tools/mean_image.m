% :Usage:
% ::
%
%     OUT = mean_image(input file names, output file name, [weights], [string command args]);
%
% Creates a weighted mean of images based on weights for each image
%
% :Inputs:
%
%   **'normlike':**
%        do in style of mean_warped_image, for norm toolbox;
%        ignore other string arguments
%
%
%   **'reweight':**
%        weight by distance from median in 2nd round of averaging
%
%   **'plot':**
%        plot output
%
%   **'sharpen':**
%        segments and smooths within tissue classes
%
% :Examples:
% ::
%
%    w = filenames('w*T1.img', 'char')
%    w = mean_image(w, 'final_mean.img', [], 'normlike'); % trimmed, weighted best 50%
%
%    worig = filenames('hr*/wT1.img', 'char')           % no weights
%    w = mean_image(worig, 'orig_mean_noweights.img', []);
%
%    mean_image(DB.PP, 'Activation.img', DB.studyweight);
%
%    mean_image(P, 'mean_wmeanw1007.img', ones(size(P, 1), 1), 'sharpen', 'plot', 'reweight');
%
%    P = filenames('w*mr.img', 1);
%
%    % Enter empty inputs for default settings.
%    OUT.P = filenames
%    OUT.dist

function OUT = mean_image(VP, Pout, w, varargin)
    % --------------------------------------
    % inputs
    % --------------------------------------
    dosharp = 0;
    dodist = 0; 
    doplot = 0; 
    doreweight = 0; 
    donormlike = 0;
    
    dozscore = 0;

    for i = 1:length(varargin)
        if ischar(varargin{i})
            if strcmp(varargin{i}, 'sharpen'), dosharp = 1; end
            %if strcmp(varargin{i}, 'dist'), dodist = 1; end
            if strcmp(varargin{i}, 'plot'), doplot = 1; end
            if strcmp(varargin{i}, 'reweight'), doreweight = 1; end
            if strcmp(varargin{i}, 'normlike'), donormlike = 1; end
            
            if strcmp(varargin{i}, 'zscore'), dozscore = 1; end
        end
    end

    if isempty(VP)
        VP = spm_vol(spm_get(1, '*', 'Select images'));
    else
        if(iscellstr(VP)), 
            VP = char(VP); 
        end
        P = VP;
        VP = spm_vol(VP);
    end
    n = length(VP);

    if ~exist('Pout', 'var') || isempty(Pout), Pout = 'mean.nii'; end
    if ~exist('w', 'var') || isempty(w), w = ones(n, 1); end

    OUT.P = VP;

    % --------------------------------------
    % do mean in style of mean_warped_image
    % --------------------------------------
    % output is weights
    if donormlike
        OUT = do_mean_like_mean_warped_image(P, Pout);
        return
    end

    % --------------------------------------
    % make mean image
    % --------------------------------------

    status_string = sprintf('%3d', 0);
    fprintf(['Percent done: ' status_string]);

    
    switch spm('Ver')
        case 'SPM2'
            % spm_defaults is a script
            disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');
            
        case 'SPM5'
            % spm_defaults is a function
            spm_defaults()
        case 'SPM8'
            spm_defaults()
    end

    % initialize output
    V = VP(1);
    V.fname = Pout;
    V.n(1) = 1;    % write to first volume

    dat = zeros(V.dim(1:3));

    % read and weight

    for i = 1:n
        v = iimg_read_vols(VP(i));
        
        if dozscore
            v = center_scale(v);
        end
        
        dat = dat + v .* w(i);
        
        erase_string(status_string);
        status_string = sprintf('%3d', round(i/n*100));
        fprintf(status_string);
    end
    
    % divide by N
    dat = dat ./ n;
    
    fprintf('\n');

    % write output
    %strrep required to make file extension correct
    V.fname = strrep(V.fname, '.nii,1', '.nii');
    V.fname = strrep(V.fname, '.img,1', '.img');
    spm_write_vol(V, dat);
    
    if doplot
        z = round(V.dim(3)./2);
        f1 = figure;
        imagesc(dat(:,:,z));
        axis image;
        axis off;
        title('Mean');  drawnow
    end





    % --------------------------------------
    % Reweight
    % --------------------------------------

    if doreweight
        fprintf(1, ' Reweighting.000');

        for i = 1:n
            fitness(i) = norm_eval_fitness(tempdat, normdat, 0);
        end

        imgMedian = zeros(n, 1);
        imgMean = imgMedian; 
        imgStd = imgMedian;

        medianDev = zeros(n, 5);
        newdat = zeros(V.dim(1:3));

        % subtract median from mean image
        [dat, datMedian] = subtract_median(dat);

        % get medians and subtract
        for i = 1:n

            fprintf(1, '\b\b\b%03d', i);

            v = spm_read_vols(spm_vol(deblank(P(i,:))));

            [v, imgMedian(i)] = subtract_median(v);

            %imgMean(i) = mean(vi);
            %imgStd(i) = std(vi);

            % distance from median-centered old mean image
            vi = v - dat;  vi = abs(vi(:));
            vi(isnan(vi) | vi==0) = [];
            if isempty(vi), warning('Image is EXACTLY the median image.'); end

            % get median deviation, and 5th and 95th prctiles
            medianDev(i,:) = prctile(vi, [50 5 20 80 95]);

            % get new weights based on distance from mean image
            w(i) = datMedian./medianDev(i, 1);

            newdat = newdat + v .* w(i);

        end

        newdat = newdat .* sum(w);          % normalize to keep in same scale

        whnull = isnan(newdat) | newdat==0;
        %newdat = newdat + mean(imgMedian);  % add median back in to average
        newdat(whnull) = 0;

        dat = newdat;
        OUT.medianDev = medianDev;
        OUT.prctiles = [50 5 20 80 95];

        if doplot
            tor_fig; 
            plot(medianDev(:,1), 'o-');
            ylabel('Med. Abs. Dev from median image');
            xlabel('Image number'); grid on;

            figure(f1); subplot(4, 1, 2);
            imagesc(dat(:,:,z));
            axis image;
            axis off;
            title('Reweighted');drawnow
        end

        % write output
        spm_write_vol(V, dat);
    end


    % --------------------------------------
    % Sharpen
    % --------------------------------------

    if dosharp
        fprintf(1, ' Sharpening: segmenting.');

        defaults.segment.estimate.reg = .1;% was .01;
        VO = spm_segment(Pout, which('T1.mnc'), defaults.segment);

        spm_write_vol(VO(1), VO(1).dat);
        spm_write_vol(VO(2), VO(2).dat);
        spm_write_vol(VO(3), VO(3).dat);

        fprintf(1, ' Lightening white matter.');

        thr = prctile(VO(2).dat(:), 90);
        whiteMask = double(VO(2).dat > thr);

        mystd = .2 .* nanstd(dat(:));
        dat = dat + whiteMask .* mystd;

        % write output
        name = [Pout(1:end-4) '_sharp.img'];
        V.fname = name;
        spm_write_vol(V, dat);

        % make brain mask
        brainMask = double(VO(1).dat>0 | VO(2).dat>0);

        name = [Pout(1:end-4) '_brain.img'];
        V.fname = name;
        spm_write_vol(V, brainMask);

        ovl = dat .* brainMask;
        name = [Pout(1:end-4) '_overlay.img'];
        V.fname = name;
        spm_write_vol(V, ovl);

        if doplot
            figure(f1); subplot(4, 1, 3);
            imagesc(dat(:,:,z));  
            axis off;
            axis image;
            title('Sharpened');drawnow

            subplot(4, 1, 4);
            imagesc(ovl(:,:,z));
            axis off; 
            axis image;
            title('Overlay');drawnow
        end
    end
end



function [dat, datMedian] = subtract_median(dat)
    vi = dat(:); 
    vi(isnan(vi) | vi==0) = [];
    datMedian = median(vi);

    whnull = isnan(dat) | dat==0;

    dat = dat - datMedian;
    dat(whnull) = 0;
end



function v = center_scale(v)
    vi = v(:); 
    vi(isnan(vi) | vi==0) = [];
    v = (v - mean(vi)) ./ std(vi);
end




function weights = do_mean_like_mean_warped_image(P, outname)
    spm_defaults;
    defaults.analyze.flip = 0;
    [volTemplate, normdat] = iimg_read_img(P);

    t1 = clock;
    n = size(P, 1);

    % Provisional mean
    % ------------------------------------
    fprintf(1, 'Initial. ');

    normdat(isnan(normdat)) = 0;

    % convert in-mask voxels to mean = gm, so averaging is not distorted by
    % overall intensity differences among images
    % and fitness weights determine relative contribution of images
    means = mean(normdat);      % mean of each image
    gm = mean(means);           % global mean
    scalef = gm ./ means;       % scaling factor to norm each image to grand mean

    for i = 1:n
        normdat(:,i) = normdat(:,i) .* scalef(i);
    end

    meandat = mean(normdat, 2);
    vox = size(meandat, 1);


    fprintf(1, 'Weighted mean of best 50%%. ');
    % Weights based on closeness to mean
    % ------------------------------------
    for i = 1:n
        fitness(i) = norm_eval_fitness(meandat, normdat(:,i), 0);
    end

    mfit = nanmedian(fitness);
    weights = fitness;
    weights(weights < mfit) = 0;

    weights = weights ./ sum(weights);

    if all(weights == 0)
        warning('All weights are zero!')
        weights = ones(1, n);
    end

    % Weighted mean
    % ------------------------------------
    meandat = zeros(vox, 1);
    for i = 1:n

        if weights(i) > 0   % just to save time
            meandat = meandat + weights(i) .* normdat(:,i);
        end

    end

    fprintf(1, 'Done: Writing\n');
    % Write output image
    % This will be the new normalization template
    % ------------------------------------
    meandata = iimg_reconstruct_3dvol(meandat, volTemplate, 'outname', outname);

    fprintf(1, 'Mean warped image completed: %3.0f s\n', etime(clock, t1));

    fprintf(1, '-----------------------------\n');
end
