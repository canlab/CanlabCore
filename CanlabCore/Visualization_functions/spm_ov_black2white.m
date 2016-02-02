function spm_ov_black2white(varargin)
% Plugin for spm_orthviews to change the black background and dark edges to
% a white background with softer gray edges, for pub. quality figures
%
% if st (a global variable) st.plugins has 'black2white' added, and the
% st.vols{i} struct has a field called 'black2white', then this will be
% called.
%
% To initialize, add this code to the calling function:
% ::
%
%    st.vols{1}.black2white = 1;
%    bwexist = strfind(st.plugins, 'black2white')
%    bwexist = any(cat(2, bwexist{:}))
%    if ~bwexist
%        st.plugins{end+1} = 'black2white';
%    end
%
% ..
%    Sept 2010, Tor Wager
% ..


    % soften edges; 0 is no softening, a range is more softening
    softenvals = [.01:.01:.35]; % higher ending values of k will be more 'softening'; lower = less
    softenvals = [.05]; 

    global st;
    if isempty(st)
        error('movie: This routine can only be called as a plugin for spm_orthviews!');
    end;

    if nargin < 2
        error('movie: Wrong number of arguments. Usage: spm_orthviews(''movie'', cmd, volhandle, varargin)');
    end;

    cmd = lower(varargin{1});
    i = varargin{2};  % volume handle, set by default in SPM


    if ~isstruct(st.vols{i}), return, end % not a valid volume with data...

    for a = 1:3 % each axis

        cd = get(st.vols{i}.ax{a}.d, 'CData');

        % turn black to white

        isblack = double(all(cd == 0, 3));

        if size(cd, 3) == 3 % is true color matrix
            cd = cd + repmat(isblack, [1 1 3]);
        else
            cd = cd + isblack;
        end

        set(st.vols{i}.ax{a}.d, 'CData', cd);

        % lighten the dark edges

%         for j = 1:3
%             cdslice = cd(:,:,j);
%         isdark = all(cdslice < softenvals, 3);
%         cdslice(isdark) = 1 - (tanh(1.*(1-cdslice(isdark)))./1);
%         cd(:, :, j) = cdslice;
%         end
        
        %cd2 = smooth3(cd, 'gaussian', 3);
        
        % higher ending values of k will be more 'softening'; lower =
        % more dark edges
        for k = softenvals % for each of these darkness values
            w = 1 - k; % weights, for wtd average of orig and inverted, 1 is all inverted
            % (bilinear weighting)

            isdark = all(cd < k, 3);
            for j = 1:size(cd, 3)

                %invert
                colorslice = cd(:, :, j);
                colorslice(isdark) = w .* (1 - colorslice(isdark)) + (1 - w) .* colorslice(isdark);

                cd(:, :, j) = colorslice;
            end

        end % the weight value loop

        set(st.vols{i}.ax{a}.d, 'CData', cd);


    end % axis


    % set the colormap, in case of color-mapped image (old-style blobs)
    % -------------------------------------------------------------------------
    fh = findobj('Type', 'Figure', 'Tag', 'Graphics'); % spm fig
    cm = get(fh, 'Colormap');

    for k = softenvals % for each of these darkness values
        w = 1 - k; % weights, for wtd average of orig and inverted, 1 is all inverted
        % (bilinear weighting)

        isdark = all(cm < k, 2);
        for j = 1:size(cm, 2)

            %invert
            colorslice = cm(:, j);
            colorslice(isdark) = w .* (1 - colorslice(isdark)) + (1 - w) .* colorslice(isdark);

            cm(:, j) = colorslice;
        end

    end % the weight value loop

    set(fh, 'Colormap', cm);


end % function

