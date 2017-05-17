function sl = timeseries_extract_slice(V, sliceno, orientation)
% For a given set of image names or memory mapped volumes (V)
% extracts data from slice # sliceno and returns an X x Y x time
% matrix of data.
%
% :Usage:
% ::
%
%     function sl = timeseries_extract_slice(V,sliceno)
%

    global defaults
    
    % defaults
    switch spm('Ver')
        case 'SPM2'
            % spm_defaults is a script
            disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');

        case {'SPM5', 'SPM8', 'SPM12'}
            % spm_defaults is a function
            spm_defaults()

        otherwise
            % unknown SPM
            disp('Unknown version of SPM!');
            spm_defaults()
    end

    
    if ischar(V)
        V = spm_vol(V);
    end

    if(~exist('orientation', 'var') || isempty(orientation))
        orientation = 'axial';
    end

    switch(orientation)
        case 'axial'
            get_slice = @get_ax_slice;
        case 'sagittal'
            get_slice = @get_sag_slice;
        case 'coronal'
            get_slice = @get_cor_slice;
        otherwise
            error('Unknown orientation: %s\n', orientation);
    end

    for i = 1:length(V)
        sl(:,:,i) = get_slice(V(i), sliceno);
    end
end
