function HRF_struct = makeHRFstruct(HRF_OBJ, CondNames, varargin)

    isSigFound = false;
    isAtlasFound = false;
    isRegsFound = false;
    r = [];
    s = [];

    if isa(HRF_OBJ, 'fmri_data')
        HRF_OBJ={HRF_OBJ};
    end

    % Initialize the parallel pool if it's not already running
    % if isempty(gcp('nocreate'))
    %     parpool;
    % end

    for k = 1:length(varargin)
        if strcmpi(varargin{k}, 'atlas')
            if isa(varargin{k+1}, 'atlas')
                at=varargin{k+1};
            else
                error('Passed in atlas not an atlas.');
            end

            isAtlasFound = true;
        end

        if strcmpi(varargin{k}, 'regions')
            if isAtlasFound
                if ischar(varargin{k+1})
                    r={varargin{k+1}};
                    isRegsFound = true;
                elseif iscell(varargin{k+1})
                    r=varargin{k+1};
                    isRegsFound = true;
                end
            else
                error('Cannot extract HRF of regions without atlas')
            end
        end

        if isAtlasFound && ~isRegsFound
            % Extract all atlas regions.
            r=at.labels;
        end

        if strcmpi(varargin{k}, 'sig') % Expected input: 'nps', 'siips', or 'all'
            if ischar(varargin{k+1})
                s={varargin{k+1}};
                isSigFound = true;
            elseif iscell(varargin{k+1})
                s=varargin{k+1};
                isSigFound = true;
            else
                error(['Input argument for sig is unknown.']);
            end

        end
    end

    % ROIs will consist of regions and signatures
    rois = [r, s];

    CondNames=matlab.lang.makeValidName(CondNames);



    HRF=extractHRF(HRF_OBJ, CondNames, varargin{:});

    HRF_struct=struct;

    if isAtlasFound
        HRF_struct.atlas=at;
    end

    HRF_struct.region=rois;
    HRF_struct.CondNames=CondNames;
    HRF_struct.fit=HRF;

end