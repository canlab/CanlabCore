function [sfmL,sfm,names] = spm2delta(SPM,varargin)
% [d_hires,d_atTR,names] = spm2delta(SPM,varargin)
% or 
% [onsets,delta_atTR,names] = spm2delta(SPM,varargin)
% 
% Example: for SPM99
% [dL,d] = spm2delta([],Sess);

if length(varargin) > 0, 
    Sess = varargin{1};,
    SPM.Sess = Sess;
else
    % SPM2
    Sess = SPM.Sess;
end

if iscell(Sess)
    
    ons = {};
    
    for i = 1:length(Sess)
    
        [sfc,sfm{i}] = downsample_delta(Sess{i}.sf,16);
    
        for j = 1:size(sfm{i},2)                % save onsets
            ons{end+1} = find(sfm{i}(:,j)) - 1; % in TRs, starting with 0
        end
        
        %sfmL{i} = cat(2,Sess{i}.sf{:});        % delta hires
    end

    sfm = cat(1,sfm{:});
    sfmL = ons;         % save onsets
    %sfmL = cat(1,sfmL{:}); % save delta_hires
    names = SPM.Sess{1}.name;
else

    % new spm2 way?
    ons = [];
    wh = [];, for i = 1:length(Sess), if isempty(Sess(i).U), wh(i) = 1;,end,end
        Sess(find(wh)) = [];
        
    delta = cell(1,length(Sess(1).U));
    
    for i = 1:length(Sess)
        sessons = [];
        for j = 1:length(Sess(i).U)
            ons{end+1} = Sess(i).U(j).ons;      % all onsets
            sessons{end+1} = Sess(i).U(j).ons;  % onsets this session
        end
        
        % add delta functions for this session to overall
        [X,sessdelt] = onsets2delta(sessons,1,size(SPM.xX.X,1) ./ length(Sess));
        for k = 1:length(sessdelt)
            
            % just in case some sessions do not have some
                    % regressors
            if k > length(delta)
            	delta{k} = zeros(size(delta{1}));
            end
                    
            try
                delta{k} = cat(1,[delta{k}; sessdelt{k}]);
            catch
                 delta{k} = cat(1,[delta{k}; sessdelt{k}']);
            end
        end
    end
    sfmL = ons;
    sfm = delta;
    
    names = cat(2,SPM.Sess(1).U(:).name);
end
            
    

return