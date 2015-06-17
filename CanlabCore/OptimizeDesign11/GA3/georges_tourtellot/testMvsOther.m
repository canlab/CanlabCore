function SIM = testMvsOther(M,SIM,varargin)
% function SIM = testMvsOther(M,SIM,varargin)
% 	M is output from GA.
%	SIM is the output of a random design SIM;
%	more cells can be added on by repeated passing
%	in of the same SIM.
%
%	SIM(x).params contains the ga/simulation params
%	for simulation x.
%
%	SIM(1).params contains the original GA M params.
%	
%	Arguments are fields in SIM.params to be modified,
%	followed by the value.
%	
% 	Example: 
%	SIM = testMvsOther(M,SIM,'ISI',4,'TR',1.5);
%
%	SIM = testMvsOther(M,[],'ISI',4,'TR',1.5);
%	for new simulation set.
%
%	Made 5/24/01 by Tor Wager

if isempty(SIM),
	firsttime = 1;
	SIM(1).params = M.ga; 
	SIM(1).name = 'GA results'; 
	SIM(1).params.type = 'ga';
	SIM(1).params.stimlist = M.stimlist;
	SIM(1).results = [];

    %gst added lines below
SIM(1).params.noise_var = 1		% noise variance
SIM(1).params.niterations = 25  %gst   was 1000		% number of iterations in sampling distribution
SIM(1).params.beta = 1;			% slope of true (simulated) effect

    %gst added below
   	SIM(2).params = M.ga; 
	SIM(2).name = 'GA results'; 
	SIM(2).params.type = 'rnd';
%	SIM(2).params.stimlist = M.stimlist;  %for testing only, gst
	SIM(2).results = [];

    %gst added lines below
SIM(2).params.noise_var = 1		% noise variance
SIM(2).params.niterations = 25		% number of iterations in sampling distribution
SIM(2).params.beta = 1;			% slope of true (simulated) effect
    
    %gst added below  %gst says remove lines below
%     SIM(3).params = M.ga; 
% 	SIM(3).name = 'GA results'; 
% 	SIM(3).params.type = 'ga';
% 	SIM(3).results = [];
% 
%     load('mymseq_stimlist.mat');
%     stimlist=a;
%     SIM(3).params.stimlist = stimlist
%     
% %gst added lines below
% SIM(3).params.noise_var = 1		% noise variance
% SIM(3).params.niterations = 25		% number of iterations in sampling distribution
% SIM(3).params.beta = 1;			% slope of true (simulated) effect

        
%gst says to put these lines back    
% else
% 	firsttime = 0;
end
simindex = size(SIM,2) + 1;
% SIM(simindex).params = M.ga; 
% SIM(simindex).name = 'Random simulation'; 
% SIM(simindex).results = [];
 

% ----------------------------------------------------------------
% * set up defaults, fixed parameters, and from-GUI parameters
% ----------------------------------------------------------------
HRF = spm_hrf(.1);
HRF = HRF/ max(HRF);
format compact
for i = 1:size(SIM,2)
%    if ~isfield(SIM(i).params,'powerTvoxels'),SIM(i).params.powerTvoxels =    30000;,end   %gst
	if ~isfield(SIM(i).params,'powerTvoxels'),SIM(i).params.powerTvoxels = 1;,end  %gst

    if ~isfield(SIM(i).params,'HRF_estlen'),SIM(i).params.HRF_estlen = 12;,end  %gst (FEATURE)
    if ~isfield(SIM(i).params,'NULL_flag'),SIM(i).params.NULL_flag = 0;,end  %gst (FEATURE)
    
    if ~isfield(SIM(i).params,'powerTdf'),SIM(i).params.powerTdf = 19;,end
	if ~isfield(SIM(i).params,'HRF'),SIM(i).params.HRF = HRF;,end
	if ~isfield(SIM(i).params,'type'),SIM(i).params.type = 'rnd';,end
	if ~isfield(SIM(i).params,'maxrestthresh'),SIM(i).params.maxrestthresh = 2;,end
end


% H = findobj('Tag', 'E12')  
% if firsttime,SIM(1).params.noise_var = str2num(get(H(1), 'String'));,end
% SIM(simindex).params.noise_var = str2num(get(H(1), 'String')); 
% 
% H = findobj('Tag', 'E14');
% if firsttime,SIM(1).params.niterations = str2num(get(H(1), 'String'));,end  
% SIM(simindex).params.niterations = str2num(get(H(1), 'String')); 
% 
% H = findobj('Tag', 'E20');
% if firsttime,SIM(1).params.beta = str2num(get(H(1), 'String'));,end   
% SIM(simindex).params.beta = str2num(get(H(1), 'String')); 


% ----------------------------------------------------------------
% * set up user arguments - these override GUI values and defaults
% ----------------------------------------------------------------

for i = 1:nargin-2
    if isstr(varargin{i})
        for jj=1:size(SIM,2) %gst added this loop, so that all params put on command line are set for SIM(jj) for all jj
            %all jj's below were simindex's

            str = (['SIM(jj).params.' varargin{i} ' = [' num2str(varargin{i+1}) '];']);
			eval(str);
			eval(['test = SIM(jj).params.' varargin{i} ';'])
			disp(['	Cmd line entered: SIM(' num2str(jj) ').params.' varargin{i} ' = ' num2str(test)])
    	end
    end
end



% ----------------------------------------------------------------
% * run the simulation for each empty structure in SIM
% ----------------------------------------------------------------

for i = 1:size(SIM,2)  %gst
    if isempty(SIM(i).results)
		disp(['Running simulation on SIM(' num2str(i) ')'])
        SIM(i).results = runsim(SIM(i).params);
	end

	% ----------------------------------------------------------------
	% * save sim
	% ----------------------------------------------------------------
	a = clock;
	b = [num2str(a(2)) '_' num2str(a(3)) '_' num2str(a(4))];
	eval(['save SIM_' b ' SIM'])   


end	% loop thru sims


return




