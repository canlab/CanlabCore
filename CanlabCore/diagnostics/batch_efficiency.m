function batch_efficiency(dwcard)
% Start in directory above individual model/results directories
%
% :Usage:
% ::
%
%     function batch_efficiency(dwcard)
% 
% :Inputs:
%
%   **dwcard:**
%        is a wildcard for directories to probe, e.g., 'subject*'
%
% ..
%    Tor Wager, 10/8/02
% ..

d = dir(dwcard);

mypwd = pwd;
eval(['diary ' mypwd filesep 'model_efficiency_stats.out'])
 
for i = 1:length(d) 

	try
		spm_efficiency([d(i).name filesep 'SPM.mat']);
	catch
		disp(['Problem with subject ' d(i).name])
	end

end

diary off

return
