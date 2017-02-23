function glm_table(stat, nms, varargin)
% glm_table(stat,nms, [betas])
%
% simple utility for printing a table from a glmfit output, intercept first
% do not include intercept name in nms (names) cell array.
%
% Can work for glmfit output or robustfit output.  If robustfit, input betas
% as 3rd argument

if nargin < 2 || isempty(nms)
    for i = 1:length(stat.beta)-1
        nms{i} = ['V' num2str(i)];
    end
end
nms = [{'Intercept'} nms];

if ~isfield(stat, 'beta')
    stat.beta = varargin{1};
end

fprintf(1,'%s\t%s\t%s\t%s\t%s\t\n', ...
    'Name','Beta','SE','t','p');


for i = 1:length(stat.beta)
    
    fprintf(1,'%s\t%3.3f\t%3.3f\t%3.3f\t%3.4f\t\n', ...
        nms{i},stat.beta(i),stat.se(i),stat.t(i),stat.p(i))
    
end

end % function


