function mytable = glm_table(stat, nms, varargin)
% mytable = glm_table(stat,nms, [betas], [any arg to suppress print])
%
% simple utility for printing a table from a glmfit output, intercept first
% do not include intercept name in nms (names) cell array.
%
% Can work for glmfit output or robustfit output.  If robustfit, input betas
% as 3rd argument

if nargin > 3
    doprint = false;
else
    doprint = true;
end

if nargin < 2 || isempty(nms)
    for i = 1:length(stat.beta)-1
        nms{i} = ['V' num2str(i)];
    end
end
nms = [{'Intercept'} nms];

if ~isfield(stat, 'beta')
    stat.beta = varargin{1};
end

if isrow(nms), nms = nms'; end

Name = nms;
Beta = stat.beta;
SE = stat.se;
t = stat.t;
p = stat.p;

sigstr = cell(size(nms));
sigstr(:) = {' '};
sigstr(p < .1) = {'+'};
sigstr(p < .05) = {'*'};
sigstr(p < .01) = {'**'};
sigstr(p < .001) = {'***'};
Sig = sigstr;

mytable = table(Name, Beta, SE, t, p, Sig);

if doprint
    disp(mytable);
end

% legacy
% mytable = [];
% 
% fprintf(1,'%s\t%s\t%s\t%s\t%s\t\n', ...
%     'Name','Beta','SE','t','p');
% 
% 
% for i = 1:length(stat.beta)
%     
%     fprintf(1,'%s\t%3.3f\t%3.3f\t%3.3f\t%3.4f\t\n', ...
%         nms{i},stat.beta(i),stat.se(i),stat.t(i),stat.p(i))
%     
% end

end % function


