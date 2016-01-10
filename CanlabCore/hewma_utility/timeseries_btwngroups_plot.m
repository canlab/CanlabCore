function h = timeseries_btwngroups_plot(dat,cov,varargin)
% :Usage:
% ::
%
%     h = timeseries_btwngroups_plot(dat,cov,[baseperiod],[dedetrend],[x],[legend text])
%
% :Inputs:
%
%   **imdat:**
%        subjects x time matrix of estimates
%
%   **baseperiod:**
%        integer; removes mean of 1:baseperiod for each subject
%
%   **cov:**
%        subjects x 1 contrast vector (individual/group differences)
%
%   **dodetrend:**
%        linear detrending of each subject's estimates
%
%   **legend text:**
%        {'name1' 'name2' ...}
%
% :Output:
%
%   **h:**
%        line handles

x = 1:size(dat,2);

baseperiod = []; dodetrend = 0; legtxt = [];

if length(varargin) > 0, baseperiod = varargin{1}; end
if length(varargin) > 1, dodetrend = varargin{2}; end
if length(varargin) > 2 && ~isempty(varargin{3}), x = varargin{3}; end
if length(varargin) > 3, legtxt = varargin{4}; end

if dodetrend
    dat = linear_detrending(dat);
end

if ~isempty(baseperiod)
    m = mean(dat(:,1:baseperiod),2);
    dat = dat - repmat(m,1,size(dat,2));
end

%cov = mediansplit(cov);

d1 = dat(cov<0,:);
m1 = nanmean(d1);
sterr = ste(d1);
d2 = dat(cov>0,:);
m2 = nanmean(d2);
sterr2 = ste(d2);
cla; hold on;
color = {'b' 'r'};

fill_around_line(m1,sterr,color{1},x);
fill_around_line(m2,sterr2,color{2},x);

h(1) = plot(x,m1,color{1},'LineWidth',2);
h(2) = plot(x,m2,color{2},'LineWidth',2);

% plot significant timepoints with black bar up top
v1 = var(d1);
v2 = var(d2);
vpool = v1 + v2;
n = size(dat,1);
vste = sqrt(vpool) ./ sqrt(n);
df = n - 2;
vt = (m1 - m2) ./ vste;
vp = 2 * (1 - tcdf(abs(vt),df));  % 2-tailed

yval = max([m1 m2]) + .12 * max([m1 m2]);
y = NaN * zeros(size(m1));
y(vp < .05) = yval;
plot(x,y,'k','LineWidth',3);

if isempty(legtxt)
    legend(h,{'Low' 'High'});
else
    legend(h,legtxt);
end

return

