function plot(r, varargin)

t = r.test_on_the_fly;
ngroups = [];
res = [];
if isfield(struct(r), 'test_on_the_fly')
	t = r.test_on_the_fly;
	if isstruct(t)
		if isfield(t, 'ngroups'), ngroups = t.ngroups(:); end
		if isfield(t, 'result'), res = t.result; end
	end
end
if isempty(res), error('no results to plot: need to set test_on_the_fly.data subfield before training'), end
if isempty(ngroups), error('found no data for number-of-features in test_on_the_fly field'), end

err = extract(res, 'Y', 1);
plot(ngroups, err, varargin{:})
