function obj = inv_fisher_z(obj)
% inv_fisher_z  Apply tanh to .dat, the inverse of fisher_z.

if numel(obj) > 1, obj = arrayfun(@inv_fisher_z, obj); return; end

obj.dat = tanh(obj.dat);
obj.history{end+1} = sprintf('%s: inv_fisher_z (tanh)', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

end
