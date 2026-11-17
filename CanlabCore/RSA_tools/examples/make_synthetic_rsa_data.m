function dat = make_synthetic_rsa_data(varargin)
% make_synthetic_rsa_data  Generate a planted-structure dataset for the RSA examples.
%
% Builds an fmri_data object with metadata columns {sub, sesno, condition,
% bodysite} and a planted representational structure: images sharing a
% condition are similar (strong), images sharing a bodysite are similar
% (weaker), plus per-subject and per-session offsets and noise.
%
% Usage
%   dat = make_synthetic_rsa_data();                  % defaults
%   dat = make_synthetic_rsa_data('n_sub', 9, 'n_ses', 5);
%
% Optional name-value: n_vox, n_sub, n_ses, cond_weight, bs_weight, noise, seed

p = inputParser;
p.addParameter('n_vox', 200);
p.addParameter('n_sub', 9);
p.addParameter('n_ses', 5);
p.addParameter('cond_weight', 0.7);
p.addParameter('bs_weight', 0.5);
p.addParameter('noise', 0.4);
p.addParameter('seed', 42);
p.parse(varargin{:});
o = p.Results;

rng(o.seed);
conditions = {'hot','warm','imagine'};
bodysites  = {'leftface','rightface','leftarm','rightarm','leftleg','rightleg','chest','abdomen'};
n_cond = numel(conditions); n_bs = numel(bodysites);

cond_sig = randn(n_cond, o.n_vox);
bs_sig   = randn(n_bs,  o.n_vox);
P = zeros(n_cond*n_bs, o.n_vox); row = 0;
for c = 1:n_cond
    for b = 1:n_bs
        row = row + 1;
        P(row, :) = o.cond_weight*cond_sig(c,:) + o.bs_weight*bs_sig(b,:) + 0.05*randn(1,o.n_vox);
    end
end

X = []; sub_v = {}; ses_v = []; cond_v = {}; bs_v = {}; idx = 0;
for s = 1:o.n_sub
    sub_offset = 0.15*randn(1, o.n_vox);
    for se = 1:o.n_ses
        ses_offset = 0.1*randn(1, o.n_vox);
        for c = 1:n_cond
            for b = 1:n_bs
                idx = idx + 1;
                X(:, idx) = P((c-1)*n_bs+b, :)' + sub_offset' + ses_offset' + o.noise*randn(o.n_vox, 1); %#ok<AGROW>
                sub_v{idx,1}  = sprintf('sub-%02d', s); %#ok<AGROW>
                ses_v(idx,1)  = se; %#ok<AGROW>
                cond_v{idx,1} = conditions{c}; %#ok<AGROW>
                bs_v{idx,1}   = bodysites{b}; %#ok<AGROW>
            end
        end
    end
end

dat = fmri_data;
dat.dat = X;
dat.metadata_table = table(sub_v, ses_v, cond_v, bs_v, ...
    'VariableNames', {'sub','sesno','condition','bodysite'});
end
