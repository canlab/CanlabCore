function [cntrngmat, scaleFactor] = get_cntrng_mat(sid)
    uniq_sid = unique(sid);
    nsubj = length(uniq_sid);
    mats = cell(nsubj,1);
    scaleFactor = zeros(length(sid),1);
    for i = 1:nsubj
        this_sid = uniq_sid(i);
        sid_idx = (this_sid == sid);
        n = sum(sid_idx);
        mats{i} = eye(n) - 1/n*ones(n);
        scaleFactor(sid_idx) = sqrt(n);
    end
    cntrngmat = blkdiag(mats{:});
end