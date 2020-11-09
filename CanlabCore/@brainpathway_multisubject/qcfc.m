% Compute QCFC correlations and plot them (Parkes et al 2018, Ciric et al
% 2018, Power et al 2012). Computes edgewise correlations between
% connectivity estimates (stored in bs.connectivity.regions.r) and a QC metric (typically head
% motion) passed in by user
function qcfc(bs, qc_metric)

    flat_conn_mat = flatten_conn_matrices(bs);
    corrs = corr(flat_conn_mat, qc_metric);
    figure; histogram(corrs) 
    fprintf('Mean (SD) corr = %3.2f (%3.2f)\n', mean(corrs), std(corrs));
end