function Vo = percent_sig_image(imgs, baseimg, outname)
% Creates a percent signal change image saved in outname by dividing each
% image by baseimg
%
% :Usage:
% ::
%
%     function Vo = percent_sig_image(imgs, baseimg, outname)
%
% ..
%    by Tor Wager
% ..



V = spm_vol(imgs);
vols = spm_read_vols(V);

Vb = spm_vol(baseimg);
basei = spm_read_vols(Vb);

basei(basei <= 10*eps) = NaN;

if size(basei, 4) > 1
    error('Baseline image cannot be 4-D')
end

for i = 1:size(vols, 4)  % for each volume if 4-D

    vols(:, :, :, i) = vols(:, :, :, i) ./ basei;

end


% -------------------------
% write
% -------------------------
fprintf('Writing:')
for i = 1:size(vols, 4)
    [d,f,e] = fileparts(V(i).fname);
    
    Q = fullfile(d, outname);
    fprintf('  %s\n', Q)
    Vo = V(i);
    Vo.fname = Q;
    Vo.descrip = ['% Change from ' Vb.fname];

    spm_write_vol(Vo, vols(:, :, :, i));
end
