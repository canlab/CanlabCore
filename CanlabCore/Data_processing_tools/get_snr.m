function snr = get_snr(data)
% snr = get_snr(data)
%
% data is a matrix whos columns index voxels, and rows index subjects (or trials, etc.)
%
% Tor Wager

mystd = nanstd(data);
mystd(mystd == 0) = NaN;
snr = nanmean(data) ./ mystd;

return

