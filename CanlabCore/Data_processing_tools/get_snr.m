function snr = get_snr(data)
% Data is a matrix whos columns index voxels, and rows index subjects (or trials, etc.)
%
% :Usage:
% ::
%
%     snr = get_snr(data)
%
% ..
%    Tor Wager
% ..

mystd = nanstd(data);
mystd(mystd == 0) = NaN;
snr = nanmean(data) ./ mystd;

return

