function iimg_weighted_ttest(image_names,varargin)
% fast weighted test
%
% :Examples:
% ::
%
%    iimg_weighted_ttest(image_names,'mask',maskname)
%
% ..
%    Tor Wager
% ..

% ..
%    Set up arguments
% ..

maskname = deblank(image_names(1,:)); % defaults
w = []; varY = []; bcon = [];

% inputs
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch lower(arg)
            case 'mask', maskname = varargin{i+1};
            case 'w', w = varargin{i+1};
            case 'btwn', bcon = contrast_code(varargin{i+1});
            case 'varY', varY = varargin{i+1};
        end
    end
end


% --------------------------------------
% * Load data
% --------------------------------------

str = display_string('Reading data.');

[dat,maskInfo] = iimg_get_data(maskname,image_names);

erase_string(str);


% --------------------------------------
% * Computation
% --------------------------------------

str = display_string('Statistics.');

[means,stats] = weighted_reg(dat,'w',w,'varY',varY,'btwn',bcon,'uni');

erase_string(str);


% --------------------------------------
% * Write output
% --------------------------------------
descrip = sprintf('dfe = %3.2f',stats.dfe);
meanvol = iimg_reconstruct_3dvol(means.Ymean',maskInfo, 'outname','mean.img', 'descrip',['mean ' descrip]);
tvol = iimg_reconstruct_3dvol(stats.t',maskInfo, 'outname','t_intercept.img', 'descrip',['t-image of intercept ' descrip]);
pvol = iimg_reconstruct_3dvol(stats.p',maskInfo, 'outname','p_intercept.img', 'descrip',['p-image of intercept ' descrip]);


return





function str = display_string(str)
str = sprintf(str); fprintf(1,'%s',str);
return


function erase_string(str)

len = length(str);
str2 = repmat('\b',1,len);

fprintf(1,str2);

return

