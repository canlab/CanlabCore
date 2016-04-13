function [int_dat,mask_vol,outname] = iimg_intersection(varargin)
% Make a mask of the intersection of n image files (pos or neg values)
% If 'name' followed by an image file name is entered, writes output to
% that image
%
% :Usage:
% ::
%
%     [int_dat,mask_vol,outname] = iimg_intersection(name1, name2, etc.)
%     [int_dat,mask_vol,outname] = iimg_intersection(dat1, dat2, etc.)
%
% If special string 'posneg' is entered, separates first two images into
% combinations of pos/neg values in each image.
% The order returned in columns of int_dat is pospos, posneg, negpos, and
% negneg
%
% fastest if one output requested.
%


% ..
%    setup input read images and special strings
% ..

dat = [];   % data
doposneg = 0;  % separate results for pos and neg combos

if nargin == 0
    imname = spm_get(1);
    [volInfo,dat] = iimg_read_img(imname);
else
    for i = 1:length(varargin)
        v = varargin{i};
        if isstruct(v), volInfo = v;

            %  put control strings (special strings) here
        elseif ischar(v) && strcmp(v,'posneg')
            doposneg = 1;   % separate positive and negative overlap

        elseif ischar(v) && strcmp(v,'name')
            outname = varargin{i+1};

        elseif i > 1 && ischar(varargin{i-1}) && strcmp(varargin{i-1},'name')
            % do nothing; this is the output name, ignore it

        else  % treat this input as an input image file or data vector
            
            if isempty(dat)  %~exist('volInfo','var') % first image, use as mask space
                imname = v;
                [volInfo,d] = iimg_read_img(imname, 2);
                
            else
                d = scn_map_image(v, imname);
                d = d(:);
            end

            if size(dat,1) ~= 0 && size(d,1) ~= size(dat,1)
                error('Not all images are the same size.');
            else
                dat = [dat d];
            end
        end
    end
end

if ~exist('volInfo','var') || isempty(volInfo)
    error('You must enter at least one image file name or enter a volInfo structure as input.')
end

% ------------------------------------------------------
% setup input read images and special strings
% ------------------------------------------------------
if doposneg
    fprintf(1,'Separating combinations of pos and neg values for first two input images.\n')
    int_dat = dat(:,1) > 0 & dat(:,2) > 0;  % pos pos
    int_dat = [int_dat dat(:,1) > 0 & dat(:,2) < 0];  % pos neg
    int_dat = [int_dat dat(:,1) < 0 & dat(:,2) > 0];  % neg pos
    int_dat = [int_dat dat(:,1) < 0 & dat(:,2) < 0];  % neg neg
else
    % 'normal' way; abs value
    int_dat = all(abs(dat),2);
end

if nargout >= 2 && ~exist('outname','var')
    outname = 'iimg_intersection_output_image.img';
    
    if doposneg
        for i = 1:4
            mask_vol(:,:,:,i) = iimg_reconstruct_3dvol(int_dat(:,i),volInfo);
        end
    else
        mask_vol = iimg_reconstruct_3dvol(int_dat,volInfo);
    end

    % ------------------------------------------------------
    % write output images
    % ------------------------------------------------------

elseif exist('outname','var')

    %name = input('Enter output filename: ','s');
    if doposneg
        name1 = append_to_name(outname,'pospos');
        mask_vol(:,:,:,1) = iimg_reconstruct_3dvol(int_dat(:,1),volInfo,'outname',name1,'descrip','Positive values for both imgs');

        name1 = append_to_name(outname,'posneg');
        mask_vol(:,:,:,2) = iimg_reconstruct_3dvol(int_dat(:,2),volInfo,'outname',name1,'descrip','Positive values for 1st, neg for 2nd img');

        name1 = append_to_name(outname,'negpos');
        mask_vol(:,:,:,3) = iimg_reconstruct_3dvol(int_dat(:,3),volInfo,'outname',name1,'descrip','Positive values for 2nd, neg for 1st img');

        name1 = append_to_name(outname,'negneg');
        mask_vol(:,:,:,4) = iimg_reconstruct_3dvol(int_dat(:,4),volInfo,'outname',name1,'descrip','Negative values for both imgs');

    else

        mask_vol = iimg_reconstruct_3dvol(int_dat,volInfo,'outname',outname,'descrip','intersection image created with iimg_intersection');
    end
end

return




function name = append_to_name(outname,str)

[d,f,e] = fileparts(outname);
name = fullfile(d,[f str e]);

return
