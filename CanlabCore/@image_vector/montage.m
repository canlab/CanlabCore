function fig_handle = montage(image_obj, varargin)
% Create a montage of an image_vector (or statistic_image or fmri_data)
% object, on top of a standard anatomical underlay.
% - Takes optional inputs to fmridisplay.addblobs (see help for options)
% - Max 4 images in object. Use get_wh_image to select images from object if needed.
%
% *Usage:
% ::
%
%    [fig_handle or o2 fmridisp object] = montage(image_obj, [optional arguments])
%
% :Optional inputs:
%   **fmridisplay:**
%        for fmridisplay object style montage [default]
%
%   **scnmontage:**
%        for circa 2008-style SCN lab montage for each image vector
%
% :Examples:
% ::
%
%    o2 = montage(mask);
%
%    o2 = montage(dat, 'trans');
%    o2 = montage(dat, 'trans', 'color', [1 0 0]);
%    o2 = montage(dat, 'trans', 'color', [1 0 0], 'transvalue', .4);
%    o2 = montage(dat, 'trans', 'maxcolor', [1 .3 0], 'mincolor', [.5 0 1]);
%    o2 = montage(dat, 'trans', 'mincolor', [.5 0 1], 'transvalue', .5);
%    
%    Options for canlab_results_fmridisplay are also passed forward, so that 
%    different named montage types can be used.  E.g.,:
% 
%    o2 = montage(t, 'trans', 'full');
% 
% Set all color maps to the same range:
% o2 = montage(dat, 'trans', 'mincolor', [.5 0 1], 'transvalue', .7, 'cmaprange', [0 3]);

meth = 'fmridisplay';
% montagetype = 'compact2';  % default. This is passed into canlab_results_fmridisplay

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'scnmontage', meth = 'scnmontage';
                
            case {'trans', 'maxcolor', 'mincolor', 'transvalue', 'cmaprange', 'full', 'compact2'}
                % ignore these - passed through
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% number of images (cols)
n = size(image_obj.dat, 2);



switch meth
    case 'fmridisplay'
                    
        if n > 4
            disp('Warning: Showing first 4 images in data object only.');
            n = 4;
        end
        
        % r = region(image_obj);
        
        if n > 1
            do_multirow = true;
            o2 = canlab_results_fmridisplay([], 'noblobs', 'nooutline', 'multirow', n);
            
        else
            do_multirow = false;
            o2 = canlab_results_fmridisplay([], 'noblobs', 'nooutline', varargin{:}); % pass forward args.
        end
        
        montage_indx = 1:2:2*n;
        
        for i = 1:n
            
            obj = get_wh_image(image_obj, i); 
            
            if do_multirow % plot only on montages for this image
                
                o2 = addblobs(o2, region(obj), 'nooutline', varargin{:}, 'wh_montages', [montage_indx(i) montage_indx(i)+1]);
                
            else % just one image, plot on all montages
                
                o2 = addblobs(o2, region(obj), 'nooutline', varargin{:});
                
            end
            
            drawnow
            
            %o2 = addblobs(o2, region(obj), 'splitcolor', {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]});
            
        end
        
        if ~do_multirow
            % Plot legend only for single-row, otherwise obscures some
            % images
            o2 = legend(o2);
        end
        
        fig_handle = o2;
        clear o2
        
    case 'scnmontage'
        
        overlay = which('SPM8_colin27T1_seg.img');
        
        for i = 1:n
            
            % data from this image
            dat = image_obj.dat(:, i);
            
            % top and bottom 10%
            %dat(dat > prctile(dat, 10) & dat < prctile(dat, 90)) = 0;
            
            cl{i} = iimg_indx2clusters(dat, image_obj.volInfo);
            
            fig_handle(i) = montage_clusters(overlay, cl{i}, [2 2]);
            
            set(fig_handle, 'Name', sprintf('Montage %3.0f', i), 'Tag', sprintf('Montage %3.0f', i))
            
        end
        
    otherwise, warning(['Unknown input string option:' varargin{i}]);
        
end % switch

end % function
