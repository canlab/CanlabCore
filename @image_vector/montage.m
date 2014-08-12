function fig_handle = montage(image_obj, varargin)
% Create a montage of an image_vector (or statistic_image or fmri_data) object
%
% Usage:
% -------------------------------------------------------------------------
% [fig_handle or o2 fmridisp object] = montage(image_obj, [optional arguments])
%
% Optional inputs:
% 'fmridisplay' for fmridisplay object style montage [default]
% 'scnmontage' for circa 2008-style SCN lab montage for each image vector
%
% Examples:
% o2 = montage(mask);

meth = 'fmridisplay';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch vararagin{i}
            case 'scnmontage', meth = 'scnmontage';
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% number of images (cols)
n = size(image_obj.dat, 2);

switch meth
    case 'fmridisplay'
        
        r = region(image_obj);
        
        o2 = canlab_results_fmridisplay(r, 'noblobs', 'nooutline');
        
        if n > 8
            disp('Warning: Showing first 8 images in data object only.');
            n = 8;
        end
        
        for i = 1:n
            
            obj = image_obj;
            obj.dat = obj.dat(:, i);
            
            if isa(image_obj, 'statistic_image')
                obj.sig = obj.sig(:, i);
                obj.dat = obj.dat .* obj.sig;
            end
                
            if i == 1
            o2 = canlab_results_fmridisplay(region(obj), o2, 'nooutline');
            
            else
                o2 = canlab_results_fmridisplay(region(obj), o2, 'nooutline', 'addmontages');
            end
            
            drawnow

            %o2 = addblobs(o2, region(obj), 'splitcolor', {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]});
   
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

end % fnction
