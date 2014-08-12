function [imgs, sess_images, spersess] = scn_get_image_names(nruns, ndisdaqs)

    if nargin < 2, ndisdaqs = 0; 
    else
        fprintf('Removing first %03d image names from each session.\n', ndisdaqs); 
    end
    
 % get image names
    % --------------------------------------------------------------
    filestr = ['r[1-' num2str(nruns) ']' filesep 'vols.V*.img'];
    disp(['Getting images: ' filestr])
    imgs = filenames(filestr, 'char', 'absolute');

    clear sess_images

    for i = 1:nruns
        filestr = ['r[' num2str(i) ']' filesep 'vols.V*.img'];
        sess_images{i} = filenames(filestr, 'char', 'absolute');
        
        % disdaqs
        if ~isempty(sess_images{i}) && size(sess_images{i}, 1) > ndisdaqs && ndisdaqs
            sess_images{i}(1:ndisdaqs, :) = [];
        end
            
    end

    wh_empty = false(1, nruns);
    spersess = zeros(1, nruns);
    for i = 1:nruns
        if isempty(sess_images{i}) || size(sess_images{1}, 1) == 0  % do seem to need 2nd expression
            wh_empty(i) = true;
        else
            spersess(i) = size(sess_images{1}, 1);
        end
    end
    
    if ndisdaqs
        % re-get all image names, w.o disdaqs
        imgs = char(sess_images{:});
    end
    
    fprintf('Sessions (runs) with images: %03d\tSessions without images:%03d\n', sum(~wh_empty), sum(wh_empty));
    spersess(wh_empty) = [];
    fprintf('Images per Session (spersess):\n');
    fprintf('%3.0f\t', spersess);
    fprintf('\n');
    
end
