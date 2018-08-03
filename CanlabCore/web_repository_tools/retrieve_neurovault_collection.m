function [files_on_disk, url_on_neurovault, mycollection, myimages] = retrieve_neurovault_collection(collection_number)
% [files_on_disk, url_on_neurovault, mycollection] = retrieve_neurovault_collection(collection_number)
%
% Downloads and saves neuroimages from the Neurovault website for a
% specified collection number. This uses the Neurovault.org API from Chris
% Gorgolewski.  
%
% Gorgolewski KJ, Varoquaux G, Rivera G, Schwartz Y, Ghosh SS, Maumet C, Sochat VV, Nichols TE, Poldrack RA,
% Poline J-B, Yarkoni T and Margulies DS (2015) NeuroVault.org: a web-based repository for collecting and 
% sharing unthresholded statistical maps of the brain.Front. Neuroinform. 9:8. doi: 10.3389/fninf.2015.00008
%
% Example: 
% -------------------------------------------------------------------------
% Download images from Kragel et al. 2018 Nature Neuroscience
% Kragel, P. A., Kano, M., Van Oudenhove, L., Ly, H. G., Dupont, P., Rubio, A., ? Wager, T. D. (2018). 
% Generalizable representations of pain, cognitive control, and negative emotion in medial frontal cortex. 
% Nature Neuroscience, 21(2), 283?289. doi:10.1038/s41593-017-0051-7 
% 
% 270 subject-level images systematically sampled from 18 studies across 3 domains
% [files_on_disk, url_on_neurovault, mycollection, myimages] = retrieve_neurovault_collection(3324);
% data_obj = fmri_data(files_on_disk)
%
% Single-subject pain data contributed by CANlab (PI, Tor Wager), from
% Chang, Krishnan et al. 2015 PLoS Biology (BMRK4 study)
% [files_on_disk, url_on_neurovault, mycollection, myimages] = retrieve_neurovault_collection(504);
% data_obj = fmri_data(files_on_disk)

collstr = num2str(collection_number); 

options=weboptions;
options.CertificateFilename=('');
options.Timeout = 10;
% Get info structure for a collection from Neurovault
mycollection = webread(['http://neurovault.org/api/collections/' collstr],options);

% Print
fprintf('Downloading: %s\nOwner: %s, Images: %d\n%s\n', ...
mycollection.name, ...
mycollection.owner_name, ...
mycollection.number_of_images, ...
mycollection.url);

% Get images
myimages = webread(sprintf('http://neurovault.org/api/collections/%s/images', collstr),options);

while length(myimages.results)~=myimages.count
    if ~exist('moreimages','var')
    moreimages = webread(myimages.next,options);
    else
    moreimages = webread(moreimages.next,options);        
    end
    myimages.results=[myimages.results; moreimages.results];
end

clear moreimages;

% Get save location
saveroot = fullfile(pwd, num2str(mycollection.id));
if ~exist(saveroot, 'dir')
    
    mkdir(saveroot)
    
end



% Pull and save all images
nimgs = length(myimages.results);

[file_to_save, url_on_neurovault, files_on_disk] = deal(cell(nimgs, 1));

fprintf('Loading: %4.0f', 0);

for i = 1:nimgs
    
    fprintf('\b\b\b\b%4.0f', i);

    [~, img, extn] = fileparts(myimages.results(i).file);
    file_to_save{i} = fullfile(saveroot, [img extn]);
    
    url_on_neurovault{i} = myimages.results(i).file;
    
    files_on_disk{i, 1} = websave(file_to_save{i}, url_on_neurovault{i},options);
    
end

% All done.
fprintf('\b\b\b\b\b\b\b\b\b\n')

end % function
