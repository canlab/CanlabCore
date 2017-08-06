function [files_on_disk, url_on_neurovault, mycollection, myimages] = retrieve_neurovault_collection(collection_number)
% [files_on_disk, url_on_neurovault, mycollection] = retrieve_neurovault_collection(collection_number)
%
% Downloads and saves neuroimages from the Neurovault website for a
% specified collection number

collstr = num2str(collection_number); 

% Get info structure for a collection from Neurovault
mycollection = webread(['http://neurovault.org/api/collections/' collstr]);

% Print
fprintf('Downloading: %s\nOwner: %s, Images: %d\n%s\n', ...
mycollection.name, ...
mycollection.owner_name, ...
mycollection.number_of_images, ...
mycollection.url);

% Get images
myimages = webread(sprintf('http://neurovault.org/api/collections/%s/images', collstr));

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
    
    files_on_disk{i, 1} = websave(file_to_save{i}, url_on_neurovault{i});
    
end

% All done.
fprintf('\b\b\b\b\b\b\b\b\b\n')

end % function
