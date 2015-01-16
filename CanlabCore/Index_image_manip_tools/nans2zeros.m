function nans2zeros(imgs)
   if(iscellstr(imgs))
       imgs = char(imgs);
   end
   fprintf('Loading meta-info...\n');
   Vimgs = spm_vol(imgs);
   
   fprintf('Zeroing... ');
   for i=1:length(Vimgs)
       status_string = sprintf('%d/%d', i, length(Vimgs));
       fprintf(status_string);
       
       data = spm_read_vols(Vimgs(i));
       data(isnan(data)) = 0;
       spm_write_vol(Vimgs(i), data);
       
       erase_string(status_string);
   end
   
   fprintf('\nFinished converting %d images.\n', length(Vimgs));
end