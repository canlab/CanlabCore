%% METHODS for canlab object-oriented tools
% The list below prints the first 200 characters of help for each method.
% You can get additional information by typing >>help objecttype.methodname
% in the Matlab command window.
%%
z = '_________________________________________';

%objname = 'fmri_data'; % for example

fprintf('%%%% %s\n', objname); 

m = methods(objname);
%%

for i = 1:length(m)
    %%% Method
    % New method begins here
    
    % *Method name*
    mname = [objname '.' m{i}];
    h = help(mname);
    %h = help(m{i});
    
    len = min(200, length(h));
    
    fprintf('\n%s\n%s\n%s\n', z, mname, z);
    
    disp(h(1:len));
    
    fprintf('\n')
    
end



