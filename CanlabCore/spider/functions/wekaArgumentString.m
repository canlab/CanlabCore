function arguments=wekaArgumentString(carray,appString)
% This function builds a weka String array used for Argument passing
% Argument has to be of the type : CELLARRAY
% where the following format is chosen for argument passing
% { '-A1', A1, '-A2', A2, ...}
%
% the odd content must be a string which will be parsed by weka. 
% Most of the weka argument strings are of the form " - CHARACTER"
%
% All even arguments must be of type number or string. 

if(nargin ==1 )
arguments=javaArray('java.lang.String',length(carray));


for i=1:2:length(carray)
  arguments(i)=java.lang.String(carray{i});
end
% 
for i=2:2:length(carray)
    if ( isnumeric(carray{i}))
        arguments(i)=java.lang.String(num2str(carray{i}));        
    elseif ( ischar(carray{i}))
        arguments(i)=java.lang.String(carray{i});            
    end
end

else

    arguments=javaArray('java.lang.String',length(carray)+length(appString));
    for i=1:length(appString)
        arguments(i)=appString(i);
    end
    s=wekaArgumentString(carray);
    for i=1:length(s)
        arguments(i+length(appString))=s(i);
    end
end