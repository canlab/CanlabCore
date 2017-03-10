function str=get_name(a)
  str=' ';
  
  str=['kernel '];
  str=[str a.ker];
  
  if  str2num(version('-release'))<=13
    writingPerm=0;
    if iscell(a.kerparam),
      str = [str ' ('];
      for i=1:length(a.kerparam),
        temp = a.kerparam{i};
        if isa(temp,'char'),
          str = [str, a.kerparam{i} ', '];
          writingPerm = 1;
        else
          if length(temp)==1&~iscell(temp),
            str = [str,num2str(temp) ', '];
            writingPerm=1;
          end
        end
      end
      if writingPerm,   
        str = [str(1:length(str)-2) ') '];
      else
        str = str(1:length(str)-2);
      end
    else
      if ~isempty(a.kerparam) 
        if length(a.kerparam)==1
          str = [str '=' num2str(a.kerparam)  ];
        end
      end
    end
    
  else %% if version
    
    kern_par = a.kerparam;
    
    if ~isempty( kern_par) 
      switch class( kern_par)
        
       case 'cell'
        %             disp( 'type = cell');
        str = [ str ' ('];
        temp = evalc( 'disp( kern_par)');
        temp = regexprep( temp, '([\]}''])\s*([\[{''])', '$1, $2'); % remove spaces between elements
        temp = regexprep( temp, '[\[\]{}]', '');                    % remove brackets and braces
        temp = regexprep( temp, '^\s*', '');                        % remove leading ...
        temp = regexprep( temp, '\s*$', '');                        % ... and trailing whitespaces
        str = [ str temp ')'];
        
       case 'struct'
        %             disp( 'type = struct');
        str = [ str ' ('];
        temp = evalc( 'disp( kern_par)');
        temp = regexprep( temp, '^\s*', '');   % remove leading ...
        temp = regexprep( temp, '\s*$', '');   % ... and trailing whitespaces
        temp = regexprep( temp, '\n', ';');    % replace linebreak with ';'
        temp = regexprep( temp, '\s*', ' ');   % leave only one whitespace in between
        str = [ str temp ')'];
        
       case 'char'
        %             disp( 'type = char');
        str = [ str '=' kern_par ];
        
       case 'double'
        %             disp( 'type = double');
        if length( kern_par) == 1
          str = [ str '=' num2str( kern_par)];
        else
          temp = evalc( 'disp( {kern_par})');
          temp = regexprep( temp, '^\s*', ''); % remove leading ...
          temp = regexprep( temp, '\s*$', ''); % ... and trailing whitespaces
          str = [ str '(' temp ')'];
        end
        
       otherwise
        str = [ str '(type = ' class( kern_par) ')'];
      end
    end
  end
