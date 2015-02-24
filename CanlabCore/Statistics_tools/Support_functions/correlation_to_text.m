function [str,sigmat] = correlation_to_text(a,cm,varargin)
    %str = correlation_to_text(matrix,value to mark with *,[cell array of names])
    %
    % input:
    % a) correlation matrix and critical value
    % b) or correlation matrix and 1/0 indicator matrix of significant values
    %   (sigmat)
    % c) or raw data and empty cm value
    %   (determines threshold based on n)
    %
    % see also print_correlation.m
    %
    % Examples:
    % correlation_to_text(corrcoef(dat), [], names)
    % correlation_to_text(dat, [], names)
    % correlation_to_text(partialrs, p < .05, cat(2, {clneg_data.shorttitle}));

    if length(varargin) > 0, nms = varargin{1}; else nms = {[]}; end

    if iscorrmtx(a)
        disp('Matrix already appears to be a correlation matrix.')
        n = [];
    else
        disp('Matrix is not a correlation matrix; Calculating corrcoef.');
        n = size(a,1);
        a = corrcoef(a);
    end


    if isempty(cm)
        % choose default value
        if isempty(n), n = input('Enter sample size: '); end
        [rci,sig,z,p,cm] = r2z(.5,n,.05);

    end

    % names

    str=sprintf('Crit=%3.2f\t',cm(1));
    for j=1:size(a,2),

        if length(nms) < j, nms{j} = ['V' num2str(j)]; end

        str=[str sprintf('%s\t',nms{j})];
    end
    str=[str sprintf('\n')];


    % table
    if ~isvector(cm);
        % already a sig matrix
        sigmat = cm;

        for i = 1:size(a,1),

            str=[str sprintf('%s\t',nms{i})];

            for j=1:size(a,2),

                if sigmat(i,j)
                    t='*';
                else t='';
                end

                str=[str sprintf('%3.3f%s\t',a(i,j),t)];
            end
            str=[str sprintf('\n')];

        end





    else
        % we have a critical value in cm

        sigmat=a*0;
        for i = 1:size(a,1),

            str=[str sprintf('%s\t',nms{i})];

            for j=1:size(a,2),

                if abs(a(i,j))>cm && a(i,j) ~= 1,t='*';
                    sigmat(i,j)=1;
                else t='';

                end,

                str=[str sprintf('%3.3f%s\t',a(i,j),t)];
            end
            str=[str sprintf('\n')];
        end
        sigmat=sigmat.*sign(a);

    end % input sigmat in cm or critical value



    disp(str)

end



function value = iscorrmtx(a)

    value = all(diag(a) == 1) && size(a,1) == size(a,2) && ~any(a(:) < -1) && ~any(a(:) > 1);

end
