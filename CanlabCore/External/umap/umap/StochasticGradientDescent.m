%   AUTHORSHIP
%   Primary Developers: cm.gautham@yahoo.in, Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef StochasticGradientDescent < handle
    properties(Constant)
        FINDING_ISLANDS='Finding data islands';
    end
    
    methods(Static)
        function cmd=GetCmd
            curPath=fileparts(mfilename('fullpath'));
            if ispc
                exe=String.ToSystem(...
                    fullfile(curPath, '/StochasticGradientDescent.exe'));
            else
                exe=String.ToSystem(fullfile(curPath,'/StochasticGradientDescent'));
                system(['xattr -r -d com.apple.quarantine ' exe]);
            end
            cmd=[ exe ' ' ];
        end
        
        function ok=IsAvailable
            cmd=StochasticGradientDescent.GetCmd;
            [status, ~]=system(cmd);
            ok=status==10;
        end
        
        function [out, cmdOut]=Go(inFile, outFile, head_embedding, ...
                tail_embedding,head, tail, n_epochs, n_vertices, ...
                epochs_per_sample, a, b, gamma, initial_alpha, ...
                negative_sample_rate, rand, dataDims, progress_callback)
            hasCallback=isequal('function_handle', class(progress_callback));
            if ~exist('head_embedding', 'var') || isempty(head_embedding) ...
               || ~exist('tail_embedding', 'var') || isempty(tail_embedding) ...
               || ~exist('head', 'var') || isempty(head) ...
               || ~exist('tail', 'var') || isempty(tail) ...
               || ~exist('n_epochs', 'var') || isempty(n_epochs) ...
               || ~exist('n_vertices', 'var') || isempty(n_vertices) ...
               || ~exist('epochs_per_sample', 'var') || isempty(epochs_per_sample) ...
               || ~exist('n_vertices', 'var') || isempty(n_vertices) ...
               || ~exist('a', 'var') || isempty(a) ...
               || ~exist('b', 'var') || isempty(b) ...
               || ~exist('gamma', 'var') || isempty(gamma) ...
               || ~exist('initial_alpha', 'var') || isempty(initial_alpha) ...
               || ~exist('negative_sample_rate', 'var') || isempty(negative_sample_rate)
                    out=[];
                    cmdOut = 1;
                    error('Not enough input parameters');
            end 
            StochasticGradientDescent.WriteText(inFile, head_embedding, ...
                tail_embedding,head, tail, n_epochs, ...
                n_vertices, epochs_per_sample, a, b, gamma, ...
                initial_alpha, negative_sample_rate);
            cmd=[ StochasticGradientDescent.GetCmd ' ' ...
                String.ToSystem(inFile) ' ' ...
                String.ToSystem(outFile) ' ' num2str(rand) ' 1'];
            fileSpec=[outFile '.*'];
            progressObj.getEpochsDone=1;
            progressObj.getEpochsToDo=n_epochs;
            progressObj.getEmbedding=head_embedding;
            delete(fileSpec);
            if hasCallback
                if ~feval(progress_callback, progressObj)
                    out=[];
                    cmdOut=1;
                    return;
                end
            else
                fprintf('0/%d epochs done\n', n_epochs);
            end
            if ismac 
                cmd=[cmd ' &'];
            elseif ispc
                cmd=[cmd ' < nul'];
            end
            delete(outFile);
            if ispc
                [flag, cmdOut]=system(['start "" /min ' cmd]);
            else
                [flag, cmdOut]=system(cmd);
            end
            if flag~=0
                error(cmdOut);
            end
            pu=[];
            pause(1);
            if dataDims>0
                sizeDsc=['<font color="blue">' ...
                    String.encodeInteger(size(head_embedding, 1)) ' x '...
                    num2str(dataDims) ' values</font><br>'];
            else
                sizeDsc='';
            end
            if ~hasCallback
                pu=PopUp(Html.WrapC...
                    (['Finding data islands (' ...
                    num2str(n_epochs) ' epochs)<hr><br>' sizeDsc...
                    '(<i>Click cancel to halt C++ executable</i>)...']), ...
                    'center', 'Stochastic gradient descent',  ...
                    false, true);
            end
            while true
                if exist(outFile, 'file')
                    break;
                end
                allFiles=dir(fileSpec);
                if ~isempty(allFiles)
                    allFiles=sortStructs(allFiles, 'datenum', 'descend');
                    fl=allFiles(1);
                    [~,~, epochs]=fileparts(fl.name);
                    progressObj.getEpochsDone=str2double(epochs(2:end));
                    epochFile=fullfile(fl.folder, fl.name);
                    if hasCallback
                        progressObj.getEmbedding=...
                            StochasticGradientDescent.ReadEmbedding(epochFile);
                        ok=feval(progress_callback, progressObj);
                    else
                        fprintf('%d/%d epochs done\n', ...
                            progressObj.getEpochsDone, n_epochs);
                        drawnow;
                        ok=~pu.cancelled;
                    end
                    delete(fileSpec);
                elseif hasCallback
                    ok=feval(progress_callback, ...
                        StochasticGradientDescent.FINDING_ISLANDS);
                else
                    drawnow;
                    ok=~pu.cancelled;
                end
                if ~ok
                    h = fopen([outFile '.STOP'], 'wb');
                    fprintf(h, "%s\n", 'STOP!');
                    fclose(h);
                    out=[];
                    cmdOut=1;
                    delete(inFile);
                    if ~hasCallback
                        pu.close;
                    end
                    return;
                end
                pause(3);
            end
            delete(fileSpec);
            if ~hasCallback
                pu.close;
            end
            if exist(outFile, 'file')
                out=StochasticGradientDescent.GetResult(outFile, inFile);
                if hasCallback
                    progressObj.getEpochsDone=n_epochs+1;
                    progressObj.getEmbedding=out;
                    feval(progress_callback, progressObj);
                else
                    fprintf('%d/%d epochs done\n', n_epochs, n_epochs);
                end
            else
                out=[];
            end
        end
        
        function [head_embedding]=GetResult(outFile, inFile)
            [head_embedding] = StochasticGradientDescent.ReadEmbedding(outFile);
            delete(outFile);
            delete(inFile);
        end
       
        function WriteEmbedding(outFile, head_embedding)
            h = fopen(outFile, 'wb');
            d=size(head_embedding);
            fprintf(h, "%d\n", d);
            fprintf(h, "%f\n",head_embedding);
        end
              
        function head_embedding=ReadEmbedding(inFile)
            h = fopen(inFile, 'rb');
            rows = fscanf(h, "%d\n", 1);
            cols=fscanf(h, "%d\n", 1);
            temp = fscanf(h, "%f\n", rows*cols);
            head_embedding = reshape(temp,rows,[]);
            %head_embedding = reshape(temp, [cols rows])';
            fclose(h);
        end
        
        
        function WriteText(inFile, head_embedding, tail_embedding,head, tail, n_epochs, ...
            n_vertices, epochs_per_sample, a, b, gamma, ...
            initial_alpha, negative_sample_rate)
            h = fopen(inFile, 'wb');
            
            d=size(head_embedding);
            fprintf(h, "%d\n", d(2));
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%f\n",head_embedding');
            if isequal(head_embedding, tail_embedding)
                fprintf(h, "0\n");
            else
                d=size(tail_embedding);
                fprintf(h, "%d\n", d(1));
                fprintf(h, "%f\n", tail_embedding');
            end
            d=size(head);
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%d\n", head);
            
            d=size(tail);
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%d\n", tail);
            
            d=size(n_epochs);
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%d\n", n_epochs);
            
            d=size(n_vertices);
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%d\n", n_vertices);
            
            d=size(epochs_per_sample);
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%f\n", epochs_per_sample);
            
            d=size(a);
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%f\n", a);
            
            d=size(b);
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%f\n", b);
            
            d=size(gamma);
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%f\n", gamma);
            
            d=size(initial_alpha);
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%f\n", initial_alpha);
            
            d=size(negative_sample_rate);
            fprintf(h, "%d\n", d(1));
            fprintf(h, "%d\n", negative_sample_rate);
                       
            fclose(h);
        end
        
        function [head_embedding, tail_embedding,head, tail, n_epochs, ...
            n_vertices, epochs_per_sample, a, b, gamma, ...
            initial_alpha, negative_sample_rate, randis] = ReadText(outFile)
            h = fopen(outFile, 'rb');
            
            d = fscanf(h, "%d\n", 1);  
            temp = fscanf(h, "%f\n", d*2);
            head_embedding = reshape(temp,d,[]);
            
            d = fscanf(h, "%d\n", 1);  
            temp = fscanf(h, "%f\n", d*2);
            tail_embedding = reshape(temp,d,[]);
            
            d = fscanf(h, "%d\n", 1);  
            head = fscanf(h, "%f\n", d);
            
            d = fscanf(h, "%d\n", 1);  
            tail = fscanf(h, "%f\n", d);
            
            d = fscanf(h, "%d\n", 1);  
            n_epochs = fscanf(h, "%f\n", d);
            
            d = fscanf(h, "%d\n", 1);  
            n_vertices = fscanf(h, "%f\n", d);
            
            d = fscanf(h, "%d\n", 1);  
            epochs_per_sample = fscanf(h, "%f\n", d);
            
            d = fscanf(h, "%d\n", 1);  
            a = fscanf(h, "%f\n", d);
            
            d = fscanf(h, "%d\n", 1);  
            b = fscanf(h, "%f\n", d);
           
            d = fscanf(h, "%d\n", 1);  
            gamma = fscanf(h, "%f\n", d);
            
            d = fscanf(h, "%d\n", 1);  
            initial_alpha = fscanf(h, "%f\n", d);
            
            d = fscanf(h, "%d\n", 1);  
            negative_sample_rate = fscanf(h, "%f\n", d);
            
            d = fscanf(h, "%d\n", 1);  
            randis = fscanf(h, "%f\n", d);
           
            fclose(h);
        end
 
    end
end
