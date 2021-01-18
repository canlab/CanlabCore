function [embedding, method] = choose_optimize_layout(head_embedding,...
    tail_embedding, head, tail, n_epochs, n_vertices,epochs_per_sample,...
    a, b, gamma, initial_alpha, negative_sample_rate, verbose, method,...
    progress_callback, epoch_reports, random_state, min_dist, ...
    move_point, dataDims)
%CHOOSE_OPTIMIZE_LAYOUT Given all the data necessary to perform stochastic
% gradient descent, use the "method" variable to decide whether to use
% Java, MatLab C coder, C++ executable, mex C++ or MATLAB to perform SGD.
%
% See also: OPTIMIZE_LAYOUT

%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
TEST_CROSS_ENTROPY=false;
doingJoinedTransform=nargin>=19 && ~isempty(move_point);
if nargin<19
    dataDims=[];
    if nargin < 16
        epoch_reports = 0;
        if nargin<15
            progress_callback=[];
            if nargin < 14
                method = 'Java';
                if nargin < 13
                    verbose = false;
                end
            end
        end
    end
end
if strcmpi(method, 'C')
    if ~exist(['optimize_layout_mex.' mexext], 'file')
        if initJava
            method='Mex';
        end
    end
end
if strcmpi(method, 'C vectorized')
    if ~exist(['optimize_layout2_mex.' mexext], 'file')
        if initJava
            method='Mex';
        end
    end
end
mexCancelled=false;
if strcmpi(method, 'Mex')
    try
        curPath=fileparts(mfilename('fullpath'));
        
        if islogical(random_state)
            if random_state
                rand = 0;
            else
                rand = -1;
            end
        end
        [~,exe]=UmapUtil.LocateMex('sgd');
        if ismac
            system(['xattr -r -d com.apple.quarantine ' String.ToSystem(exe)]);
        end
        if ~exist(exe, 'file')
            method='Java';
        elseif verbose
            if reportMexProgress(head_embedding, [1 n_epochs])==0
                embedding=[];
                return;
            end
            
            if isequal(head_embedding, tail_embedding)
                %move_other [] for tail_embedding
                embedding=mexStochasticGradientDescent(head_embedding, [],...
                    uint32(head), uint32(tail), uint32(n_epochs), ...
                    uint32(n_vertices), epochs_per_sample, a, b, gamma,...
                    initial_alpha, negative_sample_rate, rand, ...
                    @(data, epochs)reportMexProgress(data, epochs));
            else
                embedding=mexStochasticGradientDescent(head_embedding, ...
                    tail_embedding, uint32(head), uint32(tail), ...
                    uint32(n_epochs),  uint32(n_vertices), epochs_per_sample, ...
                    a, b, gamma, initial_alpha, negative_sample_rate, rand, ...
                    @(data, epochs)reportMexProgress(data, epochs));
            end

        else
            if isequal(head_embedding, tail_embedding)
                %move_other [] for tail_embedding
                embedding=mexStochasticGradientDescent(head_embedding, [],...
                    uint32(head), uint32(tail), uint32(n_epochs), ...
                    uint32(n_vertices), epochs_per_sample, a, b, gamma,...
                    initial_alpha, negative_sample_rate, rand);
            else
                embedding=mexStochasticGradientDescent(head_embedding, ...
                    tail_embedding, uint32(head), uint32(tail), ...
                    uint32(n_epochs),  uint32(n_vertices), epochs_per_sample, ...
                    a, b, gamma, initial_alpha, negative_sample_rate, rand);
            end
        end
        if mexCancelled
            embedding=[];
        end
        if strcmpi(method, 'Mex')
            return;
        end
     catch ex
        ex.getReport
        msg('Mex version unavailable using slower Java');
        method='Java';
    end
end
if strcmpi(method, 'C++')
    try
        t=datetime;
        t.Format='yyMMdd_HHmmss';
        s=char(t);
        fldr=fullfile(File.Home, '.umap');
        File.mkDir(fldr);
        inFile=fullfile(fldr, ['in_' s '.txt']);
        outFile=fullfile(fldr, ['out_' s '.txt']);
        rand=0;
        if islogical(random_state)
            if ~random_state
                rand=-1;
            end
        end
        embedding=StochasticGradientDescent.Go(inFile, ...
            outFile, head_embedding, tail_embedding, ...
            head, tail, n_epochs, n_vertices, epochs_per_sample, a, b, gamma, ...
            initial_alpha, negative_sample_rate, rand, dataDims, ...
            progress_callback);
        return;
     catch ex
        ex.getReport
        method='Java';
    end 
end
if strcmpi(method, 'Java')
    initJava;
    try
        N=size(epochs_per_sample, 1);
        if TEST_CROSS_ENTROPY
            weights = ones(N,1)./epochs_per_sample; %We probably should have passed in weights to this instead...
        end 
        javaObject=edu.stanford.facs.swing.StochasticGradientDescent(...
            head_embedding, tail_embedding,head, tail, n_epochs, ...
            n_vertices, epochs_per_sample, a, b, gamma, ...
            initial_alpha, negative_sample_rate);
        javaObject.move_other=isequal(head_embedding, tail_embedding);
        if doingJoinedTransform
            % java code only inspects move_point if move_other is TRUE
            %javaObject.move_other=true;
            javaObject.move_point=move_point;
        else
            javaObject.move_other=isequal(head_embedding, tail_embedding);
        end
        if islogical(random_state)
            if ~random_state
                javaObject.randomize;
            end
        end
        if epoch_reports>0
            javaObject.setReports(epoch_reports);
        end
        if ~reportJavaProgress
            embedding=[];
            return;
        end
        while ~javaObject.nextEpochs
            if ~reportJavaProgress
                embedding=[];
                return;
            end
        end
        reportJavaProgress;
        embedding=javaObject.getEmbedding;
        return;
    catch ex
        ex.getReport
        if doingJoinedTransform
            warning(' JAVA jar not installed? .. using MatLab');
            method='MatLab';
        else
            warning(' JAVA jar not installed? .. using C');
            method='C';
        end
    end
end
if strcmpi(method, 'C')
    try
        embedding = optimize_layout_mex(single(head_embedding), ...
            single(tail_embedding), int32(head), int32(tail), ...
            int32(n_epochs), int32(n_vertices), ...
            single(epochs_per_sample), single(a), single(b), ...
            single(gamma), single(initial_alpha), ...
            int32(negative_sample_rate), verbose);
    catch ex
        ex.getReport
        yelp;
    end
elseif strcmpi(method, 'MatLab')
    if ~doingJoinedTransform
        embedding = optimize_layout(single(head_embedding), ...
            single(tail_embedding), int32(head), int32(tail), ...
            int32(n_epochs), int32(n_vertices), ...
            single(epochs_per_sample), single(a), single(b), ...
            single(gamma), single(initial_alpha), ...
            int32(negative_sample_rate), verbose);
    else 
        embedding=trans2_optimize_layout(head_embedding, ...
            tail_embedding, head, tail, ...
            n_epochs, n_vertices, ...
            epochs_per_sample, a, b, ...
            gamma, initial_alpha, ...
            negative_sample_rate, verbose, min_dist, move_point);
    end
elseif strcmpi(method, 'C vectorized')
    embedding = optimize_layout2_mex(single(head_embedding), ...
        single(tail_embedding), int32(head), int32(tail), ...
        int32(n_epochs), int32(n_vertices), ...
        single(epochs_per_sample), single(a), single(b), ...
        single(gamma), single(initial_alpha), ...
        int32(negative_sample_rate), verbose);
elseif strcmpi(method, 'MatLab experimental')
    embedding = optimize_layout4(single(head_embedding), ...
        single(tail_embedding), int32(head), int32(tail), ...
        int32(n_epochs), int32(n_vertices), ...
        single(epochs_per_sample), single(a), single(b), ...
        single(gamma), single(initial_alpha), ...
        int32(negative_sample_rate), verbose);
elseif strcmpi(method, 'MatLab experimental 2')
    embedding = optimize_layout5(single(head_embedding), ...
        single(tail_embedding), int32(head), int32(tail), ...
        int32(n_epochs), int32(n_vertices), ...
        single(epochs_per_sample), single(a), single(b), ...
        single(gamma), single(initial_alpha), ...
        int32(negative_sample_rate), verbose);
else %method is MatLab vectorized
    embedding = optimize_layout2(single(head_embedding), ...
        single(tail_embedding), int32(head), int32(tail), ...
        int32(n_epochs), int32(n_vertices), ...
        single(epochs_per_sample), single(a), single(b), ...
        single(gamma), single(initial_alpha), ...
        int32(negative_sample_rate), verbose);
end
    function ok=reportMexProgress(data, epochs)
        if isequal('function_handle', class(progress_callback))
            progressObj.getEpochsDone=double(epochs(1));
            progressObj.getEpochsToDo=double(epochs(2));
            progressObj.getEmbedding=data;
            ok=feval(progress_callback, progressObj);
            if ok
                ok=int32(1);
            else
                mexCancelled=true;
                ok=int32(0);
            end
        else
            fprintf('%d of %d epochs\n', epochs(1)-1, epochs(2))
            ok=int32(1);  
            drawnow;
        end
    end

    function yelp
        UmapUtil.OfferFullDistribution(false);
    end

    function wantsToContinue=reportJavaProgress        
        if isequal('function_handle', class(progress_callback))
            wantsToContinue=feval(progress_callback, javaObject);
        else
            wantsToContinue=true;
            if verbose
                done=javaObject.getEpochsDone-1;
                toDo=javaObject.getEpochsToDo;
                fprintf('%d/%d epochs done\n', done, toDo);
            end
        end
        if TEST_CROSS_ENTROPY
            if size(javaObject.head_embedding, 1)*size(javaObject.tail_embedding) < 1e7
                CE = cross_entropy(javaObject.head_embedding, javaObject.tail_embedding, head, tail, weights, a, b, javaObject.move_other);
                fprintf('The cross entropy is %s\n', ...
                    String.encodeRounded(CE,1));
            end

        end
    end
end