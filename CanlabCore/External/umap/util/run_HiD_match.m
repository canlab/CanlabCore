function [result, done]=run_HiD_match(trainingSet, trainingIds, testSet, ...
    testIds, varargin)
%%RUN_QF_MATCH runs the QF match algorithm 
%   on the subsets found within a training test set of data.
%
%   The publication that introduces the algorithm is 
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5818510/
%
%   A publication that further elaborates the algorithm is
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/
%
%   [reduction,umap] = RUN_QF_MATCH(trainingSet, trainingIds, 
%       testSet, testIds, 'NAME1',VALUE1,..., 'NAMEN',VALUEN) 
%   
%
%   RETURN VALUES
%   Invoking run_umap produces 2 return values:
%   1)result; an instance of the QFHiDM class describing match results
%       in the instance variable result.matches, result.matrixHtml
%       is a web page description of the result
%
%   2)done indicating success
%
%
%   REQUIRED INPUT ARGUMENT
%   trainingSet row/col matrix of data for training set
%   trainingIds numeric identifiers of training set subsets.  There is 
%       more than 1 column when subsets are overlapping,  run_qf_match
%       asserts same number of rows in trainingSet and trainingIds.  
%   testSet row/col matrix of data for test set
%   testIds numeric identifiers of test set subsets.  There is 
%       more than 1 column when subsets are overlapping,  run_qf_match
%       asserts same number of rows in testSet and testIds.  
%
%   OPTIONAL INPUT ARGUMENTS
%   The optional argument name/value pairs are:
%
%   Name                    Value
%   pu                      instance of PopUp.  run_qf_match uses 
%                           this to describe computing progress.  
%                           if 'none' is provided then no progress
%                           reporting occurs, otherwise the default is
%                           is an internal PopUp for reporting progress.
%   trainingNames           names of subsets in the test set.  Used
%                           for result.matrixHtml and is retained 
%                            in result.tNames.
%   testNames               names of subsets in the training set.  Used
%                           to annotate result.matrixHtml and is retained 
%                            in result.sNames.
%   trainingSetComp         the compensated data if appicable for testing
%                           standard deviation unit distance.  Must be
%                           same size as trainingSet. Default is empty.
%   testSetComp             the compensated data if appicable for testing
%                           standard deviation unit distance.  Must be
%                           same size as testSet. Default is empty.
%   log10                   if true then data columns where max > 1 
%                           is converted to log10
%   mergeLimit              0 if unlimited else 2-12 for 7-17 merge
%                           maximum merge candidates.  Default is 5 for
%                           maximum of 10 merge candidates.
%
%   matchStrategy           1 if quadratic form (QF), 2 if F-measure and
%                           3 if QF + F-measure optimizing for merging.
%                           1 is default
%
%   mergeStrategy           1 if best QF matches, 2 - 8  if
%                           percrent of top matches is 150% to 500% in 
%                           steps of 50%
%
%   AUTHORSHIP
%   Bioinformatics Lead:  Darya Orlova <dyorlova@gmail.com>
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
ml=BasicMap.Global.getNumeric(QfHiDM.PROP_MERGE_LIMIT, 1);
p=parseArguments();
parse(p,varargin{:});
[args, hadPu]=checkArgs(p.Results);
result=QfHiDM(trainingSet, args.trainingSetComp, trainingIds, ...
    testSet, args.testSetComp, testIds, args.bins, args.binStrategy);
result.tNames=args.trainingNames;
result.sNames=args.testNames;
result.matchStrategy=args.matchStrategy;
result.mergeStrategy=args.mergeStrategy;
done=result.compute(args.pu, true, 4, 1);
if ~hadPu
    args.pu.close;
end
BasicMap.Global.setNumeric(QfHiDM.PROP_MERGE_LIMIT, num2str(ml));

    function [args, hadPu]=checkArgs(args)
        hadPu=~isempty(args.pu);
        [trainers, trainingIds]=checkSet(trainingSet, trainingIds, args.trainingNames);
        [testees, testIds]=checkSet(testSet, testIds, args.testNames);
        if args.log10
            trainingSet=QfHiDM.Log10(trainingSet);
            testSet=QfHiDM.Log10(testSet);
        end
        if isempty(args.trainingSetComp)
            args.trainingSetComp=trainingSet;
        end
        if isempty(args.testSetComp)
            args.testSetComp=testSet;
        end
        if isempty(args.pu)
            args.pu=PopUp(['Matching with ' num2str(length(trainers)) ...
                ' training sets & ' num2str(length(testees)) ...
                ' test sets'], 'north', 'Running QF Match', false, true);
            args.pu.setTimeSpentTic(tic);
        elseif isequal('none', args.pu)
            args.pu=[];
        end
        BasicMap.Global.setNumeric(QfHiDM.PROP_MERGE_LIMIT, ...
            num2str(args.mergeLimit));
    end

    function [u, ids]=checkSet(dataSet, ids, names)
        [R1,C1]=size(dataSet);
        [R2, C2]=size(ids);
        assert((R1==R2) || (R1==C2 ), 'Need same # of data rows and id rows');
        if R1==C2 && R2==1
            ids=ids';
        end
        u=unique(ids(ids ~= 0));
        if ~isempty(names)
            R3=length(names);
            assert(length(u)==R3, 'Need same # of names as non zero IDs');
        end
    end

    function p=parseArguments(varargin)
        p = inputParser;
        addParameter(p,'bins', 0, @isnumeric);
        addParameter(p,'binStrategy', 0, @isnumeric);
        addParameter(p,'log10', true, @islogical);
        addParameter(p,'trainingSetComp', [], @isnumeric);
        addParameter(p,'testSetComp', [], @isnumeric);
        addParameter(p,'trainingNames', {}, @(x)StringArray.IsOk(x));
        addParameter(p,'testNames', {}, @(x)StringArray.IsOk(x));
        addParameter(p,'pu',[], @(x)isempty(x) || isequal('none', x) || isa(x,'PopUp'));
        addParameter(p,'mergeLimit',3, @(x)validMergeLimit(x));
        addParameter(p,'matchStrategy',1, @(x)validMatchStrategy(x));
        addParameter(p,'mergeStrategy', 1, @(x)validMergeStrategy(x));
        
        function ok=validMergeLimit(x)
            ok=x>=1 && x <=12;
        end

        function ok=validMatchStrategy(x)
            ok=x>=1 && x <=3;
        end
        function ok=validMergeStrategy(x)
            ok=x>=1 && x <=8;
        end

    end

end