function Test_Mahal(file)
if nargin<1
    file='sample10k.csv';
end
UmapUtil.Initialize;
SAMPLE_FILE = 'sample10k.csv';
file = UmapUtil.RelocateExamples(file);
data= File.ReadCsv(file);

cov = nancov(data);
disp(['The covariance matrix for the data from ' SAMPLE_FILE ' is:']);
disp(cov);

invCov = inv(cov);
disp(['The inverse of the covariance matrix for the data from ' SAMPLE_FILE ' is:']);
disp(invCov);

X1 = knnsearch(data, data, 'K', 15, 'Distance', 'mahalanobis');
X2 = knnsearch(data, data, 'K', 15, 'Distance', 'mahalanobis', 'Cov', cov);



if isequal(X1, X2)
    disp('MATLAB''s knnsearch.m confirms that we correctly computed the covariance matrix!');
else
    disp('Hmm... did we calculate the covariance matrix incorrectly?');
end
X3 = KnnFind.Approximate(data, 15, 'mahalanobis', invCov, false, 3);
X4 = KnnFind.Approximate(data, 15, 'mahalanobis', [], false, 3);

acc = KnnFind.AssessApproximation(X1, X3);
disp(['nn_descent found nearest neighbors with ' num2str(100*acc) ' percent accuracy!']);
