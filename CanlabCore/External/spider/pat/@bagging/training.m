function [d,a] = training(a,d)

l=get_dim(d);
m=a.m;
if m<1 m=l*m; end; % m is a fraction of the whole dataset
as=[];rs=[];

% Bagging - "bootstrap aggregation", multiple versions of the training set,
% each created by drawing n' < n samples from the training set, with
% replacement.  (See Duda and Hart)
for i=1:a.bags
	r=ceil(l.*rand(min(l,m),1));
	dbag=get(d,r);

	[rss,ass]=train(a.child,dbag);
	rs{i}=rss;
	as{i}=ass;
end

a.child=group(as);

if ~a.algorithm.do_not_evaluate_training_error
	d=test(a,d);
end

