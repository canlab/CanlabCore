addpath(genpath('/dartfs-hpc/rc/home/m/f0042vm/software/fdaM'))

% create a basis set

dt = 0.0288;
order = 3;
l = 25;
degree = 8;

basis = create_bspline_basis([0,l/dt], degree+4, order);    
bf = full(eval_basis((1:l/dt),basis));
bf = bf(:,3:end-2);

obf = spm_orth(bf); % orthogonalize bf

figure(1);
subplot(3,1,1);
cla
h1 = plot(bf,'r');
title('Spline Basis Set');
subplot(3,1,2);
h2 = plot(obf,'b');
title('SPM orthogonalized splines');

% compute inverse orthogonal transform

iO = (obf'*obf)\obf'*bf;
iO(iO < 0) = 0;

% generate a random walk;

rand_walk = sum(triu(repmat(randn(length(bf),1),1,length(bf))));

figure(1);
subplot(3,1,3);
cla
h1=plot(rand_walk,'bla');
hold on;

% fit both basis fnuctions to the random walk

B = (bf'*bf)\bf'*rand_walk';
orth_B = (obf'*obf)\obf'*rand_walk';

h2=plot(obf*orth_B,'*b');
h3=plot(bf*B,'r');
h4=plot(bf.*repmat(B',size(bf,1),1),'g')

legend([h1(1),h2(1),h3(1),h4(1)],'Random Walk','spm spline model','spline model','weighted splines','location','southeast');
title('Spline model can be obtained with either basis set');

% so now we have the following
% splineModel = bf*B
% splineModel = obf*orth_B
% obf*iO = bf
% so we can rearrange
% obf*iO*b = obf*orth_B
% iO*b = orth_B
% b = inv(iO)*orth_B

i0 = (bf'*obf)'\obf'*obf;
iO(iO<0) = 0;

sum(inv(iO)*orth_B)
sum(B)
