function d =  generate(a)

if a.seed>0 rand('seed',a.seed);   end;
if a.seed>0 randn('seed',a.seed);   end;

switch a.dist
case 'gaussians'
	d=bayes({gauss([-1 0]) gauss([1 0])});
	d.l=a.l;
	d=gen(d);
otherwise
	d=make_data(a);
end
end	

