function d =  generate(a)

if a.seed>0 rand('seed',a.seed); randn('seed',a.seed);   end;

X=randn(a.l,a.n);

X=X+a.input_noiselevel*randn(a.l,a.n);

if a.nonlinear==0

    if(~isempty(a.W))    
%         if(size(a.W,1)~=a.o & size(a.W,2)~=a.n)
        a.W=randn(a.o,a.n);
    end

    Y=(a.W*X')'+a.noiselevel*randn(a.l,a.o);
else
   Y=[];
   
   Y=get_x(test(a.outpmap,data(X)));
%    if(size(Y,2)~=a.o & ~isempty(a.o))
%        if(a.algorithm.verbosity>0)
%            warning('output of map does not fit a.o member');
%        end
%    end
   
end



d=data([get_name(a) ' '] ,X,Y);