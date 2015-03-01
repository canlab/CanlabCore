
function [d,a] =  training(a,d)


%if alpha empty in child | alpha full in rss model, train child (first train/retrain).
%if alpha full in child, empty in rss, only train rss.

if isempty(a.child.alpha) | ~isempty(a.alpha) % train algorithm first
   [r a.child]=train(a.child,d);
end    

alpha=a.child.alpha; 
svs=find(sum(abs(alpha)',1)>1e-5); 

if(isfield(struct(a.child),'X'))
   a.Xsv=get(a.child.X,svs);
else
    a.Xsv=get(a.child.Xsv,svs);
end

K=calc(a.child.child,a.Xsv);
newalpha=alpha(svs,:)*0; alpha=alpha(svs,:); origalpha=alpha;
w2=[]; for i=1:size(alpha,2) w2(i)=alpha(:,i)'*K*alpha(:,i); end;
worig=w2; origsvs=svs;
loops=1;
disp('compressing..');
svs=[];
 
while max(w2)>a.tolerance & loops<a.max_loops
   %w is determined by alpha ONLY
  
   wx=alpha(:,1)*0;
   for i=1:size(alpha,2)  
    f=find(abs(alpha(:,i))>0); 
    tmp= w2(i) - (((K(:,f)*alpha(f,i)).^2) ./ diag(K)) ;  % calc all (w.x_i)^2/||x_i||^2 for each w
    if a.bal_w
     wx=wx+ (tmp/w2(i));  % was just tmp, attempt to normalize problems
    else
     wx=wx+tmp;
    end
   end       
   if a.dont_revisit 
     wx(svs)=Inf;    % ensure do not revisit existing basis function    
   end
   [val,I]=min(wx); % find arg max w.x_i    
   if ~isempty(intersect(I,svs)) & a.dont_revisit 
        break; % seen everything 
    end;
   oldsvs=svs;
   if length(oldsvs)<length(union(svs,I)) 
    svs=[oldsvs;I];
   end   
   
   if loops>0 & (mod(loops,a.backfit)==0  | (loops<a.backfit_at_start) | ...
		 (sum(loops==a.test_on)==1 & a.backfit>0) )
       % find optimal alphas by backfitting
       
       %replace this by finvupdate
      if ~strcmp(a.optimizer,'iterative')
       %inv=pinv(K(svs,svs)); % only need to compute inverse once
       minv=inv(K(svs,svs)+eye(length(svs))*1e-6);
                             % only need to compute inverse once          
       else
                newsv=I;
                 if length(oldsvs)==0  
                   R=1/sqrt(K(newsv,newsv)); minv=1/(K(newsv,newsv));
               else
                  newsv=I;
                  [R,minv]=finvupdate(a,R,K(oldsvs,newsv),K(newsv,newsv));
         end
       end   

     alpha=origalpha;
     for i=1:size(alpha,2)
      beta=minv*K(svs,:)*origalpha(:,i);
      newalpha(svs,i)=beta;  alpha(svs,i)=alpha(svs,i)-beta;
     end

	if ~isempty(a.dtst) & sum(loops==a.test_on)==1
	                 % run on separate test / validation set 
                         %  classification only right now

	  a.alpha=newalpha;a.Xsv=get(a.child.Xsv,origsvs);
	  fin=find(sum(abs(a.alpha)',1)>a.alpha_cutoff);
	  a.alpha=a.alpha(fin,:); a.Xsv=get(a.Xsv,fin);
          a.b0=a.child.b0;

         if a.reoptimize_b % reoptimize b // only for pat. rec. right now
         a.b0=0; 
	   r=test(a,d); a.b0=[];
	  for i=1:size(newalpha,2)  
	     rr=r.X; rr=rr(:,i); [x s2]=sort(rr); y=r.Y(s2,i);
            xs=cumsum(y(end:-1:1)==1);
            [m1 m2]=max(cumsum(y==-1)+xs(end:-1:1));
            a.b0=[a.b0 -x(m2)];  
            %if m2+1<=size(x,1) a.b0=[a.b0 -(x(m2)+x(m2+1))/2]; end;
         end 
         end 

        [r2]=test(a,a.dtst);
        [m1 c]=max(a.dtst.Y'); 
        [m1 c1]=max(r2.X');
        err=sum(c1~=c)/length(c); 
        a.res=[a.res; loops err sum(sum(abs(newalpha)'>0)>0) ];
	 a.res(:,1)' 
	 a.res(:,2)' 
        disp(sprintf('svs: %d    err: %f ',sum(sum(abs(newalpha)'>0)>0),err));
     end
     
   else % only update one alpha
    for i=1:size(alpha,2) 
     ai= (K(I,:)*alpha(:,i)) / K(I,I);   %   a_i = (x_i,w) / (x_i,x_i)  , optimal alpha
     alpha(I,i)=alpha(I,i)-ai;          %   w <- w - a_i x_i, project out learnt part of w 
     newalpha(I,i)=newalpha(I,i)+ai;    %   wnew <- wnew + a_i x_i, add to new w 
    end 
   end
   
   
    
    w2=[]; for i=1:size(alpha,2) f=find(abs(alpha(:,i))>0); w2(i)=alpha(f,i)'*K(f,f)*alpha(f,i); end;
    a.w2=w2; % store w2
    %% this line could be sped up

   if mod(loops,25)==0 
   txt='iteration %d : max||w_orig-w_new||^2=%1.3f svs=%d';
   disp(sprintf(txt,[loops max(w2) sum(sum(abs(newalpha)',1)>0)]));   
   end;
   loops=loops+1;
end




a.alpha=newalpha;

if(isfield(struct(a.child),'X'))
    a.Xsv=get(a.child.X,origsvs);
else
    a.Xsv=get(a.child.Xsv,origsvs);
end

fin=find(sum(abs(a.alpha)',1)>a.alpha_cutoff);
a.alpha=a.alpha(fin,:);
a.Xsv=get(a.Xsv,fin);
a.b0=a.child.b0;

if a.reoptimize_b % reoptimize b // only for pat. rec. right now
   a.b0=0; r=test(a,d); a.b0=[];
   for i=1:size(newalpha,2)  
    rr=r.X; rr=rr(:,i); [x s2]=sort(rr); y=r.Y(s2,i);
    xs=cumsum(y(end:-1:1)==1);
    [m1 m2]=max(cumsum(y==-1)+xs(end:-1:1));
    a.b0=[a.b0 -x(m2)];   
    %if m2+1<=size(x,1) a.b0=[a.b0 -(x(m2)+x(m2+1))/2]; end;
   end 
end

if a.algorithm.do_not_evaluate_training_error==1   
  d=set_x(d,get_y(d));
else
  d=test(a,d);
end
 
