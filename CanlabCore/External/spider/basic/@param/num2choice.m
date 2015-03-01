function ch =  num2choice(a,num,sz)
 
  num=num-1;
  
  for i=1:length(sz)
    ch(i)=mod(num,sz(i));
    num=floor(num/sz(i));
  end
  
  ch=ch+1; 