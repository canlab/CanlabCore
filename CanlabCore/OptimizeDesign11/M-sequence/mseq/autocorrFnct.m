function [acor]=autocorrFnct(b,x)

%[acor]=autocorrFnct(b,x)

n=length(x);
acor(1)=1;

acor(2:n)=(1-b(1))*b(2).^x(2:end);

acor=acor(:);
