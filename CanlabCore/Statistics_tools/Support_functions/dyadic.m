%adds another dimension to a matrix (the 'dyadic product').
%for example, if you want to mean subtract eeg data you might have one
%matrix 'eeg' (channels * frequency bands * time)
%   eeg=rand(65,5,100);
%say the first 10 timepoints are the baseline :
%   base=squeeze(mean(eeg(:,:,1:10),3));
%   size(base)
%   ans =
%    65     5
%you want to subtract the baseline from the eegdata. you can do this with a
%loop, but it's much easier if you make 'base' into a 65*5*100 matrix and
%simply subtract by matrix multiplication. typing
%   bigbase=dyadic(base,100);
%   size(bigbase)
%this does not entirely solve your problem, as the dimensions are in different 
%places.  however, using shiftdim you can align them:
%   bigbase=shiftdim(bigbase,1);
%   size(bigbase)
%   baseline_corrected_eeg=eeg-bigbase;
%
%....and there you are.  the format is dyadic(x,num) where x is your matrix
%and num is the number of entries in the dimension you want to have.
%dyadic works with up to 5 dimension matrices.
%© chris summerfield 2003.  summerfd@paradox.columbia.edu



function varb3=dyadic(varb,numb);

numdim=length(size(varb));
varb2=ones(numb,1)*(reshape(varb,1,prod(size(varb))));
if numdim==2;
    varb3=reshape(varb2,numb,size(varb,1),size(varb,2));
elseif numdim==3;
    varb3=reshape(varb2,numb,size(varb,1),size(varb,2),size(varb,3));
elseif numdim==4;
    varb3=reshape(varb2,numb,size(varb,1),size(varb,2),size(varb,3),size(varb,4));
elseif numdim==5;
    varb3=reshape(varb2,numb,size(varb,1),size(varb,2),size(varb,3),size(varb,4),size(varb,5));
end

