%function [meanmap,pmap]=randperm(srp_data,compmat,numperm,corp)

function [realcomp]=randperm(srp_data,compmat,corp)

if nargin <2
    error('need to input raw data, comparison matrix,and p level');
end

if nargin <3
    corp=0.05;
end

a=size(srp_data);
if length(a)>2;
    error('for 2d only. use anovacont() or permtest()');
end

numsub=size(srp_data,1);
numcond=size(srp_data,2);

%disp('NOT correcting across the third dimension');

srp_data=srp_data(:,:,:,:);
%disp([(numsub*factorial(numcond)),' of your permutations will be the same as point estimate']);


%%calculate degrees of freedom
df_tot=(numsub*numcond)-1;
df_group=numcond-1;
df_subs=numsub-1;
df_error=df_tot-df_group-df_subs;

%    disp('calculating critical F....');

    smtot=squeeze(sum(srp_data,1));               %sum of each condition
    smtot1=squeeze(sum(srp_data,2));
    smmean=squeeze(mean(srp_data,1));          %mean of each condition
    ss_srp_data=srp_data.^2;                     %matrix of squared entries
    g_srp_data=squeeze(sum(smtot,2));        % this is G, the sum of all entries
    sum_ss_srp_data=sum(ss_srp_data(:)); %sum of squares
    sstot=sum_ss_srp_data-((g_srp_data.^2)/((numcond*numsub)));  %ss total
    ssgroup=((squeeze(sum(smtot.^2)))/(numsub))-((g_srp_data.^2)/(numcond*numsub));   %ss group
    sssubs=((squeeze(sum(smtot1.^2)))/numcond)-((g_srp_data.^2)/(numcond*numsub));   %ss group
    sserror=sstot-ssgroup-sssubs;
    msgroup=ssgroup/df_group;
    mserror=sserror/df_error;  %df error=36: 56(dftot)-2(dfgroup)-18(dfsubjects)
    realF=msgroup./(mserror+0.00001);         %this is the matrix of critical F values for frequency*time
    for c=1:size(compmat,1);
        for w=1:numcond;
        smtot2(w)=smmean(w)*compmat(c,w);
        end;
        comp1=squeeze(sum(smtot2));
        estvar1=(mserror./numsub)*(sum(compmat(c,:).^2));
        realcomp=comp1./sqrt(estvar1+0.000001);    %critical F values for frequency*time for each planned comparison
    end
