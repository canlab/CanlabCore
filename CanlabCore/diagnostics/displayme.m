function [subjM,Mtotalv] = displayme(mm,txtlab,tlab2)
% Used in img_hist2 - included as internal function there.
% this function is for indepenent re-display after img_hist2 is finished.
%
% :Usage:
% ::
%
%     function [subjM,Mtotalv] = displayme(mm,txtlab,tlab2)
%
% :Example:
% ::
%
%    % TO run:
%    [O.subjM,O.Mtotalv] = displayme(O.m,txtlab,'MEANS');
%    [O.subjS,O.Stotalv] = displayme(O.s,txtlab,'STD');
%    [O.subjW,O.Wtotalv] = displayme(O.w,txtlab,'SKEWNESS');
%    [O.subjK,O.Ktotalv] = displayme(O.k,txtlab,'KURTOSIS');
%
% ..
%    Tor Wager
% ..

gm = mean(mean(mean(mm)));
Mtotalv = sum((mm(:) - gm).^2);

disp(['Variability across ' tlab2 ' within mask: ' txtlab])
disp(['Grand mean (over the experiment) is ' num2str(gm) ', SSt = ' num2str(Mtotalv)])

subjM = getmeans(mm,2,3);
Mexp = (prod(size(mm)) ./ size(mm,1)) * sum((subjM - mean(subjM)).^2);
dfb = size(mm,1) - 1; dfw = prod(size(mm)) - (dfb+1);
F = Mexp./(Mtotalv-Mexp);
fprintf(1,'SESSION\tm %3.2f\tstd %3.2f\t SSb %3.2f\t%%Var explained %3.2f\tF(%3.0f,%3.0f) = %3.2f\tp = %3.5f\n', ...
	mean(subjM),std(subjM),Mexp,100*(Mexp./Mtotalv),dfb,dfw,F,1-fcdf(F,dfb,dfw));

subjM = getmeans(mm,1,3);
Mexp = (prod(size(mm)) ./ size(mm,2)) * sum((subjM - mean(subjM)).^2);
dfb = size(mm,2) - 1; dfw = prod(size(mm)) - (dfb+1);
F = Mexp./(Mtotalv-Mexp);
fprintf(1,'RUN\tm %3.2f\tstd %3.2f\t SSb %3.2f\t%%Var explained %3.2f\tF(%3.0f,%3.0f) = %3.2f\tp = %3.5f\n', ...
	mean(subjM),std(subjM),Mexp,100*(Mexp./Mtotalv),dfb,dfw,F,1-fcdf(F,dfb,dfw));

subjM = getmeans(mm,1,2);
Mexp = (prod(size(mm)) ./ size(mm,3)) * sum((subjM - mean(subjM)).^2);
dfb = size(mm,3) - 1; dfw = prod(size(mm)) - (dfb+1);
F = Mexp./(Mtotalv-Mexp);
fprintf(1,'CONDITION\tm %3.2f\tstd %3.2f\t SSb %3.2f\t%%Var explained %3.2f\tF(%3.0f,%3.0f) = %3.2f\tp = %3.5f\n', ...
	mean(subjM),std(subjM),Mexp,100*(Mexp./Mtotalv),dfb,dfw,F,1-fcdf(F,dfb,dfw));


return

gm = mean(mean(mean(mm)));
Mtotalv = sum((mm(:) - gm).^2);

disp(['Variability across ' tlab2 ' within mask: ' txtlab])
disp(['Grand mean (overall bias in experiment) is ' num2str(gm)])

subjM = getmeans(mm,2,3);
Mexp = (prod(size(mm)) ./ size(mm,1)) * sum((subjM - mean(subjM)).^2);
dfb = size(mm,1) - 1; dfw = prod(size(mm)) - (dfb+1);
F = Mexp./(Mtotalv-Mexp);
fprintf(1,'SESSION\tm %3.2f\tstd %3.2f\t SSb %3.2f\t%%Var explained %3.2f\tF(%3.0f,%3.0f) = %3.2f\tp = %3.5f\n', ...
	mean(subjM),std(subjM),Mexp,100*(Mexp./Mtotalv),dfb,dfw,F,1-fcdf(F,dfb,dfw));

subjM = getmeans(mm,1,3);
Mexp = (prod(size(mm)) ./ size(mm,2)) * sum((subjM - mean(subjM)).^2);
dfb = size(mm,2) - 1; dfw = prod(size(mm)) - (dfb+1);
F = Mexp./(Mtotalv-Mexp);
fprintf(1,'RUN\tm %3.2f\tstd %3.2f\t SSb %3.2f\t%%Var explained %3.2f\tF(%3.0f,%3.0f) = %3.2f\tp = %3.5f\n', ...
	mean(subjM),std(subjM),Mexp,100*(Mexp./Mtotalv),dfb,dfw,F,1-fcdf(F,dfb,dfw));

subjM = getmeans(mm,1,2);
Mexp = (prod(size(mm)) ./ size(mm,3)) * sum((subjM - mean(subjM)).^2);
dfb = size(mm,3) - 1; dfw = prod(size(mm)) - (dfb+1);
F = Mexp./(Mtotalv-Mexp);
fprintf(1,'CONDITION\tm %3.2f\tstd %3.2f\t SSb %3.2f\t%%Var explained %3.2f\tF(%3.0f,%3.0f) = %3.2f\tp = %3.5f\n', ...
	mean(subjM),std(subjM),Mexp,100*(Mexp./Mtotalv),dfb,dfw,F,1-fcdf(F,dfb,dfw));


return




function a = getmeans(a,d1,d2)
% input mat, dim to avg over, 2nd dim to avg over

a = squeeze(mean(mean(a,d1),d2));
if size(a,2) > size(a,1), a = a';, end

return
