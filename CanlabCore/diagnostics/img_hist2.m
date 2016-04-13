function O = img_hist2(subdir)
% A general function for plotting histograms of any image
% For each subject, comparing across subjects
%
% :Inputs:
%
%   **imgname:**
%        name of image file to make intensity histograms from
%
%   **subdir:**
%        cell array of text strings containing names of individual subject
%        directories (wherein are contained the file specified in imgname 
%        or each subject)
%
% Performs the histogram plot a number of times, without plotting
% and reports the variance in pdf moments as a function of subject, 
% run, and condition (beta img within run).
%
% Start in directory above individual subject results
%
% :Examples:
% ::
%
%    img_hist2(EXPT.subjects)
%    img_hist2({'020827mk' '020829jh' '020903lb'})
%
% ..
% Tor Wager
% ..

% ..
%    defaults
% ..
dwcard = '02*';									% wildcard defining ind subject directories, e.g., '02*'
%csfname = 'rnnhet1spgr_seg3.img';	% reslice first to space of functionals
%grname = 'rnnhet1spgr_seg1.img';
%csfpath = '/data/placebo/';		% before ind subject
%csfp2 = 'anatomy';							% after ind subject

Pc = getcsfname(subdir);
PP = getimgname;
 
[dummy,txtlab] = fileparts(Pc(1,:));


% rows = subjects
% cols = moment
% 
% rows = subjects
% cols = runs
% 3d   = beta (condition)

for condition = 1:size(PP,1)
   for run = 1:size(PP,2) 
    
       fprintf(1,'%s\n ',PP{condition,run})
       % hP = get_filename(dwcard,PP{condition,run});
	hP = tor_list_files(subdir,PP{condition,run});
 	hP = str2mat(hP(:));
       csfM = spm_general_hist(hP,Pc,txtlab,0);
       
       % Mean, Std, skeW, Kurtosis
       O.m(:,run,condition) = csfM(:,1);
       O.s(:,run,condition) = csfM(:,2);
       O.w(:,run,condition) = csfM(:,3);
       O.k(:,run,condition) = csfM(:,4);
    
   end
end
    

fprintf(1,'\n ') 
 
diary img_hist2.out

disp(['EXPERIMENT: ' pwd]);disp('')

[O.subjM,O.Mtotalv] = displayme(O.m,txtlab,'MEANS');disp('')
[O.subjS,O.Stotalv] = displayme(O.s,txtlab,'STD');disp('')
[O.subjW,O.Wtotalv] = displayme(O.w,txtlab,'SKEWNESS');disp('')
[O.subjK,O.Ktotalv] = displayme(O.k,txtlab,'KURTOSIS');disp('')


diary off

keyboard

return





    

function Pc = getcsfname(subdir)


doseg = input('Do anatomical images need to be segmented? (1/0)');
if doseg
    % get the list of csf images for each subject
    myp = spm_get(1,'*img',['Choose anatomical img for ' subdir{1} '-  Others should have similar path.']);
    [d f e] = fileparts(myp);
    disp([myp])
    %subc = input('Enter subject code part of this path (replaced for all subjects): ','s');
    x = findstr(d,subdir{1});
    [d2 f2] = fileparts(d(x:end));
    d = d(1:x-1);
    % use d, subdir, and f2 as path
    Pc = fullfile(d,subdir{1},f2,[f e]);

    for i = 2:length(subdir)
        Pc = str2mat(Pc,fullfile(d,subdir{i},f2,[f e]));
    end  
    
    try
        V = spm_vol(Pc);
    catch
        Pc
        disp('Error loading images!')
        keyboard
    end
    
    disp('Images must be resliced to space of functionals.')
    fname = spm_get(1,'*img','Choose a functional image from this study.');
    reslice_imgs(fname,Pc,1);
    
    canonicalT1 = spm_get(1,'*img','Choose Template Image.');
    spm_segment(canonicalT1,'C')
end


% get the list of csf images for each subject
myp = spm_get(1,'*img',['Choose csf mask (in space of functionals) for ' subdir{1} '-  Others should have similar path.']);
[d f e] = fileparts(myp);
disp([myp])
%subc = input('Enter subject code part of this path (replaced for all subjects): ','s');
x = findstr(d,subdir{1});
[d2 f2] = fileparts(d(x:end));
d = d(1:x-1);
% use d, subdir, and f2 as path
Pc = fullfile(d,subdir{1},f2,[f e]);

for i = 2:length(subdir)
    Pc = str2mat(Pc,fullfile(d,subdir{i},f2,[f e]));
end



return





function PP = getimgname

% get the list of betas within runs
mypwd = pwd;
mydir = spm_get(-1,'*','Choose an individual subject SPM results dir');
eval(['cd ' mydir])
d = dir(['beta*img']);
load SPM
P = str2mat(d.name);
ind = 1;
for i = 1:length(Sess)
    for j = 1:length(Sess{i}.col)
        PP{j,i} = deblank(P(ind,:));
        ind = ind + 1;
    end
end
eval(['cd ' mypwd])

return





function a = getmeans(a,d1,d2)
% input mat, dim to avg over, 2nd dim to avg over

a = squeeze(mean(mean(a,d1),d2));
if size(a,2) > size(a,1), a = a';, end

return






function [subjM,Mtotalv] = displayme(mm,txtlab,tlab2)

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


