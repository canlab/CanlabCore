function view_betas(cwd)

eval(['cd ' cwd])
d = dir('beta*img');
P = str2mat(d.name);

for i = 1:11:length(d)

	spm_check_registration(P(i:i+10,:))
	spm_orthviews('window',1:11,[-2 2],'interp',0)

	disp(['Session ' num2str(i)])

	keyboard

end