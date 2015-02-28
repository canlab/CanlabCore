function subdirs=spider_subdirs

% find the fragment in the path string

%d = dbstack('-completenames');
if datenum(version('-date')) >= datenum('May 6 2004')
 % look at this lovely mutual incompatibility between R14+ and previous releases
 % well done MathWorks, you've made our lives needlessly difficult AGAIN!
 d = dbstack('-completenames');
 rootdir = fileparts(d(1).file);
else
 d = dbstack; 
 rootdir = fileparts(d(1).name);
end   


if 0
    %%might be a relative path: convert to absolute
    olddir = pawd; cd(rootdir),
    rootdir = pawd;
    cd(olddir)
end

subdirs0 = {
	'Optimization'
	'functions'
	'basic'
	'clust'
	'demos'
	{'demos', 'data'}
	'density'
	'feat_sel'
	'mclass'
	'redset'  
 	'mod_sel'
	'pat'
	'reg'
	'trans'
};

t=1;
for i=1:length(subdirs0)    
	s = cellstr(subdirs0{i});
    if(isdir(fullfile(rootdir, s{:})))
        subdirs{t}=subdirs0{i};
        t=t+1;
    end
end
for i = 1:length(subdirs)
	s = cellstr(subdirs{i});
	subdirs{i} = fullfile(rootdir, s{:});
end
