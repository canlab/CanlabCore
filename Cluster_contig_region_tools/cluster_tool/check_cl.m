function clout=check_cl(clin)

% Usage:
% clout=check_cl(clin)
% 
% Checks for the existence of important fields in the cl structure and
% attempts to resolve any problems by asking for using input at the command
% line. Checks some other possible problems in the cl structure as well.
% 
% 
fdir=pwd;

if length(size(clin))>2||min(size(clin))>1
    clin=clin(:);
end

if ~isfield(clin,'M')||~isfield(clin,'dim')
    if ~isfield(clin,'mat')
        beep
        disp('Input structure does not contain either an affine transformations matrix or image dimension information.\n');
        space_src=input('Would you like to load the correct information from an image or a cluster [I/c]?\n','s');
        if strcmp(space_src,'c')||strcmp(space_src,'C')
            handles.pwd=pwd;cd(fdir);
            [FileName,PathName]=uigetfile('*.mat','Select .mat file containing cl structure');
            cd(handles.pwd);
            fdir=PathName;
            try load([PathName FileName]);M=cl(1).M;dim=cl(1).dim;catch,error('Could not retrieve affine matrix from cl structure in specified file. Possibly no cl structure in the .mat file'),end
        else
            handles.pwd=pwd;cd(fdir);
            [FileName,PathName]=uigetfile('*.img','Select .img file');
            cd(handles.pwd);
            fdir=PathName;
            try space=iimg_read_img([PathName FileName]);M=space.mat;dim=space.dim;catch,error('Could not retrieve affine matrix from specified file'),end
        end
        for k=1:length(clin)
            clin(k).M=M;
            clin(k).dim=dim;
        end
    else
        for k=1:length(clin)
            clin(k).M=cl(k).mat;
        end
        clin=rmfield(clin,'mat');
    end
end

if ~isfield(clin,'voxSize')
    for k=1:length(clin)
        clin(k).voxSize=diag(clin(k).M)';
    end
end

poss{1}=clin(1).M;
for k=1:length(clin)
    equal=0;
    for j=1:length(poss)
        if isequal(poss{j},clin(k).M)
            equal=1;
        end
    end
    if ~equal
        poss{end+1}=clin(k).M;
    end
end
if length(poss)>1
    disp('The input structure contains at least two non-identical affine matrices.\nFunctions that make reference to the image space of the structure typically use only the affine matrix stored in cl(1).M.\n');
    fix=input('Would you like to try to remedy this [Y/n]?\n','s');
    if ~strcmp(fix,'n')&&~strcmp(fix,'N')
        for k=1:length(poss)
            disp(['Matrix #' num2str(k) ':'])
            fprintf(1,'%4.2f\t%4.2f\t%4.2f\t%4.2f\n%4.2f\t%4.2f\t%4.2f\t%4.2f\n%4.2f\t%4.2f\t%4.2f\t%4.2f\n%4.2f\t%4.2f\t%4.2f\t%4.2f\n\n\n',poss{k})    
        end
        toUse=input('Please enter the Matrix # for the affine matrix you would like to use.\nIt is highly recommended that you select the matrix that corresponds to the lowest spatial resolution.\n');
        M=poss{toUse};
        for k=1:length(clin)
            if isequal(M,clin(k).M)
                dim=clin(k).dim;
                break
            end
        end
        for k=1:length(clin)
            if ~isequal(M,clin(k).M)
                [clin(k).XYZ,point_list]=mmToVoxel(clin(k).XYZmm,M,'valid');
                clin(k).XYZmm=voxelToMm(clin(k).XYZ,M);
                clin(k).M=M;
                for m=length(point_list)
                    ind=find(point_list(:)==point_list(m));
                    clin(k).Z(ind)=mean(clin(k).Z(ind));
                end
                clin(k).Z=clin(k).Z(unique(point_list));                    
            end
            clin(k).voxSize=diag(M(1:3,1:3))';
            clin(k).dim=dim;
            clin(k).numVox=length(clin(k).Z);
        end        
    end
end

if ~isfield(clin,'title')
    clin(1).title=input('cl structure has no title. Please enter one now:\n','s');
    for k=2:length(clin)
        clin(k).title=clin(1).title;
    end
end

if ~isfield(clin,'threshold')
    clin(1).threshold(1)=input('cl structure has no threshold. Please enter lower threshold:\n');
    clin(1).threshold(2)=input('Please enter upper threshold:\n');
    for k=2:length(clin)
        clin(k).threshold=clin(1).threshold;
    end
end

if ~isfield(clin,'name')
    for k=1:length(clin)
        clin(k).name='';
    end
end

if ~isfield(clin,'numVox')
    for k=1:length(clin)
        clin(k).numVox=length(clin(k).Z);
    end
end

for k=1:length(clin)
    clin(k).mm_center=center_of_mass(clin(k).XYZmm,clin(k).Z);
end

for k=1:length(clin)
    if size(clin(k).XYZ,1)~=3
        clin(k).XYZ=clin(k).XYZ';
    end
    if size(clin(k).XYZmm,1)~=3
        clin(k).XYZmm=clcin(k).XYZmm';
    end
end

clout=clin;

disp('Done checking cluster integrity.')