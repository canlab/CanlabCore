function [XYZmm_out]=space_match(XYZmm,voxsize,varargin)

% Usage:
% XYZmm_out=space_match(XYZmm,voxsize,['force'])
% 
% Takes a series of points in mm coordinates, XYZmm, and matches each point
% to the nearest point in a space defined by voxsize and the coordinate
% point 0,0,0.
% 
% Including the string 'force' will force a matched point that is
% equidistant from at least 2 points in the space to be assigned to the
% point with the highest absolute value in each dimension.
% 

if nargin>2
    if strcmp(varargin{1},'force')
        force=1;
    else
        force=0;
    end
else
    force=0;
end

if size(XYZmm,1)~=3
    XYZmm=XYZmm';
end

if length(XYZmm)>3
    l=length(XYZmm);
else
    l=min(size(XYZmm));
end

XYZmm_out=zeros(size(XYZmm));

for k=1:l
    hi=XYZmm(:,k)+voxsize(:);
    lo=XYZmm(:,k)-voxsize(:);
    lobound(1:3)=0;hibound(1:3)=0;
    for i=1:3
        if lo(i)>0
            while lobound(i)<=lo(i)
                lobound(i)=lobound(i)+voxsize(i);
            end
            lobound(i)=lobound(i)-voxsize(i);
        else
            while lobound(i)>=lo(i)
                lobound(i)=lobound(i)-voxsize(i);
            end
        end
        if hi(i)>0
            while hibound(i)<=hi(i)
                hibound(i)=hibound(i)+voxsize(i);
            end
        else
            while hibound(i)>=hi(i)
                hibound(i)=hibound(i)-voxsize(i);
            end
            hibound(i)=hibound(i)+voxsize(i);
        end
    end
    space=zeros([3 length(lobound(1):voxsize(1):hibound(1))*length(lobound(2):voxsize(2):hibound(2))*length(lobound(3):voxsize(3):hibound(3))]);
    count=0;
    for m=lobound(1):voxsize(1):hibound(1)
        for n=lobound(2):voxsize(2):hibound(2)
            for p=lobound(3):voxsize(3):hibound(3)
                count=count+1;
                space(:,count)=[m;n;p];
            end
        end
    end
    dist=zeros([1 size(space,2)]);
    for i=1:size(space,2)
        dist(i)=sqrt((space(1,i)-XYZmm(1,k))^2+(space(2,i)-XYZmm(2,k))^2+(space(3,i)-XYZmm(3,k))^2);
    end
    a=find(dist==min(dist));
    if isscalar(a)
        XYZmm_out(:,k)=space(:,a);
    elseif force
        b=find(sum(abs(space(:,a)))==max(sum(abs(space(:,a)))));
        XYZmm_out(:,k)=space(:,a(b));
    else
        r=randperm(length(a));
        XYZmm_out(:,k)=space(:,a(r(1)));
        disp(['Warning: ' num2str(XYZmm(1,k)) ',' num2str(XYZmm(2,k)) ',' num2str(XYZmm(3,k)) 'is equidistant from at least 2 points. It has been set to ' num2str(XYZmm_out(1,k)) ',' num2str(XYZmm_out(2,k)) ',' num2str(XYZmm_out(3,l)) '.'])
    end
end
% XYZmm_out=unique(XYZmm_out','rows')'; commented because it can cause
% XYZmm_out to be of different length than XYZmm, which can create problem
% outside of the function, such as when there are Z values corresponding to
% the elements  of XYZmm that also need to be matched to XYZmm_out. This
% should be handled in calling functions.
    