%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

% this wrapper around MatLab internal map is meant for 
% polymorphic function signature compatability with other
% classdef and JAVA class in CytoGenie AutoGate
classdef BasicMap < Map
    properties
        pu;
        parentCmpForPopup;
        needToAskForUmap=true;
        noDbscan=false;
        urlMap;
    end
    
    properties(SetAccess=private)
        appFolder;
        contentFolder;
        toolBarSize=0;
        toolBarFactor=0;
        highDef=false;
        supStart='<sup>';
        supEnd='</sup>';
        subStart='<sub>';
        subEnd='</sub>';
        smallStart='<small>';
        smallEnd='</small>';
        h3Start='<h3>';
        h3End='</h3>';
        h2Start='<h2>';
        h2End='</h2>';
        h1Start='<h1>';
        h1End='</h1>';
        whereMsgOrQuestion='center';
    end
    
    methods(Static)

        
        function this=Global(closeIfTrueOrMap)
            persistent singleton;
            if nargin>0 && islogical(closeIfTrueOrMap) && closeIfTrueOrMap
                clear singleton;
                singleton=[];
                disp('Resetting global BasicMap');
                this=[];
            else
                if nargin==1
                    try
                        priorMap=closeIfTrueOrMap;
                        %test method compatibility
                        priorMap.size
                        prop='IsBasicMap';
                        if ~priorMap.has(prop)
                            priorMap.set(prop, 'false')
                            priorMap.get(prop, 'true')
                            priorMap.remove(prop);
                        else
                            priorMap.get(prop, 'true')
                        end
                        singleton=priorMap;
                        this=singleton;
                        return;
                    catch ex
                        ex.getReport
                    end
                end
                if isempty(singleton) 
                    singleton=BasicMap;
                    singleton.highDef=Gui.hasPcHighDefinitionAProblem(2000, 2500, false);
                    BasicMap.SetHighDef(singleton, singleton.highDef);
                end     
                this=singleton;
            end
        end
        
        function path=Path
            path=BasicMap.Global.contentFolder;
        end
        
        function obj=SetHighDef(obj, hasHighDef)
            factor=0;
            NORMAL_FONT_SIZE=12;
            SMALL_FONT_SIZE=2;
            H3_FONT_SIZE=3;
            H2_FONT_SIZE=3.5;
            H1_FONT_SIZE=4;
            if hasHighDef
                obj.highDef=true;
                factor=javax.swing.UIManager.getFont('Label.font').getSize...
                    /NORMAL_FONT_SIZE;
            else
                if ismac
                    %factor=1.6;
                end
            end
            if factor>0
                obj.toolBarFactor=factor;
                obj.toolBarSize=floor(16*factor);
                smallSize=floor(SMALL_FONT_SIZE*factor);
                if ispc
                    smallSize=smallSize+1;
                end
                obj.smallStart=['<font size="' num2str(smallSize) '">'];
                obj.smallEnd='</font>';
                obj.subStart=obj.smallStart;
                obj.supStart=obj.smallStart;
                obj.subEnd=obj.smallEnd;
                obj.supEnd=obj.smallEnd;
                h1Size=floor(H1_FONT_SIZE *factor);
                if ispc
                    h1Size=h1Size+1;
                end
                obj.h1Start=['<center><font size="' num2str(h1Size) ...
                    '" color="blue"><b>'];
                obj.h1End='</b></font></center><br>';

                h2Size=floor(H2_FONT_SIZE *factor);
                if ispc
                    h2Size=h2Size+1;
                end
                obj.h2Start=['<center><font size="' num2str(h2Size) ...
                    '" color="blue"><b>'];
                obj.h2End='</b></font></center><br>';

                h3Size=floor(H3_FONT_SIZE *factor);
                if ispc
                    h3Size=h3Size+1;
                end
                obj.h3Start='<h1>';
                obj.h3End='</h1>';
            else
                obj.toolBarSize=0;
                obj.toolBarFactor=0;
                obj.highDef=false;
        
                obj.smallStart='<small>';
                obj.smallEnd='</small>';
                obj.h3Start='<h3>';
                obj.h3End='</h3>';
                obj.h2Start='<h2>';
                obj.h2End='</h2>';
                obj.h1Start='<h1>';
                obj.h1End='</h1>';

            end
        end
        
        function nums=GetNumbers(props, name)
            nums=props.get(name, []);
            if ~isempty(nums)
                nums=str2num(nums);
            else
                nums=[];
            end
        end
        
    end
    
    methods
        function this=BasicMap(keysOrFileName, values)
            if nargin<2
                values={};
                if nargin<1
                    keysOrFileName=[];
                end
            end
            this=this@Map(keysOrFileName, values);
            this.contentFolder=fileparts(mfilename('fullpath'));
            this.appFolder=fullfile(File.Home, '.run_umap');
            File.mkDir(this.appFolder);
            this.urlMap=Map;
            %'https://1drv.ms/u/s!AkbNI8Wap-7_jNJYg4RNDTKR4mkYOg?e=pfwWfO'
        end
        
        function closeToolTip(this)
        end
        
                
         
        function retry=reportProblem(this, exception, retryId)
            retry=false;
            exception.getReport;
        end
        
        function [name, found]=getMatchName(this, name, showEppImg, b1, b2)
            if nargin<4
                if this.highDef
                    b1='<b>';
                    b2='</b>';
                else
                    b1='<sup>';
                    b2='</sup>';
                end
                if nargin<3
                    if nargin<2
                        showEppImg=true;
                    end
                end
            end
            idx=strfind(name, ',\bf');
            if ~isempty(idx)
                sIdx=4;
                idx=idx(1);
            else
                idx=strfind(name, '\bf');
                sIdx=3;
            end
            if ~isempty(idx)
                matchName=strtrim(name(idx(1)+sIdx:end));
                eppName=name(1:idx(1)-1);
                if showEppImg
                    eppName =[ eppName ' ' this.eppImg];
                end
                if idx==1
                    name=[matchName eppName];
                else
                    name=[matchName ' ' b1 eppName b2];
                end
                found=true;
            else
                found=false;
            end
        end
        
    end
end