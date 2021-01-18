%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Funded by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef UmapPython
    properties(Constant)
        PROP_CMD='pythonCmdV3';
        PROP_DESIRED='umapDesired';
        TRY_PYTHON=false;
        PYTHON_TEMPLATE='lastUmap.python';
    end
    
    methods(Static)
        
        function [outData, cmdOut]=Go(data, inFile, this, metric, ...
            neighbors, minDist, nComponents, lbls, template, verbose)
            cmdOut=[];
            app=BasicMap.Global;
            if ~UmapPython.IsAvailable(app)
                outData=[];
                msg('Python not available...', 8, 'north');
                return;
            end
            if nargin<9
                template=[];
                if nargin<8
                    if nargin<7
                        nComponents=2;
                        lbls=[];
                        if nargin<6
                            minDist=.3;
                            if nargin<5
                                neighbors=30;
                                if nargin<4
                                    metric='manhattan';
                                    if nargin<3
                                        this=[];
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if strcmpi('cityblock', metric)
                metric='manhattan';
            end
            File.SaveMatrix(inFile, data);
            [inPath, inFileName]=fileparts(inFile);
            if ~isempty(template)
                args=sprintf(['--n_neighbors=%d --metric=%s --min_dist %d '...
                    '--firstRow 0 --verbose '...
                    ' --output_dimensions %d'], neighbors, metric, minDist, ...
                    nComponents);
                
                args=[args ' --useTemplate ' String.ToSystem(template)];
            else
                tempTemplate=String.ToSystem(...
                    fullfile(inPath, UmapPython.PYTHON_TEMPLATE));
                
                args=sprintf(['--n_neighbors=%d --metric=%s --min_dist %d '...
                    '--firstRow 0 --verbose --saveTemplate %s'...
                    ' --output_dimensions %d'], neighbors, metric, minDist, ...
                    tempTemplate, nComponents);
                
                if ~isempty(lbls)
                    lblFile=fullfile(inPath, [inFileName '.labels']);
                    if size(lbls, 2)>1
                        lbls=lbls';
                    end
                    File.SaveMatrix(lblFile, lbls);
                    args=[args ' --labels ' String.ToSystem(lblFile)];
                end
            end
            cmd=app.get(UmapPython.PROP_CMD, 'python');
            cmdline=String.ToSystem(cmd);
            pyFilePath=fileparts(mfilename('fullpath'));
            pythonScript=String.ToSystem(fullfile(pyFilePath, ...
                'doUmap.py'));
            pathArg=String.ToSystem(inFile);
            fullCmd=[cmdline ' ' pythonScript ' ' pathArg ' ' args];
            [p,f]=fileparts(inFile);
            outFile=fullfile(p,[f '.umap.csv']);
            if exist(outFile, 'file')
                delete(outFile);
            end
            fldr=fileparts(inFile);
            runInBackground=nargin<10 || ~strcmp(verbose, 'none');
            if ~runInBackground
                [RR, CC]=size(data);
                terminalName=['  >>> Stand by ... the Python UMAP of '...
                    'Leland McInnes is running on ' ...
                    String.encodeInteger(RR) ' X ' ...
                    String.encodeInteger(CC) ' values'];
            else
                terminalName=['AutoGate UMAP ' datestr(datetime)];
            end
            script=fullfile(fldr, 'umap.cmd');
            [status, stdout]=File.Spawn(fullCmd, script, terminalName, ...
                runInBackground);
            if runInBackground
                UmapPython.Wait(this, outFile);
            end
            cmdOut=strtrim(stdout);
            outData=File.ReadMatrix(outFile);
            delete(inFile);
            delete(outFile);
            if ~isempty(lbls)
                delete(lblFile);
            end
            if status==0
            end
        end
        
        function ok=IsDesired(app, interactWithUser)
            if nargin<1
                app=BasicMap.Global;
            end
            ok=app.is(UmapPython.PROP_DESIRED, false);
            if ~ok
                if app.needToAskForUmap && interactWithUser
                    [yes, cancelled]=askYesOrNo(Html.WrapHr(...
                        'Use original python implementation of UMAP?'));
                    if yes
                        app.set(UmapPython.PROP_DESIRED, 'true');
                        ok=true;
                    elseif ~cancelled
                        app.needToAskForUmap=false;
                    end
                end
            end
            %ok=false;
        end
        
       function [ok, version]=IsAvailable(app, interactWithUser)
            if nargin<2
                interactWithUser=true;
                if nargin<1
                    app=BasicMap.Global;
                end
            end
            if ~UmapPython.IsDesired(app, interactWithUser)
                ok=false;
                version=[];
                return;
            end
            cmd=app.get(UmapPython.PROP_CMD, 'python');
            if isempty(cmd)
                cmd='Python';
            end
            cmdline=String.ToSystem(cmd);
            [status, version]=system([cmdline ' -V']);
            ok=status==0;
            firstCmdOk=ok;
            if ok
                ok=String.StartsWithI(version, 'python 3');
                if ok
                    curPath=fileparts(mfilename('fullpath'));
                    pythonScript=String.ToSystem(...
                        fullfile(curPath, 'testImports.py'));
                    [status,stdout]=system([cmdline ' ' pythonScript]);
                    ok=status==0;
                    if ~ok
                        app.needToAskForUmap=true;
                        if interactWithUser
                            if ispc
                                cmdApp='Windows "cmd" window';
                            else
                                cmdApp='Mac''s "terminal"';
                            end
                            html=Html.WrapHr([...
                                'UMAP requires the python packages:'...
                                '  <b>pandas</b> & <b>umap-learn</b>']);
                            choice=Gui.Ask(html,...
                                {'Try automatic download & install', ...
                                ['Open ' cmdApp 'to install myself']}, ...
                                'umapInstall', 'UMAP packages needed');
                            
                            [fldr,python]=fileparts(cmd);
                            if ~isempty(fldr)
                                if strcmpi('python3', python)
                                    pipCmd=fullfile(fldr, 'pip3');
                                    if ~exist(pipCmd, 'file')
                                        pipCmd=fullfile(fldr, 'pip');
                                    end
                                else
                                    pipCmd=fullfile(fldr, 'pip');
                                end
                                if ~exist(pipCmd, 'file')
                                    pipCmd='pip';
                                end
                            else
                                pipCmd='pip';
                            end
                            if choice==1
                                cmds={[pipCmd ' install pandas '], ...
                                    [pipCmd ' install umap-learn ']};
                                File.Spawn(cmds, ...
                                    fullfile(app.appFolder, 'pipUmap.cmd'),...
                                    ['AutoGate is installing UMAP ' ...
                                    datestr(datetime)], false, true);
                                [ok, version]=UmapPython.IsAvailable(app, true);
                            elseif choice==2
                                msg(Html.Wrap(['The python packages '...
                                    'which UMAP needs are:'...
                                    Html.ToList({'pandas', 'umap-learn'}) ...
                                    'From the the command line type <br>'...
                                    '<br><i> ' pipCmd ' install "'...
                                    '<b>package name</b>"'...
                                    '</i><br><br>for <b>each</b> of the '...
                                    'packages in the list above.']));
                                if ispc
                                    system('cmd');
                                else
                                    system('open -b com.apple.terminal');
                                end
                            end
                        end
                    end
                    return;
                end
            end
            if interactWithUser
                bp=Gui.BorderPanel;
                if isempty(version)
                    problem='';
                else
                    if firstCmdOk
                        problem=['<br>(Incorrect version "' version '" was found)'];
                    else
                        problem=['<br>("' cmdline '" returned "'...
                            version '")'];
                    end
                end
                if ispc
                    cmdApp='Command prompt';
                    whichCmd='where';
                else
                    
                    cmdApp='Terminal app';
                    whichCmd='which';
                end
                lbl=Gui.Label(['<html><font color="red">Python version '...
                    '<b>3.x</b> is <b>required</b></font>.' ...
                    app.smallStart problem app.smallEnd ...
                    '<br><br>Please enter the path/location '...
                    'where python <br>version <b>3.x</b> is installed '...
                    'on your computer...<br><br>' app.smallStart ...
                    '<b>NOTE:  </b>To find a python installation, open <b>' ...
                    cmdApp '</b><br>and type commands like "' whichCmd ...
                    ' python" or "' whichCmd ' python3" etc.'...
                    '<hr><br></html>']);
                btn=Gui.NewBtn(['<html>' app.smallStart 'Download python?' ...
                    app.smallEnd '</html>'], @(h,e)download());
                bp.add(lbl, 'Center');
                bp2=Gui.BorderPanel;
                bp2.add(btn, 'East');
                bp.add(bp2, 'North');
                cmd=inputDlg(struct('msg', ...
                    bp, 'where', 'North'),...
                    'UMAP needs Python version 3...', cmd);
                if ~isempty(cmd)
                    app.set(UmapPython.PROP_CMD, cmd);
                    [ok, version]=UmapPython.IsAvailable(app, true);
                else
                    if askYesOrNo(struct('msg', ...
                            Html.WrapHr(['Open the download page'...
                            ' for <br>python version 3.x in your'...
                            ' browser?']),'where', 'North'))
                        web('https://www.python.org/downloads/', ...
                            '-browser')
                    end
                end
            end
            
           function download
               web('https://www.python.org/downloads/', '-browser');
           end
        end
        
        function Wait(this, outFile)
            prefix=['<html><center><b>Running the Python UMAP <br>'...
                'implementation of Leland McInnes</b><hr><br>'];
            progress='(<i>see python progress in shell window</i>)';
            if  isempty(this) || ishandle(this) 
                fig=this;
                btn=[];
                html=[prefix progress ];
            else
                html=[prefix String.RemoveTex(this.focusTitle) '<br>'...
                    '<font color="blue">' this.sizeTitle ...
                    '</font><br><br>' progress ];
                btn=this.btn;
                fig=this.h;
            end
            html=[html '</center></html>'];
            File.Wait(outFile, fig, btn, html);
        end
        
    end
end