%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function msgModalOnTop(txt, where, javaWnd, icon, title)
if nargin<5
    title='Note';
    if nargin<4
        icon='facs.gif';
        if nargin<3
            javaWnd=Gui.JFrame;
            if nargin<2
                where='center';
            end
        end
    end
end
pane=javaObjectEDT('javax.swing.JOptionPane', txt, 1);
jd=javaMethodEDT('createDialog', pane, javaWnd, title);
jd.setResizable(true);
jd.setModal(true);
if ~isempty(javaWnd)
    javaMethodEDT( 'setLocationRelativeTo', jd, javaWnd);
end
Gui.LocateJava(jd, javaWnd, where);
jd.setAlwaysOnTop(true);
pane.setIcon(Gui.Icon(icon));
jd.pack;
Gui.SetJavaVisible(jd);
end

