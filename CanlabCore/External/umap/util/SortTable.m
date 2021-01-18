%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef SortTable < handle
    properties(SetAccess=private)
        uit;
        jscrollpane;
        jtable;
        fig;
        tips;
        uil;
    end
    methods
        function this=SortTable(fig, data, columnNames, ...
                normalizedPosition, fncSelect, tips)
            if nargin<6
                tips=[];
                if nargin<5
                    fncSelect=[];
                    if nargin<4
                        normalizedPosition=[];
                        if nargin<3
                            columnNames=[];
                            if nargin<2
                                data={...
                                    'Pepper', 225, 14, 'Golden retriever';...
                                    'Fergie-roo', 44, 4, 'Golden retriever'; ...
                                    'Killarney', 4, 100, 'Black lab'};
                                columnNames={'Top dog', 'IQ', 'Age', 'Breed'};
                                if nargin<1
                                    fig=figure;
                                end
                            end
                        end
                    end
                end
            end
            if ~isempty(tips)
                this.tips=tips;
                this.uil=uicontrol('style', 'text', 'parent', fig, ...
                    'units', 'Normalized', 'Position', [.01 .01 .97 .04]);
            end
            if isempty(normalizedPosition)
                if isempty(tips)
                    y=.04;
                    height=.94;
                else
                    y=.06;
                    height=.92;
                end
                normalizedPosition=[.03 y .94 height];
            end
            on=get(fig, 'visible');
            uit = uitable(fig, 'Data',data,'units', 'normalized', ...
                'Position', normalizedPosition, 'RowName', []);
            set(uit, 'fontName', 'arial');
            set(fig, 'visible', on);
            if ~isempty(columnNames)
                uit.ColumnName=columnNames;
            end
            this.uit=uit;
            this.fig=fig;
            [this.jtable, this.jscrollpane]=SortTable.Go(uit, fig);
            if ~isempty(fncSelect)
                set(uit, 'CellSelectionCallback', fncSelect);
                set(uit, 'RowStriping', 'off');
            end
            set(fig, 'visible', on);
            if ismac
                pause(.25);
            end
            
        end
        
        function showTip(this, col)
            if col<=length(this.tips)
                if col==0
                    set(this.uil, 'String', '');
                else
                    vi=this.jtable.convertColumnIndexToView(col-1);
                    hdr=this.jtable.getColumnName(vi);
                    hdr=hdr.replaceAll('-<br>', '');
                    hdr=hdr.replaceAll('<br>', ' ');
                    hdr=char(edu.stanford.facs.swing.Basics.RemoveXml(hdr).trim);
                    set(this.uil, 'String', ['"' hdr '": ' this.tips{col}]);
                end
            end
        end
        
        function setColumnWidth(this,col, width)
            tcm=this.jtable.getColumnModel;
            t=tcm.getColumn(col-1);
            t.setPreferredWidth(width*8);
        end
    end
    
    methods(Static)
        function [out, widths]=Convert(in, fmt)
            [R, C]=size(in);
            out=cell(R,C);
            widths=zeros(1,C);
            for row=1:R
                for col=1:C
                    col_=fmt(col,2);
                    if isnan(col_)
                        widths(col)=fmt(col, 1);
                        out{row,col}=in{row,col};
                    elseif col_<0
                        widths(col)=fmt(col,1)+(0-col_)+2;
                        out{row,col}=String.encodePadHtml(...
                            in{row,col}/1*100, fmt(col,1), 0-col_, '%');
                    elseif col_==0
                        widths(col)=fmt(col,1);
                        out{row,col}=String.encodePadHtml(in{row,col},...
                            fmt(col, 1), 0);
                    else
                        widths(col)=fmt(col,1)+col_+1;
                        out{row,col}=String.encodePadHtml(in{row,col},...
                            fmt(col, 1), col_);
                    end
                end
            end            
        end
        
        function [jtable, jscrollpane]=Go(uit, fig)
            app=BasicMap.Global;
            jscrollpane = javaObjectEDT(findjobj_fast(uit));
            jtable = javaObjectEDT(jscrollpane.getViewport.getView);
            % Now turn the JIDE sorting on
            jtable.setSortable(true);
            jtable.setAutoResort(true);
            jtable.setMultiColumnSortable(true);
            jtable.setPreserveSelectionsAfterSorting(true);
            if app.highDef
                jtable.setRowHeight(35)
                jtable.setIntercellSpacing(java.awt.Dimension(10, 6))
            else
                jtable.setRowHeight(25)
                jtable.setIntercellSpacing(java.awt.Dimension(5,3))
            end
            if false
                filter=net.coderazzi.filters.gui.TableFilterHeader(jtable);
                filter.setAutoChoices(...
                    net.coderazzi.filters.gui.AutoChoices.ENABLED);
                if app.highDef
                    f=filter.getFont;
                    filter.setFont(java.awt.Font('arial', 0, 14));
                    filter.setRowHeightDelta(16)
                end
            end
            jtable.getTableHeader.setReorderingAllowed(true);
            
        end
        
        function rowIdxs=ModelRows(jtable)
            rows=jtable.getRowCount;
            rowIdxs=zeros(1,rows);
            for row=0:rows-1
                rowIdxs(row+1)=1+jtable.getActualRowAt(jtable.convertRowIndexToModel(row));
            end
        end
        
        function colIdxs=ModelCols(jtable)
            cols=jtable.getColumnCount;
            colIdxs=zeros(1,cols);
            for row=0:cols-1
                colIdxs(row+1)=1+jtable.convertColumnIndexToModel(row);
            end
        end
        
        function html=ToHtml(jtable)
            rowClr='"#FFFFDD"';
            colClr='"#AABDFF"';
            html='<table><thead>';
            R=jtable.getRowCount;
            C=jtable.getColumnCount;
            for j=1:C
                v=Html.remove(char(jtable.getColumnName(j-1)));
                html=[html '<th bgcolor=' colClr '>' v '</th>'];
            end
            html=sprintf('%s </thead>\n', html);
            for i=1:R
                html=[html '<tr>'];
                for j=1:C
                    v=jtable.getValueAt(i-1,j-1);
                    if ischar(v)
                        v=Html.remove(v);
                    else
                        v=num2str(v);
                    end
                    if j>1
                        html=[html '<td>' v '</td>'];
                    else
                        html=[html '<td bgcolor=' rowClr '>' v '</td>'];
                    end
                end
                html=sprintf('%s</tr>\n', html);
            end
            html=[html '</table>'];
        end

        function [order, widths, changed]=GetColumnOrder(jt)
            C=jt.getColumnCount;
            order=zeros(1,C);
            widths=zeros(1,C);
            tcm=jt.getColumnModel;
            for c=0:C-1
                mi=jt.convertColumnIndexToModel(c);
                order(c+1)=mi;
                widths(c+1)=tcm.getColumn(c).getPreferredWidth;
            end
            naturalOrder=1:C;
            naturalOrder=naturalOrder-1;
            changed=~isequal(order, naturalOrder);
        end
        
        
        function SetColumnOrder(jt, order, widths, force)
            if nargin<4
                force=false;
            end
            C=length(order);
            naturalOrder=1:C;
            naturalOrder=naturalOrder-1;
            tcm=jt.getColumnModel;
            if force || ~isequal(order, naturalOrder)
                for c=0:C-1
                    o=order(c+1);
                    v=jt.convertColumnIndexToView(o);
                    jt.moveColumn(v, c);
                end
            end
            if ~isempty(widths)
                for c=0:C-1
                    tcm.getColumn(c).setPreferredWidth(...
                        widths(c+1));
                end
            end

        end
        
        function order=GetRowOrder(jt)
            m=jt.getModel;
            C=jt.getColumnCount;
            order=[];
            for c=1:C
                if m.isColumnSorted(c)
                    order(end+1,:)=[c ...
                        m.isColumnAscending(c) m.getColumnSortRank(c)];
                end
            end
        end
        
        function [ok, C]=SetRowOrder(jt, order)
            [C,~]=size(order);
            if C<1
                ok=false;
                return;
            end
            drawnow;
            jt.unsort;
            ok=true;
            m=jt.getModel;
            [~,I]=sort(order(:,3));
            for c=1:C
                idx=I(c);
                jt.sortColumn(order(idx,1), c==1, order(idx,2));
            end
        end
        
        function ok=MoveSortColumnsLeft(jt)
            order=SortTable.GetRowOrder(jt);
            [C,~]=size(order);
            if C<1
                ok=false;
                return;
            end
            ok=true;
            m=jt.getModel;
            [~,I]=sort(order(:,3));
            for c=C:-1:1
                idx=I(c);
                o=order(idx, 1);
                v=jt.convertColumnIndexToView(o);
                jt.moveColumn(v, 0);
            end
        end

        function idxs=SelectedRows(jt, orUnselectedToo)
            if nargin<2
                orUnselectedToo=false;
            end
            rows=jt.getSelectedRows;
            N=length(rows);
            if N==0
                if orUnselectedToo
                    N=jt.getRowCount;
                    rows=0:N-1;
                else
                    idxs=[];
                    return;
                end
            end
            idxs=zeros(1, N);
            for i=1:N
                idx=rows(i);
                idxs(i)=jt.getActualRowAt(jt.convertRowIndexToModel(idx))+1;
            end
        end
        
    end
    
end