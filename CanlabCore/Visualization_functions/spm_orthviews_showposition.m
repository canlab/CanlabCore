function h = spm_orthviews_showposition()
% :Usage:
% ::
%
%    handles = spm_orthviews_showposition;
%
% plots x, y, and z coordinates on SPM orthviews figure
%
% bring context structure variable into this script.
% contains handles, etc.
%
% ..
%    tor wager, august 2006
%    updated: june 2007
% ..

    global st


    % get current position in mm
    pos = round(spm_orthviews('Pos'));  % x y z

    for i = 1:length(st.vols)

        if ~isempty(st.vols{i}) && ~isempty(st.vols{i}.ax) && ~isempty(st.vols{i}.ax{1})

            axish(1) = st.vols{i}.ax{1}.ax;
            axish(2) = st.vols{i}.ax{2}.ax;
            axish(3) = st.vols{i}.ax{3}.ax;

            % axial
            h(1) = put_text(axish(1),['z = ' num2str(pos(3))]);

            % coronal
            h(2) = put_text(axish(2),['y = ' num2str(pos(2))]);

            % saggital
            h(3) = put_text(axish(3),['x = ' num2str(pos(1))]);

        end

        % %
        % % % Get figure and axis handles
        % % fh = findobj('Tag','Graphics');
        % % ch = get(fh,'Children');
        % % for i= 1:length(ch),mytype = get(ch(i),'Type'); wh(i)=strcmp(mytype,'axes'); end
        % %
        % % if ~exist('wh','var')
        % %     disp('spm_orthviews_showposition: No axis handles found. Exiting.');
        % %     return
        % % end
        % %
        % % axish = ch(find(wh));
        % %
        % % if isempty(axish), return, end

        % % % get which axis is which
        % % for i = 1:length(axish)
        % % poss(i,:) = get(axish(i),'Position');
        % % end
        % %
        % % % get rid of extra axes we may have created in the 4th quadrant
        % % other_axes = find(any(poss(:,1) > .45 & poss(:,2) < .2,2));
        % % axish(other_axes) = [];
        % % poss(other_axes,:) = [];
        % %
        % % % sort into order:  axial, coronal, saggital
        % % ssum = sum(poss(:,1:2),2);
        % % [ssum,ind] = sort(ssum);
        % % axish = axish(ind);
        % %
        % %
        % % % axial
        % % h(1) = put_text(axish(1),['z = ' num2str(pos(3))]);
        % %
        % % % coronal
        % % h(2) = put_text(axish(2),['y = ' num2str(pos(2))]);
        % %
        % % % saggital
        % % h(3) = put_text(axish(3),['x = ' num2str(pos(1))]);



    end  % end loop thru vols



    return



function h = put_text(axish,str)

    axes(axish)

    % remove old text
    oldh = findobj(gca,'Type','text');
    delete(oldh)

    xdelta = diff(get(axish,'XLim'));
    ydelta = diff(get(axish,'YLim'));

    % add new text
    xpos = xdelta .* .42;
    ypos = ydelta .* .1;
    h1 = text(xpos - xdelta * .01,ypos + ydelta * .005,str,'Color','k','FontSize',24, 'FontWeight', 'b');
    h = text(xpos,ypos,str,'Color','g','FontSize',24, 'FontWeight', 'demi');

    return
