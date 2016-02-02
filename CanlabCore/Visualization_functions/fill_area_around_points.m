function h = fill_area_around_points(x, y, borderscale, color)
% :Usage:
% ::
%
%    h = fill_area_around_points(x, y, borderscale, color)
%
% fills area around list of coordinates in a color
% using spline interpolation and other stuff.
% designed for cluster imaging in nmdsfig figures.
%
% :Example:
% ::
%
%    x = randn(3,1);
%    y = randn(3,1);
%    figure; plot(x,y,'k.');
%    h = fill_area_around_points(x, y, .2, 'r');
%
% ..
%    tor wager, feb 07
%    tor recommends .2 for borderscale...
% ..

if length(x) == 1 & length(y) == 1
    % only one point; exit
    h = [];
    return
end

    % bootstrap means of x and y to get values in area
    % ---------------------------------------------------
    c = [x y];

    cm = bootstrp(2000,@mean,c);

    border = borderscale * mean(range(c));  % x% of range

    % add random noise to get coverage
    c1 = c(:,1) + border .* (rand(size(c,1),1) - .5);
    c2 = c(:,2) + border .* (rand(size(c,1),1) - .5);
    c = [c; [c1 c2]];

    cm = [cm; bootstrp(500,@mean,c)];
    cm = [cm; c];

    % plot(cm(:,1),cm(:,2),'b.');

    sampleres = border ./ 5;
    r1 = min(cm(:,1)) - 2 * border:sampleres:max(cm(:,1)) + 2 * border;
    r2 =  min(cm(:,2)) - 2 * border:sampleres:max(cm(:,2)) + 2 * border;

    % make grid
    % ---------------------------------------------------
    [X,Y] = meshgrid(r1, r2);

    xy = combvec(r1,r2);  % across x and then up y

    len = length(cm);

    Z = zeros(size(X)); % values for sig.

    % Find values within border
    % ---------------------------------------------------
    b2 = border ^ 2;
    for i = 1:len % for each point, find xy combos that are within radius
        dsquared = (cm(i,1) - X) .^2 + (cm(i,2) - Y) .^2;
        Z(dsquared <= b2) = 1;
    end



    % Plot initial contours (rough)
    % ---------------------------------------------------
    hold on;
    % %     [cout,hand] = contourf(X,Y,Z,[0 1]);
    % %     delete(hand);
    % %     cout = cout(:, 5:end);

    cout = contourc(r1,r2,Z,1);
    cout = cout(:,2:end);

    % % hand = hand(2);
    % % set(hand,'FaceColor','y')

    % Reduce and use spline interpolation
    % ---------------------------------------------------
    desiredpts = 10;
    
    wh = round(linspace(1,size(cout,2),desiredpts)); % which to save for spline
    
    cout = cout(:,wh);
    
% %     skipevery = ceil((size(cout,2) - 4) ./ desiredpts);
% % 
% %     % %    cout = cout(:,all(~isnan(cout)));
% %     cout = cout(:,1:skipevery:end);
    
    % pad
    % cool effects but artsy
    % cout = [cout(:,[end-3:end]) cout cout(:,1:3)]; cout = [cout cout(:,1)];
    
    % pad
    cout = [cout(:,end-2:end-1) cout cout(:,[2 3])];

    
    n = size(cout,2);
    t = 1:n;
    ts = 3:1/desiredpts:n-2;    % get rid of padding
    xs = spline(t,cout(1,:),ts);
    ys = spline(t,cout(2,:),ts);
    hold on
    hh = plot(xs,ys,'k');
    hold on

    h = fill(xs,ys,color,'FaceAlpha',.2);
    h = [h hh];

        % figure; imagesc(Z);
end
