function make_predictive_model_codemap(outpng)
% make_predictive_model_codemap  Render a CANlab-style code map for the
% sklearn-style @predictive_model API and save it as a PNG.
%
% Visual style follows docs/canlab_template_codemap.pptx: a title + short
% description, then a flow diagram of object types (rounded boxes) wired by
% arrows, with a colour legend (input/data, construct, fit & CV, inference
% & tuning, visualize).
%
% Usage:  make_predictive_model_codemap            % writes next to this file
%         make_predictive_model_codemap('foo.png')

    if nargin < 1 || isempty(outpng)
        here = fileparts(mfilename('fullpath'));
        outpng = fullfile(here, 'pngs', 'predictive_model_api_codemap.png');
    end

    % --- palette (by role) ---
    C.data    = [0.86 0.90 0.98];   % input / key variables
    C.constr  = [0.82 0.92 0.82];   % construct
    C.fit     = [0.99 0.92 0.74];   % fit & cross-validate
    C.infer   = [0.97 0.85 0.80];   % inference & tuning
    C.viz     = [0.90 0.85 0.96];   % visualize
    edge      = [0.30 0.30 0.30];

    f = figure('Color','w','Position',[80 80 1280 760],'Visible','on');
    ax = axes('Position',[0 0 1 1]); axis(ax,[0 100 0 100]); axis(ax,'off'); hold(ax,'on');

    % --- title + description ---
    text(3, 96, 'CANlab code map: the sklearn-style {\bf@predictive\_model} API', ...
        'FontSize', 17, 'FontWeight','bold');
    text(3, 92.0, ['Construct a model from hyperparameters, cross-validate it on (X, Y, groups), ' ...
        'add inference/tuning, then read out predictions, a weight map, and figures.'], ...
        'FontSize', 11);
    text(3, 88.8, 'Value semantics: every method returns a NEW object — write  pm = crossval(pm, X, Y, ...).', ...
        'FontSize', 10, 'FontAngle','italic', 'Color',[0.35 0.35 0.35]);

    % ---------- Lane 1: construct -> cross-validate ----------
    bx = @(x,y,w,h,s,c) draw_box(ax, x,y,w,h,s,c,edge);
    ar = @(x1,y1,x2,y2) draw_arrow(ax, x1,y1,x2,y2, edge);

    b_img = bx( 3.0, 74, 19, 10, {'fmri\_data / image\_vector','.dat  .Y  metadata\_table'}, C.data);
    b_xy  = bx(26.5, 74.5, 12, 9, {'X, Y, groups','(numeric)'}, C.data);
    b_pm  = bx(42.5, 74, 24, 10, {'predictive\_model(', '  ''algorithm'',''svm'', ''task'',...)'}, C.constr);
    b_cv  = bx(70.5, 74, 26.5, 10, {'crossval(pm, X, Y,', '  ''groups'',id, ''cv'',cv\_splitter...)'}, C.fit);

    ar(22.0, 79, 26.4, 79);
    ar(38.6, 79, 42.4, 79);
    ar(66.6, 79, 70.4, 79);
    text(33.5, 71.0, 'extract', 'FontSize',8,'Color',[0.4 0.4 0.4],'HorizontalAlignment','center');

    % fit() alternative (in-sample) under construct->cv
    b_fit = bx(42.5, 60, 24, 8, {'fit(pm, X, Y)','in-sample model'}, C.fit);
    ar(54.5, 74, 54.5, 68.2);
    text(56.0, 71.0, 'or', 'FontSize',8,'Color',[0.4 0.4 0.4]);

    % ---------- Lane 2: inference & tuning (fan out from crossval/pm) ----------
    text(3, 53.5, 'Inference & tuning  (each returns an updated pm)', 'FontSize',11,'FontWeight','bold');
    iy = 41; ih = 9;
    b_boot = bx( 3.0, iy, 22, ih, {'bootstrap(pm, X, Y)','weights.z .p .fdr\_sig'}, C.infer);
    b_perm = bx(26.5, iy, 22, ih, {'permutation\_test(pm, X, Y)','permutation\_results.p\_value'}, C.infer);
    b_grid = bx(50.0, iy, 22, ih, {'grid\_search(pm,X,Y,grid)','nested-CV tuning'}, C.infer);
    b_stab = bx(73.5, iy, 23.5, ih, {'stability\_selection(pm,X,Y)','diagnostics.stability\_selection'}, C.infer);
    b_cal  = bx(50.0, iy-12, 22, ih, {'calibrate(pm, X, Y)','-> predict\_proba(pm, Xnew)'}, C.infer);

    % connector from crossval down to inference band
    ar(83.5, 74, 83.5, 50.5);
    ar(83.5, 50.5, 14, 50.5); ar(14, 50.5, 14, iy+ih+0.2);
    ar(37.5, 50.5, 37.5, iy+ih+0.2);
    ar(61.0, 50.5, 61.0, iy+ih+0.2);
    ar(85.2, 50.5, 85.2, iy+ih+0.2);
    ar(61.0, iy, 61.0, iy-12+ih+0.2);

    % ---------- Lane 3: outputs / visualize ----------
    text(3, 21.0, 'Predict & visualize', 'FontSize',11,'FontWeight','bold');
    oy = 8; oh = 9;
    b_pred = bx( 3.0, oy, 22, oh, {'predict(pm, Xnew)','-> yhat, scores'}, C.viz);
    b_wi   = bx(26.5, oy, 27, oh, {'weight\_map\_object(pm, source)','-> statistic\_image + cache'}, C.viz);
    b_mont = bx(55.5, oy, 18, oh, {'montage(pm)','surface(pm)'}, C.viz);
    b_plot = bx(75.5, oy, 21.5, oh, {'plot(pm) / rocplot(pm)','confusionchart(pm)'}, C.viz);

    ar(14, iy, 14, oy+oh+0.2);                 % bootstrap -> predict/viz band
    ar(53.5, oy+oh*0.5, 55.4, oy+oh*0.5);      % weight_map_object -> montage
    ar(73.5, oy+oh*0.5, 75.4, oy+oh*0.5);      % montage -> plot group

    % ---------- legend (horizontal strip along the bottom) ----------
    ly = 2.2; sw = 2.2; sh = 1.6;
    items = {'input / data',C.data; 'construct',C.constr; 'fit & cross-validate',C.fit; ...
             'inference & tuning',C.infer; 'predict & visualize',C.viz};
    colx = [3 19 33 53 71];
    for i = 1:size(items,1)
        rectangle(ax,'Position',[colx(i) ly sw sh],'FaceColor',items{i,2},'EdgeColor',edge,'Curvature',0.4);
        text(colx(i)+sw+0.5, ly+sh/2, items{i,1}, 'FontSize',9,'VerticalAlignment','middle');
    end

    drawnow;
    exportgraphics(f, outpng, 'Resolution', 150);
    fprintf('wrote %s\n', outpng);
end


function h = draw_box(ax, x, y, w, h, str, facecolor, edge)
    rectangle(ax, 'Position',[x y w h], 'Curvature',0.18, ...
        'FaceColor',facecolor, 'EdgeColor',edge, 'LineWidth',1.0);
    text(ax, x+w/2, y+h/2, str, 'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', 'FontSize',9.5, 'Interpreter','tex');
end


function draw_arrow(ax, x1, y1, x2, y2, col)
    plot(ax, [x1 x2], [y1 y2], '-', 'Color',col, 'LineWidth',1.2);
    L = hypot(x2-x1, y2-y1); if L==0, return; end
    ux = (x2-x1)/L; uy = (y2-y1)/L; px = -uy; py = ux;
    a = 1.3; b = 0.7;
    xa = [x2, x2-a*ux+b*px, x2-a*ux-b*px];
    ya = [y2, y2-a*uy+b*py, y2-a*uy-b*py];
    patch(ax, xa, ya, col, 'EdgeColor',col);
end
