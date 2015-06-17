function cl = scn_check_vmpfc_signal(imgs)


    spm_image('init', imgs(1,:))
    input('Click on genu of CC and press return');

    orig_pos = spm_orthviews('Pos');
    pos = spm_orthviews('Pos');
    pos(2) = pos(2) + 20;
    pos(3) = pos(3) - 10;

    % FIRST REGION: VMPFC
    create_figure('Timeseries', 1, 1);

    cl = get_sphere(pos, imgs);
    cl(1).title = 'VMPFC';
    create_figure('Timeseries', 1, 1, 1);
    plot(cl(1).timeseries);


    % SECOND REGION: DACC

    pos = orig_pos;
    pos(3) = pos(3) + 20;
    cl(2) = get_sphere(pos, imgs);
        create_figure('Timeseries', 1, 1, 1);
    plot(cl(2).timeseries, 'r');
    
    legend({'VMPFC' 'DACC'})

end



function cl = get_sphere(pos, imgs)

    spm_orthviews('Reposition', pos);

    disp('Extracting data');
    [data, XYZvoxSphere, XYZmmSphere] = iimg_sphere_timeseries(imgs, pos, 10);

    datamean = mean(data, 2);
    glev = mean(datamean);
    snr = glev ./ std(datamean);

    % Viz sphere
    V = spm_vol(deblank(imgs(1,:)));

    cl(1).XYZ = XYZvoxSphere;
    cl(1).XYZmm = XYZmmSphere;
    cl(1).Z = ones(1, size(XYZvoxSphere, 2));
    cl(1).title = 'Sphere';
    cl(1).voxSize = diag(V.mat(1:3, 1:3))';
    cl(1).M = V.mat;
    cl(1).timeseries = datamean;
    cl(1).snr = snr;
    cl(1).global_lev = glev;

    cluster_orthviews(cl,  {[1 0 0]}, 'add', 'solid');
    drawnow

end

