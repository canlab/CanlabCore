function setFigWidth(hFig, width)
    if(~ishandle(hFig))
        error('SCANLab:invalidHandleError', 'Figure handle %d is not a valid handle', hFig);
    else
        figPos = get(hFig, 'Position');
        set(hFig, 'Position', [figPos(1:2) width figPos(4)]);
    end
end