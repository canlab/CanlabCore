function setFigHeight(hFig, height)
    if(~ishandle(hFig))
        error('SCANLab:invalidHandleError', 'Figure handle %d is not a valid handle', hFig);
    else
        figPos = get(hFig, 'Position');
        set(hFig, 'Position', [figPos(1:3) height]);
    end
end