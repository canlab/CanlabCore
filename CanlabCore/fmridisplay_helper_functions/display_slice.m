function Z = display_slice(dat, wh_slice, SPACE, varargin)
% Resample slice data in dat to SPACE and display
%
% :Usage:
% ::
%
%     Z = display_slice(dat, wh_slice, SPACE, varargin)
%

dolighten = 0;
myview = 'axial';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case 'lighten', dolighten = 1;
            case 'nolighten', dolighten = 0;
                
            case 'sagittal', myview = 'sagittal'; %disp('Warning! NOT implemented correctly yet!!!'), %pause(5)
            case 'coronal', myview = 'coronal'; %disp('Warning! NOT implemented correctly yet!!!'), pause(5)
            case 'axial', myview = 'axial';
                
            otherwise, warning('displayslice:Badinput', ['Unknown input string option:' varargin{i}]);
        end
    end
end


switch myview
    
    case 'axial'
        
        
        Zo = dat(:, :, wh_slice);
        
        Z = interp2(SPACE.Xo(:, :, wh_slice), SPACE.Yo(:, :, wh_slice), Zo, SPACE.X(:, :, wh_slice), SPACE.Y(:, :, wh_slice));
        
        % set transparent value for clear axes
        myAlphaData = double(abs(Z') > 0);
        
        % done at slice level now
        Z = lighten_underlay_edges(Z, 32);
        
        h = imagesc(SPACE.xcoords, SPACE.ycoords, Z');
        
        % If we set alphadata to clear for BG and axis color to none, we get clear
        % axes
        set(h, 'AlphaDataMapping', 'scaled', 'AlphaData', myAlphaData)
        set(gca, 'Color', 'none')
        
    case 'coronal'
        Zo = squeeze(dat(:, wh_slice, :));
        
        yo = squeeze(SPACE.Yo(:, wh_slice, :));
        zo = squeeze(SPACE.Zo(:, wh_slice, :));
        yi = squeeze(SPACE.Y(:, wh_slice, :));
        zi = squeeze(SPACE.Z(:, wh_slice, :));
        
        Z = interp2(yo', zo', Zo', yi', zi');
        
        % set transparent value for clear axes
        myAlphaData = double(abs(Z) > 0);
        
        % done at slice level now
        Z = lighten_underlay_edges(Z, 32);
        
        h = imagesc(SPACE.xcoords', SPACE.zcoords', Z); % Wani modified ycoords -> xcoords
        
        % If we set alphadata to clear for BG and axis color to none, we get clear
        % axes
        set(h, 'AlphaDataMapping', 'scaled', 'AlphaData', myAlphaData)
        set(gca, 'Color', 'none')
        
    case 'sagittal'
        Zo = squeeze(dat(wh_slice, :, :));
        
        xo = squeeze(SPACE.Xo(wh_slice, :, :));
        zo = squeeze(SPACE.Zo(wh_slice, :, :));
        xi = squeeze(SPACE.X(wh_slice, :, :));
        zi = squeeze(SPACE.Z(wh_slice, :, :));
        
        Z = interp2(xo', zo', Zo', xi', zi');
        
        % set transparent value for clear axes
        myAlphaData = double(abs(Z) > 0);
        
        % done at slice level now (slow!)
        Z = lighten_underlay_edges(Z, 32);
        
        h = imagesc(SPACE.ycoords', SPACE.zcoords', Z);
        
        % If we set alphadata to clear for BG and axis color to none, we get clear
        % axes
        set(h, 'AlphaDataMapping', 'scaled', 'AlphaData', myAlphaData)
        set(gca, 'XDir', 'reverse', 'Color', 'none')
        
    otherwise
        
        error ('Unknown view.  Check inputs...view must be axial, sagittal.');
        
end

colormap(gray);
set(gca,'YDir', 'Normal')
axis image

end


