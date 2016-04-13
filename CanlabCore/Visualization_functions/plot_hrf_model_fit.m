function plot_hrf_model_fit(m,TR,bf,varargin)
% :Usage:
% ::
%
%    plot_hrf_model_fit(m,TR,bf,[stim input function for epoch, in 1 s resolution])
%
% if m is an hrf curve sampled at 1 s,
%
% plot_model_fit(m,1,'fir');
%
% plot_model_fit(m,1,'hrf');


len = length(m); % interpolate hrf to 1 s

m = interp1((0:TR:TR*len-TR)',m',(0:1:ceil(len.*TR))');
m = m ./ max(m);

m(isnan(m)) = [];
len = length(m);




switch bf
    case 'hrf'
        hrf = spm_hrf(1)./max(spm_hrf(1));

    case 'derivatives'

        xbf = spm_get_bf(struct('dt',1,'name','hrf (with time and dispersion derivatives)'));
        %tor_fig
        %x = 0:.1:32; plot(x(1:length(hrf)),xbf.bf(:,1),'k','LineWidth',2);
        %hold on;
        %plot(x(1:length(hrf)),xbf.bf(:,2),'k--','LineWidth',2);
        %plot(x(1:length(hrf)),xbf.bf(:,3),'k.-','LineWidth',2);
        %legend({'HRF' 'Derivative (Time)' 'Derivative (Dispersion)'})
        %set(gca,'YColor',[1 1 1],'FontSize',16);
        
        hrf = xbf.bf;
        
    case 'fir'
        xbf = spm_get_bf(struct('dt',1,'name','Finite Impulse Response'));
        hrf = xbf.bf;
        
    otherwise 
        error('Enter name of basis set: ''hrf'' ''derivatives'' or ''fir''')
        
        
end

if size(hrf,1) < len
    hrf = [hrf; zeros(len - size(hrf,1),size(hrf,2))];
end

% convolve, if stim function entered
if length(varargin)>0
    stim = varargin{1};
    for i = 1:size(hrf,2)
        h = conv(stim,hrf(:,i));
        hrf(:,i) = h(1:size(hrf,1));
    end
end


X = hrf; X(:,end+1) = 1; X = X(1:len,:);

b = pinv(X) * m;
f = X*b;


tor_fig;plot(m,'k','LineWidth',2)

hold on; plot(f,'.--','Color',[.3 .3 .3],'LineWidth',2);
set(gca,'YColor',[1 1 1],'FontSize',24);
legend({'Actual' 'Fitted'})

if length(varargin) > 0
    delta = stim;    
    d = find(delta);
    z = zeros(size(d));

    plot([d d]',[z-.2 z-.1]','k','LineWidth',1);
end


