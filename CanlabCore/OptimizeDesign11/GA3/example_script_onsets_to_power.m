load leftRhtSac % contains c{1:2} = onset times in s for 2 conditions
TR = 2; 
HPlength = 100;
LPsmooth = 0;
con = [1 -1 0];
n_in_group = 12;

[X,d,dh]=onsets2delta(sac,TR);

[S,Vi,svi] = getSmoothing(HPlength,LPsmooth,TR,size(X,1),'auto');

% power for contrast

[ipow,gpow,OUT] = xpower(X,con,[],n_in_group,[],Vi,S);

% test power as a function of HP filter
% THERE's a problem with this, as it doesn't take smoothing into account
% right...maybe have to apply smoothing to sigma,
ind = 1; clear ip, clear gp

hprange = 20:10:500;
for HPlength = hprange

    [S] = getSmoothing(HPlength,LPsmooth,TR,size(X,1),'auto');
    [ip(ind),gp(ind)] = xpower(X,con,[],n_in_group,[],Vi,S);
    
    ind = ind + 1;
    
end

figure;plot(hprange,ip); hold on; plot(hprange,gp,'r'); 
xlabel('High-pass filter length'), ylabel('Power (Z)'),;legend({'Individual' 'Group'}) 


% generally useful
E = inline('inv(x'' * x) * x'' * Vi','x','Vi');
sigma = 1; t = .5; e = E(S*X,S* sigma * Vi);ee = e * e', t ./ trace(ee).^.5
sigma = 1; t = .5; e = E(X,sigma * Vi);ee = e * e', t ./ trace(ee).^.5

