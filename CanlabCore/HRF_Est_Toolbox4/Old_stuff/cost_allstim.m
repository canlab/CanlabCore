function Q = cost_allstim(V,t,tc,Run)
%
% Least-squares cost function for the IL model
%
% INPUT:
% Run = stick function
% tc = time course
% t = vector of time points
% V = parameters
%
% OUTPUT:
% Q = cost
%
% By Martin Lindquist and Tor Wager
% Edited 12/12/06
% Further edited by Christian Waugh 2/15/08 to include multiple trialtypes

numstim = length(Run);
len = length(Run{1});
h = zeros(length(t),numstim);
yhatt =zeros(len,numstim);

for k = 1:numstim
    h(:,k) = Get_Logit(V(k*7-6:k*7),t);            % Get IL model corresponding to parameters V
    yhat(:,k) = conv(Run{k}, h(:,k));              % Convolve IL model with stick function
    yhatt(:,k) = yhat(1:len,k);
end

yhat2 = sum(yhatt,2); %Sum models together to get overall estimate

Q = sum((yhat2-tc).^2);              % Calculate cost function

return