function Q = cost_allstim_2(V,tr,tc,Run,down)
%
% Least-squares cost function for the IL model
% Multi-condition case
%
% INPUT:
% down = downsampling factor
% Run = stick function
% tc = time course
% tr = repetition time
% V = parameters
%
% OUTPUT:
% Q = cost
%
% By Martin Lindquist and Tor Wager
% Edited 12/12/06
% Further edited by Christian Waugh 2/15/08 to include multiple trialtypes
% Edited by ML on 02/12/13

t = 0:(1/down):30;

numstim = length(Run);
len = length(Run{1});
h = zeros(length(t),numstim);
yhatt =zeros(len,numstim);


for k = 1:numstim

    h(:,k) = Get_Logit(V(k*7-6:k*7),t);             % Get IL model corresponding to parameters V
    yhat(:,k) = conv(Run{k}, h(:,k));               % Convolve IL model with stick function
    yhatt(:,k) = yhat(1:len,k);

end

tt = 1:(tr*down):len;

yhat2 = sum(yhatt,2); %Sum models together to get overall estimate

Q = sum((yhat2(tt)-tc).^2);              % Calculate cost function

return