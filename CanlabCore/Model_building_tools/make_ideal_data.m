function ideal_data = make_ideal_data(hrf,sf,eres,finalLength)
% function ideal_data = make_ideal_data(hrf,sf,eres,finalLength)
% 

mylen = size(sf{1},1);

% -------------------------------------------------------------------
% * make hi-res ideal data with eres timebins per TR
% -------------------------------------------------------------------

ideal_data = [];
for i = 1:length(sf)
    mydata = conv(hrf,full(sf{i}));
    mydata = mydata(1:mylen);
    ideal_data(:,i) = mydata;
end
ideal_data = sum(ideal_data,2);

% -------------------------------------------------------------------
% * downsample ideal data to TR
% -------------------------------------------------------------------
%i2 = ideal_data;

ideal_data = resample(ideal_data,1,eres);
ideal_data = ideal_data(1:finalLength);

%figure;hold on
%plot(1:1/length(i2):2-1/length(i2),i2)
%plot(1:1/length(ideal_data):2-1/length(ideal_data),ideal_data,'r')

return