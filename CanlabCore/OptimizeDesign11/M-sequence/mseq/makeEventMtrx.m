function eventMatrix=makeEventMtrx(events,HDRdur)

% eventMatrix=makeEventMtrx(events,HDRdur)
% creates a length(events) X (HDRdur) convolution matrix 
% for each row of "events" and concatenates them

eventMatrix=[];
nVals=size(events,1);
for i=1:nVals,
   eventMatrix=[eventMatrix convmtx(events(i,:)',HDRdur)];
end;

