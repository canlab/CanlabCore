function outList = transform2switches(stimList)
% assumes 2 stimtypes, 1 and 2
% classifies them into no switch (1 and 2) or switch (3 and 4)
outList(1,1) = stimList(1,1);
for i = 2:size(stimList,1)
   switch stimList(i-1,1)		% based on previous trial
   case 1
      if stimList(i,1) == 2, outList(i,1) = 4;,else outList(i,1) = stimList(i,1);,end
   case 2
      if stimList(i,1) == 1, outList(i,1) = 3;,else outList(i,1) = stimList(i,1);,end
   case 5
      if i > 2 & stimList(i,1) == 1 & stimList(i-2,1) == 2, outList(i,1) = 3;
      elseif i > 2 & stimList(i,1) == 2 & stimList(i-2,1) == 1, outList(i,1) = 4;
      else outList(i,1) = stimList(i,1);
      end
   otherwise outList(i,1) = stimList(i,1);
   end
end
return

