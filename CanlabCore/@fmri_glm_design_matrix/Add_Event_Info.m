function obj = Add_Event_Info(obj,Doc_path)
% Add event information from excel document
%    Example file
%    column1   column2   column3   column4 
%    onset     duration    Name  Customized_Name
%     ...       ...         ...     ...
% 
% 
%  %    Programmers' notes:
%    Ke Bo, 12/23/2021 Initial function
%
% :Usage:
% ::
%
%     obj = Add_Event_Info(obj,Doc_path)
%

%%%%% Add event info from example sheet %%% 
% Doc_path='C:\Users\KeBo\Documents\GitHub\CanlabCore\CanlabCore\@fmri_glm_design_matrix\Example_BehaviorSheet.xlsx'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[numericData, textData, rawData]=xlsread(Doc_path);  %%%% read excel sheet, maybe need a backup loading function

Regressor_num=length(textData(:,1)); 
for i=2:Regressor_num
    Singletrial_Regressor_name{i-1}=[textData{i,3},textData{i,4}]; %Get each individual regressor name
end
ConditionName=unique(Singletrial_Regressor_name); %Find all types of regressor

for cond=1:size(ConditionName,2)
    obj.Sess.U(cond).name=ConditionName{cond};  %%%Load name for regressor
    Condition_Index=find(strcmp(Singletrial_Regressor_name,ConditionName{cond})==1);
    
    
    obj.Sess.U(cond).ons=numericData(Condition_Index,1); %%%Load onset timing for regressor
    
    obj.Sess.U(cond).dur=numericData(Condition_Index,2); %%%Load duration for regressor
    
end



