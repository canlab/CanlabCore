function [alpha_fin,b_fin] = rsptrain(K,Y,C,T,S_min),
% Recursive Stabilization Procedure
%% Variable Description
%%
%% 
m = size(K,1);
k = ceil(log(m/S_min)/log(T/(T-1)));
current_position = [0,0];
Node_done = zeros(k,T);
Current_Path = zeros(k*T,m);
not_finished = 1;
count = 0;
Path_from_root = zeros(k,1);
while not_finished,
    
    cp1 = current_position(1);
    cp2 = current_position(2);
    
    if (cp1~=k),
        active = find(Node_done(cp1+1,1:T)==0);
    else
        active = [];
    end;%if (cp1~=k)
    
    if ~isempty(active),
        
        % The part of the execution tree is not completely explored
        if cp1==0,
            tmp=[1:m];
        else
            tmp = find(Current_Path(cp2 + (cp1-1)*T,:)==1);
        end; % if cp1==0,
        
        % Build the training set for the node
        tab_tmp = zeros(1,m);
        train_ind_tmp = tmp([1:floor((active(1)-1)*length(tmp)/T),floor(active(1)*length(tmp)/T):length(tmp)]);
        %% If all the points are from the same class then do not go further
        is_plus = find(Y(train_ind_tmp)==1);
        
        if (length(is_plus)==0|(length(is_plus)==length(train_ind_tmp))),
            Node_done(cp1+1,1:T)=1;
            Current_Path([1:T] + (cp1)*T,:) = zeros(T,m);
        else
            tab_tmp(train_ind_tmp)=1;
            current_position(1) = cp1+1;
            Path_from_root(current_position(1)) = current_position(2);
            current_position(2) = active(1);
            Current_Path(current_position(2) + (current_position(1)-1)*T,:) = tab_tmp;        
        end;                           
    else    
        if (cp1==0),
            tmp = [1:m]; 
        else
            tmp = find(Current_Path(current_position(2)+(current_position(1)-1)*T,:)==1);            
        end;
            Ylearn = Y(tmp,:);
            Clearn = (C*m)/(length(tmp));
            Klearn = K(tmp,tmp);
        
        if current_position(1)==k,
            [alpha] = quadsolve(Klearn,-ones(length(tmp),1),Ylearn',0,Clearn);
            alpha_bias=alpha;
            scale_f = sqrt(alpha'*Klearn*alpha);
            alpha = alpha/scale_f;
        else
            c = Klearn*(Current_Path(current_position(1)*T+1:current_position(1)*T+T,tmp))';
            c = mean(c,2);
            c = -1+c;
            H = 1/(T)*Klearn;
            [alpha,y] = quadsolve(H,c,Ylearn',0,Clearn);           
            alpha_bias = alpha;
            alpha = 1/T*(alpha + sum(Current_Path(current_position(1)*T+1:current_position(1)*T+T,tmp)',2));
            scale_f = sqrt(alpha'*Klearn*alpha);
            alpha = alpha/scale_f;
        end; % if current_pos...==k-1
        
        if current_position ~=0,
            Current_Path(current_position(2)+(current_position(1)-1)*T,tmp) = alpha';
            Node_done(current_position(1)+1,1:T)=0;
           
            Node_done(current_position(1),current_position(2)) = 1;
            count=count+1;
            current_position(2) = Path_from_root(current_position(1));
            current_position(1) = current_position(1)-1;
        else,
            not_finished=0;
            alpha_fin = alpha;        
            b_fin = -y/scale_f;                   
        end; % if current_pos...~=0
        
    end;% if ~isempty(active)
    
%    disp(sprintf('Node done: %d\n',count));
    
end; %while not_finished
% End function