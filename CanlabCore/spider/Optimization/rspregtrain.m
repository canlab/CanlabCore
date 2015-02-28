function [alpha_fin,b_fin] = rsptrain(K,Y,C,epsilon,T,S_min),
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
gamma = 1/T*ones(T,1);
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
        tab_tmp(train_ind_tmp)=1;
        current_position(1) = cp1+1;
        Path_from_root(current_position(1)) = current_position(2);
        current_position(2) = active(1);
        Current_Path(current_position(2) + (current_position(1)-1)*T,:) = tab_tmp;        
        
    else
        
        if (cp1==0),
            tmp = [1:m]; 
        else
            tmp = find(Current_Path(current_position(2)+(current_position(1)-1)*T,:)==1);            
        end;
            Ylearn = Y(tmp,:);
            Clearn = C/length(tmp);
            Klearn = K(tmp,tmp);
        
        if current_position(1)==k,
            % Learn a SVR with Clearn
            mtmp = length(tmp);
            K2 = [Klearn , -Klearn ; -Klearn , Klearn];
            c = zeros(2*mtmp,1);
            c(1:mtmp) = epsilon*ones(mtmp,1) - Ylearn;
            c(mtmp+1:2*mtmp) = epsilon*ones(mtmp,1) + Ylearn;
            cst = ones(2*mtmp,1); cst(mtmp+1:2*mtmp) = -1;
            [alphatmp,y] = quadsolve(K2,c,cst',0,Clearn);
            alpha=alphatmp(1:mtmp)-alphatmp(mtmp+1:2*mtmp);
        else
            % learning with stability            
            c = Klearn*(Current_Path(current_position(1)*T+1:current_position(1)*T+T,tmp))';
	    c = c*gamma;
            Ylearn = Ylearn-c;              
            mtmp = length(tmp);
            K2 = [Klearn , -Klearn ; -Klearn , Klearn];
            c = zeros(2*mtmp,1);
            c(1:mtmp) = epsilon*ones(mtmp,1) - Ylearn;
            c(mtmp+1:2*mtmp) = epsilon*ones(mtmp,1) + Ylearn;
            cst = ones(2*mtmp,1); cst(mtmp+1:2*mtmp) = -1;
            [alphatmp,y] = quadsolve(K2,c,cst',0,Clearn);
            alpha=alphatmp(1:mtmp)-alphatmp(mtmp+1:2*mtmp);                        
            %alpha = alpha + 1/T*sum(Current_Path(current_position(1)*T+1:current_position(1)*T+T,tmp)',2);
            alpha = alpha + Current_Path(current_position(1)*T+1:current_position(1)*T+T,tmp)'*gamma;
        end; % if current_pos...==k
        
        if current_position ~=0,
            Current_Path(current_position(2)+(current_position(1)-1)*T,tmp) = alpha';
            Node_done(current_position(1),current_position(2)) = 1;
            count=count+1;
            current_position(2) = Path_from_root(current_position(1));
            current_position(1) = current_position(1)-1;
        else,
            not_finished=0;
            alpha_fin = alpha;        
            b_fin = -y;                   
        end; % if current_pos...~=0
        
    end;% if ~isempty(active)
    
%    disp(sprintf('Node done: %d\n',count));
    
end; %while not_finished
% End function
