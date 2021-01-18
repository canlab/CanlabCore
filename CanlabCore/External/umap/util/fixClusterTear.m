%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [P, tears]=fixClusterTear(M, P, neighborHood, ids)
if nargin<4
    ids=unique( P(P<-1) );
end
tears=[];
fix;

    function fix        
        MM=M^2;
        N=length(ids);
        for i=1:N
          %  if i==71
          %      disp('uh')
          %  end
            id=ids(i);
            idCnt=1;
            p=find(P==id, 1, 'first');
            while ~isempty(p)
                newId=id-(idCnt*MM);
                done=[];
                toDo=p;
                while ~isempty(toDo)
                    P(toDo)=newId;
                    done=[done toDo];
                    neighbors=unique(cell2mat( neighborHood(toDo) ));
                    neighbors=neighbors(P(neighbors)==id);
                    toDo=setdiff(neighbors,done);
                end
                p=find(P==id, 1, 'first');
                idCnt=idCnt+1;
                %if sum(P==newId)<4
                %    P(P==newId)=-1;
                %end
            end
            if idCnt>2
                tears(end+1,:)=[id, idCnt-1];
            end
        end
        disp(['check sum=' sum(P)]);
    end
end