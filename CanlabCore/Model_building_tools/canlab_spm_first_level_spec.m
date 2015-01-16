bf      = SPM.xBF.bf;

V   = SPM.xBF.Volterra;
U   = SPM.Sess(s).U;
v   = length(U);

Uname     = U(i).name(1);
catch
    str       = sprintf('name for condition/trial %d ?',i);
    
    % Covariates: spm_fMRI_design 285
    %----------------------------------------------------------
    C     = SPM.Sess(s).C.C;
    Cname = SPM.Sess(s).C.name;
    
    
    % Onsets
    % spm_get_onsets.m
    %----------------------------------------------------------
    
    
    ons = U(i).ons;
    ons = ons(:);
    catch
        ons = [];
        end
        if isempty(ons)
            str      = ['vector of onsets - ' Uname{1}];
            ons      = spm_input(str,4,'r',' ',[Inf 1]);
            U(i).ons = ons(:);
        end
        
        dur = U(i).dur;
        
        xP    = U(i).P;
        Pname = xP(1).name;
        
        switch Pname
            
            case 'none'
                %----------------------------------------------------------
                xP.name  = 'none';
                xP.h     = 0;
                
        end
        
        % Parametric modulators
        %----------------------------------------------------------
        
        Pname = {'none','time','other'};
        Pname = spm_input('parametric modulation',6,'b',Pname);
        
        case 'other'
            %----------------------------------------------------------
            str   = ['# parameters (' Uname{1} ')'];
            for q = 1:spm_input(str,7,'n1',1);
                
                % get names and parametric variates
                %------------------------------------------------------
                str   = sprintf('parameter %d name',q);
                Pname = spm_input(str,7,'s');
                P     = spm_input(Pname,7,'r',[],[length(ons),1]);
                
                % order of polynomial expansion h
                %------------------------------------------------------
                h     = spm_input('polynomial order',8,'n1',1);
                
                % sub-indices and inputs
                %------------------------------------------------------
                xP(q).name = Pname;
                xP(q).P    = P(:);
                xP(q).h    = h;
                
            end
            