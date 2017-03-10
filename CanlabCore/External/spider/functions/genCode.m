function genCode(a,type,filename,funcname)
% genCode generates plain C/JAVA compilable modules of 
% a trained algorithm like SVM/SVR.  
% Example:
% 
% [r s]=train(svm(kernel('rbf',2)),gen(toy))
% genCode(s,'C','test.c','myfunc')    
% 
% will generate a  c file containing a function double myfunc(double *x)
% which implements \sum_{i=1}^{nrofsv} \alpha_i*k(x_i,x) + b0
% If you dont specify a function name 'predict' is used



if(nargin==3)
    funcname='predict';
end

fid=-1;
switch( lower(type))
    case 'java' 
        disp('generating Java MODULE');
        fid=fopen(filename,'wt');
        if(fid>=0)
	    classname = strrep(filename,'.java','');
            genJava(a,fid,funcname,classname);
            fclose(fid);
        end
        return;
    case {'perl','tcl','phyton'}
        disp('only programming languages are possible');
    case {'fortran','f77'}
        disp('..back to the future man!!');
    case 'c'
        disp('generating C MODULE');
        fid=fopen(filename,'wt');
        if(fid>=0)
            genC(a,fid,funcname);
            fclose(fid);
        end
        return;
    otherwise
        disp('unknown target language');
end






%% ================================================================
%% ================================================================
%%                                  C
%% ================================================================
%% ================================================================
function genC(a,f,funcname)
switch( class(a))
    case 'data' 
        disp('IMPLEMENT ME!');
    case {'svm','svr'} 
        disp('gen SVM/SVR MODULE');
        
        if(a.algorithm.trained==1)
            writeCHeader(a,f);
            d=a.Xsv;
            [l,n,k]=get_dim(d);
            fprintf(f,'#define _NROFSV %d\n',l);
            fprintf(f,'#define _INPUTSPACE %d\n',n);
            fprintf(f,'#define _OUTPUTSPACE %d\n',k);
            % ===========================================
            % output Support vectors
            % ===========================================
            fprintf(f,'static long double SV[_NROFSV][_INPUTSPACE] =\n');
        	fprintf(f,'{\n\t');

            for i=1:l-1
                
                if (mod(i,10)==0) 
                    fprintf('.');
                end
                if (mod(i,100)==0) 
                    fprintf('o\n');
                end
                fprintf(f,'\t{');
                for j=1:n-1
                   if( mod(j+1,10) ==0)
                        fprintf(f,'\n\t');
                    end

                      fprintf(f,'%10.8f, ',d.X(i,j));
                end
                fprintf(f,'%10.8f},\n',d.X(i,n));
            end
            fprintf(f,'\t{');
            i=l;
            for j=1:n-1
                   if( mod(j+1,10) ==0)
                        fprintf(f,'\n\t');
                    end
                  fprintf(f,'%10.8f, ',d.X(i,j));
            end
            fprintf(f,'%10.8f}\n};',d.X(l,n));
        
            % ===========================================
            % output Lagrange factors
            % ===========================================
        	fprintf(f,'\n\n\nstatic long double LAGRANGE[_NROFSV]=\n{\n');			
            fprintf(f,'\t%10.8f,\n',a.alpha(1:l-1));
            fprintf(f,'\t%10.8f};\n\n',a.alpha(l));
            
            % ===========================================
            % output dot function
            % ===========================================
            
            fprintf(f,'static double dot(double *x,double*y)\n');
            fprintf(f,'{\n');
            fprintf(f,'\t double _LResult=0;int index=0;\n');
            fprintf(f,'\t for(;index<_INPUTSPACE;index++) _LResult+=x[index]*y[index];\n');
            fprintf(f,'\t return _LResult;\n}\n\n');
            
            % ===========================================
            % output sub function
            % ===========================================
            
            
              
            % ===========================================
            % output kernel function
            % ===========================================

            if(~isempty(a.child.kerparam))
                fprintf(f,'#define _KERPARAM %d\n',a.child.kerparam);
            end


            fprintf(f,'static double kernel(double *x,double*y)\n');
            fprintf(f,'{\n');
            switch( lower( a.child.ker))
                case 'rbf'
                    fprintf(f,'\t double _LResult=0;int index=0;\n');
                    fprintf(f,'\t for(;index<_INPUTSPACE;index++) _LResult+=(x[index]-y[index])*(x[index]-y[index]);\n');
                    fprintf(f,'\t return exp( -1./(2*_KERPARAM*_KERPARAM)* _LResult);\n}');
                case 'poly'
                    fprintf(f,'\t return pow(dot(x,y)+1,_KERPARAM);\n}');
                case 'linear'
                    fprintf(f,'\t return dot(x,y);\n}');
            end
            fprintf(f,'\n\n');
            
            % ===========================================
            % output decision function
            % ===========================================
            fprintf(f,'#define _BIAS %f\n',a.b0);
          
            fprintf(f,'double %s(double *x)\n',funcname);
            fprintf(f,'{\n');
            fprintf(f,'\t double _LResult=0;int index=0;\n');
            fprintf(f,'\t for(;index<_NROFSV;index++) _LResult+=LAGRANGE[index]*kernel(SV[index],x) ;\n');
            fprintf(f,'\t return _LResult+_BIAS;\n');
            fprintf(f,'}\n\n');
                    

              
         else
            error('needed a trained svm');
        end
        
    otherwise
        disp('Unknown obj, code gen not supported');
end


function writeCHeader(a,f)

    fprintf(f,'/* ---------   +SPIDER + Generated ---------------\n');
    writeLogo(f);
	fprintf(f,'\n\n             Used File : %s          */\n',get_name(a));
	fprintf(f,'/*             Creation Date dmy: %s     */\n',date);
	fprintf(f,'#include  <math.h>\n#include  <mex.h>\n\n');

    

%% ================================================================
%% ================================================================
%%                                  Java
%% ================================================================
%% ================================================================
function genJava(a,f,funcname, classname)
switch( class(a))
    case 'data' 
        disp('IMPLEMENT ME!');
    case {'svm','svr'} 
        disp('gen SVM/SVR MODULE');
        
        if(a.algorithm.trained==1)
            writeJHeader(a,f,classname);
            d=a.Xsv;							%x support vectoren
            [l,n,k]=get_dim(d);
            fprintf(f,'protected static int NROFSV = %d;\n',l);
            fprintf(f,'protected static int INPUTSPACE = %d;\n',n);
            fprintf(f,'protected static int OUTPUTSPACE = %d;\n',k);
            % ===========================================
            % output Support vectors
            % ===========================================
            fprintf(f,'public static double SV[][] =\n');
        	fprintf(f,'{\n\t');

            for i=1:l-1
                fprintf(f,'\t{');
                for j=1:n-1
                   if( mod(j+1,10) ==0)
                        fprintf(f,'\n\t');
                    end
                      fprintf(f,'%f, ',d.X(i,j));
                end
                fprintf(f,'%f},\n',d.X(i,n));
            end
            fprintf(f,'\t{');
            for j=1:n-1
                   if( mod(j+1,10) ==0)
                        fprintf(f,'\n\t');
                    end
                  fprintf(f,'%f, ',d.X(i,j));
            end
            fprintf(f,'%f}\n};',d.X(l,n));
        
            % ===========================================
            % output Lagrange factors
            % ===========================================
            fprintf(f,'\n\n\nprotected static double LAGRANGE[]=\n{\n');			
            fprintf(f,'\t%f,\n',a.alpha(1:l-1));
            fprintf(f,'\t%f};\n\n',a.alpha(l));
            
            % ===========================================
            % output dot function
            % ===========================================
            
            fprintf(f,'protected static double dot(double[] x,double[] y){\n');
            fprintf(f,'\t double lResult=0;\n\tint index;\n');
            fprintf(f,'\t for(index = 0;index<INPUTSPACE;++index)\n\t\t lResult+=x[index]*y[index];\n');
            fprintf(f,'\t return lResult;\n}\n\n');
            
            % ===========================================
            % output sub function
            % ===========================================
            
            
              
            % ===========================================
            % output kernel function
            % ===========================================

            if(~isempty(a.child.kerparam))
                fprintf(f,'protected static int KERPARAM = %d;\n',a.child.kerparam);
            end


            fprintf(f,'protected static double kernel(double[] x,double[] y){');
            switch( lower( a.child.ker))
                case 'rbf'
                    fprintf(f,'\t double lResult=0;\n\t int index;\n');
                    fprintf(f,'\t for(index = 0;index<INPUTSPACE;++index) lResult+=(x[index]-y[index])*(x[index]-y[index]);\n');
                    fprintf(f,'\t return Math.exp( -1./(2*KERPARAM*KERPARAM)* lResult);\n}');
                case 'poly'
                    fprintf(f,'\t return Math.pow(dot(x,y)+1,KERPARAM);\n}');
                case 'linear'
                    fprintf(f,'\t return dot(x,y);\n}');
            end
            fprintf(f,'\n\n');
            
            % ===========================================
            % output decision function
            % ===========================================
            fprintf(f,'protected static double BIAS = %f;\n',a.b0);
          
            fprintf(f,'public double %s(double[] x){\n',funcname);
            fprintf(f,'\t double lResult=0;\n\tint index;\n');
            fprintf(f,'\t for(index = 0;index<NROFSV;++index) lResult+=LAGRANGE[index]*kernel(SV[index],x);\n');
            fprintf(f,'\t return lResult+BIAS;\n');
            fprintf(f,'}\n\n}');
                    

              
         else
            error('needed a trained svm');
        end
        
    otherwise
        disp('Unknown obj, code gen not supported');
end

function writeJHeader(a,f,classname)

    fprintf(f,'/* ---------   +SPIDER + Generated ---------------\n');
    writeLogo(f);
	fprintf(f,'\n\n             Used File : %s          */\n',get_name(a));
	fprintf(f,'/*             Creation Date dmy: %s     */\n',date);
	fprintf(f,'public class %s{\n',classname);

    
%% ================================================================
%% ================================================================
%%                                  geek area
%% ================================================================
%% ================================================================
function writeLogo(f)
fprintf(f,'_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|___|___|___|___|___|D)_|___|___|___|___|___|___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|___|___|___|___|___|_JMD?__|___|___|___|___|___|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|___|___|___|___|__L[<*LT___|___|___|___|___|___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|___|___|___|___|___RF_?)_]Z!___|___|___|___|___|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|___|___|___|____A+_|{+_|:E&;___|___|___|___|___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|___|___|___|___|?S{|._QA_._|<VR)~__|___|__."(__|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|___|___|___||E1`_._LCC"|___|.+&&C+/)?/ITJQ4|___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|___|___|___-=NT;___|J1_EU-_|___|_`:})))<{~.JI{_|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___[&|)<>)=L&R+~__.|_TL|`1;Z-__|___|___|___|Z.}[___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|_TV||//|<;`|__.._.EC___&_"X|___|___|___|.H/|-&_|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___*ZI_|___|__.|__?S>__|_J-|_VL-___|___|:]GW_._V~__|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|_*==C__|___|__>AL,_|___I.__|(VT?;_-"=RE|K1_._)+|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___.1_YZ){;;{|JN1;_|___|_Y:|___|;+J1LI/._]F]___`V__|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|_"C|13I[CC[?;__|___|__/LH__|___|___|__|V_V_|__I(___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___/=_+(&__|___|___|_._?L|ON,__|___|_._.U__L*__|.A_|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|_E-|={{A_._|___|___|?U_)<{U>___|___|`U<|_{L|___?J__|___|___|___|__\n');
fprintf(f,'___|___|___|___|__[I__|(_[C|___|___|_LV._>=|_JE=*.,.>]M5___|A-_|__V|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___JV_.|I!_|R>_|___|_+&JW__}[__|-=JJLJ]}VT-|__*N___|`Z.|___|___|___|__\n');
fprintf(f,'___|___|___|__.+A]_|__J;__"HOJ]+FTR[`;#<_>F|.__|___|_L}/|__|_1||___~X:_|___|___|___|\n');
fprintf(f,'_|___|___FIJJTR=.|___|&`_|}JO"}|"____|GM_Z4~_|___|._<V_?=|___|O;_|__-X>__|___|___|__\n');
fprintf(f,'___|___|_W/-~._|___|_.&|__(>{T_|___|__KO=&AO|__|___|O,_-&__|__-X___|_,NJ`__|___|___|\n');
fprintf(f,'_|___|___&D__|___|___F<__|??_LC__|___|*#K,I,VN|":|/G?|__E-___|_"Z,___|_?V|___|___|__\n');
fprintf(f,'___|___|_?1N`__|___|*E_|__]{__B)___|`}LFQ_V|_.=CF+[S"__|<J_|___|}Z.|____]JM"___|___|\n');
fprintf(f,'_|___|___[;|U~___|_!X~___|&__|=UVTJTVJ;WO.J__|___|1=[|___N"__|___{U,./V&(`E__|___|__\n');
fprintf(f,'___|___|_+||*S||-}1V___|__E|__(CV,_,___O@*F|___._]]|E__|__S-___|__}JNI-|_{C|___|___|\n');
fprintf(f,'_|___|___I{__,X3T={__|___+=_.|=*|F___._A#H4;_|__,O__C|___|;X<|__~IXQ.|_._F*__|___|__\n');
fprintf(f,'___|___|_F:|._+F&!_|.__|_O.|__+*_N>|___T#WES[..-Q)_|.O.|_._.J1<1A<]}___|_V.|___|___|\n');
fprintf(f,'_|___|___L.__|R__RE"_|_;S*___|D(_.B/{"|B#T=.[&JBV<___<O__|__:ZY)_|V,_|__~T___|___|__\n');
fprintf(f,'___|___|-T_|_.R|__{NSJ]J;__|_-L8Y.=N=F[V#N]|___L|F_|__>&:=$88QI|_-T|___|{]_|___|___|\n');
fprintf(f,'_|___|__|+___;L__|__MKE._|__.S_{R0U[||_R#J]_.|=]_E.__|{3##6Z-]{__;I__|__<=___|___|__\n');
fprintf(f,'___|___|R-_|_)/|___|I_>AL?!|A>_|.N#FT{_+#O?|_.Z._?T||L##YJ_|_T_|_"I|___|)?_|___|___|\n');
fprintf(f,'_|___|_"&|___L;__|__R.__R3JL_|_._JY#WE__%AD__F)_.|CC%#OL,F___1___.F__|__>|___|___|__\n');
fprintf(f,'___|___&*__|`N_|___,V_.||>.J&|_|FC.7#0L>8WMC&3_|__L##B_|"]_|~1_|_`R|___||[_|___|___|\n');
fprintf(f,'_|___|*V_|__)]___|_?[|__>=__-LW$V|__U#O=4WO3ND"_FB#TB|__!]__.L___|T._|__;F___|___|__\n');
fprintf(f,'___|__S-___|E__|__.U,__|[<_|__R**C{.CZWSN#M%WAGP#MO;]_.|:1_|_&`|__=/___|`R_|___|___|\n');
fprintf(f,'_|___J?__|_)F|___|R?_|__&____~&|LJX8FRW####W#Z##H"_<+|___&___C>__|_A_|___T-__|___|__\n');
fprintf(f,'___|/L_|__~Z`__|_[T|___!T__._EW#%90%###W@#@#@#Q)?__:&__|_I<|_-E|___1)__|_=/|___|___|\n');
fprintf(f,'_|_|O|___.H*_|__CJ___|_U,|__|&)"_|SQ#W######@@UEN)._C;____X._|F?_|_.U|____&__|___|__\n');
fprintf(f,'__<E,__.|E"|___FF__|_:A!`_-&Y!)[LJ#WM########@M####6CS*|__{R__.U"_.|(T_|__T"___|___|\n');
fprintf(f,'_|KLITLMAA+TTFDUVF1L&PWITL1&HJ?"`X#$%#######@@V77YM9Q##YPC|R&|_.U:_._T|__|"R_|___|__\n');
fprintf(f,'___CE._|!XR_`-.)BV.|_*VL___|_|&|"#24#@#######W#&___|)QTPEY)=VT+|!U*|__Z>_._T!__|___|\n');
fprintf(f,'_|__>S-__|/Z||___(A(_|_(S.___|(LURVE@@#####W#K0K!._]E`__`|.1L}|FTV7C_|_E{|_-A|___|__\n');
fprintf(f,'___|_.Z*____V1_|__`&C_.|"X,|___T#CJA#W###WWBAMS#3_=N__.._}Z<___|_;XMJ&]<H{_|=L_|___|\n');
fprintf(f,'_|___|~B-|.__]A._|__FC_._}T__|-EMS_"W#WWWM#&B#AP#ZT__|__(Z`__|__CO}__.>/R5T;_L<__|__\n');
fprintf(f,'___|___.X,_|__(U.__|_R?|._&:__E#B__|U#####2|T#1`A#F|___;Z__|__.A]__|___}TUPLNCG{___|\n');
fprintf(f,'_|___|__<V___|_(A|___.Z`_|>]_|APLJ___FU5LVF_<#T]UGW?_|_O-|___.U}_|___~AF~|__.*$Y_|__\n');
fprintf(f,'___|___|_F)|___|II_|__==__-L_[#|.(S;_(?|__O|L#E.`.VH__>F___|.U*|___|~Z}|__.|.V+.___|\n');
fprintf(f,'_|___|___|A`_|___E;__|_R_|.](SQ`_|~U.1._`|C:UW-__|`Z)|L-_|__1?___|_.S:___|_=U*___|__\n');
fprintf(f,'___|___|__/F___|_"R|___T.__S&R5|___?1E1&T1FHZ0_|___|V=F|___:N__|___E|__|__LR___|___|\n');
fprintf(f,'_|___|___|_V`|___|R-_|_J.|TA]PY__|__YI!__*A3HL___:=FCSP__|_J{|___|=F_|___R]__|___|__\n');
fprintf(f,'___|___|___[>__|__)|___J.E[|_`|U"_.|&!_|___ROE_.&R?*~.F&__`&___|_~S|___|&=_|___|___|\n');
fprintf(f,'_|___|___|_;L|___|.F_|`T&?_._|_:S;_.&.___._1#),X>|___|_/S..F_|___T<__|_[F|___|___|__\n');
fprintf(f,'___|___|___|R__|__`L__*BT__|___|_Z."J__|___VK_V*___|___|"T]>___|;A_|__;E___|___|___|\n');
fprintf(f,'_|___|___|__L.___|.&_*S)|&V|_|___>F/"|{))(."V>F__|___|~<>IYC_|__J}___|&;_|___|___|__\n');
fprintf(f,'___|___|___|]*_|__.C>Z~|__-&&__|__&P[EJ?(=&LUS.|___|?VT]|{"NL___V__|_"J|___|___|___|\n');
fprintf(f,'_|___|___|__=}___|~$J`___|__+V___|(4|`___|_;YM___|.U1-___|__/N`?||___C!__|___|___|__\n');
fprintf(f,'___|___|___|=}_|_`UI.__|___|_/&|__]{___|___|,E_|_*O;___|___`!ZJK___|_T.|___|___|___|\n');
fprintf(f,'_|___|___|__[|__;U<;UI.__|___|&>_|&`_|___|_._&~_-X`__|___|INC".&L`__.J___|___|___|__\n');
fprintf(f,'___|___|___|I._FA._|_/O(___|__`O_*L|_`"<(!_|_[||E;_|___|}U/.___|?U)|`F_|___|___|___|\n');
fprintf(f,'_|___|___|__R)&/_|.__|.&R.___|_C</*;LE1/|CTN)-I{L|___|_+A____|___|]&]V___|___|___|__\n');
fprintf(f,'___|___|___|FJU?"__|___|?U-|___~NLLN|__|___-CESE.__|__CV___|___|}=CUWI_|___|___|___|\n');
fprintf(f,'_|___|___|___|~?&V?__|___;X*_|__ZN>__|___|___~QY_|___|&__|___:FAJ>.__|___|___|___|__\n');
fprintf(f,'___|___|___|___|_-1A{__|__-B:__|L<_|___|___|.__V___|"S_|___~TE>.___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|_.EF___|_.Z`_,&|___|.__.___|_R.|__U.___|=U<|___|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|_CE.___|>T_?]__~+RTJJR&].__+)__C/__|_VL|___|___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|___|)A.|___&.L_.FE=~___|~?VV;,&|-A_|_,O|___|___|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|_._)A__|_;Z["Z?|___|___|_{E]R_J}__.Z<__|___|___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|___|__=J___|POV~___|___|___|_/JFT__|A>_|___|___|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|___._V>|__IX__.|___|___|___|!5!|_F[|___|___|___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|___|___,U,_"B__|__.~*{!*-._|_._T:_!N___|___|___|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|___|__(J_L:|_{[&&FC[[C1&T]!|_>]|U__|___|___|___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|___|___|_V;F,[N&(~_|___|___(TV=_[?]|___|___|___|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|___|___*$X&[,__|___|___|___.C&WY___|___|___|___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|___|___|__DF___|___|___|___|___.&[_|___|___|___|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|\n');
fprintf(f,'_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__\n');
fprintf(f,'___|___|___|___|___|___|___|___|___|___|___|___|__\n');

