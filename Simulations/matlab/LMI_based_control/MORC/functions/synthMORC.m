% MULTI-OBJECTIVE ROBUST CONTROLLER SYNTHESIS
% via LMI constraints from [Scherer97] using cvx toolbox

% Preliminaries:
%   - install cvx
%   - define system properly
%   - pick proper channel names ('z','y','w','u')

% Input:
%       P           plant in state space form
%       channels    structure for channel pair definitions with in/out distinction
%                       - channels.contr ... controller (in:u,out:y)
%                       - channels.integ ... integrator (meas:yi)
%                       - channels.H2 ... generilized H2 (in:w2,out:z2)
%                       - channels.Hinf ... Hinf norm (in:winf,out:zinf)
%                       - channels.passive ... passivity (in:wp,out:zp)
%                       - channels.delta ... uncertainties (in:wd,out:zd)
%                       - channels.nomreg ... nominal regulation (in:wnr,out:znr)
%       W           filter structure for system extension (signalwise, state space form)
%                       - W.Wof ... output signal filter
%                       - W.Wif ... intput signal filter
%                       - W.Wdf ... uncertainty filter
%       options     structure of case options
%                       - options.optHinf ... optimal or suboptimal Hinf approach
%                       - options.optGH2 ... optimal generalized H2 norm
%                       - options.optLCond ... optimal Lyapunov conditioning
%                       - options.normCond ... optimal norm conditioning
%       optin       structure for predefined optimization variables
%                       - optin.gamma ... Hinf suboptimal value
%                       - optin.alpha ... generalized H2 peak bound
%                       - optin.t ... well conditioning of lyapunov variables
%                       - optin.S ... signal generator for nominal regulation
%                       - optin.b ... conditioning of matrix variables norm


% Output:
%       Ksys        resulting controller in state space form
%       optout      structure of characterstic variables for optimization
%                   (compare with 'optin')
%       CLsys       resulting closed-loop system representation
%       OLctrl      open-loop system transfer function in SISO control channel


function [Ksys,optout,CLsys,OLctrl] = synthMORC(P,channels,W,options,optin)

    % Include filters & define systems:
    [A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,V,Pext] = extendSys(P,channels,W,optin.S);
    n = size(A,1);
    ny = size(Cy,1); nu = size(Bu,2);
    nz = size(Cz,1); nw = size(Bw,2);
    
    % Check soluability:
    if (~check_cond(A,Bu,Cy,Dyu,Dzw(channels.H2.out,channels.H2.in)))
        disp('System doesnt fulfill basic conditions')
        Ksys=ss([],[],[],[]); optout = optin; return        
    end   
    
    % Semidefinite program:
    epsilon = 1e-3;% lmi conditioning
    cvx_begin sdp
        % Variables:
        variable X(n,n) symmetric;
        variable Y(n,n) symmetric;
        variable Ah(n,n);
        variable Bh(n,ny);
        variable Ch(nu,n);  
        variable Dh(nu,ny);
        
        % Cost function:
        if(options.optHinf)
            variable gamm;
            gamm >= epsilon;
        else
            gamm = optin.gamma;
        end
        if(options.optGH2)
            variable alph;
            alph >= epsilon;
        else
            alph = optin.alpha;
        end
        if(options.optLCond)
            variable t;
            t >= epsilon;
        else
            t = optin.t;
        end
        if(options.normCond)
            variable b;
            X-b*eye(n)<=-epsilon*eye(n);
            Y-b*eye(n)<=-epsilon*eye(n);
            [b*eye(nu) Ch;Ch' b*eye(n)]>=epsilon*eye(n+nu);
            [b*eye(nu) Dh;Dh' b*eye(ny)]>=epsilon*eye(ny+nu);
        else
            b = optin.b;
        end        
              
        minimize gamm + alph - t + b
                
        % Internal stability:
        [X t*eye(n);t*eye(n) Y] >= epsilon*eye(2*n);
        
        [A*X+Bu*Ch+X*A'+Ch'*Bu' Ah'+(A+Bu*Dh*Cy);
         Ah+(A+Bu*Dh*Cy)' Y*A+Bh*Cy+A'*Y+Cy'*Bh'] <= -epsilon*eye(2*n);
     
        % Generalized H2:
        if(~isempty(channels.H2.out))
            iz = channels.H2.out;
            iw = channels.H2.in;

            B2 = Bw(:,iw);
            C2 = Cz(iz,:);
            Dy2 = Dyw(:,iw);
            D2u = Dzu(iz,:);
            Dz2 = Dzw(iz,iw);

            [A*X+Bu*Ch+X*A'+Ch'*Bu' A+Bu*Dh*Cy+Ah' B2+Bu*Dh*Dy2;
             A'+Cy'*Dh'*Bu'+Ah Y*A+A'*Y+Bh*Cy+Cy'*Bh' Y*B2+Bh*Dy2;
             B2'+Dy2'*Dh'*Bu' B2'*Y+Dy2'*Bh' -eye(length(iw))] <= -epsilon*eye(2*n+length(iw));

            [X eye(n) X*C2'+Ch'*D2u';
             eye(n) Y C2'+Cy'*Dh'*D2u';
             C2*X+D2u*Ch C2+D2u*Dh*Cy alph*eye(length(iz))] >= epsilon*eye(2*n+length(iz));            

            Dz2+D2u*Dh*Dy2 == 0;
        end
        
        % Hinf performance:
        if(~isempty(channels.Hinf.out))
            iz = channels.Hinf.out;
            iw = channels.Hinf.in;

            Binf = Bw(:,iw);
            Cinf = Cz(iz,:);
            Dyinf = Dyw(:,iw);
            Dinfu = Dzu(iz,:);
            Dzinf = Dzw(iz,iw);

            [A*X+X*A'+Bu*Ch+Ch'*Bu' Ah'+(A+Bu*Dh*Cy) Binf+Bu*Dh*Dyinf (Cinf*X+Dinfu*Ch)';
             Ah+(A+Bu*Dh*Cy)' A'*Y+Y*A+Bh*Cy+Cy'*Bh' Y*Binf+Bh*Dyinf (Cinf+Dinfu*Dh*Cy)';
             (Binf+Bu*Dh*Dyinf)' (Y*Binf+Bh*Dyinf)' -gamm*eye(length(iw)) (Dzinf+Dinfu*Dh*Dyinf)';
             Cinf*X+Dinfu*Ch Cinf+Dinfu*Dh*Cy Dzinf+Dinfu*Dh*Dyinf -gamm*eye(length(iz))] ...
                <= -epsilon*eye(2*n+length(iw)+length(iz));
        end
     
        % Passive behavior: (nz=nw)
        if(~isempty(channels.passive.out))
            iz = channels.passive.out;
            iw = channels.passive.in;
            np = length(iz);

            Bp = Bw(:,iw);
            Cp = Cz(iz,:);
            Dyp = Dyw(:,iw);
            Dpu = Dzu(iz,:);
            Dp = Dzw(iz,iw);

            ainp = 1;
            
            [A*X+Bu*Ch+X*A'+Ch'*Bu' A+Bu*Dh*Cy+Ah' (Bp+Bu*Dh*Dyp)-(Cp*X+Dpu*Ch)';
             A'+Cy'*Dh'*Bu'+Ah Y*A+A'*Y+Bh*Cy+Cy'*Bh' (Y*Bp+Bh*Dyp)-(Cp+Dpu*Dh*Cy)';
            (Bp+Bu*Dh*Dyp)'-(Cp*X+Dpu*Ch) (Y*Bp+Bh*Dyp)'-(Cp+Dpu*Dh*Cy) ainp*eye(np)-(Dp+Dpu*Dh*Dyp)-(Dp+Dpu*Dh*Dyp)'] ...
                <= -epsilon*eye(2*n+length(iw));
        end
        
        % Robust stability:
        if(~isempty(channels.delta.out))
            iz = channels.delta.out;
            iw = channels.delta.in;

            Bd = Bw(:,iw);
            Cd = Cz(iz,:);
            Dyd = Dyw(:,iw);
            Ddu = Dzu(iz,:);
            Dd = Dzw(iz,iw);

            [A*X+X*A'+Bu*Ch+Ch'*Bu' Ah'+(A+Bu*Dh*Cy) Bd+Bu*Dh*Dyd (Cd*X+Ddu*Ch)';
             Ah+(A+Bu*Dh*Cy)' A'*Y+Y*A+Bh*Cy+Cy'*Bh' Y*Bd+Bh*Dyd (Cd+Ddu*Dh*Cy)';
             (Bd+Bu*Dh*Dyd)' (Y*Bd+Bh*Dyd)' -1*eye(length(iw)) (Dd+Ddu*Dh*Dyd)';
             Cd*X+Ddu*Ch Cd+Ddu*Dh*Cy Dd+Ddu*Dh*Dyd -1*eye(length(iz))] ...
                <= -epsilon*eye(2*n+length(iw)+length(iz));
        end
        
    cvx_end
    
    optout = optin;
    optout.gamma = gamm;
    optout.alpha = alph;
    optout.t = t;
    
    % Controller computation:
    nk = n; % easy choice
    
    N = eye(n)-Y*X; % selection in case of n=nk possible
    M = eye(n);
    
    Dk = Dh;
    Ck = (Ch-Dk*Cy*X)/M;
    Bk = N\(Bh-Y*Bu*Dk);
    Ak = N\(Ah-N*Bk*Cy*X-Y*Bu*Ck*M'-Y*(A+Bu*Dk*Cy)*X)/M;
    
    if(~isempty(channels.nomreg.in)) % internal model
        Ck1 = Ck(:,1:nu); Ck2 = Ck(:,nu+1:nu+size(S,1));
        Dk1 = Dk(:,1:nu); Dk2 = Dk(:,nu+1:nu+size(S,1));
        Ak = [S Ck2;zeros(nk,size(S,2)) Ak];
        Bk = [Dk2;Bk];
        Ck = [V Ck1];
        Dk = Dk1;
    end
    
    if(~isempty(channels.integ.meas)) % integrator
        iyi = channels.integ.meas; nyi = length(iyi);
        iyn = setdiff(1:ny,iyi); nyn = length(iyn);
        
        Ak = [Ak Bk(:,iyi);zeros(nyi,nk) zeros(nyi,nyi)];
        Bk = [zeros(nk,nyi) Bk(:,iyn);eye(nyi) zeros(nyi,nyn)];
        Ck = [Ck Dk(:,nyi)];
        Dk = [zeros(nu,nyi) Dk(:,iyn)];
    end
    
    % Result & closed-loop system:
    Ksys = ss(Ak,Bk,Ck,Dk);    
    Ksys.u = 'y'; Ksys.y = 'u'; 
    CLsys = connect(Pext,Ksys,'w','z');    
    if(ny==1)&&(nu==1) % SISO test
        OLctrl = -tf(ss(P.A,P.B(:,channels.contr.in),P.C(channels.contr.out,:)...
                    ,P.D(channels.contr.out,channels.contr.in)))*tf(Ksys);
    else
        OLctrl = tf(0,1);
    end
    
end


%---------------------%
% Auxiliary functions %
%---------------------%

function [A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,V,Pext] = extendSys(P,channels,W,S)

    % Filter inclusion: (no in/out dimension change)
    if(isempty(W.Wdf)) % no uncertainty
        P = connect(P,W.Wof,{'w','u'},{'zf','y'});
        P = connect(P,W.Wif,{'wf','u'},{'zf','y'});
    else % uncertainty considered
        P = connect(P,W.Wof,{'wd','w','u'},{'zd','zf','y'});
        P = connect(P,W.Wif,{'wd','wf','u'},{'zd','zf','y'});
        P = connect(P,W.Wdf,{'wd','wf','u'},{'zdf','zf','y'});
    end
    
    % Channel selection:
    iy = channels.contr.out;
    iu = channels.contr.in;
    iz = 1:iy(1)-1;
    iw = 1:iu(1)-1;

    % System matrices:  
    A = P.A;
    Bw = P.B(:,iw);
    Bu = P.B(:,iu);
    Cz = P.C(iz,:);
    Cy = P.C(iy,:);
    Dzw = P.D(iz,iw);
    Dzu = P.D(iz,iu);
    Dyw = P.D(iy,iw);
    Dyu = P.D(iy,iu);
    
    n = size(A,1);
    nu = size(Bu,2); ny = size(Cy,1);
    nw = size(Bw,2); nz = size(Cz,1);    
    
    % Nominal regulation:
    if(~isempty(channels.nomreg.in))
        % Channel definition:
        iznr = channels.nomreg.out;
        iwnr = channels.nomreg.in;
        nnr = length(iznr);    

        Bnr = Bw(:,iwnr);
        Cnr = Cz(iznr,:);
        Dynr = Dyw(:,iwnr);
        Dnru = Dzu(iznr,:);
        Dnr = Dzw(iznr,iwnr);

        % Check requirements:
        if(check_nomreg(S,A,Bnr,Cnr,Cy,Dnr,Dynr,Dnru))
            % Solve linear condition:
            linsys = @(UV)([A Bu;Cy zeros(ny,nu)]*UV-[UV(1:n,:);zeros(ny,nnr)]*S-[Bnr;Dynr]);
            UVsol = fsolve(linsys,zeros(n+nu,nnr));
            V = UVsol(n+1:end,:);

            % Extended plant:
            A = [A Bu*V;zeros(nnr,n) S];
            Bw = [Bw;zeros(nnr,size(Bw,2))];
            Bu = [Bu zeros(n,nnr);zeros(nnr,nu) eye(nnr)];
            Cz = [Cz Dzu*V];
            Cy = [Cy zeros(ny,nnr)];
            Dzu = [Dzu zeros(nz,nnr)];
            
            n = size(A,1);
            nu = size(Bu,2); ny = size(Cy,1);
            nw = size(Bw,2); nz = size(Cz,1);
        else
            V = -1;
        end
    else
        V = -1;
    end
    
    % Integrator:
    if(~isempty(channels.integ.meas))
        iyi = channels.integ.meas; nyi = length(iyi);
        iyn = setdiff(1:ny,iyi); nyn = length(iyn);

        A = [A zeros(n,nyi);Cy(iyi,:) zeros(nyi,nyi)];  
        Bu = [Bu;zeros(nyi,nu)];
        Bw = [Bw;Dyw(iyi,:)];
        Cz = [Cz zeros(nz,nyi)];  
        Cy = [zeros(nyi,n) eye(nyi);Cy(iyn,:) zeros(nyn,nyi)];
        Dyw = [zeros(nyi,nw);Dyw(iyn,:)];
        Dyu = zeros(ny,nu);

        n = n + nyi;
    end           
    
    % Extended plant & signal labels:
    Pext = ss(A,[Bw Bu],[Cz;Cy],[Dzw Dzu;Dyw Dyu]);
    Pext.u = 'w';
    Pext.u{end} = 'u';
    Pext.y = 'z';
    if(ny==1)
        Pext.y{end} = 'y';
    elseif(ny==2)
        Pext.y{end-(size(Cy,1)-1)} = 'y(1)';
        Pext.y{end} = 'y(2)';
    end
    
end

function flag = check_nomreg(S,A,Bnr,Cnr,Cy,Dnr,Dynr,Dnru)

    % Condition list:
    c1 = Cnr == Cy;
    c2 = Dnr == Dynr;
    c3 = Dnru == zeros(size(Dnru));
    c4 = eig(S) >= 0; % non-trivial
    
    % Detectability of :
    Aaux = [A Bnr;zeros(size(S,1),size(A,2)) S];
    Caux = [Cy Dynr];
    
    tmpSys = ss(Aaux,zeros(size(Aaux,1),1),Caux,zeros(size(Caux,1),1)); % with dummy input
    [tmpS,tmpU] = stabsep(tmpSys); % seperate Stable/Unstable part
    
    O = obsv(tmpU.A,tmpU.C); % check unobservable+unstable modes
    c5 = rank(O)==size(tmpU.A,1);
    
    % Check conditions:
    flag = 1;
    if~(c1&c2&c3&c4&c5)
        flag = 0;
        disp('Nominal regulation not possible due to bad conditions.')
    end

end

function flag = check_cond(A,Bu,Cy,Dyu,Dz2)

    % Controllability and observability matrices:
    C = ctrb(A,Bu);
    O = obsv(A,Cy);
    
    % Check condition:
    flag = 1;
    epsilon = 1e-7;
    
    if(rank(C)<epsilon)
        flag = 0;
        disp('Control system not controllable!')
    end
    if(rank(O)<epsilon)
        flag = 0;
        disp('Control system not observable!')
    end
    if (max(max(abs(Dyu)))~=0) % check if element ~=0 exists
        flag = 0;
        disp('No control feedthrough allowed!')
    end
    if (max(max(abs(Dz2)))~=0)
        flag = 0;
        disp('No feedthrough from w2 to z allowed!')
    end

end
