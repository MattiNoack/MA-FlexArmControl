% MULTI-OBJECTIVE ROBUST CONTROLLER SYNTHESIS
% via LMI constraints from [Scherer97] using cvx

% Preliminaries:
%   - initiate cvx
%   - define system properly

% Input:
%       P           plant in state space form
%       channels    structure for channel pair definitions (Remark: (z,w) first)
%                       - channels.contr ... controller [y:u]
%                       - channels.H2 ... generilized H2 [z2;w2]
%                       - channels.Hinf ... Hinf norm [zinf;winf]
%                       - channels.passive ... passivity [zp;wp]
%                       - channels.delta ... uncertainties [zd;wd]
%                       - channels.nomreg ... nominal regulation [znr;wnr]
%       Wof         output filter in state space form (signalwise)
%       options     structure of case options
%                       - options.optHinf ... optimal or suboptimal Hinf approach
%                       - options.integ ... integrator in control included
%       optin       structure for predefined optimization variables
%                       - optin.gamma ... Hinf suboptimal value
%                       - optin.alpha ... generalized H2 peak bound
%                       - optin.t ... well conditioning of lyapunov variables
%                       - optin.S ... signal generator for nominal regulation


% Output:
%       Ksys        resulting controller in state space form
%       optout      structure of characterstic variables for optimization
%                   (compare with 'optin')


function [Ksys,optout] = multiobjRC(P,channels,Wof,options,optin)

    % Function input handling:
    if nargin<5
        optin = struct;
    end    
    
    % Parameters:
    alph = optin.alpha;
    t = optin.t;
    gamm = optin.gamma;
    
    % Include filters & define systems:
    [A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu] = filterSys(P,channels,Wof,options.integ,optin.S);
    n = size(A,1);
    ny = size(Cy,1);
    nu = size(Bu,2);
    
    % Check soluability:
    if (~check_cond(A,Bu,Cy,Dyu,Dzw(channels.H2(1,:),channels.H2(2,:))))
        disp('System doesnt fulfill basic conditions')
        Ksys=ss([],[],[],[]); optout = optin; return        
    end   
    
    % Semidefinite program:
    cvx_begin sdp
        % Variables:
        variable X(n,n) symmetric;
        variable Y(n,n) symmetric;
        variable Ah(n,n);
        variable Bh(n,ny);
        variable Ch(nu,n);  
        variable Dh(nu,ny)
        
        % Auxilaries:
%         variable alph
%         variable gamm
%         variable W(n,n) symmetric;

        % Cost function:   
%         minimize alph
%         alph > 0;
%         minimize gamm
%         gamm > 0;
%         minimize trace(W)
       
%         variable t
%         maximize t
%         t > 0;
        
        % Internal stability:
        [X t*eye(n);t*eye(n) Y] > 0;
        
        [A*X+Bu*Ch+X*A'+Ch'*Bu' Ah'+(A+Bu*Dh*Cy);
         Ah+(A+Bu*Dh*Cy)' Y*A+Bh*Cy+A'*Y+Cy'*Bh'] < 0;
        
        % Generalized H2:            
        iz = channels.H2(1,:);
        iw = channels.H2(2,:);
        n2 = length(iz);

        B2 = Bw(:,iw);
        C2 = Cz(iz,:);
        Dy2 = Dyw(:,iw);
        D2u = Dzu(iz,:);
        Dz2 = Dzw(iz,iw);

        [A*X+Bu*Ch+X*A'+Ch'*Bu' A+Bu*Dh*Cy+Ah' B2+Bu*Dh*Dy2;
         A'+Cy'*Dh'*Bu'+Ah Y*A+A'*Y+Bh*Cy+Cy'*Bh' Y*B2+Bh*Dy2;
         B2'+Dy2'*Dh'*Bu' B2'*Y+Dy2'*Bh' -eye(n2)] < 0;

        [X eye(n) X*C2'+Ch'*D2u';
         eye(n) Y C2'+Cy'*Dh'*D2u';
         C2*X+D2u*Ch C2+D2u*Dh*Cy alph*eye(n2)] > 0;            

        Dz2+D2u*Dh*Dy2 == 0;     
        
        % Hinf performance:            
        iz = channels.Hinf(1,:);
        iw = channels.Hinf(2,:);
        ninf = length(iz);

        Binf = Bw(:,iw);
        Cinf = Cz(iz,:);
        Dyinf = Dyw(:,iw);
        Dinfu = Dzu(iz,:);
        Dzinf = Dzw(iz,iw);

        [A*X+X*A'+Bu*Ch+Ch'*Bu' Ah'+(A+Bu*Dh*Cy) Binf+Bu*Dh*Dyinf (Cinf*X+Dinfu*Ch)';
         Ah+(A+Bu*Dh*Cy)' A'*Y+Y*A+Bh*Cy+Cy'*Bh' Y*Binf+Bh*Dyinf (Cinf+Dinfu*Dh*Cy)';
         (Binf+Bu*Dh*Dyinf)' (Y*Binf+Bh*Dyinf)' -gamm*eye(ninf) (Dzinf+Dinfu*Dh*Dyinf)';
         Cinf*X+Dinfu*Ch Cinf+Dinfu*Dh*Cy Dzinf+Dinfu*Dh*Dyinf -gamm*eye(ninf)] < 0;
        
        % Passive behavior:            
        iz = channels.passive(1,:);
        iw = channels.passive(2,:);
        np = length(iz);

        Bp = Bw(:,iw);
        Cp = Cz(iz,:);
        Dyp = Dyw(:,iw);
        Dpu = Dzu(iz,:);
        Dp = Dzw(iz,iw);

        [A*X+Bu*Ch+X*A'+Ch'*Bu' A+Bu*Dh*Cy+Ah' (Bp+Bu*Dh*Dyp)-(Cp*X+Dpu*Ch)';
         A'+Cy'*Dh'*Bu'+Ah Y*A+A'*Y+Bh*Cy+Cy'*Bh' (Y*Bp+Bh*Dyp)-(Cp+Dpu*Dh*Cy)';
        (Bp+Bu*Dh*Dyp)'-(Cp*X+Dpu*Ch) (Y*Bp+Bh*Dyp)'-(Cp+Dpu*Dh*Cy) -(Dp+Dpu*Dh*Dyp)-(Dp+Dpu*Dh*Dyp)'] < 0;        
                
        % Robust stability: [...]

    cvx_end
    
    optout = optin;
    optout.alpha = alph;
    optout.gamma = gamm;
    
    % Controller computation:
    nk = n;
    
    N = eye(n)-Y*X; % selection in case of n=nk possible
    M = eye(n);
    
    Dk = Dh;
    Ck = (Ch-Dk*Cy*X)/M;
    Bk = N\(Bh-Y*Bu*Dk);
    Ak = N\(Ah-N*Bk*Cy*X-Y*Bu*Ck*M'-Y*(A+Bu*Dk*Cy)*X)/M;
    
    if(options.integ) % integrator
        Ak = [Ak Bk;zeros(ny,nk) zeros(nu,nu)];
        Bk = [zeros(nk,ny);eye(ny)];
        Ck = [Ck Dk];
        Dk = zeros(nu,ny);
    end

    Ksys = ss(Ak,Bk,Ck,Dk);
    
end


%---------------------%
% Auxiliary functions %
%---------------------%

function [A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu] = filterSys(P,channels,Wof,integ,S)


    % Channel indices:
    iy = channels.contr(1,:);
    iu = channels.contr(2,:);
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
    nu = size(Bu,2);
    ny = size(Cy,1);
    nw = size(Bw,2);
    nz = size(Cz,2);
    
    % Nominal regulation:
    if(~isempty(channels.nomreg))
        % Channel definition:
        iznr = channels.nomreg(1,:);
        iwnr = channels.nomreg(2,:);
        nnr = length(iznr);    

        Bnr = Bw(:,iwnr);
        Cnr = Cz(iznr,:);
        Dynr = Dyw(:,iwnr);
        Dnru = Dzu(iznr,:);
        Dnr = Dzw(iznr,iwnr);

        % Check requirements:
        if(check_nomref())
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
        end
    end
    
    % Integrator:
    if(integ)
        ny = length(iy); nu = length(iu);
        A = [A zeros(n,ny);Cy zeros(ny,ny)];  
        Bu = [Bu;zeros(ny,nu)];
        Bw = [Bw;Dyw];
        Cz = [Cz zeros(length(iz),ny)];  
        Cy = [zeros(ny,n) eye(ny)];
        Dyw = zeros(ny,length(iw));        
        Dyu = zeros(ny,nu);
        
        n = n + ny;
    end
    
    % Output filter:
    if(~isempty(Wof))
        % Filter matrices:
        Aof = Wof.A;
        nof = size(Aof,1);
        Bof = Wof.B;
        Cof = Wof.C;
        Dof = Wof.D;

        % Form extended system:
        A = [A zeros(n,nof);Bof*Cz Aof];

        Bw = [Bw;Bof*Dzw];        
        Bu = [Bu;Bof*Dzu];

        Cz = [Dof*Cz Cof];
        Cy = [Cy zeros(length(iy),nof)];

        Dzw = Dof*Dzw;        
        Dzu = Dof*Dzu;
    end
        
end

function flag = check_nomref(S,A,Bw,Cz,Cy,Dzw,Dyw,Dzu)

    % Condition list:
    c1 = Cz == Cy;
    c2 = Dzw == Dyw;
    c3 = Dzu == zeros(size(Dzu));
    c4 = eig(S) <= 0;
    
    Aaux = [A Bw;zeros(size(S,1),size(A,2)) S];
    Caux = [Cy Dyw];
    O = obsv(Aaux,Caux);
    c5 = rank(O)==(size(S,1)+size(A,2));
    
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
