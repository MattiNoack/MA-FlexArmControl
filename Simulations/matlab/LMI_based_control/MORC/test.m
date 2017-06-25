% TEST: check algorithm

% Get representation: (return P)
[P,channels] = system_def(100);

% Filter selection: (return Wof,Wif,Wdf)
W = weight_filter();

% Specify optimization parameters:
[def_channels,options,optin] = initMORC();
options.integ = 1;
integ = options.integ;
optin.t = 1;
optin.alpha = 50^2;
optin.gamma = 8;

% PROGRAM

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
            
            n = size(A,1);
            nu = size(Bu,2); ny = size(Cy,1);
            nw = size(Bw,2); nz = size(Cz,1);
        end
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
    
    n = size(A,1);
    ny = size(Cy,1); nu = size(Bu,2);
    nz = size(Cz,1); nw = size(Bw,2);
    
    % Optimization parameters:
    alph = optin.alpha;
    t = optin.t;
    gamm = optin.gamma;
    
    X = 1*eye(n); Y = 1*eye(n);
    Ah = magic(n); Bh = ones(n,ny); Ch = ones(nu,n); Dh = ones(nu,ny);
    
    % CHECK: dimensions (quadratic)
    size([X t*eye(n);t*eye(n) Y])
    size([A*X+Bu*Ch+X*A'+Ch'*Bu' A+A';
         Ah+Ah' Y*A+Bh*Cy+A'*Y+Cy'*Bh'])
     
     % CHECK: symmetry
    [X t*eye(n);t*eye(n) Y] == [X t*eye(n);t*eye(n) Y]'
    [A*X+Bu*Ch+X*A'+Ch'*Bu' Ah'+(A+Bu*Dh*Cy);
         Ah+(A+Bu*Dh*Cy)' Y*A+Bh*Cy+A'*Y+Cy'*Bh']==[A*X+Bu*Ch+X*A'+Ch'*Bu' Ah'+(A+Bu*Dh*Cy);
         Ah+(A+Bu*Dh*Cy)' Y*A+Bh*Cy+A'*Y+Cy'*Bh']'
     
     % CHECK: integrator selection
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
     
    % CHECK: passivity    
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

        tmp_pas=[A*X+Bu*Ch+X*A'+Ch'*Bu' A+Bu*Dh*Cy+Ah' (Bp+Bu*Dh*Dyp)-(Cp*X+Dpu*Ch)';
         A'+Cy'*Dh'*Bu'+Ah Y*A+A'*Y+Bh*Cy+Cy'*Bh' (Y*Bp+Bh*Dyp)-(Cp+Dpu*Dh*Cy)';
        (Bp+Bu*Dh*Dyp)'-(Cp*X+Dpu*Ch) (Y*Bp+Bh*Dyp)'-(Cp+Dpu*Dh*Cy) ainp*eye(np)-(Dp+Dpu*Dh*Dyp)-(Dp+Dpu*Dh*Dyp)'];
    
        eig(tmp_pas)
    end
    
    
%% Passivity example SDP:
clear

% Test system:
A = [0 1;-0.5 0];
Bu = [0;1];
Cy = [0 1];
Dyu = 0;
n = size(A,1);
ny = size(Cy,1);
nu = size(Bu,2);

Bp = [0;1];
Cp = Cy;
Dyp = 0;
Dpu = Dyu;
Dp = 0;
np = 1;

G = tf(ss(A,Bp,Cp,Dp));
isminphase(G.num{:},G.den{:})

% Passivity test:
epsilon = 1e-3;% lmi conditioning
cvx_begin sdp
    variable P(n,n) symmetric
%     [A'*P+P*A P*Bp-Cp';Bp'*P-Cp -(Dp+Dp')]<=-epsilon*eye(n+np);    
    A'*P+P*A<=0;%-epsilon*eye(n);
    Bp'*P==Cp;
cvx_end
% eig([A'*P+P*A P*B-C';B'*P-C -(D+D')])

% Optimization:
ainp = 0;

cvx_begin sdp
    % Variables:
    variable X(n,n) symmetric;
    variable Y(n,n) symmetric;
    variable Ah(n,n);
    variable Bh(n,ny);
    variable Ch(nu,n);  
    variable Dh(nu,ny);
    
    [X eye(n);eye(n) Y] >= epsilon*eye(2*n);
        
    [A*X+Bu*Ch+X*A'+Ch'*Bu' Ah'+(A+Bu*Dh*Cy);
     Ah+(A+Bu*Dh*Cy)' Y*A+Bh*Cy+A'*Y+Cy'*Bh'] <= -epsilon*eye(2*n);
    
    [A*X+Bu*Ch+X*A'+Ch'*Bu' A+Bu*Dh*Cy+Ah' (Bp+Bu*Dh*Dyp)-(Cp*X+Dpu*Ch)';
     A'+Cy'*Dh'*Bu'+Ah Y*A+A'*Y+Bh*Cy+Cy'*Bh' (Y*Bp+Bh*Dyp)-(Cp+Dpu*Dh*Cy)';
    (Bp+Bu*Dh*Dyp)'-(Cp*X+Dpu*Ch) (Y*Bp+Bh*Dyp)'-(Cp+Dpu*Dh*Cy) ainp*eye(np)-(Dp+Dpu*Dh*Dyp)-(Dp+Dpu*Dh*Dyp)'] ...
        <= -epsilon*eye(2*n+np);
cvx_end

% Extra system output with relative deg 1:
shift = [zeros(2*N+1,1) eye(2*N+1);zeros(1,2*(N+1))]; % shift to velocities
Cz = [C1+gamm*C2;C2;zeros(1,n);(C1+gamm*C2)*shift];

