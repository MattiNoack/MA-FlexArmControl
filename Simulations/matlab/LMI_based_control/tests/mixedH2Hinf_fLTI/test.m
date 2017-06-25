% DEBUGGING

clear
clc


[P,channels] = system_def();
Wof = weight_filter();

[def_channels,options,optin] = initMORC();
options.integ = 1;
optin.t = 1;
optin.alpha = 200;
optin.gamma = 10;


% Parameters:
alph = optin.alpha;
t = optin.t;
gamm = optin.gamma;
integ = options.integ;


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


cvx_begin sdp
    % Variables:
    variable X(n,n) symmetric;
    variable Y(n,n) symmetric;
    variable Ah(n,n);
    variable Bh(n,ny);
    variable Ch(nu,n);  
    variable Dh(nu,ny)

    % Internal stability:
    [X t*eye(n);t*eye(n) Y] > 0;

    [A*X+Bu*Ch+X*A'+Ch'*Bu' A+A';
     Ah+Ah' Y*A+Bh*Cy+A'*Y+Cy'*Bh'] < 0;

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

cvx_end
