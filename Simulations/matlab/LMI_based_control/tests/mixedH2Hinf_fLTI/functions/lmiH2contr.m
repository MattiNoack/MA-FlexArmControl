% H2 CONTROLLER SYNTHESIS
% using cvx for LMI constrained optimization (SDP)

% Preliminaries:
%   - initiate cvx
%   - define system properly

% Input:
%       P           plant in state space form
%       indim       input channel dimensions (n2,nu)
%       outdim      output channel dimensions (nz2,ny)
%       Wof         output filter in state space form (signalwise)
%       options     case options
%                       - options(1)=integralON ... integrator in control
%                                                   included

% Output:
%       Ksys        resulting controller in state space form


function Ksys = lmiH2contr(P,indim,outdim,Wof,options)

    % Function input handling:
    if nargin<5
        options = 0; % no integrator
    end  

    % Include filters & define systems:
    [A,B2,Bu,Cz,Cy,Dz2,Dzu,Dy2,Dyu] = extend_system(P,indim,outdim,Wof,options(1));
    
    n = size(A,1);
    n2 = indim(1);    
    nu = indim(2);
    nz = outdim(1);
    ny = outdim(2);
    
    % Check soluability:
    if (Dyu~=0)
        disp('No control feedthrough allowed')
        Ksys=NaN; return
    elseif (Dz2~=0)
        disp('No feedthrough from w2 to z2 allowed')
        Ksys=NaN; return
    elseif (~check_cond(A,Bu,Cy))
        disp('System doesnt fulfill basic conditions')
        Ksys=NaN; return
    end

    % Semidefinite program:
    cvx_begin sdp
        % Variables:
        variable X(n,n) symmetric;
        variable Y(n,n) symmetric;
%         variable W(nz,nz) symmetric;
        variable alph;
%         alph = 200;
%         variable t;
        t = 40;
        variable Ah(n,n);
        variable Bh(n,ny);
        variable Ch(nu,n);                

        % Cost function:
%          minimize trace(W)        
%         maximize t
       
        % Internal stability:
        [X t*eye(n);t*eye(n) Y] > 0;
        
        [A*X+Bu*Ch+X*A'+Ch'*Bu' A+A';
         Ah+Ah' Y*A+Bh*Cy+A'*Y+Cy'*Bh'] < 0;
         
        % H2 conditions:
        [A*X+Bu*Ch+X*A'+Ch'*Bu' A+Ah' B2;
            A'+Ah Y*A+A'*Y+Bh*Cy+Cy'*Bh' Y*B2+Bh*Dy2;
            B2' B2'*Y + Dy2'*Bh' -eye(n2,n2)] < 0;
        
        [X eye(n) X*Cz'+Ch'*Dzu';
            eye(n) Y Cz';
            Cz*X+Dzu*Ch Cz alph*eye(nz)] > 0;

    cvx_end

    % Calculate controller parameterization:
    nk = n; % control order    

    N = eye(n)-Y*X; % selection in case of n=nk possible
    M = eye(n);

    Bk = N\Bh; % calculation in case of n=nk possible
    Ck = Ch/(M');
    Ak = N\(Ah-Y*A*X-N*Bk*Cy*X-Y*Bu*Ck*M')/(M');
    
    if(options(1))
        Ak = [Ak zeros(nk,nu);Ck zeros(nu,nu)];
        Bk = [Bk;zeros(nu,ny)];
        Ck = [zeros(nu,nk) eye(nu)];
    end

    disp(''); disp('Obtained output feedback controller:')
    Ksys = ss(Ak,Bk,Ck,zeros(nu,ny));
    
end


%---------------------%
% Auxiliary functions %
%---------------------%

function [A,B2,Bu,Cz,Cy,Dz2,Dzu,Dy2,Dyu] = extend_system(P,indim,outdim,Wof,int)

    % System matrices:  
    A = P.A;        
    n = size(A,1);
    n2 = indim(1);
    nu = indim(2);
    nz = outdim(1);
    ny = outdim(2);
    
    B2 = P.B(:,1:n2);
    Bu = P.B(:,1+n2:n2+nu);
    Cz = P.C(1:nz,:);
    Cy = P.C(1+nz:nz+ny,:);
    Dz2 = P.D(1:nz,1:n2);    
    Dzu = P.D(1:nz,1+n2:n2+nu);
    Dy2 = P.D(1+nz:nz+ny,1:n2);
    Dyu = P.D(1+nz:nz+ny,1+n2:n2+nu);
    
    % Integrator:
    if(int)
        A = [A Bu;zeros(nu,nu+n)];  
        Bu = [zeros(n,nu);eye(nu)];
        B2 = [B2;zeros(nu,n2)];
        Cz = [Cz Dzu];
        Dzu = zeros(nz,nu);
        Cy = [Cy Dyu];
        Dyu = zeros(ny,nu);
        
        n = n + nu;
    end
    
    % Filter:
    if(~isempty(Wof))
        % Filter matrices:
        Aof = Wof.A;
        nof = size(Aof,1);
        Bof = Wof.B;
        Cof = Wof.C;
        Dof = Wof.D;

        % Form extended system:
        A = [A zeros(n,nof);Bof*Cz Aof];

        B2 = [B2;Bof*Dz2];
        Bu = [Bu;Bof*Dzu];

        Cz = [Dof*Cz Cof];
        Cy = [Cy zeros(ny,nof)];

        Dz2 = Dof*Dz2;
        Dzu = Dof*Dzu;
    end

end

function flag = check_cond(A,Bu,Cy)

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

end


