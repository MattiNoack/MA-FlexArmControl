% MIXED H2/Hinf CONTROLLER SYNTHESIS
% using cvx for LMI constrained optimization (SDP)

% Preliminaries:
%   - initiate cvx
%   - define system properly

% Input:
%       P           plant in state space form
%       indim       input channel dimensions (n2,ninf,nu)
%       outdim      output channel dimensions (nz,ny)
%       Wof         output filter in state space form (signalwise)
%       options     case options
%                       - options(1)=optimalON ... optimal or suboptimal
%                                                  Hinf approach
%                       - options(2)=integralON ... integrator in control
%                                                   included
%       mu          optimization constant given in suboptimal case

% Output:
%       Ksys        resulting controller in state space form
%       eta         Hinf bound optimal value mu^2


function [Ksys,eta] = mixH2Hinf(P,indim,outdim,Wof,options,mu)

    % Function input handling:
    if nargin<5
        options = [0 0]; % subobtimal, no integrator
        mu = 10;
    elseif nargin<6
        mu = 10; % subobtimal, no integrator
    end  
    
    % Include filters & define systems:
    [A,B2,Binf,Bu,Cz,Cy,Dz2,Dzinf,Dzu,Dy2,Dyinf,Dyu] = filterSys(P,indim,outdim,Wof);
    
    n = size(A,1);
    n2 = indim(1);
    ninf = indim(2);
    nu = indim(3);
    nz = outdim(1);
    ny = outdim(2);
    
    % Check soluability:
    if (Dyu~=0)
        disp('No control feedthrough allowed')
        Ksys=NaN; eta=-1; return
    elseif (Dz2~=0)
        disp('No feedthrough from w2 to z allowed')
        Ksys=NaN; eta=-1; return
    elseif (~check_cond(A,Bu,Cy))
        disp('System doesnt fulfill basic conditions')
        Ksys=NaN; return        
    end   
    
    % Semidefinite program:
    cvx_begin sdp %quiet
        % Variables:
        variable X(n,n) symmetric;
        variable Y(n,n) symmetric;
        variable W(nz,nz) symmetric;
        variable Ah(n,n);
        variable Bh(n,ny);
        variable Ch(nu,n);        
        
        
        % Integral state extension:
        if (options(2)) % integral action
            
        end

        % Cost function:
        if (options(1)) % optimal case
            variable eta; % mu^2
            minimize eta+trace(W)
        else % suboptimal case
            eta = mu^2;
            minimize trace(W)
        end        
       
        % H2 conditions:
        [W, Cz*X + Dzu*Ch, Cz;
        X*Cz' + Ch'*Dzu', X, eye(n,n);
        Cz', eye(n,n), Y] > 0;
        
        [A*X + Bu*Ch + X*A' + Ch'*Bu', A+Ah', B2;
        A'+Ah, Y*A + A'*Y + Bh*Cy + Cy'*Bh', Y*B2 + Bh*Dy2;
        B2', B2'*Y + Dy2'*Bh', -eye(n2,n2)] < 0;

        % Hinf conditions:    
        [X, eye(n,n);
        eye(n,n), Y] > 0;
    
        [A*X + Bu*Ch + X*A' + Ch'*Bu', A+Ah', Binf, X*Cz' + Ch'*Dzu';
        A'+Ah, Y*A + A'*Y + Bh*Cy + Cy'*Bh', Y*Binf + Bh*Dyinf, Cz';
        Binf', Binf'*Y + Dyinf'*Bh', -eye(ninf,ninf), Dzinf';
        Cz*X + Dzu*Ch, Cz, Dzinf, -eta*eye(nz,nz)] < 0;
    cvx_end

    % Calculate controller parameterization:
    nk = n; % control order    

    N = eye(n)-Y*X; % selection in case of n=nk possible
    M = eye(n);

    Bk = N\Bh; % calculation in case of n=nk possible
    Ck = Ch/(M');
    Ak = N\(Ah-Y*A*X-N*Bk*Cy*X-Y*Bu*Ck*M')/(M');

    disp(''); disp('Obtained output feedback controller:')
    Ksys = ss(Ak,Bk,Ck,zeros(nu,ny));
    
end


%---------------------%
% Auxiliary functions %
%---------------------%

function [A,B2,Binf,Bu,Cz,Cy,Dz2,Dzinf,Dzu,Dy2,Dyinf,Dyu] = filterSys(P,indim,outdim,Wof)

    % System matrices:  
    A = P.A;        
    n = size(A,1);
    n2 = indim(1);
    ninf = indim(2);
    nu = indim(3);
    nz = outdim(1);
    ny = outdim(2);
    
    B2 = P.B(:,1:n2);
    Binf = P.B(:,1+n2:n2+ninf);
    Bu = P.B(:,1+n2+ninf:n2+ninf+nu);
    Cz = P.C(1:nz,:);
    Cy = P.C(1+nz:nz+ny,:);
    Dz2 = P.D(1:nz,1:n2);
    Dzinf = P.D(1:nz,1+n2:n2+ninf);
    Dzu = P.D(1:nz,1+n2+ninf:n2+ninf+nu);
    Dy2 = P.D(1+nz:nz+ny,1:n2);
    Dyinf = P.D(1+nz:nz+ny,1+n2:n2+ninf);
    Dyu = P.D(1+nz:nz+ny,1+n2+ninf:n2+ninf+nu);
    
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
        Binf = [Binf;Bof*Dzinf];
        Bu = [Bu;Bof*Dzu];

        Cz = [Dof*Cz Cof];
        Cy = [Cy zeros(ny,nof)];

        Dz2 = Dof*Dz2;
        Dzinf = Dof*Dzinf;
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
