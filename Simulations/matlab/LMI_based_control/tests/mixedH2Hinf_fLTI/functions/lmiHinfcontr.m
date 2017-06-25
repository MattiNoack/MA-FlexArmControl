% Hinf CONTROLLER SYNTHESIS
% using cvx for LMI constrained optimization (SDP)

% Preliminaries:
%   - initiate cvx
%   - define system properly

% Input:
%       P           plant in state space form
%       indim       input channel dimensions (ninf,nu)
%       outdim      output channel dimensions (nzinf,ny)
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


function [Ksys,eta] = lmiHinfcontr(P,indim,outdim,Wof,options,mu)

    % Function input handling:
    if nargin<5
        options = [0 0]; % subobtimal, no integrator
        mu = 10;
    elseif nargin<6
        mu = 10; % subobtimal, no integrator
    end  

    % Include filters & define systems:
    [A,Binf,Bu,Cz,Cy,Dzinf,Dzu,Dyinf,Dyu] = filterSys(P,indim,outdim,Wof);
    
    n = size(A,1);
    ninf = indim(1);    
    nu = indim(2);
    nz = outdim(1);
    ny = outdim(2);
    
    % Check soluability:
    if (Dyu~=0)
        disp('No control feedthrough allowed')
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

        % Cost function:
        if (options(1)) % optimal case
            variable eta; % mu^2
            minimize eta
        else % suboptimal case
            eta = mu^2;
        end            

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

function [A,Binf,Bu,Cz,Cy,Dzinf,Dzu,Dyinf,Dyu] = filterSys(P,indim,outdim,Wof)

    % System matrices:  
    A = P.A;        
    n = size(A,1);
    ninf = indim(1);
    nu = indim(2);
    nz = outdim(1);
    ny = outdim(2);
    
    Binf = P.B(:,1:ninf);
    Bu = P.B(:,1+ninf:ninf+nu);
    Cz = P.C(1:nz,:);
    Cy = P.C(1+nz:nz+ny,:);
    Dzinf = P.D(1:nz,1:ninf);    
    Dzu = P.D(1:nz,1+ninf:ninf+nu);
    Dyinf = P.D(1+nz:nz+ny,1:ninf);
    Dyu = P.D(1+nz:nz+ny,1+ninf:ninf+nu);
    
    if(~isempty(Wof))
        % Filter matrices:
        Aof = Wof.A;
        nof = size(Aof,1);
        Bof = Wof.B;
        Cof = Wof.C;
        Dof = Wof.D;

        % Form extended system:
        A = [A zeros(n,nof);Bof*Cz Aof];

        Binf = [Binf;Bof*Dzinf];
        Bu = [Bu;Bof*Dzu];

        Cz = [Dof*Cz Cof];
        Cy = [Cy zeros(ny,nof)];

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
