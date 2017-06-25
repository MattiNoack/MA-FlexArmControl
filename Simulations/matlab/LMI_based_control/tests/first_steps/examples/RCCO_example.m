%EXAMPLE FROM Robust Control & Convex Optimization
% Hinf control with state feedback

clear all
% close all


%% System definition
% Input:    (w,u)
% Output:   (e,y)

% State space representation:
A = [0 1;-1 -2];
Bu = [0;1];
Bw = [0;1];
Ce = [1 1];
Cy = eye(2,2); % state FB
Dew = zeros(1,1);
Deu = 0;

% System parameters:
n = size(A,1);
nw = size(Bw,2);
nu = size(Bu,2);

% Total system transfer behaviour
P = ss(A,[Bw Bu],[Ce;Cy],[Dew Deu; zeros(n,nw+nu)]);


%% Control synthesis

% Convex optimization:
cvx_begin sdp
    variable Q(n,n) symmetric;
    variable F(nu,n);
    variable eta;
    
    minimize(eta);
    Q > 0;
    [Q*A'+F'*Bu'+A*Q+Bu*F Bw Q*Ce'+F'*Deu';
        Bw' -eye(nu,nu) Dew';
        Ce*Q+Deu*F Dew -eta*eye(ne,ne)] < 0;
cvx_end

% Controller & closed loop:
K = F/Q;

Acl= A + Bu*K;
disp(eig(Acl));