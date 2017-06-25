% CVX EXAMPLES
% from [lmi-cvx.pdf]

clear all

addpath('c:\program files\matlab\cvx')
addpath('c:\program files\matlab\cvx\structures')
addpath('c:\program files\matlab\cvx\lib')
addpath('c:\program files\matlab\cvx\functions')
addpath('c:\program files\matlab\cvx\commands')
addpath('c:\program files\matlab\cvx\builtins')


%% Standard Lyapunov equation 
% no strict inequalities in cvx!
% Problem: A'P+PA<0, P>0

A = [0 1;-1 -2];
n = size(A,1);

% Reformulation:
cvx_begin sdp
    variable P(n,n) symmetric
    A'*P + P*A <= -eye(n)
    P >= eye(n)
cvx_end

disp(P)

% -> newer version: strict inequalities seem possible


%% Adding linear objective function

% Miniize trace:
cvx_begin sdp
    variable P(n,n) symmetric
    minimize(trace(P))
    A'*P + P*A <= -eye(n)
    P >= eye(n)
cvx_end

disp(P)


%% Diagonal solution

cvx_begin sdp
    variable D(n,n) diagonal
    A'*D + D*A <= -eye(n)
    D >= eye(n)
cvx_end

disp(D)


%% Multiple Lyapunov equations

A1 = [0 1;-1 -3];
A2 = [0 1;-1 -2];
n = size(A1,1);

cvx_begin sdp
    variable P(n,n) symmetric
    A1'*P + P*A1 <= -eye(n)
    A2'*P + P*A2 <= -eye(n)
    P >= eye(n)
cvx_end

disp(P)


%% Bounded real lemma

B = [0;1];
C = [1 0];
m = size(B,2);

gamma = 2;

cvx_begin sdp
    variable P(n,n) symmetric
    P >= 0
    [ A'*P + P*A + C'*C P*B; ...
    B'*P -gamma^2*eye(m)] <= 0
cvx_end

disp(P)


%% Hinf gain minimization

cvx_begin sdp
    variable rho
    variable P(n,n) symmetric
    minimize (rho)
    P >= 0
    [ A'*P + P*A + C'*C P*B; ...
    B'*P -rho*eye(m)] <= 0
cvx_end

gamma = sqrt(rho);
