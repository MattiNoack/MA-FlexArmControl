% SIMPLE 2nd ORDER SYSTEM
% connected to the finite dimensional flexible arm approximated model

clc
clear all
close all

addpath('simulations')
addpath('functions')
addpath('aux_scripts')
cvx_init % inlcude cvx library

cc = 0; % figure counter

% System form:  qdd + w*q = k*r
% Inputs:   (r,d,u)
% Outputs:  (e,u,y)
% H2:       r(,n)
% Hinf:     d


%% System definition - simplified case
% System parameters:
omega = -1; % unstable
k = 1;
x0 = [1;0];

% State-space representation:
A = [0 1;-omega 0];
B2 = [0;0];
Binf = [0;1];
Bu = [0;k];
Cz = [-1 0;0 0]; % error
Dz2 = [1;0];
Dzu = [0;1];
Cy = eye(2); % state FB
Dy2 = [-1;0]; % no noise case

n = size(A,1);
nz = size(Cz,1);
ny = size(Cy,1);
n2 = size(B2,2);
ninf = size(Binf,2);
nu = size(Bu,2);

% Total system:
P1 = ss(A,[B2 Binf Bu],[Cz;Cy],[Dz2 zeros(nz,ninf) Dzu; Dy2 zeros(ny,ninf+nu)]);
Ptf1 = tf(P1);

% Filter design:
apol=poly([-1e-3 -1e-2]); % char pol for output
czer=poly([]);
Au = 10; % weighting factor for control
weighting_filter % create filter functions

% Frequency domain test:
G_of = tf(sys_of);
G_cf = tf(sys_cf);
G = Ptf1(3,3); % r->q
cc = cc + 1; figure(cc)
hold on
bode(G_of)
bode(G_cf)
bode(G)
legend('Gof','Gcf','G'); grid
hold off


%% Nominal controller synthesis - state feedback
% Paramters:
gamma = 0.4; % Hinf suboptimum

% SDP optimization:
cvx_begin sdp
    variable X1(n,n) symmetric
    variable X2(nf,nf) symmetric
    variable Z(nu,n)
    variable W(nz,nz) symmetric
%     variable gamma
    minimize (trace(W))%+gamma)
%     gamma > 0
    X1 > 0
    X2 > 0
    [A*X1+X1*A'+Bu*Z+Z'*Bu' X1*Cz'*Bf'+Z'*Dzu'*Bf' zeros(n,n2); ...
        Bf*Cz*X1+Bf*Dzu*Z Af*X2+X2*Af' Bf*Dz2; ...
        zeros(n2,n) Dz2'*Bf' -eye(n2)] < 0
    [W zeros(nz,n) Cf*X2; ...
        zeros(n,nz) X1 zeros(n,nf)
        X2*Cf' zeros(nf,n) X2] > 0
    [A*X1+X1*A'+Bu*Z+Z'*Bu' X1*Cz'*Bf'+Z'*Dzu'*Bf' Binf zeros(n,nz); ...
        Bf*Cz*X1+Bf*Dzu*Z Af*X2+X2*Af' zeros(nf,ninf) X2*Cf'; ...
        Binf' zeros(ninf,nf) -gamma*eye(ninf) zeros(ninf,nz); ...
        zeros(nz,n) Cf*X2 zeros(nz,ninf) -gamma*eye(nz)] < 0
cvx_end

% Resulting control law & closed-loop:
K = Z/X1;
disp(' '); disp('Resulting gain:')
disp(K)
disp('Closed-loop eigenvalues:')
disp(eig(A+Bu*K))

% Simulation test:
sim_illustrate_stateFB % run simulink model & generate plots


%% Include disturbances
% state feedback required

% % Include noise: w2=(qref,n)
% Dy2 = [-1 1;0 2];
% Dz2 = [1 0;0 0];
% B2 = zeros(2,2);
% n2 = size(B2,2);

% % Reduce output: ny = 1
% Cy = [1 0]; % state FB
% Dy2 = [-1 1]; % no noise case
% ny = size(Cy,1);

P2 = ss(A,[B2 Binf Bu],[Cz;Cy],[Dz2 zeros(nz,ninf) Dzu; Dy2 zeros(ny,ninf+nu)]);
Ptf2 = tf(P2);

% Optimization with auxiliary variables:
% gamma =10;

cvx_begin sdp
    variable W(nz,nz) symmetric
    variable P(n,n) symmetric
    variable Q(n,n) symmetric
    variable Pof(nf,nf) symmetric
    variable Ah(n,n)
    variable Bh(n,ny)
    variable Ch(nu,n)
    variable gamma
    
    minimize (trace(W)+gamma)
    
    gamma > 0
    [Q eye(n) zeros(n,nf); ...
        eye(n) P zeros(n,nf); ...
        zeros(nf,2*n) Pof] > 0
    
    [W zeros(nz,2*n) Cf*Pof; ...
        zeros(n,nz) Q eye(n) zeros(n,nf); ...
        zeros(n,nz) eye(n) P zeros(n,nf); ...
        Pof*Cf' zeros(nf,2*n) Pof] > 0
    [Q*A+Bh*Cy+A'*Q+Cy'*Bh' Ah+A' Cz'*Bf' Bh*Dy2; ...
        A+Ah' A*P+Bu*Ch+P*A'+Ch'*Bu' P*Cz'*Bf'+Ch'*Dzu'*Bf' zeros(n,n2); ...
        Bf*Cz Bf*Cz*P+Bf*Dzu*Ch Af*Pof Bf*Dz2; ...
        Dy2'*Bh' zeros(n2,n) Dz2'*Bf' -eye(n2)] < 0
    
    [Q*A+Bh*Cy+A'*Q+Cy'*Bh' Ah+A' Cz'*Bf' Q*Binf zeros(n,nz); ...
        A+Ah' A*P+Bu*Ch+P*A'+Ch'*Bu' P*Cz'*Bf'+Ch'*Dzu'*Bf' Binf zeros(n,nz); ...
        Bf*Cz Bf*Cz*P+Bf*Dzu*Ch Af*Pof zeros(nf,ninf) Pof'*Cf'; ...
        Binf'*Q' Binf' zeros(ninf,nf) -gamma*eye(ninf) zeros(ninf,nz); ...
        zeros(nz,2*n) Cf*Pof zeros(nz,ninf) -gamma*eye(nz)] < 0
   
cvx_end

% Calculate controller parameterization:
%TODO: generate function
nk = n; % control order

M = sym('M%d%d',[n,nk],'real');
N = sym('N%d%d',[n,nk],'real');
eqn = M*N' == eye(n) - P*Q;

sol=solve(eqn,[M,N]);
M1 = double([sol.M11(1) sol.M12(1);sol.M21(1) sol.M22(1)]);
M2 = double([sol.M11(2) sol.M12(2);sol.M21(2) sol.M22(2)]);
M3 = double([sol.M11(3) sol.M12(3);sol.M21(3) sol.M22(3)]);
N1 = double([sol.N11(1) sol.N12(1);sol.N21(1) sol.N22(1)]);
N2 = double([sol.N11(2) sol.N12(2);sol.N21(2) sol.N22(2)]);
N3 = double([sol.N11(3) sol.N12(3);sol.N21(3) sol.N22(3)]);
% Problem: how to get solution pairs?!

Msol = eye(n)-P*Q; % selection in case of n=nk possible
Nsol = eye(n);

Ak = sym('Ak%d%d',[nk,nk],'real');
Bk = sym('Bk%d%d',[nk,ny],'real');
Ck = sym('Ck%d%d',[nu,nk],'real');
eqn1 = Ah == Q*A*P+Msol*Bk*Cy*P + Q*Bu*Ck*Nsol' + Msol*Ak*Nsol';
eqn2 = Bh == Msol*Bk;
eqn3 = Ch == Ck*Nsol';

sol_total=solve([eqn1,eqn2,eqn3],[Ak,Bk,Ck']);
Ak1 = double([sol_total.Ak11(1) sol_total.Ak12(1);sol_total.Ak21(1) sol_total.Ak22(1)]);
Bk1 = double([sol_total.Bk11(1) sol_total.Bk12(1);sol_total.Bk21(1) sol_total.Bk22(1)]);
Ck1 = double([sol_total.Ck11(1) sol_total.Ck12(1)]);

disp(''); disp('Obtained output feedback controller:')
Ksys = ss(Ak1,Bk1,Ck1,zeros(nu,ny))

% Simulation test:
xK0 = zeros(nk,1);
% sim_illustrate_outputFB % run simulink model & generate plots


%% Hinf approach without filters
% System redefinition: (without noise)
A = A;
Bu = Bu;
Bw = [Binf [0;0]]; % w = (d,qref)'
Ce = Cz;
Cy = Cy;
Dew = [0 1;0 0];
Deu = [0;1];
Dyw = [0 -1;0 0];

n = size(A,1);
ny = size(Cy,1);
nu = size(Bu,2);
ne = size(Ce,1);
nw = size(Bw,2);

P2 = ss(A,[Bw Bu],[Ce;Cy],[Dew Deu;Dyw zeros(ny,nu)]);
Ptf2 = tf(P2);

% Optimization: [ATIC lecture 8]
eta = (3)^2; % eta = mu^2

cvx_begin sdp
    variable X(n,n) symmetric;
    variable Y(n,n) symmetric;
    variable Ah(n,n);
    variable Bh(n,ny);
    variable Ch(nu,n);
%     variable eta;
    
%     minimize eta;
    
    [X, eye(n,n);
    eye(n,n), Y] > 0;

    [A*X + Bu*Ch + X*A' + Ch'*Bu', A+Ah', Bw, X*Ce' + Ch'*Deu';
    A'+Ah, Y*A + A'*Y + Bh*Cy + Cy'*Bh', Y*Bw + Bh*Dyw, Ce';
    Bw', Bw'*Y + Dyw'*Bh', -eye(nw,nw), Dew';
    Ce*X + Deu*Ch, Ce, Dew, -eta*eye(ne,ne)] < 0;
cvx_end

% Calculate controller parameterization:
nk = n; % control order

M = sym('M%d%d',[n,nk],'real');
N = sym('N%d%d',[n,nk],'real');
eqn = N*M' == eye(n) - Y*X;

sol=solve(eqn,[M,N]);
M1 = double([sol.M11(1) sol.M12(1);sol.M21(1) sol.M22(1)]);
M2 = double([sol.M11(2) sol.M12(2);sol.M21(2) sol.M22(2)]);
M3 = double([sol.M11(3) sol.M12(3);sol.M21(3) sol.M22(3)]);
N1 = double([sol.N11(1) sol.N12(1);sol.N21(1) sol.N22(1)]);
N2 = double([sol.N11(2) sol.N12(2);sol.N21(2) sol.N22(2)]);
N3 = double([sol.N11(3) sol.N12(3);sol.N21(3) sol.N22(3)]);
% Problem: how to get solution pairs?!

Nsol = eye(n)-Y*X; % selection in case of n=nk possible
Msol = eye(n);

Ak = sym('Ak%d%d',[nk,nk],'real');
Bk = sym('Bk%d%d',[nk,ny],'real');
Ck = sym('Ck%d%d',[nu,nk],'real');
eqn1 = Ah == Y*A*X+Nsol*Bk*Cy*X + Y*Bu*Ck*Msol' + Nsol*Ak*Msol';
eqn2 = Bh == Nsol*Bk;
eqn3 = Ch == Ck*Msol';

sol_total=solve([eqn1,eqn2,eqn3],[Ak,Bk,Ck']);
Ak1 = double([sol_total.Ak11(1) sol_total.Ak12(1);sol_total.Ak21(1) sol_total.Ak22(1)]);
Bk1 = double([sol_total.Bk11(1) sol_total.Bk12(1);sol_total.Bk21(1) sol_total.Bk22(1)]);
Ck1 = double([sol_total.Ck11(1) sol_total.Ck12(1)]);

Bk1 = Nsol\Bh; % calculation in case of n=nk possible
Ck1 = Ch/(Msol');
Ak1 = Nsol\(Ah-Y*A*X-Nsol*Bk1*Cy*X-Y*Bu*Ck1*Msol')/(Msol');

disp(''); disp('Obtained output feedback controller:')
Ksys = ss(Ak1,Bk1,Ck1,zeros(nu,ny))

% Simulation test:
xK0 = zeros(nk,1);
sim_illustrate_outputFB % run simulink model & generate plots


%% Consider uncertainty


%TODO: write symbolic function for creating SDPs automatically by giving structure