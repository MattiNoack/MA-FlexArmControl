% SIMPLE 2nd ORDER SYSTEM
% system definition

% System form:  qdd + w*q = k*r
% Inputs:   (r,d,u)
% Outputs:  (e,u,y)
% H2:       r(,n)
% Hinf:     d

function [P,indim,outdim] = system_def()

    % System parameters:
    omega = -1; % unstable
    k = 1;

    % State-space representation:
    A = [0 1;-omega 0];
    B2 = [[];[]];%[0;0];
    Binf = [0 0;0 1];%[0;1];
    Bu = [0;k];
    Cz = [-1 0;0 0]; % error
    Dz2 = [[];[]];%[1;0];
    Dzinf = [1 0;0 0];%[0;0]
    Dzu = [0;1];
    Cy = eye(2); % state FB
    Dy2 = [[];[]];%[-1;0]; % no noise case
    Dyinf = [-1 0;0 0];%[0;0];

    n = size(A,1);
    nz = size(Cz,1);
    ny = size(Cy,1);
    n2 = size(B2,2);
    ninf = size(Binf,2);
    nu = size(Bu,2);

    % Dimension vectors:
    indim = [n2 ninf nu];
    outdim = [nz ny];

    % Total system:
    P = ss(A,[B2 Binf Bu],[Cz;Cy],[Dz2 Dzinf Dzu; Dy2 Dyinf zeros(ny,nu)]);

end