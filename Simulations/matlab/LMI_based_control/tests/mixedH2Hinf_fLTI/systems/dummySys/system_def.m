% DUMMY SYSTEM OF 1st ORDER
% system definition

% System form:  dy + a*y = b*u + d
% Inputs:       (r,d,u)
% Outputs:      (e,u,y)
% Disturbance:  d

function [P,indim,outdim] = system_def()

    % System parameters:
    a = 1; % stable
    b = 1;

    % State-space representation:
    A = -a;
    Bw = [1 0];
    Bu = b;
    Cz = [1;0]; % error+input
    Dzw = [0 -1;0 0];
    Dzu = [0;1];
    Cy = 1; % state FB
    Dyw = [0 -1];    

    n = size(A,1);
    nz = size(Cz,1);
    ny = size(Cy,1);
    nw = size(Bw,2);    
    nu = size(Bu,2);

    % Dimension vectors:
    indim = [nw nu];
    outdim = [nz ny];

    % Total system:
    P = ss(A,[Bw Bu],[Cz;Cy],[Dzw Dzu; Dyw zeros(ny,nu)]);

end