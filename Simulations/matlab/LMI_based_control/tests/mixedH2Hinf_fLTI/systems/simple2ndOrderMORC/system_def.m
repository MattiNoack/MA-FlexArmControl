% SIMPLE 2nd ORDER SYSTEM
% system definition

% System form:  xdd + w*x = k*u
% Inputs:          (r,d,u)
% Outputs:         (e,u,y)
% Reference:       r
% Disturbance:     d

function [P,channels] = system_def()

    % System parameters:
    omega = -1; % stable
    k = 1;

    % State-space representation:
    A = [0 1;-omega 0];
    Bw = [0 0;0 1];
    Bu = [0;k];
    Cz = [1 0;0 0]; % error
    Dzw = [-1 0;0 0];
    Dzu = [0;1];
    Cy = [1 0];
    Dyw = [-1 0];

    % Channel structure:
    channels = initMORC(); % init
    channels.contr = [3;3];
    channels.H2 = [1 2;
                   2 1];
%     channels.Hinf = [1;
%                      2];               
%     channels.passive = [1;2];

    % Total system:
    P = ss(A,[Bw Bu],[Cz;Cy],[Dzw Dzu; Dyw zeros(1,1)]);

end