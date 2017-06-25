% DUMMY SYSTEM OF 1st ORDER
% system definition

% System form:  dy + a*y = b*u + d
% Inputs:       (r,d,u)
% Outputs:      (e,u,y)
% Disturbance:  d

function [P,channels] = system_def()

    % System parameters:
    a = 1; % stable
    b = 1;

    % State-space representation:
    A = -a;
    Bw = [0 1];
    Bu = b;
    Cz = [-1;0]; % error+input
    Dzw = [1 0;0 0];
    Dzu = [0;1];
    Cy = 1; % state FB
    Dyw = [-1 0];    

    % Channel structure:
    channels = initMORC(); % init
    channels.contr = [3;3];
%     channels.H2 = [1;
%                    2];
%     channels.Hinf = [1 2;
%                      2 1];               
%     channels.passive = [1;2];

    % Total system:
    P = ss(A,[Bw Bu],[Cz;Cy],[Dzw Dzu; Dyw zeros(1,1)]);

end