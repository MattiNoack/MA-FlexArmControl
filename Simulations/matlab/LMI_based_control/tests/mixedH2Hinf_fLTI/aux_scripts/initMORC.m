% INITIALIZE STRUCTURES FOR MULTI-OBJECTIVE ROBUST CONTROL
% setting of required fields to default values

% Output:
%       channels    structure for channel pair definitions (Remark: (z,w) first)
%                       - channels.contr ... controller [y:u]
%                       - channels.H2 ... generilized H2 [z2;w2]
%                       - channels.Hinf ... Hinf norm [zinf;winf]
%                       - channels.passive ... passivity [zp;wp]
%                       - channels.delta ... uncertainties [zd;wd]
%                       - channels.nomreg ... nominal regulation [znr;wnr]
%       options     structure of case options
%                       - options.optHinf ... optimal or suboptimal Hinf approach
%                       - options.integ ... integrator in control included
%       optin       structure for predefined optimization variables
%                       - optin.gamma ... Hinf suboptimal value
%                       - optin.alpha ... generalized H2 peak bound
%                       - optin.t ... well conditioning of lyapunov variables


function [channels,options,optin] = initMORC()

    % Default channel definitions:
    channels.contr = zeros(2,0);
    channels.H2 = zeros(2,0);
    channels.Hinf = zeros(2,0);
    channels.passive = zeros(2,0);
    channels.delta = zeros(2,0);
    channels.nomreg = zeros(2,0);
    
    % General options:
    options.optHinf = 0; % suboptimal
    options.integ = 0; % no integral action
    
    % Optimization predefinitions:
    optin.gamma = 10;
    optin.alpha = 1e3;
    optin.t = 1;
    optin.S = 0; % single input constant signal

end
