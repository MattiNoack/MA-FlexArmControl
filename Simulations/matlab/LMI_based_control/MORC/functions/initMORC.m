% INITIALIZE STRUCTURES FOR MULTI-OBJECTIVE ROBUST CONTROL
% setting of required fields to default values

% Output:
%       channels    structure for channel pair definitions with in/out distinction
%                       - channels.contr ... controller (in:u,out:y)
%                       - channels.integ ... integrator (meas:yi)
%                       - channels.H2 ... generilized H2 (in:w2,out:z2)
%                       - channels.Hinf ... Hinf norm (in:winf,out:zinf)
%                       - channels.passive ... passivity (in:wp,out:zp)
%                       - channels.delta ... uncertainties (in:wd,out:zd)
%                       - channels.nomreg ... nominal regulation (in:wnr,out:znr)
%       options     structure of case options
%                       - options.optHinf ... optimal or suboptimal Hinf approach
%                       - options.optGH2 ... optimal generalized H2 norm
%                       - options.optLCond ... optimal Lyapunov conditioning
%                       - options.normCond ... optimal norm conditioning
%       optin       structure for predefined optimization variables
%                       - optin.gamma ... Hinf suboptimal value
%                       - optin.alpha ... generalized H2 peak bound
%                       - optin.t ... well conditioning of lyapunov variables
%                       - optin.S ... signal generator for nominal regulation
%                       - optin.b ... conditioning of matrix variables norm


function [channels,options,optin] = initMORC()

    % Default channel definitions:
    channels.contr.in = zeros(1,0); channels.contr.out = zeros(1,0);
    channels.integ.meas = zeros(1,0);
    channels.H2.in = zeros(1,0); channels.H2.out = zeros(1,0);
    channels.Hinf.in = zeros(1,0); channels.Hinf.out = zeros(1,0);
    channels.passive.in = zeros(1,0); channels.passive.out = zeros(1,0);
    channels.delta.in = zeros(1,0); channels.delta.out = zeros(1,0);
    channels.nomreg.in = zeros(1,0); channels.nomreg.out = zeros(1,0);
    
    % General options:
    options.optHinf = 0; % suboptimal
    options.optGH2 = 0; % no optimization
    options.optLCond = 0; % no optimization
    options.normCond = 0; % no matrix variable norm conditioning
    
    % Optimization predefinitions:
    optin.gamma = 10;
    optin.alpha = 1e3;
    optin.t = 1;
    optin.S = 0; % single input constant signal
    optin.b = 0;

end
