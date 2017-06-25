% MIXED H2/Hinf CONTROL PREPARATION SCRIPT

clc
clear
% close all

addpath('functions')
addpath('aux_scripts')

cc = 0; % figure counter


%% System definition
% Select system:
% sys_path = 'systems/simple2ndOrder';
% sys_path = 'systems/dummySys';
% sys_path = 'systems/dummySysMORC';
% sys_path = 'systems/simple2ndOrderMORC';
sys_path = 'systems/flexible_robot_arm';
addpath(sys_path)

% Get representation: (return P)
[P,channels] = system_def();

% Filter selection: (return Wof)
Wof = weight_filter();
Gof = tf(Wof);

% cc=cc+1; figure(cc)
% bode(Gof(1,1)); grid; legend('We')


%% Multi-objective robust control
% Specify optimization parameters:
[def_channels,options,optin] = initMORC();
options.integ = 1;
optin.t = 1;
optin.alpha = 50^2;
optin.gamma = 4;

% Controller synthesis:
[Ksys,optout] = multiobjRC(P,channels,Wof,options,optin);
disp(''); disp('Obtained controller:'); Ksys


%% Simulation and evaluation
% Simulation:
run('run_sim')


%% Clean up
% rmpath(sys_path)
