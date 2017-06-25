% EXECUTION OF MULTI-OBJECTIVE ROBUST CONTROL PROCEDURE

clc
clear
% close all

addpath('functions')
cc = 0; % figure counter


%% System definition
% Select system:
sys_path = 'systems/flexible_arm';
% sys_path = 'systems/flexible_arm_damped';
% sys_path = 'systems/flexible_arm_q0';
addpath(sys_path)

% Get representation & filter:
[P,channels] = system_def();
W = weight_filter();


%% Multi-objective robust control
% Specify optimization parameters:
[def_channels,options,optin] = initMORC();
%optin.S = 0; 
optin.t = 1;
optin.alpha = (11)^2;%(3)^2
optin.gamma = 4;%20

% Test step calculations:
% run('tuning')

% Controller synthesis:
[Ksys,optout,CLsys,OLctrl] = synthMORC(P,channels,W,options,optin);
disp(''); disp('Obtained controller:'); Ksys

% Simulation:
run('run_sim')


%% Superordinated optimization
% % Strain position choice:
% load('eigfun.mat') % load spatial parameters
% samp = 2; % sampling distance
% Jopt = inf; % optimal cost value
% iopt = 1; % optimal position index
% normarr = [];
% tmpopt = options; tmpopt.optGH2 = 1;
% 
% for i=1:samp:length(X)-2
%     [Piter,chan_iter] = system_def(i);
%     [Kiter,optparam] = synthMORC(Piter,chan_iter,W,tmpopt,optin);
%     
%     J = optparam.alpha;
%     if(J<Jopt)
%         Jopt = J;
%         iopt = i;
%     end
%     normarr = [normarr J];
% end
% disp(''); disp('Optimal sensor position:'); disp(X(iopt))
% cc=cc+1; figure(cc)
% plot(X(1:samp:length(X)-2),normarr,'--s'); grid
% xlabel('x'); ylabel('J'); title('Sensor location optimization')


%% Performance analysis
% Performance:
%   - sensitivity functions (bode plot)
CLplot = CLsys([1 3],:);
CLplot.u{1}='r'; CLplot.u{2}='d';
CLplot.y{1}='e'; CLplot.y{2}='u';
cc=cc+1; figure(cc)
bodemag(CLplot); grid; title('Closed-loop sensitivities')

% Fragile controller test:
%   - analyse phase & gain margin of ol system
%   - characteristic values (rise time, etc.)
[Gm,Pm,Wcg,Wcp]=margin(OLctrl);
disp(''); disp('Open loop margins:');
disp(strcat('Pm = ',num2str(Pm),', Gm = ',num2str(Gm)))


%% Comparison controllers
% Controller synthesis:
% [K1,K2,infoPID,gammH2] = cmpCtrl(P,channels,W);

% Simulation test:
% run('run_comparison')


%% Clean up
% rmpath(sys_path)
