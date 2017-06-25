% TEST SIMULATION FOR VALIDATION OF ROBUST CONTROL
% using simulink

% Inputs & perturbations:
ref = 0;
Tref = 2;
dist = 1;
Tdist = 8;

% Simulation:
Tsim = 30; % simulation time
Ts = 1e-2; % sample time

x0 = [0;zeros(size(P.A,1)-1,1)];
xK0 = zeros(size(Ksys.A,1),1);

sim('sim_test')

% Visualization:
cc = 0;
cc = cc + 1;
figure(cc); hold on
plot(simulation_test.time,simulation_test.signals(1).values, 'LineWidth', 1.5)
title('Test simulation results')
ylabel('z'); xlabel('t'); grid
% text(20,-0.8,'\frac{1}{2}','FontSize',14)
cc = cc + 1;
figure(cc); hold on
plot(simulation_test.time,control_test.signals(1).values, 'LineWidth', 1.5)
title('Test simulation results')
ylabel('u'); xlabel('t'); grid
