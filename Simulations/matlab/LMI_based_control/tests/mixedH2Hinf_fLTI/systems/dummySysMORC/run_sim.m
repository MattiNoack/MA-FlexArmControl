% TEST SIMULATION FOR VALIDATION OF ROBUST CONTROL
% using simulink

% Inputs & perturbations:
ref = 1;
Tref = 2;
dist = 0.5e0;
Tdist = 10;

% Simulation:
Tsim = 30; % simulation time
Ts = 1e-2; % sample time

x0 = [0;zeros(size(P.A,1)-1,1)];
xK0 = zeros(size(P.A,1),1);%zeros(size(P.A,1)+size(Wof.A,1),1);

sim('sim_test')

% Visualization:
cc = cc + 1;
figure(cc); hold on
plot(simulation_test.time,simulation_test.signals(1).values, 'LineWidth', 2)
title('Test simulation results')
ylabel('z'); xlabel('t'); grid
cc = cc + 1;
figure(cc); hold on
plot(simulation_test.time,control_test.signals(1).values, 'LineWidth', 2)
title('Test simulation results')
ylabel('u'); xlabel('t'); grid
