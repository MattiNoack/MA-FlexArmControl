% TEST SIMULATION FOR VALIDATION OF ROBUST CONTROL
% using simulink

% Inputs & perturbations:
ref = 1;
Tref = 2;
dist = 0.5e0;
Tdist = 8;

% Simulation:
Tsim = 15; % simulation time
Ts = 1e-2; % sample time

x0 = [0;zeros(size(P.A,1)-1,1)];
xK0 = zeros(size(P.A,1),1);%zeros(size(P.A,1)+size(Wof.A,1),1);

sim('sim_test')

% Visualization:
cc = cc + 1;
figure(cc)
% plot(simulation_test.time,simulation_test.signals(1).values, 'LineWidth', 2)
% title('Test simulation results')
% ylabel('z'); xlabel('t'); grid
% cc = cc + 1;
% figure(cc)
% plot(simulation_test.time,control_test.signals(1).values, 'LineWidth', 2)
% title('Test simulation results')
% ylabel('u'); xlabel('t'); grid
subplot(3,1,1)
hold on
plot(simulation_test.time,simulation_test.signals(1).values)
plot(simulation_test.time,simulation_test.signals(2).values)
legend(simulation_test.signals(1:2).label)
ylabel('position'); grid
title('Test simulation results - state feedback')
hold off
subplot(3,1,2)
hold on
plot(simulation_test.time,simulation_test.signals(3).values)
legend(simulation_test.signals(3).label)
ylabel('velocity'); grid
hold off
subplot(3,1,3)
hold on
plot(simulation_test.time,control_test.signals.values)
legend(control_test.signals.label)
ylabel('control'); xlabel('time t'); grid
hold off
