% TEST SIMULATION FOR VALIDATION OF ROBUST CONTROL
% using simulink

% Inputs & perturbations:
Tref = 0;
ref = 0;
Tdist = 4;
dist = 2;
Tutest = 1;
utest = 0;

% Simulation:
Tsim = 3; % simulation time
Ts = 1e-3; % sample time
sim('sim_test_stateFB')

% Visualization:
cc = cc + 1;
figure(cc)
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