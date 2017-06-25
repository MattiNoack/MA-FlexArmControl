% TEST SIMULATION FOR VALIDATION OF ROBUST CONTROL
% using simulink

% Inputs & perturbations:
ref = 1;
Tref = 2;
dist = 0.5e0;
Tdist = 25;

% Simulation:
Tsim = 80; % simulation time
Ts = 1e-2; % sample time
x0 = [0;zeros(size(P.A,1)-1,1)];
xK0 = zeros(size(Ksys.A,1),1);%zeros(size(P.A,1)+size(Wof.A,1),1);

sim('sim_test')

% Visualization:
cc = cc + 1; figure(cc); 
subplot(3,1,1); hold on
plot(simulation_test.time,simulation_test.signals(1).values, 'LineWidth', 2)
plot(simulation_test.time,simulation_test.signals(3).values,'--r','LineWidth', 2)
% legend('z_1','y')
title('Rotation angle error'); ylabel('\Delta\theta'); grid
subplot(3,1,2); hold on
plot(simulation_test.time,simulation_test.signals(2).values, 'LineWidth', 2)
title('Mode output'); ylabel('q'); grid
subplot(3,1,3); hold on
plot(control_test.time,control_test.signals.values, 'LineWidth', 2)
title('Control action'); ylabel('u'); xlabel('t'); grid

% Deformation plot:
load('eigfun.mat')
time = simulation_test.time;

nappr = 80; % grid approximation
lx = length(X);
indx = 1:round((lx-1)/nappr):lx;
Xpl = X(indx);
lt = length(time);
indt = 1:round((lt-1)/nappr):lt;
Tpl = time(indt)';

q = simulation_test.signals(2).values(indt)';
phi = phi(indx);
wpl = phi*q;

cc = cc + 1; figure(cc); surfl(Tpl,Xpl,wpl)
xlabel('t'); ylabel('x'); zlabel('w(x,t)'); title('Ideal bending deformation')

