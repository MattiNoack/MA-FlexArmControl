% TEST SIMULATION FOR VALIDATION OF ROBUST CONTROL
% using simulink

% Compute comparison controller:
[K1,K2] = cmpCtrl(P,channels,W);
Kcomp = ss(K1,'minimal');

% Inputs & perturbations:
ref = 1;
Tref = 0;
dist = 0.5e0;
Tdist = 25;

% Simulation:
Tsim = 80; % simulation time
Ts = 1e-2; % sample time
x0 = [0;zeros(size(P.A,1)-1,1)];
xKcomp0 = zeros(size(Kcomp.A,1),1);%zeros(size(P.A,1)+size(Wof.A,1),1);

% Output selection:
ny = length(channels.contr.out);
C_scope = [eye(2) zeros(2,1+ny)]; % (Dtheta,q)
C_meas = -[zeros(ny,3) eye(ny)]; % (y1,y2)

sim('sim_comparison')

% Visualization:
cc = cc + 1; figure(cc);
subplot(3,1,1); hold on
plot(simulation_test.time,simulation_test.signals(1).values, 'LineWidth', 2)
% plot(simulation_test.time,simulation_test.signals(3).values(:,1),'--r','LineWidth', 2)
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

% cc = cc + 1; figure(cc); surfl(Tpl,Xpl,wpl)
% xlabel('t'); ylabel('x'); zlabel('w(x,t)'); title('Ideal bending deformation')

