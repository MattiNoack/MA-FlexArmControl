% EVALUATION OF UNCONSTRAINED MODE FLEXIBLE ARM SIMULATION

close all


%% Deformation plot
% Select rougher grid:
nappr = 50; % grid approximation

lx = length(X);
indx = 1:round((lx-1)/nappr):lx;
Xpl = X(indx);

lt = length(time);
indt = 1:round((lt-1)/nappr):lt;
Tpl = time(indt)';

% 3D plot:
wpl = deform_data.signals.values(indt,indx);
figure; surfl(Xpl,Tpl,wpl)
xlabel('x'); ylabel('t'); zlabel('w(x,t)'); title('Bending deformation')


%% Rotation angle
% Get data:
theta = theta_data.signals.values;
eta = eta_data.signals.values;

% Comparison plot:
figure; hold on
plot(time,theta,'LineWidth',1.5); plot(time,eta,'LineWidth',1.5)
grid; xlabel('time t'); ylabel('\theta(t)'); title('Rotation angle')
legend('including deformation','undisturbed')


%% Approximation error
% Get data:
delta_theta = error_data.signals.values;

% Plot:
figure; plot(time,delta_theta,'LineWidth',1.5)
grid; xlabel('time t'); ylabel('\Delta\theta(t)');
titlestr = strcat('Approximation error for order (N,Na)=(',num2str(N),',',num2str(Na),')');
title(titlestr)
