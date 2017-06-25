% INITIALIZATION OF UNCONSTRAINED MODE FLEXIBLE ARM
% solve BVP for eigenfunctions and construct LTI system for damped system

clear
close all
addpath functions


%% Parameters
% Physical system:
E = 1.96e11;     % [N/m^2]
I = 2.08e-12;    % [m^4]
Al = 25e-6;      % [m^2]
L = 0.25;        % [m]
rho = 10.667e3;  % [kg/m^3]
m = 0.2;         % [kg]
J0 = 1e-3;       % [kg*m^2]
Jp = 1e-3;       % [kg*m^2]
J = J0+Al*rho*L^3/3+Jp+m*L^2;
Th = 1e-3;       % [m]

% E = 0.1; I = 1;
% Al = 1; L = 5;
% rho = 1; m = 1;
% J0 = 1; Jp = 1;
% J = 1; %J0+Al*rho*L^3/3+Jp+m*L^2;
% Th = 1e-2;

% Extra damping:
D = 5e-1; %1e-4;
zeta = 0.5*D;

% Simulation:
Tf = 15;
Ts = 1e-4;
time = 0:Ts:Tf;     % time scale

N = 9;              % approbetamation order
h = 5e-3;           % spatial step size
X = 0:h:L;          % spatial domain
epsilon = 1e-2;     % numeric accuracy


%% Eigenfunctions
% Characteristic matrix:
syms beta real
M = charMat([1 0 0 0;0 0 1 0;0 beta^4/(Al*rho) 1/(J0+Jp) 0;m/(Al*rho)*beta^4 0 0 1],[0;L;0;L],beta);

% Transcendental equation roots:
char = matlabFunction(det(M));
bsmp = 1:1:120; % samples (200)
bsol = 0; % storage
for i=1:length(bsmp)
    tmp = fzero(char,bsmp(i));
    if (tmp>bsol(end)+epsilon)
        bsol = [bsol tmp];      
    end
end

% Auxiliary eigenfunction construction:
Mf = matlabFunction(M); % A(beta)
Phiarr = zeros(N,length(X));
cnr = zeros(1,N);
for i=1:N
    be = bsol(i+1); % exclude trivial solution beta=0
    P = null(Mf(be));
    if (size(P,2)==1)
        tmpfcn = @(x)(P(1)*sin(be*x)+P(2)*cos(be*x)+P(3)*sinh(be*x)+P(4)*cosh(be*x));
        tmpnorm = sqrt(integral(@(x)(tmpfcn(x).^2),0,L)); % auxiliary normilization
        Phiarr(i,:) = tmpfcn(X)/tmpnorm;
    else        
        disp(strcat('Ill conditioned characteristic matrix (rank<3) at beta=',...
        num2str(be),', i=',num2str(i+1),'.'))
    end
end
figure; hold on
for i=1:5
    plot(X,Phiarr(i,:))
end
legend('i=1','i=2','i=3','i=4','i=5')
grid; xlabel('x'); ylabel('\Phi_i(x)'); title('Auxiliary eigenmodes'); hold off

% Original eigenfunction:
omeg = sqrt(E*I/(Al*rho))*bsol(2:N+1)'.^2;
D2Phiarr = (diff(Phiarr',2)/h^2)';
gamma = -E*I/(J0+Jp).*D2Phiarr(:,1)./(omeg.^2);
phiarr = Phiarr - gamma.*X;
figure; hold on
for i=1:5
    plot(X,phiarr(i,:))
end
legend('i=1','i=2','i=3','i=4','i=5')
grid; xlabel('x'); ylabel('\phi_i(x)'); title('Eigenmodes'); hold off


%% Dynamical system
% State space representation:
karr = -(E*I*D2Phiarr(:,1)./(omeg.^2)+m*L*Phiarr(:,end))/(J*Al*rho);
B = [0;1/J;kron(karr,[0;1])];
A = blkdiag([0 1;0 -D/J],kron(-2*diag(omeg.*zeta),[0 0;0 1])+kron(-diag(omeg.^2),[0 0;1 0])...
        + kron(eye(N),[0 1;0 0]));

% Initial state & input:
z0 = zeros(2*(N+1),1);
u = @(t)(0.05*sin(t));
% z0 = [0 0 kron(exp(-[1:N]),[1 0])]'; % excitation by IC
% u = @(t)(0*t);
% u = @(t)(0.05+0*t);

input_data.time = time;
input_data.signals.values = u(time)';
input_data.signals.dimensions = 1;


%% Lower order approximation
% Parameter selection:
Na = 1;
Aa = A(1:2*(Na+1),1:2*(Na+1));
Ba = B(1:2*(Na+1));
z0a = z0(1:2*(Na+1));

% Signal routing:
C1a = [1 zeros(1,2*Na+1)];
if (Na~=0)
    ga = gamma(1:Na);
    phia = phiarr(1:Na,:);
    D2phia = (diff(phia',2)/h^2)';
    C2a = [zeros(Na,2) kron(eye(Na),[1 0])];
else
    ga = 0;
    phia = 0;
    D2phia = 0;
    C2a = zeros(1,2);
end

% Simulation:
% sim('sim_damped')


%% Save data
% ctrl_path = ('C:\Users\user\Dropbox\projekte\Lehrveranstaltungen\MA\Simulations\matlab\LMI_based_control\mixedH2Hinf_fLTI\systems\flexible_robot_arm\');
ctrl_path = ('C:\Users\user\Dropbox\projekte\Lehrveranstaltungen\MA\Simulations\matlab\LMI_based_control\MORC\systems\flexible_arm_damped\');

% Control design:
Adyn = Aa; Bdyn = Ba; gamm = ga;
save(strcat(ctrl_path,'dyn_param.mat'),'Adyn','Bdyn','gamm','Th')
save('dyn_param.mat','Adyn','Bdyn','gamm','Th')

% Eigenfunctions:
phi = phia';
D2phi = D2phia';
save(strcat(ctrl_path,'eigfun.mat'),'X','phi','D2phi')
save('eigfun.mat','X','phi','D2phi')

% Higher order dynamics:
D2phiarr = (diff(phiarr',2)/h^2)';
% save(strcat(ctrl_path,'high_order.mat'),'omeg','karr','gamma','D2phiarr')
save('high_order.mat','omeg','karr','gamma','D2phiarr')

