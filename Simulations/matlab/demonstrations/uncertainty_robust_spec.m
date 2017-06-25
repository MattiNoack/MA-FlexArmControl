% UNCERTAINTY MOTIVIATION & ROBUST PERFORMANCE SPECIFCATION
% by weight bound on sensitivity function

clear all
close all
clc


%% Additive uncertainty motivation
% System definitions:
G = tf(8,conv([10 1],conv([0.05 1],[0.05 1]))); % system model
K = 1.5*tf([0.1 1],conv([0.65 1],[0.03 1])); % controller

% W = tf([],[]); % weight filter in characteristic bandwidth
% del = ultidyn('Delta',[1 1],'Bound',1); % uncertainty with gain<1
% delS = usample(del,1); % select one delta
load('generated_real_plant.mat') 
delComp = usample(del,25);

% H = G + delS; % real system
H = Htf;

% Uncertain magnitudes:
% preference input: ctrlpref
figure(1)
hold on
bodemag(H,'r')
bodemag(G,'b')
title('System magnitudes')
grid; legend('Plant','Model')
hold off

figure(2)
bodemag(del,'b')
title('Uncertainty magnitude'); grid
figure(3); hold on
bodemag(G+delComp,'r')
bodemag(G,'b')
title('Uncertain system magnitude')
grid; legend('Possible plants','Model')
hold off

% Controller test:
check = minreal(K/(1+K*G));

figure(4)
nyquist(G*K)
figure(5); hold on
bodemag(check)
bodemag(tf(1,1))
title('Sensitivity magnitude'); grid


%% Multiplicative uncertainty
% Contour example: [http://stackoverflow.com/questions/11345838/how-to-plot-inequalities]
% t = (0:999) / 100;
% [d1, d2] = meshgrid((0:999) / 1000, (0:999) / 1000);
% d = rand(2, 1);                            %# In this example p = [0.1, 0.2]
% ineq1 = d2 < d(2) * (1 - d(1));             %# First inequation
% ineq2 = d1 < d(1) * (1 - (d2 / (1 - d(1)))); %# Second inequation
% both = ineq1 & ineq2;                      %# Intersection of both inequations
% 
% figure, hold on
% c = 1:3;                                   %# Contour levels
% % contourf(c(1) * ineq1, [c(1), c(1)], 'b')  %# Fill area for first inequation
% % contourf(c(2) * ineq2, [c(2), c(2)], 'g')  %# Fill area for second inequation
% contourf(c(3) * both, [c(3), c(3)], 'r')   %# Fill area for both inequations
% legend('First', 'Second', 'Both')
% set(gca, ...                               %# Fixing axes ticks
%     'XTickLabel', {t(get(gca, 'XTick'))}, 'YTickLabel', {t(get(gca, 'YTick'))})
% 
% colors = zeros(size(X))+ineq1+ineq2;
% scatter(X(:),Y(:),3,colors(:),'filled')

% Condition region:
t = (-9999:9999) / 10000; % scaling
[d1, d2] = meshgrid((-9999:9999) / 10000, (-9999:9999) / 10000);
c = 1:2; % contour index

figure(5); hold on
a = 10;
ineq1 = 2 + d1 + d2 <= 0;
ineq2 = (1+d1+d2) + (a^2+1).*d1.*d2 <= 0;
both = ineq1 | ineq2;
contourf(both, [c(1), c(1)], 'b') 

a = 5;
ineq1 = 2 + d1 + d2 <= 0;
ineq2 = (1+d1+d2) + (a^2+1).*d1.*d2 <= 0;
both = ineq1 | ineq2;
contourf(both, [c(1), c(1)], 'r--')
a = 1;
ineq1 = 2 + d1 + d2 <= 0;
ineq2 = (1+d1+d2) + (a^2+1).*d1.*d2 <= 0;
both = ineq1 | ineq2;
contourf(both, [c(1), c(1)], 'm-.')

plot([900 1100 1100 900 900],[1100 1100 900 900 1100],'black')
set(gca,'XTickLabel', {t(get(gca, 'XTick'))}, 'YTickLabel', {t(get(gca, 'YTick'))})
grid; xlabel('\delta_1'); ylabel('\delta_2'); legend('a=10','a=5','a=1')


