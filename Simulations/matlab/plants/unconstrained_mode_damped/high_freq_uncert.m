% UNCERTAINTY ANALYSIS
% additive uncertainty due to high frequency component on output

function [W1,W2] = high_freq_uncert(zeta)

    % Exception handling:
    if(nargin<1)
        zeta = 0.01;
    end

    % Load parameter:    
    load('dyn_param.mat')
    N = length(gamm);
    load('high_order.mat')
    Ntot = length(omeg);
    
    % Select high order:
    xi = 10; % measurement point index x(i)
    
    g = gamma(N+1:Ntot);
    k = karr(N+1:Ntot);
    w = omeg(N+1:Ntot);
    p = D2phiarr(N+1:Ntot,xi);

    % Build transfer funtions:
    s = tf('s');
    D1=0; D2=0;
    for i=1:(Ntot-N)
        D1 = D1 + g(i)*k(i)/(s^2+2*zeta*w(i)*s+w(i)^2);
        D2 = D2 - Th/2*p(i)*k(i)/(s^2+2*zeta*w(i)*s+w(i)^2);
    end
    
    % Analyze magnitude:
    figure
    subplot(1,2,1); bodemag(D1); grid
    subplot(1,2,2); bodemag(D2); grid
    
    % Construct adequate filters:    
    K1 = (1+0.1)*getPeakGain(D1);
    W1 = K1/(1/w(4)*s+1); % low pass for medium/high damping
    
    K2 = (1+0.1)*getPeakGain(D2);
    W2 = K2/(1/w(4)*s+1);

end

% COMMENT:
% - filter type depends highly on damping factor (high or low pass)
% - cover uncertainty class wider because still missing modes
