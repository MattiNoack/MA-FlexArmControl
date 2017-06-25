% UNCERTAINTY ANALYSIS
% additive uncertainty due to high frequency component on output

function [D1,D2] = high_freq_uncert()

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
        D1 = D1 + g(i)*k(i)/(s^2+w(i)^2);
        D2 = D2 - Th/2*p(i)*k(i)/(s^2+w(i)^2);
    end
    
    % Analyze magnitude:
    figure
    subplot(1,2,1); bodemag(D1); grid; title('Uncertainty \Delta_1(j\omega)')
    subplot(1,2,2); bodemag(D2); grid; title('Uncertainty \Delta_2(j\omega)')

end

% COMMENT:
% - without damping there are resonance effects
% - cannot be covered with filter for guaranteering norm bound
% - in the idealistic case RS cannot be implented that way!
% => include damping

