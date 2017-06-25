% COMPARISON CONTROLLER SYNTHESIS

function [K1,K2,infoPID,gammH2,G] = cmpCtrl(P,channels,W)

    % Compound system setup:
    G = extractSysPID(P,channels,W);
    P2 = extractSysH2(P,channels,W);

    % PID compensator:
%     K1= constructPID(G);
    opt = pidtuneOptions('PhaseMargin',60); % not possible!
    [K1,infoPID] = pidtune(G,'PIDF',opt);
    
    % Optimal H2 performance:
    [K2ss,CL,gammH2,infoH2] = h2syn(P2,1,1);
    K2 = -tf(K2ss);
    
    % Evaluate:
    L1=minreal(G*K1); L2=minreal(G*K2);
    T1=minreal(L1/(1+L1)); T2=minreal(L2/(1+L2));
    figure; margin(L1)
    figure; margin(L2)
    figure; hold on
    step(T1); step(T2)
    figure; hold on
    nyquist(T1); nyquist(T2)

end

%---------------------%
% Auxiliary functions %
%---------------------%

function G = extractSysPID(P,channels,W)
    
    % Channel selection:
    iy = channels.contr.out;
    iu = channels.contr.in;
    iz = 1:iy(1)-1;
    iw = 1:iu(1)-1;

    % System matrices:  
    A = P.A;
    Bw = P.B(:,iw);
    Bu = P.B(:,iu);
    Cz = P.C(iz,:);
    Cy = P.C(iy,:);
    Dzw = P.D(iz,iw);
    Dzu = P.D(iz,iu);
    Dyw = P.D(iy,iw);
    Dyu = P.D(iy,iu);
    
    % Transfer fuction:
    G = tf(ss(A,Bu,Cy(1,:),Dyu(1,:)));
    
end

function [P2,ny,nu] = extractSysH2(P,channels,W)
    
    % Filter inclusion: (no in/out dimension change)
    if(isempty(W.Wdf)) % no uncertainty
        P = connect(P,W.Wof,{'w','u'},{'zf','y'});
        P = connect(P,W.Wif,{'wf','u'},{'zf','y'});
    else % uncertainty considered
        P = connect(P,W.Wof,{'wd','w','u'},{'zd','zf','y'});
        P = connect(P,W.Wif,{'wd','wf','u'},{'zd','zf','y'});
        P = connect(P,W.Wdf,{'wd','wf','u'},{'zdf','zf','y'});
    end
    
    % Channel defintions:
    iy = channels.contr.out;
    iu = channels.contr.in;
    iz = 1:iy(1)-1;
    iw = 1:iu(1)-1;

    % Channel selection:
    iw = [1 2];
    iz = [1 3];
    
    % System matrices:  
    A = P.A;
    Bw = P.B(:,iw);
    Bu = P.B(:,iu);
    Cz = P.C(iz,:);
    Cy = P.C(iy,:);
    Dzw = P.D(iz,iw);
    Dzu = P.D(iz,iu);
    Dyw = P.D(iy,iw);
    Dyu = P.D(iy,iu);
    
    n = size(A,1);
    nu = size(Bu,2); ny = size(Cy,1);
    nw = size(Bw,2); nz = size(Cz,1); 
    
%     % Integrator:
%     if(~isempty(channels.integ.meas))
%         iyi = channels.integ.meas; nyi = length(iyi);
%         iyn = setdiff(1:ny,iyi); nyn = length(iyn);
% 
%         A = [A zeros(n,nyi);Cy(iyi,:) zeros(nyi,nyi)];  
%         Bu = [Bu;zeros(nyi,nu)];
%         Bw = [Bw;Dyw(iyi,:)];
%         Cz = [Cz zeros(nz,nyi)];  
%         Cy = [zeros(nyi,n) eye(nyi);Cy(iyn,:) zeros(nyn,nyi)];
%         Dyw = [zeros(nyi,nw);Dyw(iyn,:)];
%         Dyu = zeros(ny,nu);
% 
%         n = n + nyi;
%     end 
    
    % Transfer fuction:
    ny = length(iy);
    nu = length(iu);
    P2 = ss(A,[Bw Bu],[Cz;Cy],[Dzw Dzu;Dyw Dyu]);
    
    % COMMENTS:
    %   - no Integrator possible because of D21=Dyw not full row rank
    
end

function Cpid = constructPID(G)

    % Select slowest poles:
    poles = 1./eig(G);
%     tau1 = -poles(end); % slowest pole
%     tau2 = -poles(end-1); % 2nd slowest pole
    tau1 = -poles(end); % BEWARE OF COMPLEX
    tau2 = -poles(end);
    
    % Derivative approximation:
    s = tf('s'); % transfer function construction
    TN = 3e-2; % choose relatively small
    K1hat = ((tau1*s+1).*(tau2*s+1)) / (s*(TN*s+1));
    
    % Calculate gain for desired phase margin:
    phi = -122; % PM=60°
    [mag,phase,wout] = bode(G*K1hat);
    mag_phi = interp1( squeeze(phase), squeeze(mag), phi);
    
    % Calculate coefficients:
    PI = 1./mag_phi;
    PP = (tau1+tau2-TN).*PI;
    PD = PI.*tau1.*tau2-PP*TN;
    Cpid = PP + PI/s + PD*s/(TN*s+1); 
    
    
    % COMMENTS:
    %   - PM=60° for auto pid tune not feasible (why?)
    %   - produces lots of oscillations
    %   - but really fast closed-loop system
    %   => NOT USED (matlab function instead)
end
