% TUNING PROCEDURE FOR MORC

%% Optimal Hinf gain
% Test-step for obtaining optimal Hinf gain:
tmpopt = options; tmpopt.optHinf = 1;
[Kopt,optparam] = synthMORC(P,channels,W,tmpopt,optin);
disp(''); disp('Optimal Hinf gain:'); disp(optparam.gamma) 

%% Optimal gH2 gain
% Test-step for obtaining optimal Hinf gain:
tmpopt = options; tmpopt.optGH2 = 1;
[Kopt,optparam] = synthMORC(P,channels,W,tmpopt,optin);
disp(''); disp('Optimal generalized H2 gain:'); disp(optparam.alpha) 

%% Numerical conditioning
% Test-step for obtaining optimal Lyapunov conditioning:
tmpopt = options; tmpopt.optLCond = 1;
[Kopt,optparam] = synthMORC(P,channels,W,tmpopt,optin);
disp(''); disp('Optimal condition coefficient:'); disp(optparam.t) 

%% Combined performance
% Test-step for combined Hinf/cond optimization:
tmpopt = options; tmpopt.optHinf = 1; tmpopt.optLCond = 1;
[Kopt,optparam] = synthMORC(P,channels,W,tmpopt,optin);
disp(''); disp('Optimal Hinf gain:'); disp(optparam.gamma)
disp(''); disp('Optimal condition coefficient:'); disp(optparam.t)

% Test-step for combined Hinf/H2 optimization:
tmpopt = options; tmpopt.optHinf = 1; tmpopt.optGH2 = 1;
[Kopt,optparam] = synthMORC(P,channels,W,tmpopt,optin);
disp(''); disp('Optimal Hinf gain:'); disp(optparam.gamma)
disp(''); disp('Optimal gH2 gain:'); disp(optparam.alpha)

%% Parameter selection