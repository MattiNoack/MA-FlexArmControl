% WEIGHTING FILTER DESIGN

% Output filter:
% apol=poly([-1e-3 -1e-2]);
% czer=poly([]);
fac_of = 1; % gain factor 
Kof = czer(end)/apol(end)/fac_of; % gain->1
nfilt = length(apol)-1;
mfilt = length(czer);

Aof = [zeros(nfilt-1,1) eye(nfilt-1);-flip(apol(2:end))];
Bof = [zeros(nfilt-1,1);1];
Cof = [flip(czer)/Kof zeros(1,nfilt-mfilt)];
Dof = zeros(1,1);
sys_of = ss(Aof,Bof,Cof,Dof);

% Control filter:
apolc=poly([-5e-1]); % constant
czerc=poly([]);
% Au = 10;
fac_cf = 1/Au; % gain factor
Kcf = czerc(end)/apolc(end)/fac_cf; % gain->1
nfiltc = length(apolc)-1;
mfiltc = length(czerc);

Acf = [zeros(nfiltc-1,1) eye(nfiltc-1);-flip(apolc(2:end))];
Bcf = [zeros(nfiltc-1,1);1];
Ccf = [flip(czerc)/Kcf zeros(1,nfiltc-mfiltc)];
Dcf = zeros(nu,nu);
sys_cf = ss(Acf,Bcf,Ccf,Dcf);

% Resulting objective filter:
nf = nfilt+nfiltc;
Af = [Aof zeros(nfilt,nfiltc);zeros(nfiltc,nfilt) Acf];
Bf = [Bof zeros(nfilt,nu);zeros(nfiltc,1) Bcf];
Cf = [Cof zeros(1,nfiltc);zeros(nu,nfilt) Ccf];
Df = [Dof zeros(1,nu);zeros(nu,1) Dcf];
sys_f = ss(Af,Bf,Cf,Df);