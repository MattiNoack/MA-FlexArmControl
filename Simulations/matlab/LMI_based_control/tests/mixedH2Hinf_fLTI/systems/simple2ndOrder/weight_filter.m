% WEIGHTING FILTER SELECTION
% for simple 2nd order system

function Wof = weight_filter()

    % Output filter: z1=e
    Ae = 2;
    Tfe1 = 1e-2;
    Tfe2 = 2e-2;
    tne = 1.5e-2;

    Ge = Ae*tf([tne 1],[Tfe1*Tfe2 Tfe1+Tfe2 1]); % low pass
    We = ss(Ge,'minimal');
    nfe = size(We.A,1);

    % Control filter: z2=u
    Au = 1e-1;
    Tfu = 1e-1;

%     Gu = Au*tf([1 0],[Tfu 1]); % high pass
    Gu = Au*tf([1],[1]); % static weight
    Wu = ss(Gu,'minimal');
    nfu = size(Wu.A,1);

    % Resulting objective filter:
    nof = nfe+nfu;
    Aof = [We.A zeros(nfe,nfu);zeros(nfu,nfe) Wu.A];
    Bof = [We.B zeros(nfe,1);zeros(nfu,1) Wu.B];
    Cof = [We.C zeros(1,nfu);zeros(1,nfe) Wu.C];
    Dof = [We.D 0;0 Wu.D];
    Wof = ss(Aof,Bof,Cof,Dof);

end

% Remark: channelwise filter (leads to 1s) -> SISO filters