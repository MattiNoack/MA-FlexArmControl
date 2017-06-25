% WEIGHTING FILTER SELECTION
% for simple 2nd order system

function Wof = weight_filter(tau)

    % Case distinction:
    if(nargin<1)
        % No weighting:
        Aof = [];
        Bof = [];
        Cof = [];
        Dof = [];
    else
        % Load precalculated & auxiliary structure:
        load('dyn_param.mat')
        N = length(gamm);

        % Output filter: z1=e
        if(tau==-1)            
            Tfe1 = 1e-2;
            Tfe2 = 1.5e-1;            
        else % select time constants         
            Tfe1 = tau(1);
            Tfe2 = tau(2);
        end
        Ae = 2;
        tne = 1.8e-0;

        Ge = Ae*tf([tne 1],[Tfe1*Tfe2 Tfe1+Tfe2 1]); % low pass
        We = ss(Ge,'minimal');
        nfe = size(We.A,1);

        % Resulting objective filter:
        Aof = We.A;
        Bof = [We.B zeros(nfe,N)];
        Cof = [We.C;zeros(N,nfe)];
        Dof = blkdiag(0,eye(N));
    end    
    
    Wof = ss(Aof,Bof,Cof,Dof);    

end