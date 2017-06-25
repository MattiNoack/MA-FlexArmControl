% WEIGHTING FILTER SELECTION
% for flexible robot arm with 1 mode considered

function W = weight_filter(tauo,taui,taud)

    % Case distinction for uncertainty:
    if(nargin<3)
        % No weighting:
        Adf = []; Bdf = []; Cdf = []; Ddf = [];            
    end

    % Case distinction for input:
    if(nargin<2)
        % No weighting:
        Aif = []; Bif = []; Cif = []; Dif = eye(2);    
    else
        % Iput filter: w1=r
        
        % Iput filter: w2=d
        
    end
        
    % Case distinction for output:
    if(nargin<1)
        % No weighting:
        Aof = []; Bof = []; Cof = []; Dof = eye(3);        
    else
        % Load precalculated & auxiliary structure:
        load('dyn_param.mat')
        N = length(gamm);

        % Output filter: z1=e
        if(tauo==-1)            
            Tfe1 = 1e-2;
            Tfe2 = 1.5e-1;            
        else % select time constants         
            Tfe1 = tauo(1);
            Tfe2 = tauo(2);
        end
        Ae = 2;
        tne = 1.8e-0;

        Ge = Ae*tf([tne 1],[Tfe1*Tfe2 Tfe1+Tfe2 1]); % low pass
        We = ss(Ge,'minimal');
        nfe = size(We.A,1);

        % Resulting objective filter:
        Aof = We.A;
        Bof = [We.B zeros(nfe,N+1)];
        Cof = [We.C;zeros(N+1,nfe)];
        Dof = blkdiag(0,eye(N+1));
    end    
    
    W.Wof = ss(Aof,Bof,Cof,Dof);
    W.Wof.u='z'; W.Wof.y='zf';
    
    W.Wif = ss(Aif,Bif,Cif,Dif);
    W.Wif.u='wf'; W.Wif.y='w';
    
    W.Wdf = ss(Adf,Bdf,Cdf,Ddf);
    W.Wdf.u='zd'; W.Wdf.y='zdf';

end