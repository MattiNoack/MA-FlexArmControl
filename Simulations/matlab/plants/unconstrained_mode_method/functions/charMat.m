% CONDITIONED CHARACTERISTIC MATRIX CALCULATION FOR BC
% For 4th order solution:
%       phi(x) = A*sin(beta*x)+B*cos(beta*x)+C*sinh(beta*x)+D*cosh(beta*x)

% Input:
%       coeff   array with BC information of the form (col->equ.):
%                    c1*phi(X)+...+c4*D3phi(X) = 0
%       X       row vector containing the respective spatial point in coeff
%       beta    symbolic eigenvalue variable

% Output:
%       M       characteristic matrix for respective BVP (symbolic)


function M = charMat(coeff,X,beta)

    % Check input:
    if (size(coeff,1)~=4)||(size(coeff,1)~=4)||(length(X)~=4)
        disp('@charMat -> wrong input argument sizes')
        M = 0; return
    end

    % Definitions:
    syms x A B C D real
    phi = A*sin(beta*x)+B*cos(beta*x)+C*sinh(beta*x)+D*cosh(beta*x);
    Dphi = diff(phi,x);
    D2phi = diff(phi,x,2);
    D3phi = diff(phi,x,3);
    
    % Gradient calculation:
    Dvec = [phi Dphi D2phi D3phi]';
    param = [A B C D]';
    for i=1:4
        f(i) = subs(coeff(i,:)*Dvec,x,X(i));
    end
    M0 = jacobian(f,param);
    
    % Row normalization/conditioning:
    M(1,:) = M0(1,:)/norm(M0(1,:));
    M(2,:) = M0(2,:)/norm(M0(2,:));
    M(3,:) = M0(3,:)/norm(M0(3,:));
    M(4,:) = M0(4,:)/norm(M0(4,:));

end