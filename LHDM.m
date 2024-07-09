
function [x,resnorm,exitflag, outeriter] = LHDM(C,d,options,verbose)

% LHDM solves underdetermined linear least squares with nonnegativity constraints
% using classic Lawson-Hanson algorithm accelerated by Deviation Maximization.

% INPUT
% C: underdetermined matrix
% d: column array of set points coordinates
% gefftol: G-efficiency threshold
% maxit: maximum number of iterations of the multiplicative algorithm
% options: structure containing the values of optimization parameters,
%    possible field are:
%    init  true if ULS initialization of Passive set is desired, false
%          otherwise
%    tol   tolerance on the projected residual, stop criterion
%    k     maximum number of indices to be added to the Passive set at each
%          iteration
%    thres threshold on the cosine of the angle between each pair of
%          columns cuncurrently added to the Passive set, value between
%          0 and 1

% OUTPUT
% x: sparse vector that minimizes NORM(d-C*x)
% resnorm: squared 2-norm of the residual: norm(d-C*X)^2
% exitflag: exit condition, possible value are:
%    1  LHDM converged with a solution X
%    0  Iteration count was exceeded, increasing the tolerance
%       (OPTIONS.Tol) may lead to a solution
% momerr: moment reconstruction error by the compressed measure

% 11/06/2020
% M. Dessole, F. Marcuzzi, M. Vianello

if nargin < 2
    error('MATLAB:LHDM:NotEnoughInputs',...
        getString(message('MATLAB:LHDM:NotEnoughInputs')));
end

if nargin < 3
    options=[];
end

if nargin < 4
    verbose=0;
end

[m,n] = size(C);

% Initialize vector of n zeros and Infs (to be used later)
nZeros = zeros(n,1);
wz = nZeros;

itmax = m*2;

% Initialize set of non-active columns to null
P = false(n,1);
cardP = 0;
% Initialize set of active columns to all and the initial point to zeros
Z = true(n,1);
x = nZeros;

% Check if options was created with optimoptions
if ~isempty(options)
    if ~isa(options,'struct')
        error('MATLAB:LDHM:ArgNotStruct',...
            getString(message('MATLAB:LHDM:commonMessages:ArgNotStruct')));
    end
    if isfield(options,'thres')
        thres = options.thres;
        if (thres <= eps)
            LHDMflag = 0;
        end
    else
        thres = 0.2222;
    end
    if isfield(options,'thres_w')
        thres_w = options.thres_w;
    else
        thres_w = 0.8;
    end
    if isfield(options,'k')
        k = options.k;
        if (k == 1)
            LHDMflag = 0;
        else
            LHDMflag = 1;
        end
    else
        k = ceil(m/20);
    end
    if isfield(options,'tol')
        tol = options.tol;
    else
        tol = 10*eps*norm(C,1)*length(C);
    end
    if isfield(options,'init')
        if options.init
            xtmp = C\d;
            Idx = find(xtmp>0);
            % Initialize set of non-active columns
            P(Idx) = true;
            % Initialize set of active columns
            Z(Idx) = false;
            % Initialize starting point
            x(P) = C(:,P)\d;
            tmp = find(x<0);
            if(size(tmp,1)>0)
                x(tmp) = 0;
                P(tmp) = false;
                Z(tmp) = true;
            end
            if verbose
                fprintf('LHDM(%d) with ULS inizialization\n', k);
            end
            cardP = sum(P);
        else
            if verbose
                fprintf('LHDM(%d) \n', k);
            end
        end
    else
        if verbose
            fprintf('LHDM(%d) \n', k);
        end
    end
else
    thres   = 0.2222;
    thres_w = 0.8;
    k = ceil(m/20);
    tol = 10*eps*norm(C,1)*length(C);
    LHDMflag = 1;
    if verbose
        fprintf('LHDM(%d) \n', k);
    end
end


if LHDMflag
    Cnorm = zeros(size(C));
    for j=1:n
        Cnorm(:,j) = C(:,j)/norm(C(:,j));
    end
end


resid = d - C*x;
w = C'*resid;

% Set up iteration criterion
outeriter = 0;
totiter = 0;


tol=tol*10^(-3);

% Outer loop to put variables into set to hold positive coefficients
while (any(Z) && (any(w(Z) > tol) || any(x(P) < 0) ) && (totiter < itmax) )
    outeriter = outeriter + 1;
    totiter = totiter+1;
    % Create wz, a Lagrange multiplier vector of variables in the zero set.
    % wz must have the same size as w to preserve the correct indices, so
    % set multipliers to -Inf for variables outside of the zero set.
    wz(P) = -Inf;
    wz(Z) = w(Z);
    
    if ((outeriter <= 1) || (~LHDMflag))
        % in first iteration we drop DM to ensure w_T > 0
        [~,t] = max(wz);
    elseif LHDMflag && (removedP ~= addedP)
        t = DM(Cnorm, wz, k, thres, thres_w);
        if (length(t) + cardP > m)
            t = t(1:m-cardP);
        end
        
    else
        [~,t] = max(wz);
    end
    % number of indices added to P
    addedP = length(t);
    
    % Reset intermediate solution z
    z = zeros(size(x));
    % Move variable t from zero set to positive set
    P(t) = true;
    Z(t) = false;
    % Compute intermediate solution using only variables in positive set
    z(P) = C(:,P)\d;
    
    % inner loop to remove elements from the positive set which no longer belong
    iter = 0;
    removedP = 0;
    while (any(z(P) <= 0) && (totiter < itmax))
        totiter = totiter +1;
        iter = iter+1;
        % Find indices where intermediate solution z is approximately negative
        Q = (z <= 0) & P;
        % Choose new x subject to keeping new x nonnegative
        b = x(Q)./(x(Q) - z(Q));
        alpha = min(b);
        x = x + alpha*(z - x);
        % number of indices removed from P
        t = find(((abs(x) < tol) & P))';
        removedP = removedP + length(t);
        % Reset Z and P given intermediate values of x
        Z = ((abs(x) < tol) & P) | Z;
        P = ~Z;
        % Reset z
        z = nZeros;
        % Compute intermediate solution using only variables in positive set
        z(P) = C(:,P)\d;
    end
    % update cardinality of P
    cardP = cardP + addedP-removedP;
    x=z;
    
    resid = d - C*x;
    w = C'*resid;
    
end

if (outeriter < itmax)
    exitflag = 1;
else
    exitflag = 0;
end

resnorm = resid'*resid;



function [p] = DM(Cnorm, wz, k, thres,thres_w)
% Deviation Maximization

[wzI, I] = sort(wz, 'descend');
t = I(1);
p = t;
thres_wloc = thres_w*wzI(1);
C = find(wzI>thres_wloc);
n = size(C,1);
add = 1;
for i = 2:n
    c = C(i);
    if (max(abs(Cnorm(:,I(c))'*Cnorm(:,p))) < thres)
        p = [I(c) p];
        add = add+1;
    end
    if (add >= k)
        break;
    end
end
