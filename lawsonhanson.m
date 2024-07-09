%==========================================================================
% COMPRESSION OF A CUBATURE RULE BY FAST LAWSON-HANSON ALGORITHM
%==========================================================================
% INPUT
%==========================================================================
% A: cubature Vandermonde matrix,
% b: moments vector,
%==========================================================================
% OUTPUT
%==========================================================================
% weights: cbature weights.
%==========================================================================

function [weights] = lawsonhanson(A, b)


%==========================================================================
% OBJECT.
% Lawson-Hanson algorithm.
%==========================================================================
% INPUT
%==========================================================================
% A --- coefficient matrix,
% b --- observation vector,
% options_general --- specification of options as returned by 'initopt_general',
% options_specific --- specification of options as returned by 'opt_lawsonhanson'.
%==========================================================================
% OUTPUT
%==========================================================================
% out:
%==========================================================================

% IMPORTANT: 
% 1. SET 'stopcrit', 1. OTHERWISE MOMENTS ERRORS MAY BECOME LARGE.
% 2. SET 'perfmeas', {1}. ALSO 'perfmeas', {1,2} IS OK BUT SEEMS SLOWER 
%    WITH NO ADVANTAGE.
% 3. WE TRIED ALL THE subcase SETTINGS. ONLY subcase=1 SEEMS SLOWER. ALL THE
%    OTHER  SUBCASES LOOKS EQUAL.
opt_gen1 = initopt_general('perfmeas', {1}, 'maxcpu', 50, 'stopcrit', 1, 'tol', 1E-16);

subcase=0;

switch subcase
    case 1
        weights_struct = lawsonhanson_main(A,b, opt_gen1, opt_lawsonhanson('up',0,'down',0));
    case 2
        weights_struct = lawsonhanson_main(A,b, opt_gen1, opt_lawsonhanson('up',0,'down',1));
    case 3
        weights_struct = lawsonhanson_main(A,b, opt_gen1, opt_lawsonhanson('up',1,'down',0));
    case 4
        weights_struct = lawsonhanson_main(A,b, opt_gen1, opt_lawsonhanson('up',1,'down',1));
    otherwise
        weights_struct = lawsonhanson_main(A,b, opt_gen1, opt_lawsonhanson());
end

weights=weights_struct.xopt;






function [out] = lawsonhanson_main(A, b, options_general, options_specific)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if max(cellfun(@(z) length(z), options_general.perfmeas)) > size(A, 2)
    error('Dimension mismatch of A and the minimizer x')
end

if ~all(size(options_general.stopcrit) == size(options_general.tol))
    error('Structure of stopcrit not compatible with structure of tol')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options_general.usegram
    
    if options_specific.up
        
        switch options_specific.down
            case 0
                [out] = fnnls_up(A, b, options_general, options_specific);
            case 1
                [out] = fnnls_up_down(A, b, options_general, options_specific);
            case 2
                [out] = fnnls_up_down_c(A, b, options_general, options_specific);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        [out] = fnnls_pure(A, b, options_general, options_specific);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    if options_specific.up
        
        switch options_specific.down
            case 0
                [out] = lsqnnls_up(A, b, options_general, options_specific);
            case 1
                [out] = lsqnnls_up_down(A, b, options_general, options_specific);
            case 2
                [out] = lsqnnls_up_down_c(A, b, options_general, options_specific);
        end
    else
        [out] =  lsqnnls_pure(A,b, options_general, options_specific);
    end
    
end






%==========================================================================
%
% ADDITIONAL ROUTINES.
%
% 0. function options = initopt_general(varargin)
% 1. function options = opt_lawsonhanson(varargin)
% 2. function [out] = fnnls_pure(A, b, options_general, options_specific)
% 3. function [out] = fnnls_up_down_c(A, b, options_general, options_specific)
% 4. function [out] = fnnls_up_down(A, b, options_general, options_specific)
% 5. function [out] = fnnls_up(A, b, options_general, options_specific)
% 6. function [out] = lsqnnls_pure(A, b, options_general, options_specific)
% 7. function [out] = fnnls_up_down_c(A, b, options_general, options_specific)
% 8. function [out] = lsqnnls_up_down(A, b, options_general, options_specific)
% 9. function [out] = lsqnnls_up(A, b, options_general, options_specific)
%==========================================================================








%==========================================================================
%
% 1. function options = opt_lawsonhanson(varargin)
%
%==========================================================================

% Function to set-up options specific to
% the Lawson-Hanson algorithm
%
% Inputs
%
% 'up'  --- Updating of the Cholesky factorization
%           is done ('up' = 1) or not done ('up' = 0).
%           In the second case, a Cholesky factorization
%           is re-computed from scratch each time
%           a new variable is added into the active set.
%
% 'down' --- Downdating of the Cholesky factorization
%            is done using MATLAB code ('down' = 1),
%            C-code ('down' = 2) or not done ('down' = 0).
%
% 'outeritmax' --- maximum number of outer iterations (integer).
%                  default: 'inf'
%
% 'inneritmax' --- maximum number of inner iterations (integer).
%                  default: 'inf'


function options = opt_lawsonhanson(varargin)

p = inputParser;

validup = @(x)validateattributes(x,{'numeric'},{'scalar','binary'});
p.addParamValue('up',1,validup);

validdown = @(x)validateattributes(x,{'numeric'},{'scalar','integer','nonnegative','<=',2});
p.addParamValue('down',1,validdown);

validouteritmax = @(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'});
p.addParamValue('outeritmax', inf, validouteritmax);

validinneritmax = @(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'});
p.addParamValue('inneritmax', inf, validinneritmax);

p.parse(varargin{:});
options = p.Results;

if options.down~=0 && options.up == 0
    options.up = 1;
end







%==========================================================================
%
% 2. function [out] = fnnls_pure(A, b, options_general, options_specific)
%
%==========================================================================

% *internal code not to be called directly by the user*
function [out] = fnnls_pure(A, b, options_general, options_specific)
% solve NNLS problem by Lawson-Hanson's algorithm
% Gram matrix stored

if options_general.mode
    AtA = A;
    Atb = b;
end

%%% Initialize CPU time
t0 = cputime;
%%%
[n, d] = size(A);

%%%
if options_general.mv == 0
    if 1.1 * 2 * n < d
        mv = 2;
    else
        mv = 1;
    end
else
    mv = options_general.mv;
end

%%% un-pack Gram matrix if available; otherwise, compute Gram matrix.
if ~options_general.mode
    if ~isempty(options_general.gram)
        AtA = options_general.gram;
        if ~all(size(AtA) == [d d])
            error('wrong input for the Gram matrix');
        end
        options_general.gram = [];
    else
        AtA = A'*A;
    end
end
%%%

%initialize variables
if ~options_general.mode
    Atb = A' * b;
end

%%% Initialize objective -- global [possibly unused]
global f;
if options_general.mode
    f = @(x) x' * (AtA * x - 2 * Atb);
else
    f = @(x) norm(A * x - b)^2;
end

%%%
P = false(d,1); Z = true(d,1); % Initialize passive set P (xi>0) to null, Z (xi=0) to all
x = zeros(d,1); % final solution, unordered solution x and new temporary solution after an iteration

w = Atb;
wz = zeros(d,1);% -f'(x)


outeriter = 1; outeritmax = options_specific.outeritmax;
iter = 0; itmax = options_specific.inneritmax;

% grad = -w_
perf = getL(options_general.perfmeas, -w, x);
check = check_termination(t0, options_general, perf, -w, x, inf, x + inf);
perf = perf(~isnan(perf));
ConvSpeed =  [0 perf];
% fprev = inf; xprev = x + inf;
% check
if check.stop
    out.xopt = x;
    out.err = f(x);
    out.ConvSpeed = ConvSpeed;
    out.check = check;
end

% outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(Z) > eps) && (~check.stop)
    if outeriter == outeritmax
        disp('Exiting: Outer iteration count is exceeded.');
        break;
    end
    
    outeriter = outeriter + 1;
    % Reset intermediate solution z
    z = zeros(d,1);
    wz(P) = -Inf;
    wz(Z) = w(Z);
    % Find variable with largest Lagrange mutliplier
    [unused, t] = max(wz);
    % Move variable t from zero set to positive set
    P(t) = true;
    Z(t) = false;
    % Compute intermediate solution using only variables in positive set
    z(P) = A(:,P) \ b;
    % inner loop to remove elements from the positive set which no longer belong
    while any(z(P) <= eps)
        iter = iter + 1;
        if iter > itmax
            disp('Exiting: Inner iteration count is exceeded.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        
        % Find indices where intermediate solution z is approximately negative
        Q = (z <= eps) & P;
        % Choose new x subject to keeping new x nonnegative
        alpha = min(x(Q)./(x(Q) - z(Q)));
        if ~isfinite(alpha)
            disp('Warning: Inner iteration left. No strictly positive stepsize found.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        x = x + alpha*(z - x);
        % Reset Z and P given intermediate values of x
        Z = ((abs(x) < eps) & P) | Z; % remove indice st < tol from P to Z
        P = ~Z;
        z = zeros(d,1);           % Reset z
        z(P) = A(:,P)\ b;      % Re-solve for z
    end
    x = z;
    if mv == 1
        w = Atb - AtA(:,P) * x(P);
    else
        w = A' * (b - A(:,P) * x(P));
    end
    
    perf = getL(options_general.perfmeas, -w, x);
    check = check_termination(t0, options_general, perf, -w, x, check.fprev, check.xprev);
    ConvSpeed =  [ConvSpeed; [(cputime - t0) perf]];
    perf = perf(~isnan(perf));
    
end

out.xopt = x;
out.err = f(x);
out.ConvSpeed = ConvSpeed;
out.check = check;

function[perf] = getL(perfmeas, grad, x)
global f;

perf = NaN(1, 3);
%%% entry 1: --- KKT 
%%% entry 2: --- obj. 
%%% entry 3: --- error x - x^{\ast} 

for i = 1:length(perfmeas)
 
    if length(perfmeas{i}) > 1    
        perf(3) = getd(x, perfmeas{i});
    else
        if length(perfmeas{i} == 1)
            switch perfmeas{i}
                case 1
                perf(1) = getkktopt(grad, x);
                case 2
                perf(2) = f(x);
            end
        else
        continue 
        end
    end
end


function[kktopt] = getkktopt(grad, x)

num = abs(grad) .* (grad > eps);
score = num;

Pc = num < eps;

if sum(Pc) > 0
    score(Pc) = num(Pc);
end

if sum(~Pc) > 0
   score(~Pc) = num(~Pc) ./ x(~Pc);
end

Pc = num < eps | score < 1;
P  = ~Pc;

if sum(P) > 0
   kktP = max(abs(x(P)));
else
   kktP = 0;
end

if sum(Pc) > 0
   kktPc = max(abs(grad(Pc)));
else
   kktPc = 0; 
end

N = x < -eps;

if sum(N) > 0
   kktN = max(abs(x(N)));
else
   kktN = 0; 
end

kktopt = max([kktP kktPc kktN]);




% t0 -- cputime at the beginning
% options_general
% perf -- after calling function for performance measures
% grad --- gradient
% x --- current iterate
% fprev --- objective at previous iteration
% xprev --- iterate at previous iteration

function[check] = check_termination(t0, options_general, perf, grad, x, fprev, xprev)
global f;

reason = [];
tolreturn = [];
fcur = inf;

stop = 0;
if toc(t0) > options_general.maxcpu
    stop = 1;
    reason = 0; % code: reason=0 (maximum cpu time exceeded) 
end

stopcrit = options_general.stopcrit;
tol = options_general.tol;

[nr, nc] = size(stopcrit);

for i =1:nr
   tempstop = zeros(nc, 1); 
   
    for j = 1:nc
        
        critnumber = stopcrit(i,j);
        tolval = tol(i,j);
        
        if ~isnan(perf(critnumber))
            switch critnumber
                case 1
                    if perf(critnumber) < tolval
                        tolreturn = [tolreturn perf(critnumber)];
                        tempstop(j) = 1;
                        reason = [reason 1];
                    end
                case 2
                    fcur = perf(critnumber);
                    relobj = fprev - fcur;
                    if relobj < tolval
                        tempstop(j) = 1;
                        reason = [reason 2];
                        tolreturn = [tolreturn relobj];
                    end
                case 3
                    %%% * this block appears twice * %%%
                    relx = norm(x - xprev)/norm(0.5 * (x + xprev));
                    if relx < tolval
                        tempstop(j) = 1;
                        reason = [reason 3];
                        tolreturn = [tolreturn relx];
                    end
                    %%%
            end
        else
            switch critnumber
                case 1
                    kktopt = getkktopt(grad, x);
                    if kktopt < tolval
                       tempstop(j) = 1;
                       reason = [reason 1];
                       tolreturn = [tolreturn kktopt];
                    end 
                case 2
                     fcur = f(x);
                     relobj =  fprev - fcur;
                     if relobj < tolval
                       tempstop(j) = 1; 
                       reason = [reason 2];
                       tolreturn = [tolreturn relobj];
                     end      
                case 3
                    %%% * this block appears twice * %%%
                    relx = norm(x - xprev)/norm(0.5 * (x + xprev));
                    if relx < tolval
                        tempstop(j) = 1;
                        reason = [reason 3];
                        tolreturn = [tolreturn relx];
                    end
                    %%%
            end
        end
    end
   
    
    if(all(tempstop == 1)) 
        stop = 1;
    end
end



check.stop = stop;
%%% order
[sorted, ix] = sort(reason);
check.reason = reason(ix);
if ~isempty(tolreturn)
    if min(sorted) == 0
        ix = ix - 1;
        ix(sorted == 0) = [];
        check.tol = tolreturn(ix);
    else
        check.tol = tolreturn(ix);
    end
else
    check.tol = [];
end
%%% to be re-used at next call
check.xprev = x;
check.fprev = fcur;












%==========================================================================
%
% * 3. function [out] = fnnls_up_down_c(A, b, options_general, options_specific)
%
%==========================================================================

% *internal code not to be called directly by the user*
function [out] = fnnls_up_down_c(A, b, options_general, options_specific)
% solve NNLS problem by Lawson-Hanson's algorithm
% Gram matrix stored
% Cholesky updating/downdating used, downdating is implemented in C

if options_general.mode
    AtA = A;
    Atb = b;
end


%%% Initialize CPU time
t0 = tic;
%%%
[n, d] = size(A);

%%%
if options_general.mv == 0
    if 1.1 * 2 * n < d
        mv = 2;
    else
        mv = 1;
    end
else
    mv = options_general.mv;
end

%%% un-pack Gram matrix if available; otherwise, compute Gram matrix.
if ~options_general.mode
    %%% un-pack Gram matrix if available; otherwise, compute Gram matrix.
    if ~isempty(options_general.gram)
        AtA = options_general.gram;
        if ~all(size(AtA) == [d d])
            error('wrong input for the Gram matrix');
        end
        options_general.gram = [];
    else
        AtA = A'*A;
    end
end
%%%

%initialize variables
if ~options_general.mode
    Atb = A'*b;
end

%%% Initialize objective -- global [possibly unused]
global f;
if options_general.mode
    f = @(x) x' * (AtA * x - 2 * Atb);
else
    f = @(x) norm(A * x - b)^2;
end

%%%
P = false(d,1); Z = true(d,1); % Initialize passive set P (xi>0) to null, Z (xi=0) to all
Pnew = P; sizeP = 0; sizePnew = 0; % number of points all ready added into P
ltb = zeros(d,1); % look-up table: the points added into passive set are not ordered
Lt = zeros(d,d); % upper triangular matrix for Cholesky Decomposition
y = zeros(d,1); % intermediate solution: L * Lt * x = V >> L * y = V || Lt* x = y
x = zeros(d,1); % final solution, unordered solution x and new temporary solution after an iteration
w = Atb;  wz = zeros(d,1);% -f'(x)


outeriter = 1; outeritmax = options_specific.outeritmax;
iter = 0; itmax = options_specific.inneritmax;

% grad = -w_
perf = getL(options_general.perfmeas, -w, x);
check = check_termination(t0, options_general, perf, -w, x, inf, x + inf);
perf = perf(~isnan(perf));
ConvSpeed =  [0 perf];
% fprev = inf; xprev = x + inf;

% check
if check.stop
    out.xopt = x;
    out.err = f(x);
    out.ConvSpeed = ConvSpeed;
    out.check = check;
end

% outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(Z) > eps) && (~check.stop)
    if outeriter == outeritmax
        disp('Exiting: Outer iteration count is exceeded.');
        break;
    end
    
    outeriter = outeriter + 1;
    % Move k from active set to positive set
    wz(P) = -Inf;wz(Z) = w(Z);
    [unused,k] = max(wz);
    Pnew(k) = true; Z(k) = false; sizePnew = sizePnew + 1; ltb(sizePnew) = k;
    % Cholesky Updating to get xnew
    xnew = zeros(d,1);
    if sizePnew == 1
        Lt(1,1)= sqrt(AtA(k,k));
        y(1) = Lt(1,1) \ Atb(k);
        xnew(k) = Lt(1,1) \ y(1);
    else
        l = Lt(1:sizeP,1:sizeP)' \ AtA(ltb(1:sizeP),k) ;
        lkk = sqrt(AtA(k,k)- dot(l,l));
        Lt(1:sizeP,sizePnew) = l;
        Lt(sizePnew,sizePnew) = lkk;
        y(sizePnew) = (Atb(k) - dot(l,y(1:sizeP))) / lkk;
        xnew(ltb(1:sizePnew)) = Lt(1:sizePnew, 1:sizePnew) \ y(1:sizePnew);
    end
    P = Pnew;
    sizeP = sizePnew;
    
    % inner loop to remove elements from the positive set which no longer belong
    while any( xnew(P) <= eps )
        iter = iter + 1;
        if iter > itmax
            disp('Exiting: Inner iteration count is exceeded.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        
        % Find indices where temporary solution xnew is approximately
        % negative, backtrack x and reset Z/P(move some indice out from P),
        % and recompute x by Cholesky downdating
        Q = (xnew <= eps) & P;
        alpha = min(x(Q)./(x(Q) - xnew(Q)));
        if ~isfinite(alpha)
            disp('Warning: Inner iteration left. No strictly positive stepsize found.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        x = x + alpha*(xnew - x);
        
        
        tmp = find(abs(x)<eps & P);
        moveout = flipud(find(ismember(ltb, tmp)));
        
        % Cholesky Downdating
        for i = 1:length(moveout)
            k = moveout(i); ltb(k) = 0;
            Lt(1:sizeP-1, 1:sizeP-1) = choldown(Lt(1:sizeP,1:sizeP),k);
            Lt(:,sizeP) = 0;
            sizeP = sizeP - 1;
        end
        
        % update Z/P/look up table, resolve x/y
        Z = (abs(x)<eps & P) | Z;  P = ~Z; Pnew = P; sizePnew = sizeP;
        ltb1 = nonzeros(ltb); ltb = zeros(d,1); ltb(1:sizeP) = ltb1;
        y = zeros(d,1); xnew = zeros(d,1);
        y(1:sizeP) = Lt(1:sizeP,1:sizeP)'\ Atb(ltb(1:sizeP));
        xnew(ltb(1:sizeP)) = Lt(1:sizeP,1:sizeP) \ y(1:sizeP);
    end
    x = xnew;
    if mv == 1
        w = Atb - AtA(:,P) * x(P);
    else
        w = A' * (b - A(:,P) * x(P));
    end
    perf = getL(options_general.perfmeas, -w, x);
    check = check_termination(t0, options_general, perf, -w, x, check.fprev, check.xprev);
    perf = perf(~isnan(perf));
    ConvSpeed =  [ConvSpeed; [toc(t0) perf]];
    
    
    
    
end

out.xopt = x;
out.err = f(x);
out.ConvSpeed = ConvSpeed;
out.check = check;










%==========================================================================
%
% * 4. function [out] = fnnls_up_down(A, b, options_general, options_specific)
%
%==========================================================================

% *internal code not to be called directly by the user*
function [out] = fnnls_up_down(A, b, options_general, options_specific)
% solve NNLS problem by Lawson-Hanson's algorithm
% Gram matrix stored
% Cholesky updating/downdating used

if options_general.mode
    AtA = A;
    Atb = b;
end

%%% Initialize CPU time
t0 = tic;
%%%
[n, d] = size(A);

%%%
if options_general.mv == 0
    if 1.1 * 2 * n < d
        mv = 2;
    else
        mv = 1;
    end
else
    mv = options_general.mv;
end



%%% un-pack Gram matrix if available; otherwise, compute Gram matrix.
if ~options_general.mode
    %%% un-pack Gram matrix if available; otherwise, compute Gram matrix.
    if ~isempty(options_general.gram)
        AtA = options_general.gram;
        if ~all(size(AtA) == [d d])
            error('wrong input for the Gram matrix');
        end
        options_general.gram = [];
    else
        AtA = A'*A;
    end
end
%%%

%initialize variables
if ~options_general.mode
    Atb = A'*b;
end

%%% Initialize objective -- global [possibly unused]
global f;
if options_general.mode
    f = @(x) x' * (AtA * x - 2 * Atb);
else
    f = @(x) norm(A * x - b)^2;
end

P = false(d,1); Z = true(d,1); % Initialize passive set P (xi>0) to null, Z (xi=0) to all
Pnew = P; sizeP = 0; sizePnew = 0; % number of points all ready added into P
ltb = zeros(d,1); % look-up table: the points added into passive set are not ordered
Lt = zeros(d,d); % upper triangular matrix for Cholesky Decomposition
y = zeros(d,1); % intermediate solution: L * Lt * x = V >> L * y = V || Lt* x = y
x = zeros(d,1); % final solution, unordered solution x and new temporary solution after an iteration
w = Atb;  wz = zeros(d,1);% -f'(x)

outeriter = 1; outeritmax = options_specific.outeritmax;
iter = 0; itmax = options_specific.inneritmax;

% grad = -w_
perf = getL(options_general.perfmeas, -w, x);
check = check_termination(t0, options_general, perf, -w, x, inf, x + inf);
perf = perf(~isnan(perf));
ConvSpeed =  [0 perf];
% fprev = inf; xprev = x + inf;

% check
if check.stop
    out.xopt = x;
    out.err = f(x);
    out.ConvSpeed = ConvSpeed;
    out.check = check;
end

% outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(Z) > eps) && (~check.stop)
    if outeriter == outeritmax
        disp('Exiting: Outer iteration count is exceeded.');
        break;
    end
    
    outeriter = outeriter + 1;
    % Move k from active set to positive set
    wz(P) = -Inf;wz(Z) = w(Z);
    [unused,k] = max(wz);
    Pnew(k) = true; Z(k) = false; sizePnew = sizePnew + 1; ltb(sizePnew) = k;
    % Cholesky Updating to get xnew
    xnew = zeros(d,1);
    if sizePnew == 1
        Lt(1,1)= sqrt(AtA(k,k));
        y(1) = Lt(1,1) \ Atb(k);
        xnew(k) = Lt(1,1) \ y(1);
    else
        l = Lt(1:sizeP,1:sizeP)' \ AtA(ltb(1:sizeP),k) ;
        lkk = sqrt(AtA(k,k)- dot(l,l));
        Lt(1:sizeP,sizePnew) = l;
        Lt(sizePnew,sizePnew) = lkk;
        y(sizePnew) = (Atb(k) - dot(l,y(1:sizeP))) / lkk;
        xnew(ltb(1:sizePnew)) = Lt(1:sizePnew, 1:sizePnew) \ y(1:sizePnew);
    end
    P = Pnew;
    sizeP = sizePnew;
    
    % inner loop to remove elements from the positive set which no longer belong
    while any( xnew(P) <= eps )
        iter = iter + 1;
        if iter > itmax
            disp('Exiting: Inner iteration count is exceeded.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        
        % Find indices where temporary solution xnew is approximately
        % negative, backtrack x and reset Z/P(move some indice out from P),
        % and recompute x by Cholesky downdating
        Q = (xnew <= eps) & P;
        alpha = min(x(Q)./(x(Q) - xnew(Q)));
        if ~isfinite(alpha)
            disp('Warning: Inner iteration left. No strictly positive stepsize found.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        x = x + alpha*(xnew - x);
        
        
        tmp = find(abs(x)<eps & P);
        moveout = flipud(find(ismember(ltb, tmp)));
        
        % Cholesky Downdating
        for i = 1:length(moveout)
            k = moveout(i); ltb(k) = 0;
            Lt(1:sizeP-1, 1:sizeP-1) = choldownmatlab1(Lt(1:sizeP,1:sizeP),k);
            Lt(:,sizeP) = 0;
            sizeP = sizeP - 1;
        end
        
        % update Z/P/look up table, resolve x/y
        Z = (abs(x)<eps & P) | Z;  P = ~Z; Pnew = P; sizePnew = sizeP;
        ltb1 = nonzeros(ltb); ltb = zeros(d,1); ltb(1:sizeP) = ltb1;
        y = zeros(d,1); xnew = zeros(d,1);
        y(1:sizeP) = Lt(1:sizeP,1:sizeP)'\ Atb(ltb(1:sizeP));
        xnew(ltb(1:sizeP)) = Lt(1:sizeP,1:sizeP) \ y(1:sizeP);
    end
    x = xnew;
    if mv == 1
        w = Atb - AtA(:,P) * x(P);
    else
        w =  A' * (b - A(:,P) * x(P));
    end
    perf = getL(options_general.perfmeas, -w, x);
    check = check_termination(t0, options_general, perf, -w, x, check.fprev, check.xprev);
    perf = perf(~isnan(perf));
    ConvSpeed =  [ConvSpeed; [toc(t0) perf]];
    
    
    
    
end

out.xopt = x;
out.err = f(x);
out.ConvSpeed = ConvSpeed;
out.check = check;




function Lkt =  choldownmatlab1(Lt, k)
% This is for non-sparse matrix
% cholesky downdating
% A in R^(n,p)
% G = A'* A = L * L', where L, L' come from cholesky decomposition
% now  removes kth column from A, denoted by Ak. Gk := Ak' * Ak
% Given L' and k, choldown computes the chol. decomposition of  Gk
% i.e. Lk' * Lk = Gk, without processing of A, G

p = length(Lt);

% drop the kth clm of Lt
Temp = Lt;
Temp(:,k) = []; % Temp in R^(p,p-1)

% Givens Rotations
for i = k:p-1,
    a = Temp(i,i);
    b = Temp(i+1,i);
    r = sqrt(sum(Lt(:,i+1).^2) - sum(Temp(1:i-1,i).^2));
    c =  r * a / (a^2+b^2);
    s =  r * b / (a^2+b^2);
    % ith row of rotation matrix H
    Hrowi = zeros(1,p); Hrowi(i) = c; Hrowi(i+1) = s;
    % (i+1)th row of ration matrix H
    Hrowi1 = zeros(1,p); Hrowi1(i) = -s; Hrowi1(i+1) = c;
    % modify the ith and (i+1)th rows of Temp
    v = zeros(2,p-1);
    v(1,i:p-1) = Hrowi * Temp(:,i:p-1);
    v(2,i+1:p-1) = Hrowi1 * Temp(:,i+1:p-1);
    Temp(i:i+1,:) =  v;
end

% drop the last row
Lkt = Temp(1:p-1,:);









%==========================================================================
%
% * 5. function [out] = fnnls_up(A, b, options_general, options_specific)
%
%==========================================================================

% *internal code not to be called directly by the user*
function [out] = fnnls_up(A, b, options_general, options_specific)
% solve NNLS problem by Lawson-Hanson's algorithm
% Gram matrix stored
% Cholesky updating used

if options_general.mode
    AtA = A;
    Atb = b;
end

%%% Initialize CPU time
t0 = tic;
%%%
[n, d] = size(A);

%%%
if options_general.mv == 0
    if 1.1 * 2 * n < d
        mv = 2;
    else
        mv = 1;
    end
else
    mv = options_general.mv;
end

if ~options_general.mode
    %%% un-pack Gram matrix if available; otherwise, compute Gram matrix.
    if ~isempty(options_general.gram)
        AtA = options_general.gram;
        if ~all(size(AtA) == [d d])
            error('wrong input for the Gram matrix');
        end
        options_general.gram = [];
    else
        AtA = A'*A;
    end
end
%%%

%initialize variables
if ~options_general.mode
    Atb = A'*b;
end

%%% Initialize objective -- global [possibly unused]
global f;
if options_general.mode
    f = @(x) x' * (AtA * x - 2 * Atb);
else
    f = @(x) norm(A * x - b)^2;
end

P = false(d,1); Z = true(d,1); % Initialize passive set P (xi>0) to null, Z (xi=0) to all
Pnew = P; sizeP = 0; sizePnew = 0; % number of points all ready added into P
ltb = zeros(d,1); % look-up table: the points added into passive set are not ordered
Lt = zeros(d,d); % upper triangular matrix for Cholesky Decomposition
y = zeros(d,1); % intermediate solution: L * Lt * x = V >> L * y = V || Lt* x = y
x = zeros(d,1); % final solution, unordered solution x and new temporary solution after an iteration
w = Atb;  wz = zeros(d,1);% -f'(x)


outeriter = 1; outeritmax = options_specific.outeritmax;
iter = 0; itmax = options_specific.inneritmax;

% grad = -w_
perf = getL(options_general.perfmeas, -w, x);
check = check_termination(t0, options_general, perf, -w, x, inf, x + inf);
perf = perf(~isnan(perf));
ConvSpeed =  [0 perf];
% fprev = inf; xprev = x + inf;

% check
if check.stop
    out.xopt = x;
    out.err = f(x);
    out.ConvSpeed = ConvSpeed;
    out.check = check;
end

% outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(Z) > eps) && (~check.stop)
    if outeriter == outeritmax
        disp('Exiting: Outer iteration count is exceeded.');
        break;
    end
    
    outeriter = outeriter + 1;
    % Move k from active set to positive set
    wz(P) = -Inf;wz(Z) = w(Z);
    [unused,k] = max(wz);
    Pnew(k) = true; Z(k) = false; sizePnew = sizePnew + 1; ltb(sizePnew) = k;
    % Cholesky Updating to get xnew
    xnew = zeros(d,1);
    if sizePnew == 1
        Lt(1,1)= sqrt(AtA(k,k));
        y(1) = Lt(1,1) \ Atb(k);
        xnew(k) = Lt(1,1) \ y(1);
    else
        l = Lt(1:sizeP,1:sizeP)' \ AtA(ltb(1:sizeP),k) ;
        lkk = sqrt(AtA(k,k)- dot(l,l));
        Lt(1:sizeP,sizePnew) = l;
        Lt(sizePnew,sizePnew) = lkk;
        y(sizePnew) = (Atb(k) - dot(l,y(1:sizeP))) / lkk;
        xnew(ltb(1:sizePnew)) = Lt(1:sizePnew, 1:sizePnew) \ y(1:sizePnew);
    end
    P = Pnew;
    sizeP = sizePnew;
    
    % inner loop to remove elements from the positive set which no longer belong
    while any( xnew(P) <= eps )
        iter = iter + 1;
        if iter > itmax
            disp('Exiting: Inner iteration count is exceeded.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        
        % Find indices where temporary solution xnew is approximately
        % negative, backtrack x and reset Z/P(move some indice out from P),
        % and recompute x.
        Q = (xnew <= eps) & P;
        alpha = min(x(Q)./(x(Q) - xnew(Q)));
        if ~isfinite(alpha)
            disp('Warning: Inner iteration left. No strictly positive stepsize found.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        x = x + alpha*(xnew - x);
        
        % Reset Z and P given intermediate values of x
        Z = ((abs(x) < eps) & P) | Z;
        P = ~Z; Pnew = P;
        xnew = zeros(d,1);           % Reset xnew
        % solve for xnew
        ltb = zeros(d,1);
        y = zeros(d,1);
        [K, unused_j, unused_s] = find(P);
        sizeP = length(K); sizePnew = sizeP;
        ltb(1:sizeP) = K;
        
        Lt = zeros(d,d);
        Lt(1:sizeP,1:sizeP) = chol(AtA(K,K));
        y(1:sizeP) = Lt(1:sizeP,1:sizeP)' \ Atb(K);
        xnew(K) = Lt(1:sizeP,1:sizeP) \ y(1:sizeP);
    end
    x = xnew;
    if mv == 1
        w = Atb - AtA(:,P) * x(P);
    else
        w = A' * (b - A(:,P) * x(P));
    end
    perf = getL(options_general.perfmeas, -w, x);
    check = check_termination(t0, options_general, perf, -w, x, check.fprev, check.xprev);
    perf = perf(~isnan(perf));
    ConvSpeed =  [ConvSpeed; [toc(t0) perf]];
    
end

out.xopt = x;
out.err = f(x);
out.ConvSpeed = ConvSpeed;
out.check = check;



















%==========================================================================
%
% * 6. function [out] = lsqnnls_pure(A, b, options_general, options_specific)
%
%==========================================================================

% *internal code not to be called directly by the user*
function [out] = lsqnnls_pure(A, b, options_general, options_specific)
% solve NNLS problem by Lawson-Hanson's algorithm
% Cholesky updating/downdating used

%%% Initialize objective -- global [possibly unused]
global f;
f = @(x) norm(A * x - b)^2;
%%% Initialize CPU time
t0 = tic;
%initialize variables
[n, d] = size(A);
x = zeros(d,1);

% active/passive set
P = false(d,1);
Z = true(d,1);
wz = zeros(d,1);
w = A' * b;

outeriter = 1; outeritmax = options_specific.outeritmax;
iter = 0; itmax = options_specific.inneritmax;

% grad = -w_
perf = getL(options_general.perfmeas, -w, x);
check = check_termination(t0, options_general, perf, -w, x, inf, x + inf);
perf = perf(~isnan(perf));
ConvSpeed =  [0 perf];
% fprev = inf; xprev = x + inf;

% check
if check.stop
    out.xopt = x;
    out.err = f(x);
    out.ConvSpeed = ConvSpeed;
    out.check = check;
end


% outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(Z) > eps) && (~check.stop)
    
    if outeriter == outeritmax
        disp('Exiting: Outer iteration count is exceeded.');
        break;
    end
    
    outeriter = outeriter + 1;
    % Move k from active set to positive set
    z = zeros(d,1);
    wz(P) = -Inf;
    wz(Z) = w(Z);
    % Find variable with largest Lagrange mutliplier
    [unused,t] = max(wz);
    % Move variable t from zero set to positive set
    P(t) = true;
    Z(t) = false;
    % Compute intermediate solution using only variables in positive set
    z(P) = A(:,P) \ b;
    
    % check feasibility of z, if not reset x and compute z again
    while any(z(P) <= eps)
        iter = iter + 1;
        if iter > itmax
            disp('Exiting: Inner iteration count is exceeded.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        
        % Find indices where temporary solution xnew is approximately
        % negative, backtrack x and reset Z/P(move some indice out from P),
        % and recompute x by Cholesky downdating
        Q = (z <= eps) & P;
        alpha = min(x(Q)./(x(Q) - z(Q)));
        if ~isfinite(alpha)
            disp('Warning: Inner iteration left. No strictly positive stepsize found.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        x = x + alpha*(z - x);
        % Reset Z and P given intermediate values of x
        Z = ((abs(x) < eps) & P) | Z; % remove indice st < tol from P to Z
        P = ~Z;
        z = zeros(d,1);           % Reset z
        z(P) = A(:,P)\ b;      % Re-solve for z
    end
    x = z;
    w = A' * (b - A(:,P) * x(P));
    
    perf = getL(options_general.perfmeas, -w, x);
    check = check_termination(t0, options_general, perf, -w, x, check.fprev, check.xprev);
    perf = perf(~isnan(perf));
    ConvSpeed =  [ConvSpeed; [toc(t0) perf]];
    
end

out.xopt = x;
out.err = f(x);
out.ConvSpeed = ConvSpeed;
out.check = check;





















%==========================================================================
%
% * 7. function [out] = fnnls_up_down_c(A, b, options_general, options_specific)
%
%==========================================================================

% *internal code not to be called directly by the user*
function [out] = lsqnnls_up_down_c(A, b, options_general, options_specific)
% solve NNLS problem by Lawson-Hanson's algorithm
% Cholesky updating/downdating used

%%% Initialize objective -- global [possibly unused]
global f;
f = @(x) norm(A * x - b)^2;
%%% Initialize CPU time
t0 = tic;
%initialize variables
[n, d] = size(A);
Atb = A' * b;

% active/passive set
P = [];
Z = 1:d;  % switch from Boolean representation [Q.] to index representation [M.]
Pnew = P;
sizeP = 0;
sizePnew = 0; % number of points already added into P
%

Lt = []; % upper triangular matrix for Cholesky Decomposition
x = zeros(d, 1); % current iterate
xnew = zeros(d, 1); % iterate obtained after solving the LS problem.
y = []; % temporary vector of variable size
w = A' * (b - A*x);

outeriter = 1; outeritmax = options_specific.outeritmax;
iter = 0; itmax = options_specific.inneritmax;

% grad = -w_
perf = getL(options_general.perfmeas, -w, x);
check = check_termination(t0, options_general, perf, -w, x, inf, x + inf);
perf = perf(~isnan(perf));
ConvSpeed =  [0 perf];
% fprev = inf; xprev = x + inf;

% check
if check.stop
    out.xopt = x;
    out.err = f(x);
    out.ConvSpeed = ConvSpeed;
    out.check = check;
end


% outer loop to put variables into set to hold positive coefficients
while ~isempty(Z) && any(w(Z) > eps) && (~check.stop)
    
    if outeriter == outeritmax
        disp('Exiting: Outer iteration count is exceeded.');
        break;
    end
    
    outeriter = outeriter + 1;
    % Move k from active set to positive set
    wz = w(Z);
    [unused, k] = max(wz);
    k = Z(k);
    Pnew = [P k];
    Z = setdiff(Z, k);
    sizePnew = length(Pnew);
    
    % Cholesky Updating to get xnew (temporary)
    xnew = zeros(d,1);
    
    if sizePnew == 1
        Lt(1,1) = norm(A(:,k));
        y(1) = Lt(1,1) \ Atb(k);
        xnew(k) = Lt(1,1) \ y(1);
    else
        l = Lt' \ (transpose(A(:,P)) * A(:, k));
        lkk = sqrt(sum(A(:,k).^2) - sum(l.^2));
        Lt(1:sizeP,  sizePnew) = l;
        Lt(sizePnew, sizePnew) = lkk;
        y = [y; (Atb(k) - dot(l,y)) / lkk];
        xnew(Pnew) = Lt \ y;
    end
    P = Pnew;
    sizeP = sizePnew;
    
    % check feasibility of xnew, if not reset x and compute xnew again
    while any( xnew(P) <= eps )
        iter = iter + 1;
        if iter > itmax
            disp('Exiting: Inner iteration count is exceeded.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        
        % Find indices where temporary solution xnew is approximately
        % negative, backtrack x and reset Z/P(move some indice out from P),
        % and recompute x by Cholesky downdating
        Q = P(xnew(P) <= eps);
        alpha = min(x(Q)./(x(Q) - xnew(Q)));
        if ~isfinite(alpha)
            disp('Warning: Inner iteration left. No strictly positive stepsize found.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        x = x + alpha * (xnew - x);
        
        
        tmp = find(abs(x(P)) < eps);
        Z = [Z P(tmp)];
        P(tmp) = [];
        
        
        % Cholesky Downdating
        for i = length(tmp):-1:1
            Lt = choldown(Lt, tmp(i));
            sizeP = sizeP - 1;
        end
        
        % solve for x/y
        y = Lt' \ Atb(P);
        xnew(P) = Lt \ y;
        
    end
    x = xnew;
    w = A' * (b - A(:,P) * x(P));
    
    perf = getL(options_general.perfmeas, -w, x);
    check = check_termination(t0, options_general, perf, -w, x, check.fprev, check.xprev);
    perf = perf(~isnan(perf));
    ConvSpeed =  [ConvSpeed; [toc(t0) perf]];
    
end

out.xopt = x;
out.err = f(x);
out.ConvSpeed = ConvSpeed;
out.check = check;





















%==========================================================================
%
% * 8. function [out] = lsqnnls_up_down(A, b, options_general, options_specific)
%
%==========================================================================


% *internal code not to be called directly by the user*
function [out] = lsqnnls_up_down(A, b, options_general, options_specific)
% solve NNLS problem by Lawson-Hanson's algorithm
% Cholesky updating/downdating used

%%% Initialize objective -- global [possibly unused]
global f;
f = @(x) norm(A * x - b)^2;
%%% Initialize CPU time
t0 = tic;
%initialize variables
[n, d] = size(A);
Atb = A' * b;

% active/passive set
P = [];
Z = 1:d;  % switch from Boolean representation [Q.] to index representation [M.]
Pnew = P;
sizeP = 0;
sizePnew = 0; % number of points already added into P
%

Lt = []; % upper triangular matrix for Cholesky Decomposition
x = zeros(d, 1); % current iterate
xnew = zeros(d, 1); % iterate obtained after solving the LS problem.
y = []; % temporary vector of variable size
w = A' * (b - A*x);

outeriter = 1; outeritmax = options_specific.outeritmax;
iter = 0; itmax = options_specific.inneritmax;

% grad = -w_
perf = getL(options_general.perfmeas, -w, x);
check = check_termination(t0, options_general, perf, -w, x, inf, x + inf);
perf = perf(~isnan(perf));
ConvSpeed =  [0 perf];
% fprev = inf; xprev = x + inf;
% check
if check.stop
    out.xopt = x;
    out.err = f(x);
    out.ConvSpeed = ConvSpeed;
    out.check = check;
end


% outer loop to put variables into set to hold positive coefficients
while ~isempty(Z) && any(w(Z) > eps) && (~check.stop)
    
    if outeriter == outeritmax
        disp('Exiting: Outer iteration count is exceeded.');
        break;
    end
    
    outeriter = outeriter + 1;
    % Move k from active set to positive set
    wz = w(Z);
    [unused, k] = max(wz);
    k = Z(k);
    Pnew = [P k];
    Z = setdiff(Z, k);
    sizePnew = length(Pnew);
    
    % Cholesky Updating to get xnew (temporary)
    xnew = zeros(d,1);
    
    if sizePnew == 1
        Lt(1,1) = norm(A(:,k));
        y(1) = Lt(1,1) \ Atb(k);
        xnew(k) = Lt(1,1) \ y(1);
    else
        l = Lt' \ (transpose(A(:,P)) * A(:, k));
        lkk = sqrt(sum(A(:,k).^2) - sum(l.^2));
        Lt(1:sizeP,  sizePnew) = l;
        Lt(sizePnew, sizePnew) = lkk;
        y = [y; (Atb(k) - dot(l,y)) / lkk];
        xnew(Pnew) = Lt \ y;
    end
    P = Pnew;
    sizeP = sizePnew;
    
    % check feasibility of xnew, if not reset x and compute xnew again
    while any( xnew(P) <= eps )
        iter = iter + 1;
        if iter > itmax
            disp('Exiting: Inner iteration count is exceeded.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        
        % Find indices where temporary solution xnew is approximately
        % negative, backtrack x and reset Z/P(move some indice out from P),
        % and recompute x by Cholesky downdating
        Q = P(xnew(P) <= eps);
        alpha = min(x(Q)./(x(Q) - xnew(Q)));
        if ~isfinite(alpha)
            disp('Warning: Inner iteration left. No strictly positive stepsize found.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        x = x + alpha * (xnew - x);
        
        
        tmp = find(abs(x(P)) < eps);
        Z = [Z P(tmp)];
        P(tmp) = [];
        
        
        % Cholesky Downdating
        for i = length(tmp):-1:1
            Lt = choldownmatlab(Lt, tmp(i));
            sizeP = sizeP - 1;
        end
        
        % solve for x/y
        y = Lt' \ Atb(P);
        xnew(P) = Lt \ y;
        
    end
    x = xnew;
    w = A' * (b - A(:,P) * x(P));
    
    perf = getL(options_general.perfmeas, -w, x);
    check = check_termination(t0, options_general, perf, -w, x, check.fprev, check.xprev);
    perf = perf(~isnan(perf));
    ConvSpeed =  [ConvSpeed; [toc(t0) perf]];
    
end

out.xopt = x;
out.err = f(x);
out.ConvSpeed = ConvSpeed;
out.check = check;



function Lkt =  choldownmatlab(Lt, k)
% This is for non-sparse matrix
% cholesky downdating
% A in R^(n,p)
% G = A'* A = L * L', where L, L' come from cholesky decomposition
% now  removes kth column from A, denoted by Ak. Gk := Ak' * Ak
% Given L' and k, choldown computes the chol. decomposition of  Gk
% i.e. Lk' * Lk = Gk, without processing of A, G

p = length(Lt);

% drop the kth clm of Lt
Temp = Lt;
Temp(:,k) = []; % Temp in R^(p,p-1)

% Givens Rotations
for i = k:p-1,
    a = Temp(i,i);
    b = Temp(i+1,i);
    r = sqrt(sum(Lt(:,i+1).^2) - sum(Temp(1:i-1,i).^2));
    c =  r * a / (a^2+b^2);
    s =  r * b / (a^2+b^2);
    % ith row of rotation matrix H
    Hrowi = zeros(1,p); Hrowi(i) = c; Hrowi(i+1) = s;
    % (i+1)th row of ration matrix H
    Hrowi1 = zeros(1,p); Hrowi1(i) = -s; Hrowi1(i+1) = c;
    % modify the ith and (i+1)th rows of Temp
    v = zeros(2,p-1);
    v(1,i:p-1) = Hrowi * Temp(:,i:p-1);
    v(2,i+1:p-1) = Hrowi1 * Temp(:,i+1:p-1);
    Temp(i:i+1,:) =  v;
end

% drop the last row
Lkt = Temp(1:p-1,:);













%==========================================================================
%
% * 9. function [out] = lsqnnls_up(A, b, options_general, options_specific
%
%==========================================================================

% *internal code not to be called directly by the user*
function [out] = lsqnnls_up(A, b, options_general, options_specific)
% solve NNLS problem by Lawson-Hanson's algorithm
% Cholesky updating/downdating used

%%% Initialize objective -- global [possibly unused]
global f;
f = @(x) norm(A * x - b)^2;
%%% Initialize CPU time
t0 = tic;
%initialize variables
[n, d] = size(A);
Atb = A' * b;

% active/passive set
P = [];
Z = 1:d;  % switch from Boolean representation [Q.] to index representation [M.]
Pnew = P;
sizeP = 0;
sizePnew = 0; % number of points already added into P
%

Lt = []; % upper triangular matrix for Cholesky Decomposition
x = zeros(d, 1); % current iterate
xnew = zeros(d, 1); % iterate obtained after solving the LS problem.
y = []; % temporary vector of variable size
w = A' * (b - A*x);

outeriter = 1; outeritmax = options_specific.outeritmax;
iter = 0; itmax = options_specific.inneritmax;

% grad = -w_
perf = getL(options_general.perfmeas, -w, x);
check = check_termination(t0, options_general, perf, -w, x, inf, x + inf);
perf = perf(~isnan(perf));
ConvSpeed =  [0 perf];
% fprev = inf; xprev = x + inf;

% check
if check.stop
    out.xopt = x;
    out.err = f(x);
    out.ConvSpeed = ConvSpeed;
    out.check = check;
end


% outer loop to put variables into set to hold positive coefficients
while ~isempty(Z) && any(w(Z) > eps) && (~check.stop)
    
    if outeriter == outeritmax
        disp('Exiting: Outer iteration count is exceeded.');
        break;
    end
    
    outeriter = outeriter + 1;
    % Move k from active set to positive set
    wz = w(Z);
    [unused, k] = max(wz);
    k = Z(k);
    Pnew = [P k];
    Z = setdiff(Z, k);
    sizePnew = length(Pnew);
    
    % Cholesky Updating to get xnew (temporary)
    xnew = zeros(d,1);
    
    if sizePnew == 1
        Lt(1,1) = norm(A(:,k));
        y(1) = Lt(1,1) \ Atb(k);
        xnew(k) = Lt(1,1) \ y(1);
    else
        l = Lt' \ (transpose(A(:,P)) * A(:, k));
        lkk = sqrt(sum(A(:,k).^2) - sum(l.^2));
        Lt(1:sizeP,  sizePnew) = l;
        Lt(sizePnew, sizePnew) = lkk;
        y = [y; (Atb(k) - dot(l,y)) / lkk];
        xnew(Pnew) = Lt \ y;
    end
    P = Pnew;
    sizeP = sizePnew;
    
    % check feasibility of xnew, if not reset x and compute xnew again
    while any( xnew(P) <= eps )
        iter = iter + 1;
        if iter > itmax
            disp('Exiting: Inner iteration count is exceeded.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        
        % Find indices where temporary solution xnew is approximately
        % negative, backtrack x and reset Z/P(move some indice out from P),
        % and recompute x by Cholesky downdating
        Q = P(xnew(P) <= eps);
        alpha = min(x(Q)./(x(Q) - xnew(Q)));
        if ~isfinite(alpha)
            disp('Warning: Inner iteration left. No strictly positive stepsize found.');
            out.xopt = x;
            out.err = f(x);
            out.ConvSpeed = ConvSpeed;
            out.check = check;
            return
        end
        x = x + alpha * (xnew - x);
        
        
        tmp = find(abs(x(P)) < eps);
        Z = [Z P(tmp)];
        P(tmp) = [];
        sizeP = length(P);
        
        Lt = chol(transpose(A(:,P)) * A(:,P));
        y = Lt' \ Atb(P);
        xnew(P) = Lt \ y;
        
    end
    x = xnew;
    w = A' * (b - A(:,P) * x(P));
    
    perf = getL(options_general.perfmeas, -w, x);
    check = check_termination(t0, options_general, perf, -w, x, check.fprev, check.xprev);
    perf = perf(~isnan(perf));
    ConvSpeed =  [ConvSpeed; [toc(t0) perf]];
    
end

out.xopt = x;
out.err = f(x);
out.ConvSpeed = ConvSpeed;
out.check = check;








%==========================================================================
%
% * 0. function options = initopt_general(varargin)
%
%==========================================================================

function options = initopt_general(varargin)
%
% Function to set-up a series of options regarding operation
% and termination of the various NNLS solvers.
% The function is called prior to calling one of the solvers,
% and the output (a structure with names equal to those of the inputs listed below)
% is then passed as argument 'options_general'.
%
%
% Input
%
%   'perfmeas'      - performance measure(s) to be tracked while the solver progresses
%		      Specification is via a cell of length at most 3,
%		      the elements of which can take either one of the two values {1,2}
%                     or a non-negative vector that serves as 'ground truth' x^* of the
%                     NNLS problem. In case that such a vector is provided, norm(x^k - x^*)
%                     is tracked. The values {1,2} encode the following two
%                     performance measures:
%			1: KKT optimality
%			2: function value f(x^k).
%		      default: {1}, i.e. only KKT optimality is tracked.
%
%
%   'maxcpu'        - maximum cpu time (in seconds).
%  			The solver terminates after the iteration during which 'maxcpu'
%                       is exceeded first. Time is measured using MATLAB's 'tic' and 'toc'.
%			default: inf, i.e. the solver runs until all of the
%			stopping criteria are fulfilled (see below).
%
%    'stopcrit'      - stopping criteria. A vector of length at most 3, the elements
%                      of which take values in the set {1,2,3}.
%			1: KKT optimality
%			2: f(x^k) - f(x^(k+1)), i.e. convergence of the function values
%			3: norm(x^(k+1) - x^k)/(0.5 * norm((x^(k+1) + x^k))),i.e. convergence of the iterates.
%		      If multiple stopping criteria are specified, _all_ of them have to be
%		      satisfied upon termination of the algorithm.
%
%    'tol'          - Numerical tolerances for the stopping criteria as specified in 'stopcrit'.
%                     Must be a vector of the same length as 'stopcrit'. default: [1E-10]
%
%
%   'usegram'       - whether to compute ('usegram' = 1) or not to compute ('usegram' = 0, default) the Gram matrix A' * A
%
%   'gram'          - precomputed Gram matrix to be used by the solver. default: empty matrix [].
%
%   'mv'            - How to perform the matrix-vector multiplications for gradient evaluation
%		      0: 'mv' is determined automatically based on the dimensions of A.
% 		      1: the gradient is evaluated as (A' * A) * x - A' * b [favourable if A is a tall matrix]
%		      2: the gradient is evaluated as A' * (b - A * x) [favourable if A is a fat matrix]
%
%
%
%   'mode'          - 0: The NNLS problem with inputs A and b is solved (default).
%                   - 1: The corresponding linear complementarity problem is solved,
%                        in which the arguments (A,b) passed to the solver are interpreted
%	                 as follows:
%                        A ~ A' * A,
%                        b ~ A' * b.
%                      	 Note that the specification 'mode = 1' only makes sense if 'usegram=1'.

psr =  inputParser;

validPerfmeas = @(x)iscell(x) && ~isempty(x) && length(x)<=3;
psr.addParamValue('perfmeas', {1}, validPerfmeas);

validMaxcpu = @(x)validateattributes(x,{'numeric'},{'scalar','positive'});
psr.addParamValue('maxcpu',inf,validMaxcpu);

validStopcrit = @(x)validateattributes(x,{'numeric'},{'vector','positive','integer','<=',3});
psr.addParamValue('stopcrit',1,validStopcrit);

validTol = @(x)validateattributes(x,{'numeric'},{'vector','positive'});
psr.addParamValue('tol',10^-10,validTol);

validUsegram = @(x)validateattributes(x,{'numeric'},{'scalar','binary'});
psr.addParamValue('usegram',0,validUsegram);

validGram = @(x)validateattributes(x,{'numeric'},{'2d'});
psr.addParamValue('gram',[],validGram);

validMv = @(x)validateattributes(x,{'numeric'},{'scalar','integer','nonnegative','<=',2});
psr.addParamValue('mv',0,validMv);

validMode = @(x)validateattributes(x,{'numeric'},{'scalar','integer','nonnegative','<=',1});
psr.addParamValue('mode',0,validMode);


% Parse and validate all input arguments.
psr.parse(varargin{:});
options = psr.Results;

options.stopcrit = unique(options.stopcrit);
if length(options.tol)~= length(options.stopcrit)
    error('The length of the vector of tolerances should be equal to the number of stopping criteria.');
end

if ~isempty(options.gram)
    options.usegram = true;
end

if  options.usegram == 0
    
    
    if  options.mode == 1
        display('mode 1: automatically set usegram = 1.')
        options.usegram = 1;
    end
    
    if options.mv == 1
        display('mv 1: automatically set usegram = 1.')
        options.usegram = 1;
    end
    
end

% example
% initopt_general()
% initopt_general('usegram',true)
% initopt_general('mv',2,'tol',10^-8);