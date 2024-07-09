function [xyw,res,Z,Zin,cmom,bbox,itlevel] = cubRS(ade,structure_RS,...
    structure_RSH,extraction_type,Nbox,safe_mode)

%-------------------------------------
% OBJECT:
%-------------------------------------
% This routine computes an algebraic formula with positive weights and
% interior nodes with degree of exactness "ade" with respect to the
% Lebesgue/Legendre measure in a Jordan NURBS domain (or one more generally
% defined parametrically by piecewise rational splines) by the arrays of
% splines "SxN", "SxD", "SyN", "SyD", sharing the same checkpoints.
% This code works with boundaries defined by parametric splines or
% composite Bezier curves.
%-------------------------------------
% INPUT:
%-------------------------------------
% ade: Algebraic Degree of Exactness of the rule.
%
% structure_RS: array whose k-th element is a "structure" containing the
%    informations of the k-th portion of the boundary "dD" of the domain
%    "D".
%
%  * structure_RS.P: it is a N x 2 matrix of coordinates of control
%        points. The k-th control point is described as the k-th row of P.
%      In case it is a spline it represents the couples (X(i),Y(i)) that
%      the parametric splines x=x(t), y=y(t) must match.
%
%  * structure_RS.w: vector of weights; it is a row vector 1 x N,
%        where "N" is the number of control points.
%
%  * structure_RS.knots: vector of knots;
%     1. in case the boundary is described by NURBS, following reference
%        [1], it is a row vector 1 x (N+L) of the form
%                   [k0*ones(1,L) u k1*ones(1,L)],
%        where "L" is the order of the NURBS and "u" is a vector of non
%        decreasing values, of length "N-L" where "N" is the number of
%        control points. Observe that
%
%       max(structure_RS{k}.knots) <= min(structure_RS{k+1}.knots).
%
%     2. in case the boundary is defined by a splines the variable
%        "structure_RS.knots" represents the break points of the spline.
%
%   Note that for understanding how to define the array of structures
%   "structure_RS", one can use the pertinent demos present in this
%   package.
%
%  * structure_RS.order: NURBS or spline order (i.e. "NURBS degree +1").
%
% A. The routine can be modified to study boundaries, that are defined by
% rational piecewise functions, in which the degree of the numerator and
% denominator is the same, or alternatively the denominator is equal to 1.
% This structure is typical of common splines or Composite Bezier curves.
%
% B. The variable "structure_RS" can be easily defined (see "demo").
%
% structure_RSH: as structure_RS but for the internal boundary;
%
% extraction_type: it sets how Tchakaloff rule is extracted
%      1: lsqnonneg, 2: lawsonhanson, 3: LHDM 4: \
%      This variable is not mandatory. If not declared or equal to "[]",
%      the default is LHDM.
%
% Nbox: not mandatory parameter, useful for indomain (default: 100). It is
%     a technical variable that if correctly set allows the indomain
%     routine to be faster.
%--------------------------------
% OUTPUT:
%-------------------------------
% xyw: 3-column array of cartesian coordinates of the nodes stored in the
%      first two columns and relative weights, stored in the third column.
%      Thus, the k-th node is "(xyw(k,1),xyw(k,2))" and has weight
%      "xyw(k,3)". Only this output variable is mandatory for the cubature
%      process.
%
% res: norm 2 of compressed moments vs moments, i.e.
%                          res=norm(V'*w-cmom);
%      where "V" is the Vandermonde matrix in the nodes of a certain Cheb.
%      tensorial basis of total degree "ade" and "cmom" is the vector of
%      the Chebyshev-Vandermonde moments;
%
% Zin: matrix M x 2 of M internal points from which the nodes are extracted
%       by Lawson-Hanson algorithm;
%
% ZL  : all analysed points from which only Zin are used in compression;
%
% cmom: vector of the Chebyshev-Vandermonde moments; as basis it is
%      considered the tensor product Chebyshev polynomials of total degree
%      "ade" with a suitable ordering, scaled respect to the aforementioned
%      rectangle "R".
%
% bbox: it is a vector that defines the smaller rectangle with side
%      parallel to axis, containing the domain.
%      Such rectangle is [bbox(1), bbox(2)] x [bbox(3), bbox(4)];
%
% itlevel: iterations made be the process; the routine defines at each
%      iteration a set of points from which Lawson-Hanson algorithm
%      attempts to extract the nodes. In case of failure, it adds a denser
%      set, trying again to find Tchakaloff nodes and weights. The higher
%      "itlevel", the higher is the cardinality of the set "Z" described
%      above.
%------------------------------------
% Routines.
%-----------------------------------
% 1. chebmom (attached, as well as subroutines, except "legendre_10000.m");
% 2. inRS (external);
% 3. cvand (attached).
% 4. LHDM (external)
%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2007-2024 Chiara Zara, Alvise Sommariva, Marco Vianello.
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%
%% Author:  
%%          Chiara Zara
%%          Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Date: July 8, 2024 -
%---------------------------------
% Dates
%---------------------------------
% First version: October 31, 2021;
% Checked: July 8, 2024.
%---------------------------------


% ....................... troubleshooting ...


% Extraction method: 1: lsqnonneg, 2: lawsonhanson, 3: LHDM 4: \
if nargin < 4, extraction_type=3; end
if length(extraction_type) == 0, extraction_type=3; end

% Number of boxes inside each monotone box.
if nargin < 5, Nbox=100; end
if length(Nbox) == 0, Nbox=100; end

% Number of boxes inside each monotone box.
if nargin < 6, safe_mode=100; end
if length(safe_mode) == 0, safe_mode=1; end

% ......................... initialize .......

% Discretization parameters, defining the points of the initial grid and
% following ones in case of search of finer discretizations.
alpha=1.5; beta=1.5;

% The initial pointset is a tensorial grid in a small rectangle "R"
% containing the domain, with sides parallel to the axis. A tensorial rule
% of "k" equispaced points per direction is used. Thus the initial set will
% have "k^2=k*k" points.
k=floor(max(ade,5)^(alpha));

% Note: on the setting of restol: 10^(-14) seems fine, 10^(-15) is often
% not achieved. This quantity is normalized with the number of extracted
% nodes.
restol=(1*10^(-14))*sqrt((ade+1)*(ade+2)/2);

% The variable "Z" will contain all the points used in the small rectangle
% "R" containing the domain, defined above.
Z=[];

% The variable "Zin" will contain all the points of the set defined by "Z",
% that are internal to the domain.
Zin=[];

% The matrix "V" will store the Vandermonde matrix, w.r.t. to the tensor
% product Chebyshev polynomials of total degree "ade", scaled respect
% to the aforementioned rectangle "R".
V=[];

% The variable "Lmax" defines the maximum number of iterations of the
% indomain process to obtain a set "Zin" good enough for extraction of
% Tchakaloff points.
Lmax=6;



% ....................... Determine rectangle "R" ..

% Determine NURBS domain control points.
S=length(structure_RS); P=[];
for kk=1:S
    structure_RSL=structure_RS(kk); PL=structure_RSL.P; P=[P;PL];
end

% Establish indomain structures.
[~,on,sp_xnV,sp_xdV,sp_ynV,sp_ydV,boxVx,singular_points_X,...
    singular_points_Y]=inRS([],structure_RS,Nbox,safe_mode);

SH=length(structure_RSH); PH=[];
for kk=1:SH
    structure_RSLH=structure_RSH(kk); PLH=structure_RSLH.P; PH=[PH;PLH];
end

% Establish indomain structures.
[~,onH,sp_xnVH,sp_xdVH,sp_ynVH,sp_ydVH,boxVxH,singular_points_XH,...
    singular_points_YH]=inRS([],structure_RSH,Nbox,safe_mode);

% Compute Tens.Chebyshev-Vandermonde moments.
[cmomE,bboxE]=chebmom(ade,sp_xnV,sp_xdV,sp_ynV,sp_ydV);
[cmomH,bboxH]=chebmom(ade,sp_xnVH,sp_xdVH,sp_ynVH,sp_ydVH,bboxE);
cmom=cmomE-cmomH;



% ......................... Method iterations .....


for itlevel=1:Lmax

    % .......................... indomain ..........

    % Set test point "ZL" from which we try to extract Tchakaloff nodes.
    x=linspace(bboxE(1),bboxE(2),k); y=linspace(bboxE(3),bboxE(4),k);
    [u,v]=meshgrid(x,y); ZL=[u(:) v(:)];

    % Find points of ZL in the domain and store them in "X".
    [insideE,onE]=inRS(ZL,structure_RS,Nbox,safe_mode,...
        sp_xnV,sp_xdV,sp_ynV,sp_ydV,boxVx,singular_points_X,...
        singular_points_Y);

    [insideH,onH]=inRS(ZL,structure_RSH,Nbox,safe_mode,...
        sp_xnVH,sp_xdVH,sp_ynVH,sp_ydVH,boxVxH,singular_points_XH,...
        singular_points_YH);

    % Points in which "on" is "NaN" are not too close to the boundary,
    % i.e. with distance larger that a tolerance 10^(-12), or with abscissa
    % too close to vertical points.
    ind=find((insideE == 1) & (isnan(onE)) & ...
        not((insideH == 1) & (isnan(onH))) );
    X=ZL(ind,:);

    if length(X) == 0
        % Mesh refinement in case of failures (no points in domain)
        % fprintf('\n \t * no point in domain');
        k=floor(beta*(k+1));

    else

        % The variable "Zin" contains all the points that along the
        % iteration are determined as internal to the domain.
        Zin=[Zin; X];


        % The matrix "VL" will store the Vandermonde matrix, w.r.t. to the
        % tensor product Chebyshev polynomials of total degree "ade",
        % scaled in the respect to the aforementioned rectangle "R",
        % relatively to the pointset "X".

        VL=cvand(ade,X,bboxE);
        V=[V; VL];      % Vandermonde of all points inside the domain
        Z=[Z; ZL];      % All points: inside/outside domain
        [Q,R]=qr(V,0);  % Orthogonalize Vandermonde matrix (more stable).

        cmom1=R'\cmom;  % Moments w.r.t. the orthog. basis.

        % .......................... compression ..........

        % Perform extraction of the compressed formula from nodes in "Zin".
        switch extraction_type
            case 1
                w=lsqnonneg(Q',cmom1);
            case 2
                w = lawsonhanson(Q', cmom1);
            case  3
                w = LHDM(Q', cmom1);
            otherwise
                w=Q'\cmom1;
        end

        % The nodes of the Tchakaloff formula are those points of "Zin"
        % with strictly positive weights in case a Lawson-Hanson algorithm
        % has been used.
        % In case the "backslash" method is adopted, the nodes are those
        % whose weights are nonzero.

        if (extraction_type == 1) | (extraction_type == 1) | ...
                (extraction_type == 3)
            ind1=find(w>0); w=w(ind1);
        else
            ind1=find(abs(w)>0); w=w(ind1);
        end

        % ......................... error analysis ..........

        % Evaluate residual in norm 2, that says hw much the compressed
        % rule matches the moments.

        res=norm(V(ind1,:)'*w-cmom);

        % In case of failure determine the number of equispaced points per
        % direction that determine the rule, otherwise exit from the
        % routine.

        if res>restol | length(w) == 0
            k=floor(beta*k);
            % In case of more iterations, relax the tolerance
            %  "restol" for the norm 2 of the moment matching.
            restol=10*restol;
        else
            xyw=[Zin(ind1,:) w];
            return;
        end

    end

end


% In case of failure, provide warnings and try to give a rule with internal
% nodes but possibly negative weights. In most of the cases the formula has
% "low" cubature condition number (i.e. is almost optimally stable).

fprintf(2,'\n \t Compression not completed. Moment error: %1.3e',res);

w=Q'\cmom1; ind1=find(abs(w)>0); w=w(ind1);
res=norm(V(ind1,:)'*w-cmom);
if res>restol | length(w) == 0
    fprintf(2,'\n \t Tried backslash (PO rule): fail. Mom. error: %1.3e',...
        res);
else

    ineg=length(find(w < 0));
    fprintf(2,'\n \t Tried backslash (PO rule): OK. Mom. error: %1.3e ',...
        res);
    fprintf(2,'\n \t NEGATIVE WEIGHTS: %8.0f ',...
        ineg);
end

xyw=[Zin(ind1,:) w];



%=====================
% ADDITIONAL ROUTINES
%=====================


%=====================
% cvand
%====================
function V = cvand(n,pts,bbox)

%-----------------------------
% OBJECT:
%-----------------------------
% The routine computes by recurrence the Chebyshev-Vandermonde matrix
% a) at a arbitrarily located mesh defined by "pts",
% b) in a total-degree shifted product Chebyshev basis on a given rectangle
%   defined by the vector "bbox" via minimum/maximum of abscissas and
%   ordinates.
%-----------------------------
% INPUT:
%-----------------------------
% n: polynomial degree
% pts: 2-column array of mesh point coordinates
% bbox: [bbox(1),bbox(2)] x [bbox(3),bbox(4)] bounding box for the mesh.
%----------------------------
% OUTPUT:
%----------------------------
% V: Chebyshev-Vandermonde matrix
%----------------------------
% DATES:
%----------------------------
% built: March 2019
% check: November 13, 2019
%----------------------------

% Default rectangle containing the mesh if not passed
if isempty(bbox)
    bbox=[min(pts(:,1)) max(pts(:,1)) min(pts(:,2)) max(pts(:,2))];
end

a=bbox(1);b=bbox(2);c=bbox(3);d=bbox(4);

Vx=chebpolys(n,(2*pts(:,1)-b-a)/(b-a));
Vy=chebpolys(n,(2*pts(:,2)-d-c)/(d-c));

k=0;
for i=0:n
    for j=0:n-i
        k=k+1; V(:,k)=Vx(:,i+1).*Vy(:,j+1);
    end
end




%========================
% chebpolys
%========================

function T=chebpolys(deg,x)

%------------------------
% OBJECT:
%------------------------
% computes the Chebyshev-Vandermonde matrix on the real line by recurrence
%-----------------------
%INPUT:
%-----------------------
% deg = maximum polynomial degree
% x = 1-column array of abscissas
%----------------------
% OUTPUT:
%----------------------
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg
%----------------------
% DATES:
%----------------------
% built: March 2019
% check: May 2, 2019
%----------------------

% Inizialization.
T=zeros(length(x),deg+1);
t0=ones(length(x),1); T(:,1)=t0;
t1=x; T(:,2)=t1;

% 3-term recurrence.
for j=2:deg, t2=2*x.*t1-t0; T(:,j+1)=t2; t0=t1; t1=t2; end




%========================
% chebmom
%========================

function [cmom,bbox] = chebmom(ade,sp_xnV,sp_xdV,sp_ynV,sp_ydV,bbox)

%-------------------------
% OBJECT:
%-------------------------
% It computes the moments up to degree m of a total-degree product
% Chebyshev basis with respect to the Lebesgue measure in a Jordan domain
% whose boundary is given counterclockwise by the rational splines
% arcs (Sx,Sy) defined by "sp_xn", "sp_xd", "sp_yn", "sp_yd" as
%         Sx=sp_xn/sp_xd         Sy=sp_yn/sp_yd
%-------------------------
% INPUT:
%-------------------------
% ade: polynomial degree
% sp_xn,sp_xd,sp_yn,sp_yd: arrays of spline structures;
%   (Sx(i),Sy(i)) is the i-th arc the arcs must be counterclockwise
%   concatenated forming a Jordan curve
%              Sx(i).breaks(end)=Sx(i+1).breaks(1),
%              Sy(i).breaks(end)=Sy(i+1).breaks(1)
%   i=1,...,end, with "Sx=sp_xn/sp_xd", "Sy=sp_yn/sp_yd".
%--------------------------
% OUTPUT:
%--------------------------
% cmom: array of bivariate Chebyshev moments w.r.t. a certain tensor Cheb.
%       basis of total degree "ade".
% bbox: vector that the determines an approximation of the domain bounding
%       box.
%---------------------------
% DATES:
%---------------------------
% First version: October 31, 2021;
% Checked: November 12, 2021.
%---------------------------
% NOTE:
%---------------------------
% From this routine, we compute the approximation of a rectangle "bbox"
% that is a bounding box for the spline polygon. For low "n" such rectangle
% may not be good enough. Thus, in this case, we increase a little the
% degree so to be acceptable.
%----------------------------

xyw=lineint(ade,sp_xnV,sp_xdV,sp_ynV,sp_ydV);

if nargin < 6
    bbox=[min(xyw(:,1)) max(xyw(:,1)) min(xyw(:,2)) max(xyw(:,2))];
end

intV=intcvand(ade,xyw(:,1:2),bbox);
cmom=intV'*xyw(:,3);





%==============================
% lineint
%==============================

function xyw = lineint(ade,sp_xnV,sp_xdV,sp_ynV,sp_ydV)

%------------------------------
% OBJECT:
%------------------------------
% This routine computes weights and nodes for the line integral on Jordan
% domains whose boundary is given counterclockwise by the rational splines.
% The formula is exact on bivariate polynomials up to algebraic degree of
% precision "ade".
%-------------------------------
% INPUT:
% ade: polynomial degree to be integrated, via its x-primitive, on the
%    domain.
% sp_xn,sp_xd,sp_yn,sp_yd: arrays of spline structures;
%    (Sx(i),Sy(i)) is the i-th arc the arcs must be counterclockwise
%    concatenated forming a Jordan curve
%         Sx(i).breaks(end)=Sx(i+1).breaks(1),
%         Sy(i).breaks(end)=Sy(i+1).breaks(1)
%         i=1,...,end, with end+1:=1
%    with "Sx=sp_xn/sp_xd", "Sy=sp_yn/sp_yd".
%-----------------------------
% OUTPUT:
%----------------------------
% xyw: 3-column array of nodes coords (cols 1-2) and weights (col 3)
%----------------------------
% DATES:
%-----------------------------
% built: March 2019
% check: November 12, 2021
%----------------------------


LL=length(sp_xnV);
xyw=[];

% analysing k-th block
for kk=1:LL

    % splines that define locally the boundary
    sp_xn=sp_xnV(kk); sp_xd=sp_xdV(kk); % x(t)=sp_xn(t)/sp_xd(t)
    sp_yn=sp_ynV(kk); sp_yd=sp_ydV(kk); % y(t)=sp_yn(t)/sp_yd(t)

    intv=sp_xn.breaks;

    % analysing i-th block subinterval
    for i=1:length(intv)-1

        p1=sp_xn.coefs(i,:); p2=sp_yn.coefs(i,:);
        q1=sp_xd.coefs(i,:); q2=sp_yd.coefs(i,:);

        ord_RSnum=length(p1);
        a=intv(i); b=intv(i+1);

        % .............. determining line-integral formula  ....

        % appropriate Gauss-Legendre rule
        if (norm(q1,1) == 1 & q1(end) == 1) & ...
                (norm(q2,1) == 1 & q2(end) == 1)
            % polynomial type (ord_RSnum=ord_RSden);
            xw=gaussian_quadrature(ord_RSnum,a,b,ade);
            % xw=rational_quadrature(ord_RSnum,a,b,ade);
        else
            % rational type (assuming ord_RSnum=ord_RSden)
            xw=rational_quadrature(ord_RSnum,a,b,ade);
        end

        t=xw(:,1); w=xw(:,2); tL=t-a;

        % determination of the line formula nodes and weights

        p1_tL=polyval(p1,tL); q1_tL=polyval(q1,tL);
        p2_tL=polyval(p2,tL); q2_tL=polyval(q2,tL);

        % A. nodes.
        XL=p1_tL./q1_tL; YL=p2_tL./q2_tL;

        % B. weights.
        p1L=length(p1);
        if p1L == 1
            p2_diff=0;
        else
            derV=p1L-1:-1:1; p2_diff=p2(1:end-1).*derV;
        end

        q1L=length(q1);
        if q1L == 1
            q2_diff=0;
        else
            derV=q1L-1:-1:1; q2_diff=q2(1:end-1).*derV;
        end

        p2_diff_tL=polyval(p2_diff,tL); q2_diff_tL=polyval(q2_diff,tL);

        y1_term=(p2_diff_tL.*q2_tL-p2_tL.*q2_diff_tL)./q2_tL.^2;
        WL=w.*y1_term;

        % C. Assembly the rule.
        xyw=[xyw; XL YL WL];

    end

end









%====================
% intcvand
%====================

function intV = intcvand(n,pts,bbox)

%-----------------------------------
% OBJECT:
%-----------------------------------
% This subroutine computes by recurrence an x-primitive
% Chebyshev-Vandermonde matrix on a 2d arbitrarily located mesh, in a
% total-degree product Chebyshev basis of a given rectangle.
%------------------------------------
% INPUT:
%-----------------------------------
% n: polynomial degree
% pts: 2-column array of mesh point coordinates
% bbox: [bbox(1),bbox(2)] x [bbox(3),bbox(4)] bounding box for the mesh
%-----------------------------------
% OUTPUT:
%-----------------------------------
% intV: x-primitive of the Chebyshev-Vandermonde matrix
%-----------------------------------
% DATES:
%-----------------------------------
% built: March 2019
% check: November 12, 2021
%-----------------------------------

% default rectangle containing the mesh if not passed

if isempty(bbox)
    bbox=[min(pts(:,1)) max(pts(:,1)) min(pts(:,2)) max(pts(:,2))];
end

a=bbox(1);b=bbox(2);c=bbox(3);d=bbox(4);

Vx=chebpolys(n+1,(2*pts(:,1)-b-a)/(b-a));
Vy=chebpolys(n,(2*pts(:,2)-d-c)/(d-c));

k=0;
for i=0:n
    for j=0:n-i
        k=k+1;

        if i==0
            intV(:,k)=pts(:,1).*Vy(:,j+1);
        end

        if i==1
            intV(:,k)=0.25*(b-a)*(((2*pts(:,1)-b-a)/(b-a)).^2).*Vy(:,j+1);
        end

        if i>1
            intV(:,k)=0.25*(b-a)*(Vx(:,i+2).*Vy(:,j+1)/(i+1)-...
                Vx(:,i).*Vy(:,j+1)/(i-1));
        end

    end
end





%===========================
% rational_quadrature
%===========================

function xw=rational_quadrature(ord_RS,a,b,ade)

%---------------------------
% OBJECT:
%---------------------------
% This routine computes a quadrature rule on an arc of Jordan domains
% whose boundary is given counterclockwise by the rational splines.
%----------------------------
% INPUT:
%----------------------------
% sp_xd_coefs,sp_yd_coefs: coefficients in [a,b] of the denominators
%     defining x(t), y(t); it is supposed that in the interval [a,b] the
%     nurbs is actually a rational function (i.e. quotient of 2
%     polynomials);
% a,b: two successive break points of the splines defining "Sx", "Sy" (i.e.
%    the parametric rational splines defining the domain boundary.
% ade: algebraic degree of exactness of the rule.
%----------------------------
% OUTPUT:
%----------------------------
% xw: matrix N x 2, whose first column are the nodes "t", while the second
%     are the weights "w".
%----------------------------
% DATES:
%---------------------------
% First version: October 31, 2021;
% Checked: November 7, 2021.
%---------------------------
% NOTE:
%--------------------------
% Though we have also used "rfejer routine", we do not suggest its use
% since numerical tests do not seem to provide always good results.
%--------------------------

A=(a+b)/2; B=(b-a)/2;

xw=legendre_rules(10000);


xx=A+B*xw(:,1); ww=xw(:,2); ww=B*ww; xw=[xx ww];









%========================
% rational_quadrature
%========================

function xw=gaussian_quadrature(ord_S,a,b,ade)

%------------------------
% OBJECT:
%------------------------
% This routine computes a quadrature rule on an arc of Jordan domains
% whose boundary is given counterclockwise by the rational splines.
%------------------------
% INPUT:
%------------------------
% sp_xd_coefs,sp_yd_coefs: coefficients in [a,b] of the denominators
%     defining x(t), y(t); it is supposed that in the interval [a,b] the
%     nurbs is actually a rational function (i.e. quotient of 2
%     polynomials);
% a,b: two successive break points of the splines defining "Sx", "Sy" (i.e.
%    the parametric rational splines defining the domain boundary.
% ade: algebraic degree of exactness of the rule.
%------------------------
% OUTPUT:
%------------------------
% xw: matrix N x 2, whose first column are the nodes "t", while the second
%     are the weights "w".
%------------------------
% DATES:
%------------------------
% First version: October 31, 2021;
% Checked: November 7, 2021.
%------------------------
% NOTE:
%------------------------
% Though we have also used "rfejer routine", we do not suggest its use
% since numerical tests do not seem to provide always good results.
%------------------------

A=(a+b)/2; B=(b-a)/2;
k=ceil(((ord_S-1)*(ade+2))/2)+2;
method=0;
switch method
    case 0
        ab=r_jacobi(k,0,0); xw=gauss(k,ab);
    case 1
        xw=legendre_rules(k);
end

t=xw(:,1); w=xw(:,2);

xx=A+B*xw(:,1); ww=xw(:,2); ww=B*ww; xw=[xx ww];











%====================
% r_jacobi
%====================

function ab=r_jacobi(N,a,b)

%R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
%   AB=R_JACOBI(N,A,B) generates the Nx2 array AB of the first
%   N recurrence coefficients for the monic Jacobi polynomials
%   orthogonal on [-1,1] relative to the weight function
%   w(x)=(1-x)^A*(1+x)^B. The call AB=R_JACOBI(N,A) is the same
%   as AB=R_JACOBI(N,A,A) and AB=R_JACOBI(N) the same as
%   AB=R_JACOBI(N,0,0).
%
%   Supplied by Dirk Laurie, 6-22-1998; edited by Walter
%   Gautschi, 4-4-2002.
%   Edited by Walter Leopardi 10-22-2006.

if nargin<2, a=0; end
if nargin<3, b=a; end
if((N<=0)|(a<=-1)|(b<=-1)) error('parameter(s) out of range'), end
nu=(b-a)/(a+b+2);
if a+b+2 > 128
    mu=exp((a+b+1)*log(2)+((gammaln(a+1)+gammaln(b+1))-gammaln(a+b+2)));
else
    mu=2^(a+b+1)*((gamma(a+1)*gamma(b+1))/gamma(a+b+2));
end
if N==1, ab=[nu mu]; return, end
N=N-1; n=1:N; nab=2*n+a+b;
A=[nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
n=2:N; nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
ab=[A' [mu; B1; B']];





%=================
% gauss
%=================

function xw=gauss(N,ab)

%GAUSS Gauss quadrature rule.
%   GAUSS(N,AB) generates the Nx2 array XW of Gauss quadrature
%   nodes and weights for a given weight function W. The nodes,
%   in increasing order, are placed into the first column of XW,
%   and the corresponding weights into the second column. The
%   weight function W is specified by the Nx2 input array AB
%   of recurrence coefficients for the polynomials orthogonal
%   with respect to the weight function W.
%   Supplied by Dirk Laurie and Walter Gautschi.

N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
    J(n,n-1)=sqrt(ab(n,2));
    J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];










