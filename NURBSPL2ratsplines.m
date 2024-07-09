

function [SxN,SxD,SyN,SyD]=NURBSPL2ratsplines(structure_RS,spline_type)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine computes the boundary of a certain domain (whose frontier is
% a piecewise rational polynomial) and describes it in piecewise form (in t
% the so called "pp"-form or "B-" form).
% Such a boundary is typically given by piecewise NURBS, composite Bezier
% curves or splines defined by "structure_RS".
%
% Suppose that the boundary is defined parametrically as
%                x=x(t), y=y(t),   t in [alpha,beta]
%
% Next assume that the interval "[alpha,beta]" is subdivided in blocks
% [t(k),t(k+1)] for some k=1,...,L-1 where x(t), y(t) are rational splines.
%
% In each of those intervals [t(k),t(k+1)], we have that
%
%     x=x(t)=(SxN(k))(t)/(SxD(k))(t), y=y(t)=(SyN(k))(t)/(SyD(k))(t)
%
% where SxN(k), SxD(k), SyN(k), SyD(k) are splines of a suitable order,
% determined in view of the NURBS/splines boundary.
%
% The vectors of spline structures SxN, SxD, SyN, SyD take into
% account the functions necessary to define "x(t)", "y(t)" in each block.
% In the general case, the vectors of splines "SxN", "SxD", "SyN", "SyD".
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% structure_RS: array whose k-th element is a
%     "structure_RS array" containing the informations of the k-th
%    NURBS/Composite Bezier Curve/spline.
%
%  * structure_RS.P: it is a N x 2 matrix of coordinates of control
%        points. The k-th control point is described as the k-th row of P.
%
%  * structure_RS.w: vector of NURBS weights; it is a row vector 1 x N,
%        where "N" is the number of control points; for Composite Bezier
%        Curve/spline we set "structure_RS.w=[]".
%
%  * structure_RS.knots: vector of knots;
%     * if the curve is a NURB, following reference [1], it is
%     a row vector 1 x (N+L) of the form
%                   [k0*ones(1,L) u k1*ones(1,L)],
%     where "L" is the order of the NURBS and "u" is a vector of non
%     decreasing values, of length "N-L" where "N" is the number of control
%     points. Observe that
%
%       max(structure_RS{k}.knots) <= min(structure_RS{k+1}.knots).
%     * if the curve is a Composite Bezier Curve or a parametric spline,
%     then it is an increasing sequence of check points that discretize the
%     interval "[a,b]" of parameters "t" defining "x=x(t)", "y=y(t)".
%
%  * structure_RS.order: NURBS, Composite Bezier curve or spline order,
%     i.e. degree plus 1.
%
%  The variable "structure_RS" can be easily defined for NURBS, Composite
%  Bezier Curves, splines (see the pertinent demo).
%
% spline_type: parameter that determines how the splines SxN, SxD, SyN, SyD
%    are described (e.g. set it as "B-" if you prefer B-splines, set it as
%    "pp" if you prefer the pp-form (piecewise polynomials).
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% SxN,SxD,SyN,SyD: vectors of spline structure_RSs describing the
%     boundary, as written above, in "spline_type" form.
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Checked: November 14, 2021.
%--------------------------------------------------------------------------

% Troubleshooting.
if nargin < 2, spline_type='pp'; end

if strcmp(spline_type,'PP-') | strcmp(spline_type,'PP') ...
        | strcmp(spline_type,'pp-'),
    spline_type='pp'; 
end



% Analyse the iteratively each element of the structure.
for k=1:length(structure_RS)

    structure_RSL=structure_RS(k); structure_RSL_type=structure_RSL.type;
    
    % NURBS case.
    if strcmp(structure_RSL_type,'NURBS')
        [SxNL,SxDL,SyNL,SyDL]=NURBS2ratsplines_sub(structure_RSL,...
            spline_type);
    end

    % spline type.
    if strcmp(structure_RSL_type,'spline')
        SPLtypestring='not-a-knot';
        [SxNL,SxDL,SyNL,SyDL]=SPLINEpp(structure_RSL,SPLtypestring,...
            spline_type);
    end

    % Composite Bezier Curve type.
    if strcmp(structure_RSL_type,'bezier')
        [SxNL,SxDL,SyNL,SyDL]=BEZIERpp(structure_RSL,...
            spline_type);
    end

    % Add analysed elements.
    SxN(k)=SxNL; SxD(k)=SxDL; SyN(k)=SyNL; SyD(k)=SyDL;

end










function [SxN,SxD,SyN,SyD]=NURBS2ratsplines_sub(structure_RS,spline_type)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine computes the boundary of a certain domain (whose frontier is
% a piecewise rational polynomial) and describes it in piecewise form (in t
% the so called "pp"-form or "B-" form).
%
% Such a boundary is typically given by piecewise NURBS, defined by
% "structure_RS".
%
% Suppose that the boundary is defined parametrically as
%                x=x(t), y=y(t),   t in [alpha,beta]
% Next suppose that the interval "[alpha,beta]" is subdivided in blocks
% [t(k),t(k+1)] for some k=1,...,L-1 where x(t), y(t) are NURBS.
%
% In each of those intervals [t(k),t(k+1)], we have that
%
%     x=x(t)=(SxN(k))(t)/(SxD(k))(t), y=y(t)=(SyN(k))(t)/(SyD(k))(t)
%
% where SxN(k), SxD(k), SyN(k), SyD(k) are splines of a suitable order,
% determined in view of the NURBS boundary.
%
% The vectors of spline structures "SxN", "SxD", "SyN", "SyD" take into
% account the functions necessary to define x(t), y(t) in each block.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% structure_RS: "structure_RS" element containing the informations of the
%    generic k-th NURBS component.
%
%  * structure_RS.P: it is a N x 2 matrix of coordinates of control
%        points. The k-th control point is described as the k-th row of P.
%
%  * structure_RS.w: vector of NURBS weights; it is a row vector 1 x N,
%        where "N" is the number of control points;
%
%  * structure_RS.knots: vector of knots;
%     * if the curve is a NURB, following reference [1], it is
%     a row vector 1 x (N+L) of the form
%                   [k0*ones(1,L) u k1*ones(1,L)],
%     where "L" is the order of the NURBS and "u" is a vector of non
%     decreasing values, of length "N-L" where "N" is the number of control
%     points. Observe that
%
%       max(structure_RS{k}.knots) <= min(structure_RS{k+1}.knots).
%
%  * structure_RS.order: NURBS order (i.e. "NURBS degree +1").
%
%  The variable "structure_RS" can be easily defined for NURBS (see the
%  pertinent demo).
%
% spline_type: parameter that determines how the splines SxN, SxD, SyN, SyD
%    are described (e.g. set it as "B-" if you prefer B-splines, set it as
%    "PP" if you prefer the pp-form (piecewise polynomials).
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% SxN,SxD,SyN,SyD: vectors of spline structure_RSs describing the
%     component of the boundary relatively to "structure_RS".
%--------------------------------------------------------------------------
% Routines.
%--------------------------------------------------------------------------
% 1. compute_NURBS_boundary_sub (attached, as well as subroutines).
%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2007-2021 Alvise Sommariva, Marco Vianello.
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
%% Author:  Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Dates: October 31, 2021 -
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Checked: November 14, 2021.
%--------------------------------------------------------------------------

L=length(structure_RS);
SxNV=[]; SxDV=[]; SyNV=[]; SyDV=[];

% Piecewise splines computation.
for k=1:L
    structure_RSL=structure_RS(k);
    
    PL=structure_RSL.P;
    knotsL=structure_RSL.knots;
    wL=structure_RSL.w;
    orderL=structure_RSL.order; 
    
    [SxNL,SxDL,SyNL,SyDL]=compute_NURBS_boundary_sub(PL,knotsL,wL,...
        orderL,spline_type);
    
    SxN(k)=SxNL; SxD(k)=SxDL; SyN(k)=SyNL; SyD(k)=SyDL;
end









function [SxN,SxD,SyN,SyD]=compute_NURBS_boundary_sub(P,knots,w,order,...
    spline_type)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Given the 2D NURBS control points that are stored in "P", the knots
% "knots", the weights "w" and the order "order" we describe the NURB in
% "piecewise" form in the type described by "spline_type".
% The result is stored in "SxN", "SxD", "SyN", "SyD", so that the boundary
% is defined parametrically as (x(t),y(t)), t in [a,b] in terms of splines
% of "spline_type" form as
%       x(t)=SxN(t)/SxD(t), y(t)=SyN(t)/SyD(t)     t in [a,b]
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% P: it is a N x 2 matrix of coordinates of control points. The k-th
%    control point is described as the k-th row of P.
%
% knots: vector of knots; following reference [1] it is a row vector
%     1 x (N+L) of the form
%                  [k0*ones(1,L) u k1*ones(1,L)],
%     where "L" is the  order of the NURBS and "u" is a vector of non
%      decreasing values in
%     (k0,k1), of length "N-L" where "N" is the number of control points.
%
% w: vector of NURBS weights; it is a row vector 1 x N, where "N" is the
%    number of control points.
%
% order: NURBS order (i.e. "NURBS degree +1").
%
% spline_type: it is a string, not mandatory variable,
%    * "B-" if one needs the splines "SxN", "SxD", "SSyN", "SyD" in "B-"
%       form
%    * "pp" if one needs the splines "SxN", "SxD", "SSyN", "SyD" in "PP"
%       form
%    (default) 'pp' form is variable is not declared, 'B-' if the variable
%       is declared but s.t. "upper(spline_type)" is not equal to 'PP'.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% SxN: numerator of x(t), as spline as declared in "spline_type";
% SxD: denominator of x(t), as spline as declared in "spline_type";
% SyN: numerator of y(t), as spline as declared in "spline_type";
% SyD: denominator of y(t), as spline as declared in "spline_type";
%--------------------------------------------------------------------------
% REFERENCE:
%--------------------------------------------------------------------------
% Les Piegl, Wayne Tiller,
% The NURBS Book,
% Second edition (Springer Verlag, 1995 and 1997).
%
% Note: at page 117, in U, we set "m=N-L" where "L" is the order of the
%      NURBS, and "N" the number of control points.
%--------------------------------------------------------------------------
%% Copyright (C) 2007-2021 Alvise Sommariva, Marco Vianello.
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
%% Author:  Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Date: October 31, 2021.
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Checked: November 14, 2021.
%--------------------------------------------------------------------------

% .......................... TROUBLESHOOTING  .............................
if nargin < 4, order=length(knots)-length(w); end
if nargin < 5, spline_type='pp'; end

if strcmp(spline_type,'PP-') | strcmp(spline_type,'pp-')
    spline_type='pp';
end

% ............................ main code  .................................
N=size(P,1); L=length(knots);

% check weights
if not(length(w) == N), error('wrong dimension of the weights'); end
if min(w) < 0, error('weights must be larger than 0'); end

% check order
if not(order == L-length(w)), error('unadmissible order'); end

% BUILDING "SxN" IN "B-" FORM.
SxN.form='B-';
SxN.knots=knots;
SxN.coefs=w.*(P(:,1))';
SxN.order=order;
SxN.number=length(SxN.knots)-SxN.order;
SxN.dim=1;

% BUILDING "SxD" IN "B-" FORM.
SxD=SxN; SxD.coefs=w;

% BUILDING "SyN" IN "B-" FORM.
SyN=SxN; SyN.coefs=w.*(P(:,2))';

% BUILDING "SyD" IN "B-" FORM.
SyD=SxD; % notice that for NURBS we have the same denominators.

% SPLINE CONVERSIONS, IF REQUIRED.
if strcmp(spline_type,'pp')
    SxN=fn2fm(SxN,'pp'); SxD=fn2fm(SxD,'pp');
    SyN=fn2fm(SyN,'pp'); SyD=SxD;
end









function [SxN,SxD,SyN,SyD]=SPLINEpp(structure_BOUNDARY,SPLtypestring,...
            spline_type)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This subroutine computes the parametric spline arrays "SxN", "SxD",
% "SyN", "SyD" that are used to describe the closed parametric curves
% defining the boundary "dD" of the required domain "D", i.e.
%            x=x(t)=(SxN(k))(t)/,(SxD(k))(t), t in [t(k),t(k+1)]
%            y=y(t)=(SyN(k))(t)/,(SyD(k))(t), t in [t(k),t(k+1)]
% where [t(k),t(k+1)] is a subinterval of [a,b], defining its k-th block,
% when "structure_BOUNDARY" defines a parametric spline boundary.
%
% It is assumed that both "x=x(t)", "y=y(t)" have the same sub-blocks 
% [t(k),t(k+1)].
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% structure_BOUNDARY: structure of the boundary of the domain, it contains
%    all the parameters that determine the spline.
%
%         structure_BOUNDARY.P: spline control points,
%         structure_BOUNDARY.order: spline order,
%         structure_BOUNDARY.knots: spline knots.
%
% SPLtypestring: string with the spline type i.e.
%             'complete'   : match endslopes (as given in VALCONDS, with
%                     default as under *default*).
%             'not-a-knot' : make spline C^3 across first and last interior
%                     break (ignoring VALCONDS if given).
%             'periodic'   : match first and second derivatives at first
%                     data point with those at last data point (ignoring
%                     VALCONDS if given).
%             'second'     : match end second derivatives (as given in
%                    VALCONDS, with default [0 0], i.e., as in variational).
%             'variational': set end second derivatives equal to zero
%                     (ignoring VALCONDS if given).
%
% spline_type: parameter that determines how the splines SxN, SxD, SyN, SyD
%    are described (e.g. set it as "B-" if you prefer B-splines, set it as
%    "PP" if you prefer the pp-form (piecewise polynomials).
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% SxN: spline "x=x(t)" data in "pp" form,
% SyN: spline "y=y(t)" data in "pp" form,
% SxD, SyD: if the polynomial spline is seen more generally in rational pp
%      form, we define the relative denominators (equal to "1").
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% Written: November 7, 2021.
% Modified: November 14, 2021.
%--------------------------------------------------------------------------

if nargin < 2, SPLtypestring='not-a-knot'; end
if nargin < 3, SPLtypestring='pp'; end

P=structure_BOUNDARY.P;
s=structure_BOUNDARY.knots;
SPLINE_order=structure_BOUNDARY.order;
x=P(:,1); y=P(:,2);

switch SPLINE_order
    case 4
        % Cubic splines, using taylored routines.
        SxN=csape(s,x,SPLtypestring); SyN=csape(s,y,SPLtypestring);
    otherwise % other routines
        SxN=spapi(SPLINE_order,s,x); SyN=spapi(SPLINE_order,s,y);
end

% conversion in pp-form
if strcmp(SxN.form,spline_type) == 0
    SxN=fn2fm(SxN,spline_type);
    SyN=fn2fm(SyN,spline_type);
end

if nargout > 3
    % If the spline is seen as a rational spline, we define the
    % denominators equal to "1".
    SxD=SxN; SxD.order=1; SxD_coefs=SxD.coefs;
    SxD_coefs=ones(size(SxD_coefs,1),1); SxD.coefs=SxD_coefs;

    SyD=SxD;
end









function [SxN,SxD,SyN,SyD]=BEZIERpp(structure_BOUNDARY,spline_type)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Given the structure array "structure_BOUNDARY" of a closed Bezier
% composite curve that describes the boundary of the domain, we
% redefine the curve in piecewise form via "SxN", "SxD", "SyN", "SyD".
% In other words, if the Bezier composite curve is
%                   x=x(t), y=y(t), t in [a,b]
% then in the k-th interval [t(k),t(k+1)] we have
%         x(t)=(SxN(k))(t)/(SxD(k))(t)=(SxN(k))(t),
%         y(t)=(SyN(k))(t)/(SyD(k))(t)=(SyN(k))(t).
%
% It is assumed that both "x=x(t)", "y=y(t)" have the same sub-blocks  
% [t(k),t(k+1)].
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% structure_BOUNDARY: structure of the boundary of the domain, it contains
%    all the parameters that determine the spline.
%    More exactly:
%         structure_BOUNDARY.P: Composite Bezier control points,
%         structure_BOUNDARY.order: Composite Bezier order,
%         structure_BOUNDARY.knots: Composite Bezier knots.
%
% spline_type: parameter that determines how the splines SxN, SxD, SyN, SyD
%    are described (e.g. set it as "B-" if you prefer B-splines, set it as
%    "PP" if you prefer the pp-form (piecewise polynomials).
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% SxN: spline "x=x(t)" data in "pp" form,
% SyN: spline "y=y(t)" data in "pp" form,
% SxD, SyD: if the polynomial spline is seen more generally in rational pp
%      form, we define the relative denominators (equal to "1").
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% Written: November 13, 2021.
% Modified: November 14, 2021.
%--------------------------------------------------------------------------

if nargin < 2, spline_type='pp'; end

P=structure_BOUNDARY.P;
order=structure_BOUNDARY.order;
s=structure_BOUNDARY.knots;

% ........................ troubleshooting  ...............................
deg=order-1;
LP=size(P,1);

if not(rem(LP-1,deg) == 0)
    warning('Wrong cardinality for Bezier composite structures: deg: %3.0f #P: %6.0f', ...
    deg,LP);
    Lblocks=floor((LP-1)/deg);
else
    Lblocks=(LP-1)/deg;
end

% ........................... main code  ..................................

% Bezier curve: number of sub-blocks.
sB=s(1:deg:length(s));

SxN.form=spline_type; SxN.breaks=sB; SxN.coefs=[];
SxN.pieces=Lblocks; SxN.order=order; SxN.dim=1;

SxD.form=spline_type; SxD.breaks=sB; SxD.coefs=[];
SxD.pieces=Lblocks; SxD.order=order; SxD.dim=1;

SyN.form=spline_type; SyN.breaks=sB; SyN.coefs=[];
SyN.pieces=Lblocks; SyN.order=order; SyN.dim=1;

SyD.form=spline_type; SyD.breaks=sB; SyD.coefs=[];
SyD.pieces=Lblocks; SyD.order=order; SyD.dim=1;

M=bezier2powermatrix(order-1);

for k=1:Lblocks

    PL=P(deg*(k-1)+1:deg*k+1,:);

    % Bernstein polynomials are typically used in [0,1]. Since the block
    % may be defined in a general interval [sB(k),sB(k+1)], to convert
    % correctly in "pp" form, it is easy to see that we must scale the
    % coefficients.

    LL=sB(k+1)-sB(k); scal=LL.^(-deg:0);

    coefs=M*PL;
    coefsX=coefs(:,1);  coefsX=fliplr(coefsX').*scal;
    coefsY=coefs(:,2);  coefsY=fliplr(coefsY').*scal;

    SxN.coefs=[SxN.coefs; coefsX];
    SyN.coefs=[SyN.coefs; coefsY];

    coefsXD=zeros(size(coefsX)); coefsXD(end)=1;
    SxD.coefs=[SxD.coefs; coefsXD];
    SyD.coefs=[SyD.coefs; coefsXD];

end









function M=bezier2powermatrix(p)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Given an order "p" it provides a matrix that converts from Bernstein form
% to monomial power form.
% Example:
% [(1-t)^3 3*t*(1-t)^2 3*t^2*(1-t) t^3]=[1 t t^2 t^3]*M.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% p: Composite Bezier curve order of the block in analysis.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% M: Conversion matrix from Bernstein form to monomial power form.
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Checked: November 14, 2021.
%--------------------------------------------------------------------------


switch p
    case 1
        M=[1 0; -1 1];

    otherwise
        
        for i=0:p-1
            ii=i+1; pp=p+1;
            for j=ii+1:p
                M(ii,ii)=1; M(pp,pp)=1;
                if rem(p,2) == 0
                    M(pp,1)=-1;
                else
                    M(pp,1)=1;
                end
            end
        end

        signL=-1;
        for i=1:p
            ii=i+1;
            M(ii,ii)=nchoosek(p,i);
            M(ii,1)=signL*M(ii,ii);
            M(p+1,p+1-i)=M(ii,1);
            signL=-signL;
        end

        k1=(p+1)/2; pk=p-1;

        for k=1:k1-1
            signL=-1;
            for j=k+1:pk
                M(j+1,k+1)=signL*nchoosek(p,k)*nchoosek(p-k,j-k);
                M(pk+1,p-j+1)=M(j+1,k+1);
                signL=-signL;
            end
            pk=pk-1;
        end

end



