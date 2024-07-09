
function demo_nurbs_cubature_01(example)

%----
% OBJECT:
%----
% Demo illustrating CUBATURE:
% 1. how to define a NURBS on a composite boundary, using arcs of disks,
%    ellipses, segments, polygons and "free NURBS". The routines work on
%    piecewise NURBS of different order.
%    To this purpose see the MATLAB function "define_domain" at the
%    bottom of this demo.
% 2. how to get an algebraic cubature rule on the desired domain, via
%    "cubRS" routine.
% 3. testing the routines over random polynomials.
%---
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
%---
% Dates
%---
% First version: October 31, 2021;
% Checked: July 8, 2024.
%---

warning off;

% ------------------------------ Settings ---

% Algebraic degree of precision of the cubature rule. Suggested values are
% mild, e.g. not larger than "ade=15". Typical "ade" in the applications
% are between "3" and "5".
adeV=2:2:10;

% Number of tests to determine the median cputime. We suggest to set
% "tests=10", while we do not recommend "tests=1" because it often turns
% out to provide unreliable cputimes.
tests=100;

% Cubature extraction method: the variable "extraction_type" sets how
% Tchakaloff rule is extracted.
%      1: lsqnonneg, 2: lawsonhanson, 3: LHDM 4: \
% If not declared or equal to "[]", the default is LHDM, see
% "dCATCH: a numerical package for d-variate near G-optimal Tchakaloff
% regression via fast NNLS"
% M. Dessole, F. Marcuzzi and M. Vianello
% MDPI-Mathematics 8 (2020) - Special Issue "Numerical Methods"
extraction_type=3;

% Nbox is a technical parameter for indomain routine; in doubt set 100.

Nbox=100;

% In indomain routine one must be certain that points are inside the
% domain. In our cubature needs we just wants some points in the domain.
% Consequently, if we are not fully sure that a point is in the domain, due
% to geometrical issues, we do not take it into account.

safe_mode=1;

% The variable "domain_type" set the domain taken into account. 
if nargin < 1
  example=3;
end









% ----------------------------- Main code ---

%---
% 1. Make NURBS structure
%----

[geometry_NURBS,geometry_NURBSH,domain_str]=define_domain(example);

REV=[]; logerrV=[];
for jj=1:length(adeV)

    ade=adeV(jj);

    %---
    % 2. Compute algebraic cubature rule of the domain (several tests!).
    %---
    xyw = cubRS(ade,geometry_NURBS,geometry_NURBSH,extraction_type,Nbox,...
        safe_mode);

    %---
    % 3. Compute algebraic cubature rule of the boundary.
    %---

    [inside,in_doubts,sp_xnV,sp_xdV,sp_ynV,sp_ydV,boxVx,...
        singular_points_X,singular_points_Y]=inRS([],geometry_NURBS,...
        Nbox,safe_mode);
    xywL = lineint(ade+1,sp_xnV,sp_xdV,sp_ynV,sp_ydV);

    [insideH,in_doubtsH,sp_xnVH,sp_xdVH,sp_ynVH,sp_ydVH,boxVxH,...
        singular_points_XH,singular_points_YH]=inRS([],geometry_NURBSH,...
        Nbox,safe_mode);
    xywLH = lineint(ade+1,sp_xnVH,sp_xdVH,sp_ynVH,sp_ydVH);

    %---
    % 4. Making tests.
    %---

    for k=1:tests
        a=rand(1); b=rand(1); c=rand(1);
        f=@(x,y) (a+b*x+c*y).^ade;

        xx=xyw(:,1); yy=xyw(:,2); ww=xyw(:,3);
        I0=ww'*feval(f,xx,yy);

        xxL=xywL(:,1); yyL=xywL(:,2); wwL=xywL(:,3);
        g=@(x,y) (a+b*x+c*y).^(ade+1)/(b*(ade+1));
        I1E=wwL'*feval(g,xxL,yyL);

        xxLH=xywLH(:,1); yyLH=xywLH(:,2); wwLH=xywLH(:,3);
        I1H=wwLH'*feval(g,xxLH,yyLH);

        I1=I1E-I1H;

        RE(k,1)=abs((I0-I1)./I1);
    end

    %----
    % 5. Statistics.
    %----
    fprintf('\n \n \t ADE   : %2.0f',ade);
    fprintf('\n \t RE max: %1.3e',max(RE));
    kpos=find(RE > 0);
    logerr=10^(sum(log10(RE(kpos)))/length(kpos));
    fprintf('\n \t RE log: %2.3e',logerr);

    cell_tick{jj}=num2str(ade);

    REV=[REV RE];
    logerrV=[logerrV; logerr];
end

%---
% 6. Plots.
%---

h=figure(1);
f1=ishandle(h)&&strcmp(get(h,'type'),'figure'); if f1,clf(1);end
figure(1)
axis equal;
for jj=1:length(adeV)
    n=adeV(jj);
    reL=REV(:,jj); log_reL=logerrV(jj);
    plot_errors(jj,n,reL,log_reL);
    hold on;
end
ax=gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
xlim([adeV(1)-1,adeV(end)+1]);
xticks(adeV)
hold off
fprintf('\n \n');

h=figure(2);
f2=ishandle(h)&&strcmp(get(h,'type'),'figure'); if f2,clf(2);end
figure(2)

axis equal;
plotNURBSPL(geometry_NURBS);
hold on;
plotNURBSPL(geometry_NURBSH);
plot(xx,yy,'mo','MarkerEdgeColor','k',...
                       'MarkerFaceColor','m',...
                       'MarkerSize',4)      
hold off;

%--------------------------------------------------------------------------
% ATTACHED ROUTINES
%--------------------------------------------------------------------------















%==========================================================================
% lineint
%==========================================================================

function xyw = lineint(ade,sp_xnV,sp_xdV,sp_ynV,sp_ydV)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine computes weights and nodes for the line integral on Jordan
% domains whose boundary is given counterclockwise by the rational splines.
% The formula is exact on bivariate polynomials up to algebraic degree of
% precision "ade".
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% xyw: 3-column array of nodes coords (cols 1-2) and weights (col 3)
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% built: March 2019
% check: November 12, 2021
%--------------------------------------------------------------------------


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

        % .............. determining line-integral formula  ...............

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









%==========================================================================
% intcvand
%==========================================================================

function intV = intcvand(n,pts,bbox)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This subroutine computes by recurrence an x-primitive
% Chebyshev-Vandermonde matrix on a 2d arbitrarily located mesh, in a
% total-degree product Chebyshev basis of a given rectangle.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% n: polynomial degree
% pts: 2-column array of mesh point coordinates
% bbox: [bbox(1),bbox(2)] x [bbox(3),bbox(4)] bounding box for the mesh
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% intV: x-primitive of the Chebyshev-Vandermonde matrix
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% built: March 2019
% check: November 12, 2021
%--------------------------------------------------------------------------

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





%==========================================================================
% rational_quadrature
%==========================================================================

function xw=rational_quadrature(ord_RS,a,b,ade)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine computes a quadrature rule on an arc of Jordan domains
% whose boundary is given counterclockwise by the rational splines.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% sp_xd_coefs,sp_yd_coefs: coefficients in [a,b] of the denominators
%     defining x(t), y(t); it is supposed that in the interval [a,b] the
%     nurbs is actually a rational function (i.e. quotient of 2
%     polynomials);
% a,b: two successive break points of the splines defining "Sx", "Sy" (i.e.
%    the parametric rational splines defining the domain boundary.
% ade: algebraic degree of exactness of the rule.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% xw: matrix N x 2, whose first column are the nodes "t", while the second
%     are the weights "w".
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Checked: November 7, 2021.
%--------------------------------------------------------------------------
% NOTE:
%--------------------------------------------------------------------------
% Though we have also used "rfejer routine", we do not suggest its use
% since numerical tests do not seem to provide always good results.
%--------------------------------------------------------------------------

A=(a+b)/2; B=(b-a)/2;

xw=legendre_rules(10000);


xx=A+B*xw(:,1); ww=xw(:,2); ww=B*ww; xw=[xx ww];









%==========================================================================
% rational_quadrature
%==========================================================================

function xw=gaussian_quadrature(ord_S,a,b,ade)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine computes a quadrature rule on an arc of Jordan domains
% whose boundary is given counterclockwise by the rational splines.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% sp_xd_coefs,sp_yd_coefs: coefficients in [a,b] of the denominators
%     defining x(t), y(t); it is supposed that in the interval [a,b] the
%     nurbs is actually a rational function (i.e. quotient of 2
%     polynomials);
% a,b: two successive break points of the splines defining "Sx", "Sy" (i.e.
%    the parametric rational splines defining the domain boundary.
% ade: algebraic degree of exactness of the rule.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% xw: matrix N x 2, whose first column are the nodes "t", while the second
%     are the weights "w".
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Checked: November 7, 2021.
%--------------------------------------------------------------------------
% NOTE:
%--------------------------------------------------------------------------
% Though we have also used "rfejer routine", we do not suggest its use
% since numerical tests do not seem to provide always good results.
%--------------------------------------------------------------------------

A=(a+b)/2; B=(b-a)/2;
k=ceil(((ord_S-1)*(ade+2))/2)+2;
method=0;
        xw=legendre_rules(k);


t=xw(:,1); w=xw(:,2);

xx=A+B*xw(:,1); ww=xw(:,2); ww=B*ww; xw=[xx ww];


%================================================================
% plot_errors
%================================================================

function plot_errors(ii,n,reV,log_re)

if ii <= 5
    switch ii
        case 1
            plotstr='m+'; linestr='m-';
        case 2
            plotstr='g+'; linestr='g-';
        case 3
            plotstr='r+'; linestr='r-';
        case 4
            plotstr='b+'; linestr='b-';
        case 5
            plotstr='c+'; linestr='b-';
        otherwise
            plotstr='k.'; linestr='w.';
    end
    %     semilogy(1:number_experiments,reV,plotstr); hold on;
    %     semilogy(1:number_experiments,log_reV(ii)*ones(1,number_experiments),...
    %         linestr,'Linewidth',3);
    semilogy(n*ones(size(reV)),reV,plotstr,'LineWidth',2); hold on;
    semilogy(n,log_re, 'ko','MarkerSize',30,'MarkerEdgeColor','k',...
        'LineWidth',2);
else
    semilogy(n*ones(size(reV)),reV,'color',rand(1,3)); hold on;
    semilogy(n,log_re, 'ko','MarkerSize',29,'MarkerEdgeColor','k');
end

















