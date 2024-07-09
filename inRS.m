function [in,on,sp_xnV,sp_xdV,sp_ynV,sp_ydV,boxVxy,turn_pts_X, ...
    turn_pts_Y,pgon]=inRS(pts,structure_RS,Nbox,safe_mode,varargin)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine, once it is defined a bivariate boundary by the structured
% array "structure_RS", tests if the points "pts" are strictly inside
% (or not strictly inside) the domain whose boundary is defined
% parametrically by piecewise rational functions.
%
% Due to the numerical approximation of the boundary, some points that are
% marked as interior are intended as
%      "interior to the numerical approximation of the boundary".
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% pts: matrix N x 2 of bivariate points to be tested;
%
% structure_RS: array whose k-th element is a "structure" containing the
%    informations of the k-th portion of the boundary.
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
%        [1] it is a row vector 1 x (N+L) of the form
%                   [k0*ones(1,L) u k1*ones(1,L)],
%        where "L" is the order of the NURBS and "u" is a vector of non
%        decreasing values, of length "N-L" where "N" is the number of
%        control points. Observe that
%
%       max(structure_RS{k}.knots) <= min(structure_RS{k+1}.knots).
%
%     2. in case the boundary is defined by a splines it represents the
%        break points of the spline
%
%  * structure_RS.order: NURBS or spline order (i.e. "NURBS degree +1").
%
%  The variable "structure_RS" can be easily defined (see "demo").
%
% Nbox: not mandatory parameter, useful for indomain (default: 100). It is
%     a technical variable that if correctly set allows the indomain
%     routine to be faster.
%
% safe_mode: 0: problematic points are not analysed, 1: further
%       investigations are done. Default: "safe_mode=1"
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% in: if "in(k)=1" then the point described by the k-th row of "P" is
%     inside the domain.
%
% on: vector, of the same dimension of "in" where if the i-th component is
%      * Nan    : point is "safely" outside or inside the domain.
%      * finite : point is "close" to the boundary and an approx. of the
%                 distance "d(point,boundary)" is given.
%      * -Inf   : problems with the determination inside/on/outside.
%
% sp_xnV, sp_xdV: vectors of splines in "pp" form, whose k component
%   defines in a suitable subinterval [t(k),t(k+1)] of the domain,
%   say [alpha,beta], respectively the numerator and the denominator
%   of x=x(t), i.e. with some abuse of notation:
%
%          x(t)=Sx(t)=(sp_xnV(k))(t)/(sp_xdV(k))(t), t in [t(k),t(k+1)]
%
% sp_ynV, sp_ydV: vectors of splines in "pp" form, whose k component
%   defines in a suitable subinterval [t(k),t(k+1)] of the domain of the
%   NURBS, say [alpha,beta], respectively the numerator and the denominator
%   of y=y(t), i.e. with some abuse of notation:
%
%          y(t)=Sy(t)=(sp_ynV(k))(t)/(sp_ydV(k))(t), t in [t(k),t(k+1)]
%
% boxVxy: S x 8 matrix, whose i-th row describes the properties of the i-th
%       box of the boundary, i.e.
%
%      *           [xm,xM,ym,yM,t0,t1,ispline,iblock]
%
%      where [xm,xM,ym,yM] defines the boundary of the box, i.e. the
%      vertices
%
%                   X=[xm xM xM xm]   Y=[ym ym yM yM]
%
%      such that the box consists of the rectangle with vertices
%      (X(k),Y(k)), k=1,2,3,4.
%
%      * the parametrization of x=x(t), y=y(t) has values of "t" in the
%        interval [t0,t1];
%
%      * "ispline" says what rat. spline index is defined in the box, i.e.
%
%             x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t);
%
%      * "iblock" says what spline block index is defined in the box, i.e.
%
%             x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t)
%
%        with
%
%        "t in [t1,t2]", where "t1=(Sx(ispline)).breaks(iblock)" and
%        "t2=(Sx(ispline)).breaks(iblock+1)".
%
%      This variable is not mandatory.
%
%
% turn_pts_X, turn_pts_Y: points where the first derivatives of "Sx(k)" or
%    "Sy(k)" are null.
%    They are a M x 3 matrix, whose k-th row defines the k-th point. This
%    triple, say (x,y,t) defines the cartesian coordinates of the point P
%    as well as P as "(Sx(k)(t),Sy(k)(t))".
%
%
%
% Note: the variables "sp_xnV", "sp_xdV", "sp_ynV", "sp_ydV" are vectors,
%   whose k-th component is a structured array defining a spline.
%   Now, suppose that "L" is the dimension of these vectors.
%   For any k=1,...,L we have that the breakpoints and order of
%   "sp_xnV(k)", "sp_xdV(k)", "sp_ynV(k)", "sp_ydV(k)" coincide.
%
% pgon: polygon based on samples of the box, that approximates the domain.
%--------------------------------------------------------------------------
% SUBROUTINES (directly called):
%--------------------------------------------------------------------------
% 1. NURBSPL2ratsplines: external;
% 2. boxerRSboundary:  attached;
% 3. indomainRS: attached.
%--------------------------------------------------------------------------
% REFERENCE:
%--------------------------------------------------------------------------
% [1] Les Piegl, Wayne Tiller,
% The NURBS Book,
% Second edition (Springer Verlag, 1995 and 1997).
%
% Note: at page 117, in U, we set "m=N-L" where "L" is the order of the
%      NURBS, and "N" the number of control points.
%--------------------------------------------------------------------------
% NOTE:
%--------------------------------------------------------------------------
% There are several usages of this code. The first one is to determine
% if some points are inside the NURBS domain, in some others one just wants
% to determine a dense set of nodes in the domain and some doubtful points
% can be discarded.
% If "safe_mode=0" doubtful points are rejected for successive purpose,
% while if "safe_mode=1" some more checkings are done for these more
% problematic points.
%--------------------------------------------------------------------------
% Copyrights.
%--------------------------------------------------------------------------
%% Copyright (C) 2007-2022 Alvise Sommariva, Marco Vianello.
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
%% Date: October 31, 2021 -
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Modified/checked: January 05, 2022 (17:38).
% Note: This version uses binary search.
%--------------------------------------------------------------------------

% ................. Troubleshoouting and defaults .........................


if nargin < 3, Nbox=100; end
if nargin < 4, safe_mode=1; end
if nargin >= 5
    sp_xnV=varargin{1}; sp_xdV=varargin{2};
    sp_ynV=varargin{3}; sp_ydV=varargin{4};
    boxVxy=varargin{5};
    turn_pts_X=varargin{6}; turn_pts_Y=varargin{7};
    skip_boxing = 1;
else
    skip_boxing = 0;
end

% .........................  Main code below ..............................

% If "skip_boxing == 0" then analysis of the domain was not done before.
% We compute the spline arrays "sp_xnV", "sp_xdV", "sp_ynV", "sp_ydV" as
% well as the monotone boxes "boxVxy" and turning points "turn_pts_X",
% "turn_pts_Y".

if skip_boxing == 0
    % Compute the vectors of splines defining in some subinterval the
    % rational functions defining the functions x=x(t), y=y(t) that
    % describe parametrically the boundary of the domain.

    [sp_xnV,sp_xdV,sp_ynV,sp_ydV]=NURBSPL2ratsplines(structure_RS,'PP');

    % Compute "monotone boxes" and turning points.
    [boxVxy,turn_pts_X,turn_pts_Y,pgon]=...
        boxerRSboundary(sp_xnV,sp_xdV,sp_ynV,sp_ydV,Nbox);
    if isempty(turn_pts_X) == 0, turn_pts_X=turn_pts_X(:,1:3); end
    if isempty(turn_pts_Y) == 0, turn_pts_Y=turn_pts_Y(:,1:3); end
else
    pgon=[];
end





% .................. First checks: X direction, bottom ....................

% INDOMAIN: from bottom, x direction.
% Once the analysis of the domain is available, we have computed the spline
% arrays "sp_xnV", "sp_xdV", "sp_ynV", "sp_ydV" as well as the monotone
% boxes "boxVxy" and turning points "turn_pts_X", "turn_pts_Y".
% Using them we can perform the core of the "indomain" procedure, i.e.
% apply the subroutine "indomainRS".

if isempty(pts) == 0
    if isempty(turn_pts_X)
        turn_pts_Xarg=[];
    else
        turn_pts_Xarg=turn_pts_X(:,1);
    end

    [in,on]=indomainRS(pts,sp_xnV,sp_xdV,sp_ynV,sp_ydV,boxVxy,...
        turn_pts_Xarg,10^(-12),safe_mode);

else
    in=[]; on=[];
end





% ..................Further checks: Y direction, left  ....................

% INDOMAIN: crossing from left, y-direction (in case of doubts)

if safe_mode & (isempty(pts) == 0)

    doubts_index=find(on == -Inf);
    pts_doubts=pts(doubts_index,:);

    if isempty(pts_doubts) == 0
        if isempty(turn_pts_Y)
            turn_pts_Yarg=[];
        else
            turn_pts_Yarg=turn_pts_Y(:,2);
        end

        [inD,onD]=indomainRS(pts_doubts(:,[2 1]),sp_ynV,sp_ydV,...
            sp_xnV,sp_xdV,boxVxy(:,[3 4 1 2 5 6 7 8]),turn_pts_Yarg,...
            10^(-12),safe_mode);

        in(doubts_index)=inD; on(doubts_index)=onD;
    end
else
    pts_doubts=[];
end





% ..................Further checks: X direction, top  ...................

% INDOMAIN: crossing from top, x-direction (in case of doubts)

if safe_mode & (isempty(pts_doubts) == 0)
    doubts_index=find(on == -Inf);
    pts_doubts=pts(doubts_index,:);

    if isempty(pts_doubts) == 0
        boxVxyRX=boxVxy;
        boxVxyRX(:,3)=-boxVxy(:,4); boxVxyRX(:,4)=-boxVxy(:,3);
        ptsRX=pts_doubts; ptsRX(:,2)=-ptsRX(:,2);

        for isp=1:length(sp_ynV)
            sp_ynVL=sp_ynV(isp);
            sp_ynVL.coefs=-sp_ynVL.coefs;
            sp_ynV(isp)=sp_ynVL;
        end

        [inD,onD]=indomainRS(ptsRX,sp_xnV,sp_xdV,sp_ynV,sp_ydV,boxVxyRX,...
            turn_pts_Xarg,10^(-12),safe_mode);
        in(doubts_index)=inD; on(doubts_index)=onD;
    end
else
    pts_doubts=[];
end





% ..................Further checks: Y direction, right  ...................

% INDOMAIN: crossing from right, y-direction (in case of doubts)

if safe_mode & (isempty(pts_doubts) == 0)
    doubts_index=find(on == -Inf);
    pts_doubts=pts(doubts_index,:);

    if isempty(pts_doubts) == 0
        boxVxyRY=boxVxy;
        boxVxyRY(:,1)=boxVxy(:,3); boxVxyRY(:,2)=boxVxy(:,4);
        boxVxyRY(:,3)=boxVxy(:,1); boxVxyRY(:,4)=boxVxy(:,2);
        boxVxyRY(:,3)=-boxVxy(:,4); boxVxyRY(:,4)=-boxVxy(:,3);
        for isp=1:length(sp_xnV)
            sp_ynVL=sp_ynV(isp); sp_ynVL.coefs=-sp_ynVL.coefs;
            sp_ynV(isp)=sp_ynVL;
        end
        ptsRX=pts_doubts; ptsRX(:,1)=-ptsRX(:,1);
        [inD,onD]=indomainRS(ptsRX,sp_ynV,sp_ydV,...
            sp_xnV,sp_xdV,boxVxyRY,turn_pts_Yarg,...
            10^(-12),safe_mode);
        in(doubts_index)=inD; on(doubts_index)=onD;
    end
else
    pts_doubts=[];
end





% ..................... Further checks: winding ...........................

% INDOMAIN: applying winding algorithm

if safe_mode & (isempty(pts) == 0)
    doubts_index=find(on == -Inf);
    if length(doubts_index) > 0
        % fprintf(2,'\n \t Additional test: winding number');
        pts_doubts=pts(doubts_index,:);
        inW=windingRS(sp_xnV,sp_xdV,sp_ynV,sp_ydV,pts_doubts);
        in(doubts_index)=inW;
        on(doubts_index)=NaN*ones(size(doubts_index));
    end
end

















%--------------------------------------------------------------------------
% spline_boxer
%--------------------------------------------------------------------------

function [boxVx,turning_points_Sx,turning_points_Sy,pgon]=...
    boxerRSboundary(SxNV,SxDV,SyNV,SyDV,Nsub)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Given a sequence of splines "SxNV", "SxDV", "SyNV", "SyDV" representing
% the boundary of the domain as, i.e. with a little abuse of notation,
%
% *  x(t)=(Sx(k))(t)=(SxNV(k))(t)/(SxDV(k))(t),     t in [t(k),t(k+1)],
% *  y(t)=(Sy(k))(t)=(SYNV(k))(t)/(SYDV(k))(t),     t in [t(k),t(k+1)],
%
% this routine defines the so called "monotone boxes".
%
% The routine covers the boundary, say "dD", with union of rectangles
% "B(i)", i=1,2,..., that we call "monotone boxes".
%
% Each rectangle "B(i)" is such that:
% 1. only one rat.spline defines in each variable the portion dD(i) of the
%    boundary "dD" in "B(i)". Let these rational splines be resp. "Sx(k)",
%    "Sy(k)" and "[t(k),t(k+1)]" their definition interval such that the
%    boundary is "( (Sx(k))(t) , (Sy(k))(t) )" with "t" in "[t0,t1]" subset
%    of "[t(k),t(k+1)]".
% 2. The restrictions of "Sx(k)" and "Sy(k)" to "[t0,t1]" are monotone
%    functions.
% The code will work as follows (see "spline_boxer").
% 1. Looks for maximal intervals "[T0,T1]" for which this happens. To this
%    purposes it detects when first derivatives of "Sx(k)", "Sy(k)" are
%     null and then determines each "[T0,T1]" after some considerations,
%     concerning the variables "turning_points_Sx", "turning_points_Sy"
%     that are outputs of a "simple_boxer" routine.
% 2. Subdivides, for performance purposes, "[T0,T1]" in "Nsub" equispaced
%    subinterintervals. Our experience is that if "Nsub" is choosen too
%     small or too large, the performance will be in general worst.
%
% Let "t(k)" and "t(k+1)" be two subsequent breaks of "Sx(j)" and "Sy(j)".
% Notice that in the code we suppose that "Sx(j).breaks" is equal to
% "Sy(j).breaks".

% A box is a rectangle [x0,x1] x [y0,y1] , where
%    x0=min(Sx(t))   x1=max(Sx(t))  y0=min(Sy(t))   y1=max(Sy(t))
% It collects all this boxes in "boxVx".
%
% Turning points "turning_points_Sx" are points where the derivative of
%     "Sx" in "t(k)" and "t(k+1)" is null.
%--------------------------------------------------------------------------
% INPUT:
%
% SxNV,SxDV,SyNV,SyDV: vectors of splines, representing the boundary of the
% domain as, with a little abuse of notation,
%
% *  x(t)=(Sx(k))(t)=(SxNV(k))(t)/(SxDV(k))(t),     t in [t(k),t(k+1)],
% *  y(t)=(Sy(k))(t)=(SYNV(k))(t)/(SYDV(k))(t),     t in [t(k),t(k+1)],
%
% Nsub:  The routine covers the boundary, say "dD", with union of
%   rectangles "B(i)", that we call "monotone boxes".
%
% Each square "B(i)" is such that:
% 1. only one spline defines in each variable the portion dD(i) of the
%    boundary "dD" in "B(i)". Let these splines be respectively "Sx(k)",
%    "Sy(k)" and "[t(k),t(k+1)]" their definition interval such that the
%    boundary is "( (Sx(k))(t) , (Sy(k))(t) )" with "t" in "[t0,t1]".
% 2. The restrictions of "Sx(k)" and "Sy(k)" to "[t0,t1]" are monotone
%    functions.
%
% The code will work as follows.
%
% 1. It looks for maximal intervals "[T0,T1]" for which this happens.
%    It detects when first derivatives of "Sx(k)" and "Sy(k)" are
%     null and then determines each "[T0,T1]" after some considerations,
%     concerning the variables "turning_points_Sx", "turning_points_Sy"
%     that are outputs of "simple_boxer".
% 2. Subdivides, for performance purposes, "[T0,T1]" in "Nsub" equispaced
%    subinterintervals. The number "210" below is obtained after many
%    trials and in general provides a good performance by the routine. Our
%    experience is that if "Nsub" is choosen too small or too large, the
%    performance will be in general worst.
%
% In case of doubts, do not declare this variable (value "210" is given as
% default)
%--------------------------------------------------------------------------
% OUTPUT:
%
% boxVxy: S x 8 matrix, whose i-th row describes the properties of the i-th
%       box of the boundary, i.e.
%
%      *           [xm,xM,ym,yM,t0,t1,ispline,iblock]
%
%      where [xm,xM,ym,yM] defines the boundary of the box, i.e. the
%      vertices
%
%                   X=[xm xM xM xm]   Y=[ym ym yM yM]
%
%      such that the box consists of the rectangle with vertices
%      (X(k),Y(k)), k=1,2,3,4.
%
%      * the parametrization of x=x(t), y=y(t) has values of "t" in the
%        interval [t0,t1];
%
%      * "ispline" says what rat. spline index is defined in the box, i.e.
%
%             x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t);
%
%      * "iblock" says what spline block index is defined in the box, i.e.
%
%             x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t)
%
%        with
%
%        "t in [t1,t2]", where "t1=(Sx(ispline)).breaks(iblock)" and
%        "t2=(Sx(ispline)).breaks(iblock+1)"
%
%      This variable is not mandatory.
%
% turning_points_Sx: points when first derivatives of "Sx(k)" are null. It
%    is a M x 3 matrix, whose k-th row defines the k-th point. This triple,
%    say (x,y,t) defines the cartesian coordinates of the point P as well
%    as P as "(Sx(k)(t),Sy(k)(t))".
%
% turning_points_Sy: points when first derivatives of "Sy(k)" are null. It
%    is a M x 3 matrix, whose k-th row defines the k-th point. This triple,
%    say (x,y,t) defines the cartesian coordinates of the point P as well
%    as the variable "t" such that "P=(Sx(k)(t),Sy(k)(t))".
%
% pgon: polygon based on samples of the box, that approximates the domain.
%--------------------------------------------------------------------------
% SUBROUTINES (directly called):
% 1. turning_points_analysis.
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
%% Date: October 31, 2021 -
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Checked: December 16, 2021.
%--------------------------------------------------------------------------

if nargin < 3, Nsub=210; end

% troubleshooting

if nargin < 5, Nsub=1; end
if Nsub < 1, Nsub=1; end

L=length(SxNV);

% Initialisations.
boxVx=[]; turning_points_Sx=[]; turning_points_Sy=[];
S1xV=[]; S1yV=[];  S1tV=[]; SxV=[]; SyV=[]; pgon=[];

for j=1:L

    SxN=SxNV(j); SxD=SxDV(j); SyN=SyNV(j); SyD=SyDV(j);

    knotsV=SxN.breaks;
    n_subintervals=length(knotsV)-1;

    S1tV=[S1tV; knotsV(2:end)'];

    % Degrees of the splines defining "x(t)".
    % degL=SxN.order-1;

    % The vector "derV" will be used to compute derivatives.

    % Analysis of the k-th rat. spline subinterval.
    for k=1:n_subintervals

        % .................. retrieving data ..................

        % Extremes of the k-th subinterval.
        a=knotsV(k); b=knotsV(k+1);

        % Spline coefficients defining "Sx" and "Sy" in the k-th
        % subinterval "[a,b]".
        SxN_coeffsL=SxN.coefs(k,:); SxD_coeffsL=SxD.coefs(k,:);
        SyN_coeffsL=SyN.coefs(k,:); SyD_coeffsL=SyD.coefs(k,:);

        SxN_deg=length(SxN_coeffsL)-1;
        if SxN_deg > 0
            derV=SxN_deg:-1:1; SxN1_coeffsL=SxN_coeffsL(1:end-1).*derV;
        else
            SxN1_coeffsL=0;
        end

        SxD_deg=length(SxD_coeffsL)-1;
        if SxD_deg > 0
            derV=SxD_deg:-1:1; SxD1_coeffsL=SxD_coeffsL(1:end-1).*derV;
        else
            SxD1_coeffsL=0;
        end

        SyN_deg=length(SyN_coeffsL)-1;
        if SyN_deg > 0
            derV=SyN_deg:-1:1; SyN1_coeffsL=SyN_coeffsL(1:end-1).*derV;
        else
            SyN1_coeffsL=0;
        end

        SyD_deg=length(SyD_coeffsL)-1;
        if SyD_deg > 0
            derV=SyD_deg:-1:1; SyD1_coeffsL=SyD_coeffsL(1:end-1).*derV;
        else
            SyD1_coeffsL=0;
        end

        % Coefficients of the numerators derivatives. Note that these
        % derivatives are computed to detect zeros of the derivatives of
        % the rational splines "Sx" and "Sy" that in the interval "[a,b]"
        % are rational functions, say
        %            "SxL=SxN(k)/SxD(k)" and "SyL=SyN(k)/SyD(k)".
        % Denominators of the derivatives are equal to  "(SxD)^2" and
        % "(SyD)^2" and do not contribute to the detection of these zeros.
        % Note that they are strictly positive since "SxL", "SyL" are
        % assumed continuous functions.
        % Thus only the numerators of the derivatives of "Sx" and "Sy" are
        % taken into account.
        %
        % The code becomes weird due to different components in the next
        % convolutions, e.g. using non rational splines.

        if length(SxD_coeffsL) == 1 & SxD_coeffsL(end) == 1
            Sx1_coeffsL=SxN1_coeffsL;
        else
            % length of convolution vector
            LL(1)=length(SxN1_coeffsL)+length(SxD_coeffsL)-1;
            LL(2)=length(SxD1_coeffsL)+length(SxN_coeffsL)-1;
            MM=max(LL);

            term1=zeros(1,MM); res=conv(SxN1_coeffsL,SxD_coeffsL);
            term1(MM-LL(1)+1:MM)=res;

            term2=zeros(1,MM); res=conv(SxD1_coeffsL,SxN_coeffsL);
            term2(MM-LL(2)+1:MM)=res;

            Sx1_coeffsL=term1-term2; % numerator derivative "x" direction
        end


        % length of convolution vector
        if length(SyD_coeffsL) == 1 & SyD_coeffsL(end) == 1
            Sy1_coeffsL=SyN1_coeffsL;
        else
            LL(1)=length(SyN1_coeffsL)+length(SyD_coeffsL)-1;
            LL(2)=length(SyN_coeffsL)+length(SyD1_coeffsL)-1;
            MM=max(LL);

            term1=zeros(1,MM); res=conv(SyN1_coeffsL,SyD_coeffsL);
            term1(MM-LL(1)+1:MM)=res;

            term2=zeros(1,MM); res=conv(SyN_coeffsL,SyD1_coeffsL);
            term2(MM-LL(2)+1:MM)=res;

            Sy1_coeffsL=term1-term2;  % numerator derivative "y" direction
        end


        % .................. turning points analysis ..................

        % We determine points "turning_points_xL" and "turning_points_yL",
        % i.e. points in which the polynomials defined by "SxL1_coeffsL"
        % and "SyL1_coeffsL" are null.
        %
        % Notice that locally the spline is stored as a polynomial
        % "q(t)=p(t+a)". Thus for evaluation purposes, we shift "[a,b]"
        % into "[t1,t2]".

        t1=0; t2=b-a;

        turning_points_xL=turning_points_analysis(SxN_coeffsL,...
            SxD_coeffsL,SyN_coeffsL,SyD_coeffsL,Sx1_coeffsL,a,t1,t2);

        turning_points_yL=turning_points_analysis(SxN_coeffsL,...
            SxD_coeffsL,SyN_coeffsL,SyD_coeffsL,Sy1_coeffsL,a,t1,t2);

        turning_points_Sx=[turning_points_Sx; turning_points_xL];
        turning_points_Sy=[turning_points_Sy; turning_points_yL];

        % Computing derivatives at end points, useful to determine
        % difficult junctions.
        tV=[t1 t2];

        Sx_tV=polyval(SxN_coeffsL,t2)./polyval(SxD_coeffsL,t2);
        Sy_tV=polyval(SyN_coeffsL,t2)./polyval(SyD_coeffsL,t2);
        SxV=[SxV; Sx_tV]; SyV=[SyV; Sy_tV];

        S1x_tV=polyval(Sx1_coeffsL,tV);  S1y_tV=polyval(Sy1_coeffsL,tV);
        S1xV=[S1xV; S1x_tV]; S1yV=[S1yV; S1y_tV];

        % ..................  making monotone boxes ..................

        % Determining box (or boxes) in intervals relative to "a", "b" and
        % "turning points", so that in each box the rat. spline is a
        % monotone function.
        %
        % The k-th interval is subdivided first in "Nsub" subintervals.

        t12=linspace(t1,t2,Nsub+1); t12=t12';

        % Table "SxLvalues" with values (SxL(t),SyL(t),t) where:
        % 1. "t" is curvilinear coordinate of the extrema of the Nsub
        %     subintervals;
        % 2. "t" is curvilinear coordinate of an "x" turning point (i.e. of
        %    a point where the derivative of SxL is zero).
        % 3. "t" is curvilinear coordinate of an "y" turning point (i.e. of
        %    a point where the derivative of SyL is zero).

        if length(SxD_coeffsL) == 1 & SxD_coeffsL(end) == 1
            SxLvalues=[polyval(SxN_coeffsL,t12) ...
                polyval(SyN_coeffsL,t12) a+t12; ...
                turning_points_xL; ...
                turning_points_yL];
        else
            SxLvalues=[polyval(SxN_coeffsL,t12)./polyval(SxD_coeffsL,t12) ...
                polyval(SyN_coeffsL,t12)./polyval(SyD_coeffsL,t12) a+t12; ...
                turning_points_xL; ...
                turning_points_yL];
        end

        % Sort the curvilinear coordinates of "subinterval extrema" and "x"
        % or "y" turning points and consequently rearrange function values.


        % Sorted "t" curvilinear coordinates of relevant points for making
        % monotone boxes are saved in "tt".
        tt=SxLvalues(:,3);[tt_sort,I_unique]=unique(tt);
        boxtable_0=SxLvalues(I_unique,:);

        % Main information "preboxes" of monotone boxes, i.e. the k-th row
        % of "preboxes" is [Xm(k) XM(k) Ym(k) YM(k)] and thus the k-th box
        % is the rectangle "[Xm(k),XM(k)] x [Ym(k),YM(k)]".
        Xv=boxtable_0(:,1); Yv=boxtable_0(:,2); Tv=boxtable_0(:,3);

        Xm=min(Xv(1:end-1),Xv(2:end));
        XM=max(Xv(1:end-1),Xv(2:end));
        Ym=min(Yv(1:end-1),Yv(2:end));
        YM=max(Yv(1:end-1),Yv(2:end));

        preboxes=[Xm XM Ym YM Tv(1:end-1) Tv(2:end)];

        % Keeping track in the box additional information about the box.
        ss=ones(size(preboxes,1),1);
        boxL=[preboxes j*ss k*ss];
        boxVx=[boxVx; boxL];

        XXL=Xv(1:end-1); YYL=Yv(1:end-1); TTL=Tv(1:end-1);

        pgonL=[XXL YYL TTL];
        pgon=[pgon; pgonL];


    end

end


% Looking for "turning junctions", that is with non costant sign at
% junctions.
S1xVinit= [S1xV(2:end,1); S1xV(1,1)]; S1xVend=S1xV(:,2);
iXsing=find(S1xVinit.*S1xVend <= 0);

turning_junctionsX=[SxV(iXsing) SyV(iXsing) S1tV(iXsing)];

S1yVinit= [S1yV(2:end,1); S1yV(1,1)]; S1yVend=S1yV(:,2);
iYsing=find(S1yVinit.*S1yVend <= 0);
turning_junctionsY=[SxV(iYsing) SyV(iYsing) S1tV(iYsing)];

% Make unique set of "singular points".
turning_points_Sx=[turning_points_Sx; turning_junctionsX];
turning_points_Sy=[turning_points_Sy; turning_junctionsY];








%--------------------------------------------------------------------------
% turning_points_analysis
%--------------------------------------------------------------------------

function turning_points=turning_points_analysis(SxN_coefs,SxD_coefs,...
    SyN_coefs,SyD_coefs,Sx1_coefs,a,t1,t2)

%--------------------------------------------------------------------------
% OBJECT:
% Analysis of the x-turning points of the portion of the boundary of a
% bivariate domain, locally represented by (px(t),py(t)) with "t" in
% "[t1,t2]", where
%
%         px(t)=Sx_coefs(1)*t^n+Sx_coefs(2)*t^(n-1)+...
%         py(t)=Sy_coefs(1)*t^n+Sy_coefs(2)*t^(n-1)+...
%
% In general these coefficients "Sx_coefs", "Sy_coefs" are those of the
% k-th block with arguments belonging to [a+t1,a+t2] of two splines
% "Sx(i)", "Sy(i)" that are part of a parametric representation of the
% boundary of a spline-curvilinear domain.
% If "u" is in the interval "[a+t1,a+t2]" then
%
%                (Sx(i))(u)=px(u-a) and (Sy(i))(u)=py(u-a)
%
% Turning points are those points "v" in "[a+t1,a+t2]" such that the
% derivative of "(Sx(i))", say "(Sx1(i))", is null in "v".
%
% The y-turning points are obtaining applying:
%    turning_points=turning_points_analysis( SyN_coefs,SyD_coefs,...
%    SxN_coefs,SxD_coefs,Sy1_coefs,a,t1,t2)
%--------------------------------------------------------------------------
% INPUT:
% SxN_coefs, SxD_coefs: coefficients of the polynomial representing the
%     numerators and denominators of rational spline "Sx" in the
%    interval "[a+t1,a+t2]";
% SyN_coefs, SyD_coefs: coefficients of the polynomial representing the
%     numerators and denominators of rational spline "Sy" in the
%    interval "[a+t1,a+t2]";
% Sx1_coefs: coefficients of the polynomial at numerator of the derivative
%     of the spline "Sx" in the interval "[a+t1,a+t2]";
% a, t1, t2: parameters describing the interval "[a+t1,a+t2]".
%--------------------------------------------------------------------------
% OUTPUT:
%
% turning_points: N x 3 matrix of the turning points, where the k-th row is
%     such that P=(turning_points(k,1),turning_points(k,2)) are the
%     cartesian coordinates of the turning point P, while
%     "u=turning_points(k,3)" is the coordinate in "[a+t1,a+t2]" such that
%     P=((Sx(i))(u), (Sy(i))(u)).
%--------------------------------------------------------------------------
% SUBROUTINES (directly called):
% 1. myroots
%--------------------------------------------------------------------------
% Checked: October 21, 2021.
%--------------------------------------------------------------------------

% Determine the "t" curvilinear parameters in "[t1,t2]" of the turning
% points. Remember that the original rat.spline has variables in the interval
% [a,b] that are suitably scaled. This is due to the representation of the
% spline coefficients. More precisely, if "coefs" are the coefficients
% of the polynomial representation the spline "S" in its k-th
% block with breakpoints "a", "b" then, for "x" in "[a,b]", if "p" is the
% polynomial "p(t)=coefs(1)*t^n+coefs(2)*t^(n-1)+...", we have that
% "S(x)=p(x-a)".


[t_tp,flag]=myroots(Sx1_coefs,t1,t2);

% Below the values on "flag" can be:
%       0: finite complex solutions,
%       Inf: infinite complex solutions (problem p=[0 ... 0]).
%       NaN: no complex solutions (problem p=[0 ... 0 c] with "c" not 0).

if flag == 0
    if length(t_tp) > 0
        % Some turning points in the interval.
        t_tp=sort(unique(t_tp));
        if length(SxD_coefs) == 1 & SxD_coefs(end) == 1
            den=1;
        else
            den=polyval(SxD_coefs,t_tp);
        end
        SxLvalues_t_tp=polyval(SxN_coefs,t_tp)./den;
        SyLvalues_t_tp=polyval(SyN_coefs,t_tp)./den;
        StLvalues_t_tp=a+t_tp;
        turning_points=[SxLvalues_t_tp SyLvalues_t_tp StLvalues_t_tp];
    else
        turning_points=[];
    end
    %     else
    %         % No turning points in the interval. Consider just its extrema.
    %         t_tp=[t1; t2];
    %         SxLvalues_t_tp=polyval(SxN_coefs,t_tp)./polyval(SxD_coefs,t_tp);
    %         SyLvalues_t_tp=polyval(SyN_coefs,t_tp)./polyval(SyD_coefs,t_tp);
    %         StLvalues_t_tp=a+t_tp;
    %         turning_points=[SxLvalues_t_tp SyLvalues_t_tp StLvalues_t_tp];
    %    end
else
    if flag == Inf % Analysis: infinite turning points.
        % Geometrical view: vertical segment.
        t_tp=[t1; t2];

        if length(SxD_coefs) == 1 & SxD_coefs(end) == 1
            den=1;
        else
            den=polyval(SxD_coefs,t_tp);
        end
        SxLvalues_t_tp=polyval(SxN_coefs,t_tp)./den;
        SyLvalues_t_tp=polyval(SyN_coefs,t_tp)./den;
        StLvalues_t_tp=a+t_tp;
        turning_points=[SxLvalues_t_tp SyLvalues_t_tp StLvalues_t_tp];
    else % Analysis: no turning points.
        turning_points=[];
    end
end









%--------------------------------------------------------------------------
% myroots.
%--------------------------------------------------------------------------

function [sols,flag]=myroots(p,t1,t2)

%--------------------------------------------------------------------------
% OBJECT:
% This routine solves the polynomial equation
%                     p(1)*x^n+...+p(n+1)=0
% determining the solutions in the interval (t1,t2).
%--------------------------------------------------------------------------
% INPUT:
% p: vector determining the polynomial "p(1)*x^n+...+p(n+1)=0".
% t0,t1: the solutions must be in the interval "[t0,t1]".
%--------------------------------------------------------------------------
% OUTPUT:
% sols: vector of solutions.
% flag: 0: finite complex solutions,
%       Inf: infinite complex solutions (problem p=[0 ... 0]).
%       NaN: no complex solutions (problem p=[0 ... 0 c] with "c" not 0).
%--------------------------------------------------------------------------
% Checked: September 20, 2020.
%--------------------------------------------------------------------------

% .................... troubleshooting and defaults .......................

if nargin < 3
    t2=+inf;
end

if nargin < 2
    t1=-inf;
end

% ........................ computing roots ................................
sols=roots(p);


if isempty(sols) == 0
    % ..... In this case the equation has complex roots .....
    isols=imag(sols); rsols=real(sols);
    isols_in=find(abs(isols)< 10^(-12) & rsols >= t1 & rsols <= t2);
    sols=rsols(isols_in);
    flag=0;
else
    % ... In this case the equation has infinite or no complex roots ......
    if p(end) == 0 % infinite complex roots
        flag=Inf;
    else
        flag=NaN; % no complex roots
    end
end









function [in,on]=indomainRS(P,SxNV,SxDV,SyNV,SyDV,boxVxy,...
    turning_points_x,tol,safe_mode)

%--------------------------------------------------------------------------
% OBJECT:
% This routine checks if a point is inside a region "D" whose boundary is
% defined by the rational splines "Sx" and "Sy".
%
% The piecewise rational functions "x=Sx(t)", "y=Sy(t)", are stored via the
% splines "SxN", "SxD", "SyN", "SyD" (in "pp" form), so that the boundary
% of "D" is defined parametrically as (x(t),y(t)), t in [0,1] where
%       x(t)=Sx(t)=SxN(t)/SxD(t), y(t)=Sy(t)=SyN(t)/SyD(t)     t in [0,1].
%
%
% Note: the splines "SxN", "SxD", "SyN", "SyD" have the same order and
%    breakpoints and are in "pp" form.
%    Alternatively "SxD" and "SyD" are the constant function "1".
%--------------------------------------------------------------------------
% INPUT:
% P: points to test, defined as an N x 2 matrix.
%
% SxN: numerator of x(t), as spline in "pp" form;
%
% SxD: denominator of x(t), as spline in "pp" form;
%
% SyN: numerator of y(t), as spline in "pp" form;
%
% SyD: denominator of y(t), as spline in "pp" form;
%
% boxVxy: S x 8 matrix, whose i-th row describes the properties of the i-th
%       box of the boundary, i.e.
%
%      *           [xm,xM,ym,yM,t0,t1,ispline,iblock]
%
%      where [xm,xM,ym,yM] defines the boundary of the box, i.e. the
%      vertices
%
%                   X=[xm xM xM xm]   Y=[ym ym yM yM]
%
%      such that the box consists of the rectangle with vertices
%      (X(k),Y(k)), k=1,2,3,4.
%
%      * the parametrization of x=x(t), y=y(t) has values of "t" in the
%        interval [t0,t1];
%
%      * "ispline" says what rat. spline index is defined in the box, i.e.
%
%             x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t);
%
%      * "iblock" says what spline block index is defined in the box, i.e.
%
%             x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t)
%
%        with
%
%        "t in [t1,t2]", where "t1=(Sx(ispline)).breaks(iblock)" and
%        "t2=(Sx(ispline)).breaks(iblock+1)".
%
% turning_points_x: abscissae of singular points. This variable is not
%       mandatory. Default: turning_points_x=[];
%
% tol: tolerance in detecting if the point is in the box by a polynomial
%      equation solver. This variable is not mandatory. Default: 10^(-11)
%--------------------------------------------------------------------------
% OUTPUT:
% in: the i-th component is
%     1 if P(i,:) is inside the domain,
%     0 if not inside.
% idoubts: points close to turjning points;
% on: the i-th component is
%      * Nan    : point is "safely" outside or inside the domain.
%      * finite : point is "close" to the boundary
%      * -Inf   : problems with the determination inside/on/outside.
%--------------------------------------------------------------------------
% SUBROUTINES (directly called):
% 1. spline_boxer_RS;
% 2. box_analysis_RS_02.
%--------------------------------------------------------------------------
% REFERENCE:
%--------------------------------------------------------------------------
% [1] Les Piegl, Wayne Tiller,
% The NURBS Book,
% Second edition (Springer Verlag, 1995 and 1997).
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
%% Date: October 31, 2021.
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Checked: December 09, 2021 (15:40).
%--------------------------------------------------------------------------


% ................. troubleshoouting and defaults .........................

if nargin < 7, turning_points_x=[]; end
if nargin < 8, tol=10^(-11); end
if nargin < 9, safe_mode=1; end % additional tests for indomain

% ......................... general study .................................

% We check if:
% a) a point P(k,:) is inside some monotone box and in this case
%    we perform some further analysis. This happens if the k-th column of
%    the variable "itestL" is not equal to the vector with all the
%    components equal to zero. The components equal to "1" say to what
%    boxes such point belong.
% b) a point P(k,:) is outside any monotone box. This happens if the k-th
%    column of the variable "itestL" is equal to the vector with all
%    the components equal to zero. In this case by its crossing number
%    "crossings(k)" we check if it is inside or outside the
%    spline-curvilinear polygon by means of the crossing theorem.
%
% Note: observe that "X" and "Y" are row vectors, while "a", "b", "c", "d"
%    are column vectors.

X=(P(:,1))'; Y=(P(:,2))';
a=boxVxy(:,1); b=boxVxy(:,2); c=boxVxy(:,3); d=boxVxy(:,4);
% iboxes=(X > a & X <=b & Y >= c); % #boxes x #X
% itestL=(iboxes & Y <= d); % #boxes x #X
% iregL=(iboxes & Y > d); % #boxes x #X
% crossings=sum(iregL,1); % #boxes x #X
[crossings,inboxV,ptsinbox]=inRSpre(P,boxVxy);


% ..................... Turning points analysis ...........................

% Here we pay attention to points close to the vturning points of the box.
% Skipping this test, due to numerical problems, the code may provide wrong
% results.

if isempty(turning_points_x) == 0
    tpts_x=turning_points_x(:,1);
    mindist=min(abs(X-tpts_x),[],1);
    isclose=mindist < 10^(-13);
else
    isclose=zeros(size(X,1),1);
end

idoubts=isclose;

% The k-thvalue of "on" is "-inf" in case the k-th point is too close
% turning points (but with a least P^*P crossing the domain), "NaN"
% otherwise.

on=-idoubts./zeros(size(X));



% ..................... Vertical sides analysis ...........................

% Points belonging to vertical sides (considering numerical issues)
if safe_mode
    tol=10^(-13);
    ivert=find(b-a <= tol);
    a0=a(ivert); c0=c(ivert); d0=d(ivert);
    iboxes0=( abs(X - a0) <= tol & Y >= c0 & Y <= d0); % #boxes x #X
    ind_vert=find(sum(iboxes0,1) > 0);
    on(ind_vert)=tol;
end


% ..................... Horizontal sides analysis .........................

% Points belonging/close to horizontal sides (considering numerical issues)
if safe_mode
    tol=10^(-13);
    ihor=find(d-c <= tol);
    a0=a(ihor); b0=b(ihor); c0=c(ihor);
    distY=abs(Y - c0);
    iboxes0=( distY <= tol & X >= a0 & X <= b0); % #boxes x #X
    ind_hor=find(sum(iboxes0,1) > 0);
    on(ind_hor)=tol;
end




% ..................... Study points inside boxes .........................

% Indices of the points requiring some special analysis.
% We analyse also points with abscissae close to those of turning points.
% In case we reject them later.
% tic;
further_analysis_index=find(ptsinbox == 1 & isnan(on));

for jj=1:length(further_analysis_index)
    Pindex=further_analysis_index(jj);
    PL=P(Pindex,:);

    index_testPL=inboxV{Pindex};

    for k=1:length(index_testPL)
        a=boxVxy(index_testPL,1); b=boxVxy(index_testPL,2);
        c=boxVxy(index_testPL,3); d=boxVxy(index_testPL,4);
    end

    [final_value0,on(Pindex)]=box_analysis_RS_02(PL,index_testPL,boxVxy,...
        SxNV,SxDV,SyNV,SyDV,tol);

    crossings(Pindex)=crossings(Pindex)+final_value0;
end



in=rem(crossings,2);

if safe_mode
    iilog=not(isnan(on)) & not(isinf(on));
    ii=find(iilog == 1);
    in(ii)=0.5;
end





% ..................... Too low difficult points ..........................

% points with abscissa at turning point but too low.
if safe_mode
    iout=find(crossings == 0 & on == -Inf);
    on(iout)=NaN*ones(size(iout));
    in(iout)=zeros(size(iout));
end

















%--------------------------------------------------------------------------
% box_analysis
%--------------------------------------------------------------------------

function [final_value,flag]=box_analysis_RS_02(P,itest,boxV,SxNV,SxDV,...
    SyNV,SyDV,tol)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine checks if a point "P" that belongs to some boxes is "below",
% "above" or "on" the portion of boundary belonging to each monotone box in
% consideration, defined by "boxV(itest,:)".
% If "P" it not "on" any box, it says how many times the vertical straight
% line below "P" intersect the boundary.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% P: point to test, defined as an 1 x 2 matrix.
%
% itest: index of the box in the set described by "boxV", to be studied; if
%     this vector is not empty then the point P belongs to the boxes
%     represented by "boxV(itest,:)".
%
% boxV: S x k vector, whose i-th component defines the
%      properties of the i-th box of the boundary, i.e.
%      *           [xm xM ym yM tm tM isp ibl];
%      where [xm,xM,ym,yM] defines the boundary of the box, i.e. vertices
%                   X=[xm xM xM xm]   Y=[ym ym yM yM]
%      * [tm,tM] t variable extremes
%      * ispline says what spline index is defined in the box.
%      * iblock says what spline block index is defined in the box, i.e.
%        x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t) with
%        "t in [t1,t2]", where "t1=(Sx(ispline)).breaks(iblock)" and
%        "t2=(Sx(ispline)).breaks(iblock+1)".
%
% SxN: numerator of x(t), as spline in "pp" form;
%
% SxD: denominator of x(t), as spline in "pp" form;
%
% SyN: numerator of y(t), as spline in "pp" form;
%
% SyD: denominator of y(t), as spline in "pp" form;
%
% tol: tolerance in detecting the point in the box by polynomial equation
%    solver.
%    This variable is not mandatory. Experiments suggest to use
%     "tol=10^(-11)" as default.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% final_value: crossing value, if
%
%       0:    the point is "numerically" on the boundary
%       odd:  the point is in the spline curvilinear domain
%       even: the point is not in the spline curvilinear domain
%
% flag:
%      * Nan    : point is "safely" outside or inside the domain.
%      * finite : point is "close" to the boundary
%      * -Inf   : problems with the determination inside/on/outside.
%--------------------------------------------------------------------------
% Initial version: October 18, 2021.
%--------------------------------------------------------------------------


% ........................ Some troubleshooting ...........................

if nargin < 5, singptsX=[]; end

if nargin < 6, tol=10^(-11); end
x=P(1); y=P(2);

flag=NaN;
final_value=0;

% .................... Analysis of special cases ..........................

if length(itest) > 0

    for ii=1:length(itest)

        % "itestL": indices to test

        itestL=itest(ii); % Index of the to analyse.

        % The variable "isplL" is the index of the spline related to the
        % box, i.e. the box contains a portion of the boundary
        % parametrically defined by "SxL=Sx(isplL)" and "SyL=Sy(isplL)".

        isplL=boxV(itestL,7);
        icoefsL=boxV(itestL,8); % Coeff. of the spline related to the box.

        SxN=SxNV(isplL);
        SxD=SxDV(isplL);
        SyN=SyNV(isplL);
        SyD=SyDV(isplL);

        % In the next part we test if the point is "above", "below" or even
        % "on" the portion of the boundary of the spline-curvilinear domain
        % that is in the box.
        %
        % Due to the different analysis we separate the case of linear
        % spline (i.e. order equal to 2) from those of different order.


        % ................ SPLINES OF ORDER > 2 ................


        if SxN.order > 2

            SxNL_coefs=SxN.coefs(icoefsL,:);
            SxDL_coefs=SxD.coefs(icoefsL,:);
            SyNL_coefs=SyN.coefs(icoefsL,:);
            SyDL_coefs=SyD.coefs(icoefsL,:);

            a=boxV(itestL,5);
            b=boxV(itestL,6);

            t0=SxN.breaks(icoefsL);
            t1=SxN.breaks(icoefsL+1);

            [local_value,flag]=spline_over_RS(P,SxNL_coefs,SxDL_coefs,...
                SyNL_coefs,SyDL_coefs,a,b,t0,t1,tol);


            % Technical note.
            % At the end of the routine "flag" is a number, if on the
            % boundary, or NaN (point has been well analysed in the
            % crossing).

            if isfinite(flag) == 1
                final_value=0; % point is on the boundary and
                % not inside the domain
                return;
            end

        else

            % ................ SPLINES OF ORDER = 2 ................

            nodes=SxN.breaks;
            a=boxV(itestL,5);
            b=boxV(itestL,6);

            P1=[ppval(SxN,a)./ppval(SxD,a) ppval(SyN,a)./ppval(SyD,a)];
            P2=[ppval(SxN,b)./ppval(SxD,b) ppval(SyN,b)./ppval(SyD,b)];

            % ... ... ... sides are not vertical ... ... ...

            if abs(P1(1)-P2(1)) > 0

                % He we compute the intersection of the vertical straight
                % line containing P with the segment that is the portion of
                % boundary in the box.

                x0=P1(1); y0=P1(2);
                x1=P2(1); y1=P2(2);
                m=(y1-y0)/(x1-x0);
                yy=y0+m*(x-x0);

                if abs(yy-y) <= tol % point "numerically" on the boundary
                    final_value=0;
                    flag=abs(yy-y); % "closeness" to the boundary
                    return;
                else % point is not on the boundary
                    if yy < y
                        local_value=1; % crossing
                    else
                        local_value=0; % not crossing
                    end
                end


            else
                % In this case the point has coordinates (x,y) with
                % x=P1(1)=P2(1). It is a difficult instance and must be
                % analyzed by winding if the point is not on the side.
                % Winding requirement is communicated by "flag=-Inf".

                min_y=min(P1(2),P2(2));
                max_y=max(P1(2),P2(2));

                % point on a vertical side
                if (y >= min_y) & (y <= max_y) % Treatable case
                    flag=0;
                    final_value=0;
                    return;
                else % Suggesting the use of "winding number algorithm"
                    flag=-Inf;
                    final_value=NaN;
                    return;
                end

            end


        end

        if isnan(flag) == 1 % Crossing number can be computed
            final_value=final_value+local_value;
        else % points has some issues
            final_value=0;
            return;
        end

    end

end











%--------------------------------------------------------------------------
% spline_over
%--------------------------------------------------------------------------

function [val,flag]=spline_over_RS(P,xNtv,xDtv,yNtv,yDtv,a,b,t0,t1,tol)

%--------------------------------------------------------------------------
% INPUT:
% P   : point to test, 1 x 2 vector. It is supposed that the point is in
%       the rectangle with "SxL([t0,t1]) x SyL([t0,t1])" whose local
%       representation is provided respectively by "xtv" and "ytv".
%
% xNtv,xDtv: "local" coefficients of the rational spline "SxL" defining
%       x=x(t).
%
% yNtv,yDtv: "local" coefficients of the rational spline "SyL" defining
%       y=y(t).
%
% a, b : spline "SxL","SyL" extrema (subsequent break points of the "SxL",
%       "SyL").
%
% t0, t1: the interval [t0,t1] is a subset of [0,1] (the portion of the
%        boundary to analyse has the form (SxL(t),SyL(t)) with
%        "t in [t0,t1]").
%
% tol : tolerance of the boundary (to say how much P is close to the
%      boundary).
%--------------------------------------------------------------------------
% OUTPUT:
%
% val: this variable has as value how many times a straight line with
%      x=P(1) intersects the spline boundary provided respectively by "xtv"
%      and "ytv".
%
% flag:
%     * NaN           : normal exit, point not on the boundary.
%     * finite number : parameter that says how close is to the border; in
%                       this case the point is accepted as on the boundary
%                       since its distance from the boundary is below a
%                       certain threshold (given by "tol").
%     * -Inf          : the point P has abscissa not in the range of
%                       "SxL([t0,t1])".
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
%% Date: October 31, 2021.
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Checked: October 31, 2021.
%--------------------------------------------------------------------------

% ..................... setting missing defaults ..........................

if nargin < 6, tol=10^(-12); end

% ......................... initial settings ..............................

val=0;
flag=NaN;
x=P(1); y=P(2);
s0=a-t0; s1=b-t0;

% ...................... compute intersections ............................

LN=length(xNtv); LD=length(xDtv);

if LN == LD
    xtv0=xNtv-x*xDtv;
else
    if LD == 1
        xtv0=xNtv;
        xtv0(end)=xtv0(end)-x*xDtv(end);
    else
        if LN > LD
            xDtv=[zeros(LN-LD,1) xDtv];
        else
            xNtv=[zeros(LD-LN,1) xNtv];
        end
        xtv0=xNtv-x*xDtv;
    end
end

% ..... special case: the monotone box is a vertical line .....

% NOTE: In the present code, in view of the calls from outside, the next
% condition never holds.
if norm(xNtv(1:end-1)) == 0 & norm(xDtv(1:end-1)) == 0

    % In this case the spline "Sx" has locally numerator coefficients
    % "xNtv=[0 ... 0 c]" denoting the fact that it is actually a vertical
    % line.

    if xtv0(end) == 0
        % In this case, we study for what "t" it is "0=0".
        val=0;
        flag=0;

    else
        % In this case, we study for what "t" it is "C=0" for C not 0.
        val=NaN;
        flag=-Inf;
    end
    return;
end


% ..... normal case: the monotone box is not a vertical line .....
tR=myroots(xtv0,s0,s1);

% check relevant boxes
for j=1:length(tR)

    tRL=tR(j);
    yL=polyval(yNtv,tRL)./polyval(yDtv,tRL);

    % point "safely" over the boundary
    if yL < y+tol
        val=val+1;
    end

    % point "safely" close to the boundary
    if abs(yL-y) <= tol
        val=0;
        flag=abs(yL-y);
        return;
    end

end







%--------------------------------------------------------------------------
% winding_algorithm
%--------------------------------------------------------------------------

function in=windingRS(SxNV,SxDV,SyNV,SyDV,P)

% P: points of which we compute the winding number.

L=size(P,1);

for ii=1:L
    windn = winding_algorithm(P(ii,:),SxNV,SxDV,SyNV,SyDV);
    in(ii)=rem(windn,2);
end





%--------------------------------------------------------------------------
% winding_algorithm
%--------------------------------------------------------------------------

function windn = winding_algorithm(P,SxNV,SxDV,SyNV,SyDV)

%--------------------------------------------------------------------------
% OBJECT:
% computes the winding number of P with respect to the counterclockwise
% concatened cubic spline arcs (Sx,Sy) forming a Jordan spline polygon,
% by contour integration
%--------------------------------------------------------------------------
% INPUT:
% P: 2-column array of planar points
% Sx,Sy: arrays of spline structures; (SX(i),Sy(i)) is the i-th arc
% the arcs must be counterclockwise concatenated forming a Jordan curve
% Sx(i).breaks(end)=Sx(i+1).breaks(1),Sy(i).breaks(end)=Sy(i+1).breaks(1)
% i=1,...,end, with end+1:=1
%--------------------------------------------------------------------------
% OUTPUT:
% windn: winding numbers of P with respect to the Jordan curve (Sx,Sy)
% windn(j)==1 point is inside, windn(j)==0 point is outside
%--------------------------------------------------------------------------
% DATA:
% built: October 2021
% checked: November 12, 2021
%--------------------------------------------------------------------------

windn=zeros(1,length(P(:,1)));

KK=length(SxNV);
xw=legendre_rules(10000);

for ii=1:KK

    SxN=SxNV(ii); SxD=SxDV(ii); SyN=SyNV(ii); SyD=SyDV(ii);
    L=length(SxN.breaks);

    derV=(SxN.order-1):-1:1;

    for k=1:L-1


        a=SxN.breaks(k); b=SxN.breaks(k+1);

        % Spline coefficients defining "Sx" and "Sy" in the k-th
        % subinterval "[a,b]".
        SxN_coeffsL=SxN.coefs(k,:); SxD_coeffsL=SxD.coefs(k,:);
        SyN_coeffsL=SyN.coefs(k,:); SyD_coeffsL=SyD.coefs(k,:);

        SxN_deg=length(SxN_coeffsL)-1;
        if SxN_deg > 0
            derV=SxN_deg:-1:1; SxN1_coeffsL=SxN_coeffsL(1:end-1).*derV;
        else
            SxN1_coeffsL=0;
        end

        SxD_deg=length(SxD_coeffsL)-1;
        if SxD_deg > 0
            derV=SxD_deg:-1:1; SxD1_coeffsL=SxD_coeffsL(1:end-1).*derV;
        else
            SxD1_coeffsL=0;
        end

        SyN_deg=length(SyN_coeffsL)-1;
        if SyN_deg > 0
            derV=SyN_deg:-1:1; SyN1_coeffsL=SyN_coeffsL(1:end-1).*derV;
        else
            SyN1_coeffsL=0;
        end

        SyD_deg=length(SyD_coeffsL)-1;
        if SyD_deg > 0
            derV=SyD_deg:-1:1; SyD1_coeffsL=SyD_coeffsL(1:end-1).*derV;
        else
            SyD1_coeffsL=0;
        end

        % Coefficients of the numerators derivatives. Note that these
        % derivatives are computed to detect zeros of the derivatives of
        % the rational splines "Sx" and "Sy" that in the interval "[a,b]"
        % are rational functions, say "SxL" and "SyL".
        % Denominators of the derivatives are equal to  "(SxD)^2" and
        % "(SyD)^2" and do not contribute to the detection of these zeros.
        % So only the numerators of the derivatives of "Sx" and "Sy" are
        % taken into account.
        %
        % The code becomes weird due to different components in the next
        % convolutions, e.g. using non rational splines.

        % length of convolution vector
        LL(1)=length(SxN1_coeffsL)+length(SxD_coeffsL)-1;
        LL(2)=length(SxD1_coeffsL)+length(SxN_coeffsL)-1;
        MM=max(LL);

        term1=zeros(1,MM); res=conv(SxN1_coeffsL,SxD_coeffsL);
        term1(MM-LL(1)+1:MM)=res;

        term2=zeros(1,MM); res=conv(SxD1_coeffsL,SxN_coeffsL);
        term2(MM-LL(2)+1:MM)=res;

        Sx1_coeffsL=term1-term2; % numerator derivative "x" direction



        % length of convolution vector
        LL(1)=length(SyN1_coeffsL)+length(SyD_coeffsL)-1;
        LL(2)=length(SyN_coeffsL)+length(SyD1_coeffsL)-1;
        MM=max(LL);

        term1=zeros(1,MM); res=conv(SyN1_coeffsL,SyD_coeffsL);
        term1(MM-LL(1)+1:MM)=res;

        term2=zeros(1,MM); res=conv(SyN_coeffsL,SyD1_coeffsL);
        term2(MM-LL(2)+1:MM)=res;

        Sy1_coeffsL=term1-term2;  % numerator derivative "y" direction

        x=@(s) polyval(SxN_coeffsL,s-a)./polyval(SxD_coeffsL,s-a);
        dx=@(s) polyval(Sx1_coeffsL,s-a)./(polyval(SxD_coeffsL,s-a)).^2;
        y=@(s) polyval(SyN_coeffsL,s-a)./polyval(SyD_coeffsL,s-a);
        dy=@(s) polyval(Sy1_coeffsL,s-a)./(polyval(SyD_coeffsL,s-a)).^2;

        f=@(s,j) (dy(s).*(x(s)-P(j,1))-dx(s).*(y(s)-P(j,2)))./...
            ((x(s)-P(j,1)).^2+(y(s)-P(j,2)).^2);



        nodes=xw(:,1)*(b-a)/2+(b+a)/2;
        weights=xw(:,2)*(b-a)/2;

        [u,v]=meshgrid(nodes,1:length(P(:,1)));
        fval=f(u(:),v(:));
        FV=reshape(fval,length(P(:,1)),length(nodes));
        int=FV*weights;

        windn=windn+int';


    end
end

windn=round(windn/(2*pi));












function [crossing,inboxV,ptsinbox]=inRSpre(pts,boxV)

%--------------------------------------------------------------------------
% OBJECT:
%  Given the pointset "pts" and the monotone boxes stored in "boxV", we
%  provide the "sure" crossing values for each points, while "inboxV" is an
%  array of cells in which the k-th component contains pointers to all the
%  boxes containing "pts(k,:)".
%  If "inboxV{k}" is not empty then "ptsinbox(k)=1", otherwise
%  "ptsinbox(k)=0"
%--------------------------------------------------------------------------
% INPUT:
% pts: matrix N x 2 of cartesian coordinates of the pointset.
% boxV: M x 8 matrix that describes the monotone boxes; the i-th row
%       describes the properties of the i-th box of the boundary, i.e.
%
%      *           [xm,xM,ym,yM,t0,t1,ispline,iblock]
%
%      where [xm,xM,ym,yM] defines the boundary of the box, i.e. the
%      vertices
%
%                   X=[xm xM xM xm]   Y=[ym ym yM yM]
%
%      such that the box consists of the rectangle with vertices
%      (X(k),Y(k)), k=1,2,3,4.
%
%      * the parametrization of x=x(t), y=y(t) has values of "t" in the
%        interval [t0,t1];
%
%      * "ispline" says what rat. spline index is defined in the box, i.e.
%
%             x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t);
%
%      * "iblock" says what spline block index is defined in the box, i.e.
%
%             x(t)=[Sx(ispline)](t), y(t)=[Sy(ispline)](t)
%
%        with
%
%        "t in [t1,t2]", where "t1=(Sx(ispline)).breaks(iblock)" and
%        "t2=(Sx(ispline)).breaks(iblock+1)".
%--------------------------------------------------------------------------
% OUTPUT:
% crossing: N x 1 array containing the "sure" crossing values of each pure;
% inboxV: array of cells in which the k-th component contains pointers to
%          all the boxes containing "pts(k,:)"; in other words, if j
%          belongs to the vector inboxV{k}, then the box described by
%          "boxV(k,:)" contains "pts(k,:)" and the values described by
%          crossing(k) require further investigations.
% ptsinbox: N dimensional array, in which if "inboxV{k}" is not empty then
%        "ptsinbox(k)=1", otherwise "ptsinbox(k)=0"
%--------------------------------------------------------------------------
% DATA:
% built: January 2, 2022.
% checked: January 5, 2022.
%--------------------------------------------------------------------------

X=pts(:,1); Y=pts(:,2); Npts=size(X,1);
[Px,isort]=sort(X); Py=Y(isort);

Nbox=size(boxV,1);

crossing0=zeros(Npts,1);
inboxV0=cell(Npts,1);
ptsinbox0=zeros(Npts,1);

for j=1:Nbox
    boxL=boxV(j,:);
    aL=boxL(1); bL=boxL(2); cL=boxL(3); dL=boxL(4);
    i0=binarySearchL01(Px,aL);

    if not(isnan(i0))
        k=i0; flag=0;
        while k <= Npts & flag== 0
            if Px(k) > bL
                flag=1;
            else
                if Py(k) > dL
                    crossing0(k)=crossing0(k)+1;
                else
                    if Py(k) >= cL
                        ptsinbox0(k)=1;
                        inboxV0{k}=[inboxV0{k} j];
                    end
                end
            end
            k=k+1;
        end
    end
end


crossing(isort)=crossing0;
inboxV(isort)=inboxV0;
ptsinbox(isort)=ptsinbox0;








function index = binarySearchL01(xV,numL)

%--------------------------------------------------------------------------
% OBJECT:
%  Given a vector xV, this routine computes the first component "index"
%  such that "xV(index) > numL".
%  The vector "xV" is sorted in increasing order.
%--------------------------------------------------------------------------
% INPUT:
% xV: vector of scalars, sorted in increasing order.
% numL: scalar value.
%--------------------------------------------------------------------------
% OUTPUT:
% index: positive integer equal to the first component such that
% "xV(index) > numL".
%--------------------------------------------------------------------------
% DATA:
% built: January 2, 2022.
% checked: January 5, 2022.
%--------------------------------------------------------------------------

left = 1;
right = length(xV);
flag = 0;

while left < right-1
    mid = ceil((left + right) / 2);
    if xV(mid) <= numL
        left=mid+1;
    else
        right=mid;
    end
end

if xV(left) <= numL
    index=right;
else
    index=left;
end

if xV(index) < numL
    index = NaN;
end




