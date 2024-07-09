function structure_NURBSPL=joinNURBSPLarcs(geometry_NURBSPL)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine from structured object describing the NURBS/spline arc,
% joinNURBSPLarcs useful NURBS/spline control points, weights, knots, order.
%
% Each component of "geometry_NURBSPL" is a structured array with
%
%    * geometry_NURBSPL.P as control points,
%    * geometry_NURBSPL.w as weights,
%    * geometry_NURBSPL.knots as knots (important: with elements in [0,1]),
%    * geometry_NURBSPL.order as order.
%
% 1. Piecewise NURBS/spline of the same order are joined into one
%    NURBS/spline.
% 2. The piecewise NURBS/spline are joined so to define a "logical"
%    parametric interval (increasing along the components of the vector
%    "geometry_NURBSPL").
%--------------------------------------------------------------------------

Npieces=0; % Number of pieces of different order being saved.
L=length(geometry_NURBSPL);

for k=1:L
    geometry_NURBSPL_L=geometry_NURBSPL(k);

    % knots
    knotsL=geometry_NURBSPL_L.knots;
    m=min(knotsL); M=max(knotsL);
    knotsL=(knotsL-m)/(M-m); % normalize to [0,1]
    knotsL=((k-1)+knotsL)/L;

    geometry_NURBSPL_L.knots=knotsL;
    structure_NURBSPL(k)=geometry_NURBSPL_L;
end







