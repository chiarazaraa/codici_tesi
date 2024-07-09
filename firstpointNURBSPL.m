

function P0=firstpointNURBSPL(geometry_NURBSPL)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine computes the first point of a NURBS/SPLINE object.
%--------------------------------------------------------------------------

geometry_NURBSPL1=geometry_NURBSPL(1);
geometry_type=geometry_NURBSPL1.type;

switch geometry_type
    case 'bezier'
        P0=geometry_NURBSPL1.P; P0=P0(1,:);
    case 'spline'
        P0=geometry_NURBSPL1.P; P0=P0(1,:);
    otherwise % NURBS
        [SxN,SxD,SyN,SyD]=NURBSPL2ratsplines(geometry_NURBSPL1,'B-');
        tV=geometry_NURBSPL1.knots; t0=tV(1);
        P0=[fnval(SxN,t0)/fnval(SxD,t0) fnval(SyN,t0)/fnval(SyD,t0)];
end