function Pend=lastpointNURBSPL(geometry_NURBSPL)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine computes the last point of a NURBSPL object.
%--------------------------------------------------------------------------

% geometry_NURBSPLL=geometry_NURBSPL(end);
% [SxN,SxD,SyN,SyD]=NURBSPL2ratsplines(geometry_NURBSPLL,'B-');
% 
% 
% tV=geometry_NURBSPLL.knots; t1=tV(end);
% Pend=[fnval(SxN,t1)/fnval(SxD,t1) fnval(SyN,t1)/fnval(SyD,t1)];

geometry_NURBSPL1=geometry_NURBSPL(end);
geometry_type=geometry_NURBSPL1.type;

switch geometry_type
    case 'bezier'
        P0=geometry_NURBSPL1.P; Pend=P0(end,:);
    case 'spline'
        P0=geometry_NURBSPL1.P; Pend=P0(end,:);
    otherwise % NURBS
        [SxN,SxD,SyN,SyD]=NURBSPL2ratsplines(geometry_NURBSPL1,'B-');
        tV=geometry_NURBSPL1.knots; t1=tV(end);
        Pend=[fnval(SxN,t1)/fnval(SxD,t1) fnval(SyN,t1)/fnval(SyD,t1)];
end
