
function geometry_NURBS=makeNURBSarc(domain_str,varargin)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine, given a certain domain by its geometrical definitions
% provides the relative NURBS parameters (control points, weights, knots,
% order) storing them in a structured array "structure_arc".
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% First version: October 31, 2021;
% Checked: November 18, 2021.
%--------------------------------------------------------------------------

switch domain_str
    case 'disk_arc'
        %% Example of the outer call:
        % geometry_NURBS=addNURBSarc('disk_arc',...
        % 'center',[0 0],'angles',[0 pi/2],'radius',1);
        center=varargin{2}; theta=varargin{4}; r=varargin{6};
        [P,w,knots,order]=compute_arc_parms(center,r,theta(1),theta(2));
        geometry_NURBS.shape='disk_arc';
        
    case 'elliptical_arc'
        %% Example of the outer call:
        % geometry_NURBS=addNURBSarc('elliptical_arc',...
        % 'center',Pend-[1 0],'angles',[0 pi+pi/4],'ell_axis',[1 2],...
        % 'tilt_angle',0);
        center=varargin{2}; theta=varargin{4}; ell_axis=varargin{6};
        if nargin < 8
            tilt_angle = 0;
        else
            tilt_angle=varargin{8};
        end
        [P,w,knots,order]=compute_ellarc_parms(center,theta,ell_axis,...
            tilt_angle);
        geometry_NURBS.shape='elliptical_arc';
        
        
    case 'polygonal_arc'
        %% Example of the outer call:
        % geometry_NURBS=addNURBSarc('polygonal_arc',...
        % 'vertices',[1 0; 1 1; 0 0]);
        vertices=varargin{2};
        [P,w,knots,order]=compute_polygon_parms(vertices);
        geometry_NURBS.shape='polygonal_arc';
        
    case 'segment'
        %% Example of the outer call:
        % geometry_NURBS=addNURBSarc('polygonal_arc',...
        % 'vertices',[1 0; 1 1]);
        vertices=varargin{2};
        [P,w,knots,order]=compute_segments_parms(vertices);
        geometry_NURBS.shape='segment';
        
    case 'free'
        %% Example of the outer call:
        % geometry_NURBS=addNURBSarc('free',...
        % 'P',[-2 1; -1.9 0.3; -1.8 0.5; -1.7 5],...
        % 'knots',[0 0 0 0.5 1 1 1],'weights',[1 1 2 1],'order',3);
        P=varargin{2}; knots=varargin{4}; w=varargin{6}; order=varargin{8};
        geometry_NURBS.shape='free';
        
end

geometry_NURBS.P=P;
geometry_NURBS.w=w;
geometry_NURBS.knots=knots;
geometry_NURBS.order=order;
geometry_NURBS.type='NURBS';












%--------------------------------------------------------------------------
% ATTACHED ROUTINES:
%--------------------------------------------------------------------------
% 1. compute_arc_parms
% 2. compute_polygon_parms
% 3. compute_ellarc_parms
%--------------------------------------------------------------------------


function [P,w,knots,order]=compute_arc_parms(center,r,theta0,theta1)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine defines an arc of circle.
% Arc equation (without tilting):
% * x(t)=center(1)+r*cos(t)
% * y(t)=center(2)+r*sin(t)
% with "theta" in "[theta0,theta1]".
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% center: center of the circle containing the arc.
% r: radius of the circle containing the arc.
% theta0,theta1: angles of the extrema of the arc.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% P: control points.
% w: NURBS weights.
% knots: NURBS weights.
% order: NURBS order.
%--------------------------------------------------------------------------
% DATE:
%--------------------------------------------------------------------------
% Written: November 4, 2021.
% Checked: November 18, 2021.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Troubleshooting
%--------------------------------------------------------------------------

if nargin < 1, center=[0 0]; end
if nargin < 2, r=1; end
if nargin < 3, theta0=0; end
if nargin < 4, theta1=2*pi; end

if isempty(center), center=[0 0]; end
if isempty(r), r=1; end
if isempty(theta0), theta0=0; end
if isempty(theta1), theta1=1; end

%--------------------------------------------------------------------------
% Control Points
%--------------------------------------------------------------------------

if theta0 > theta1
    flag=0;
    [theta1,theta0]=deal(theta0,theta1);
else
    flag=1;
end

% case 1. full circle
if (theta1-theta0) == 2*pi
    R = [cos(theta0) -sin(theta0); sin(theta0) cos(theta0)];
    t=linspace(theta0,theta1,9); t=t';
    PP=r*[1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1;  1 0]';
    PP=(R*PP)';
    P=bsxfun(@plus,center,PP);
    knots=[0 0 0 .25 .25 .5 .5 .75 .75 1 1 1];
    c=1/sqrt(2); w=[1 c 1 c 1 c 1 c 1];
    order=3;

    if flag == 0
        P=flipud(P);
    end

    return;
end

% Case 2. Wide at least half a circle but not a circle.
if (theta1-theta0) >= pi
    R = [cos(theta0) -sin(theta0); sin(theta0) cos(theta0)];
    PP=r*[1 0; 1 1; 0 1; -1 1; -1 0]; PP=(R*PP')';
    P=bsxfun(@plus,center,PP);

    c=1/sqrt(2); w=[1 c 1 c 1];
    theta0=theta0+pi;
else
    P=[]; w=[];
end

% Case 3. Adding some arcs.
if (theta1-theta0) > 0
    % 1. compute initial control point.
    P0=center+r*[cos(theta0) sin(theta0)];

    % 2. compute final control point.
    P2=center+r*[cos(theta1) sin(theta1)];

    % 3. compute intermediate control point.
    bL=(theta1-theta0)/2; % t in [-bL,bL]
    r1=r/cos(bL);
    P1=center+r1*[cos((theta0+theta1)/2) sin((theta0+theta1)/2)];

    if length(P) > 0 % wider than half a circle
        P=[P; P1; P2];
        c=cos(bL); w=[w c 1];
        knots=[0 0 0 .25 .25 .75 .75 1 1 1];
    else % smaller than half a circle
        P=[P0; P1; P2];
        c=cos(bL); w=[1 c 1];
        knots=[0 0 0 1 1 1];
    end
else % theta1=theta0
    if length(w) > 0 % half circle
        knots=[0 0 0 .5 .5 1 1 1];
    else
        P0=[]; P1=[]; P2=[];
        knots=[]; w=[];
    end
end


if flag == 0
    w=fliplr(w);
    P=flipud(P);
end


%--------------------------------------------------------------------------
% Knots and order.
%--------------------------------------------------------------------------
order=3;












function [P,w,knots,order]=compute_ellarc_parms(center,angles,ell_axis,...
    tilt_angle)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine defines an arc of ellipse (possibly tilted).
% Ellipse equation (without tilting):
% * x(t)=center(1)+a*cos(t)
% * y(t)=center(2)+b*sin(t)
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% center: center of the ellipse
% angles: angles of the extrema of the elliptical arc.
% ell_axis: it is a 2 x 1 vector, where the first component describes the
%         half width "a" while the second component described the half
%         height "b".
% tilt_angle: angle of rotation of the ellipse with respect to the x-axis;
%         this variable is not mandatory (default value: 0).
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% P: control points.
% w: NURBS weights.
% knots: NURBS weights.
% order: NURBS order.
%--------------------------------------------------------------------------
% DATE:
%--------------------------------------------------------------------------
% Written: November 4, 2021.
% Modified: November 18, 2021.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Troubleshooting
%--------------------------------------------------------------------------

if nargin < 1, center=[0 0]; end
if nargin < 2, angles=[0 2*pi]; end
if nargin < 3, ell_axis=[1 1]; end
if nargin < 4, tilt_angle=0; end
if nargin < 5, Npts=100; end

if isempty(center), center=[0 0]; end
if isempty(angles), angles=[0 2*pi]; end
if isempty(ell_axis), ell_axis=[1 1]; end
if isempty(tilt_angle), tilt_angle=0; end
if isempty(Npts), Npts=100; end

theta0=angles(1); theta1=angles(2);

if theta0 > theta1
    flag=0; [theta1,theta0]=deal(theta0,theta1);
else
    flag=1;
end

%--------------------------------------------------------------------------
% Defining local variables
%--------------------------------------------------------------------------

theta0=angles(1);
theta1=angles(2);
a=ell_axis(1);
b=ell_axis(2);
r=1;

%--------------------------------------------------------------------------
% Control Points
%--------------------------------------------------------------------------

% Basis set via circle.
center0=[0 0]; r0=1;
[P,w,knots,order]=compute_arc_parms(center0,r0,theta0,theta1);

% Adapting to ellipse.
P=[a*P(:,1) b*P(:,2)];

% Tilting
if not(tilt_angle == 0)
    R=[cos(tilt_angle) -sin(tilt_angle); sin(tilt_angle) cos(tilt_angle)];
    P=P*R';
end

% Shifting to center.
P=bsxfun(@plus,center,P);












function [P,w,knots,order]=compute_polygon_parms(vertices)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine defines the sides of a polygonal arc (possibly open!).
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% vertices: vertices of the arcs of the polygon (oriented counterclockwise)
%           as N x 2, where "N" is the number of vertices, that are
%           described in cartesian coordinates.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% P: control points.
% w: NURBS weights.
% knots: NURBS weights.
% order: NURBS order.
%--------------------------------------------------------------------------
% DATE:
%--------------------------------------------------------------------------
% Written: November 4, 2021.
%--------------------------------------------------------------------------

L=size(vertices,1);
P=[]; order=3;

for k=1:L-1
    P1=vertices(k,:);
    P2=vertices(k+1,:);
    [PL,wL]=compute_segments_parms([P1; P2]);
    P=[P; PL];
end

knots=repmat(0:L-1,order,1); knots=(knots(:))'/(L-1);
w=ones(1,size(P,1));












function [P,w,knots,order]=compute_segments_parms(vertices)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This routine defines a segment with vertices, via NURBS.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% center=varargin{2};: extrema of the segment joining P0 to P2.
% It is a 2 x 2 matrix where the k-th row describes the k-th point.
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% P: control points.
% w: NURBS weights.
% knots: NURBS weights.
% order: NURBS order.
%--------------------------------------------------------------------------
% DATE:
%--------------------------------------------------------------------------
% Written: November 7, 2021.
%--------------------------------------------------------------------------

P0=vertices(1,:); P2=vertices(2,:);

method_segment=2;

switch method_segment
    case 1

        P=[P0; P2];
        w=[1 1];
        knots=[0 0 1 1];
        order=2;

    otherwise
        %------------------------------------------------------------------
        % Control Points
        %------------------------------------------------------------------

        % 1. compute intermediate control point.
        P1=(P0+P2)/2;

        % 2. Determine set of control points.
        P=[P0; P1; P2];

        %------------------------------------------------------------------
        % Weights.
        %------------------------------------------------------------------
        w=[1 1 1];

        %------------------------------------------------------------------
        % Knots and order.
        %------------------------------------------------------------------

        knots=[0 0 0 1 1 1];
        order=3;

end



