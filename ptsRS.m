
function [pts,boxXY]=ptsRS(structure_RS,N,boxXY,method)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Given a structure of a NURBS, composite Bezier curve or parametric
% spline, this routine determines a suitable mesh of tensorial points, on
% equispaced meshes along the directions.
% The points are in a bounding box of the domain.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% structure_RS: structure of the boundary.
% N: determine the points of the rule depending by the string in "method".
% boxXY: matrix relative to the monotone boxes (not mandatory).
% method: if 'grid' then N^2 points of an equispaced tensor grid are
%         generated in the bounding box of "structure_RS";
%         if 'rand' then N random points are generated in the bounding box
%         of "structure_RS";
%         if 'halton' then N Halton points are generated in the bounding
%         box of "structure_RS";
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% pts: matrix N x 2 describing the cartesian coordinates of the pointset,
%      with N=N^2;
% boxXY: matrix relative to the monotone boxes.
%--------------------------------------------------------------------------
% Dates
%--------------------------------------------------------------------------
% First version: November 13, 2021;
% Checked: November 16, 2021.
%--------------------------------------------------------------------------

if nargin < 2, N =100; end
if nargin < 3, boxXY=[]; end
if nargin < 4, method='grid'; end

if isempty(boxXY), [~,~,~,~,~,~,boxXY,~,~]=inRS([],structure_RS,1,0); end
P=boxXY(:,1:4);

aa=min(P(:,1)); bb=max(P(:,2)); cc=min(P(:,3)); dd=max(P(:,4));

if strcmp(method,'grid')
    s=linspace(aa,bb,N); t=linspace(cc,dd,N);
    [X,Y]=meshgrid(s,t);
    pts=[X(:) Y(:)];
    return;
end


if strcmp(method,'rand')
    X=aa+(bb-aa)*rand(N,1); Y=cc+(dd-cc)*rand(N,1);
    pts=[X Y];
    return
end


if strcmp(method,'halton')
    p=haltonset(2);
    X=aa+(bb-aa)*p(1:N,1);
    Y=cc+(dd-cc)*p(1:N,2);
    pts=[X Y];
    return
end