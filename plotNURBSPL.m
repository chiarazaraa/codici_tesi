
function plotNURBSPL(structure_NURBSPL,plot_mode,Nbox)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Plot the boundary described by the array of structures stored in
% "structure_NURBSPL".
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% structure_NURBSPL: array of structures that describe the closed NURBS or
%    composite Bezier curve of splines that define parametrically the
%    boundary of the domain in analysis.
% plot_mode: if 'plot_mode' is 'dynamic' then it shows a point moving on the
%    boundary and it is useful to see if a point move counterclockwise.
%--------------------------------------------------------------------------
% Dates
%--------------------------------------------------------------------------
% First version: November 13, 2021;
% Checked: November 30, 2021.
%-------------------------------------------------------------------------

if nargin < 2, plot_mode='normal'; end
if nargin < 3, Nbox=1; end



% ...........................  domain plot ................................

if length(plot_mode) == 0, plot_mode='normal'; end

Npieces=length(structure_NURBSPL);
XV=[]; YV=[];
for k=1:Npieces
    geometry_NURBSPL=structure_NURBSPL(k);
    PL=geometry_NURBSPL.P;
    wL=geometry_NURBSPL.w;
    orderL=geometry_NURBSPL.order;
    knotsL=geometry_NURBSPL.knots;

    t0=knotsL(1); t1=knotsL(end); tV=linspace(t0,t1,10000); tV=tV';
    [SxN,SxD,SyN,SyD]=NURBSPL2ratsplines(geometry_NURBSPL,'pp');
    X=fnval(SxN,tV)./fnval(SxD,tV); Y=fnval(SyN,tV)./fnval(SyD,tV);
    XV=[XV; X]; YV=[YV; Y];
end


% ..................... dynamic domain plot ...............................

if strcmp(plot_mode,'dynamic')
    hold on;
    LL=length(XV);
    pause_time=0.25/LL;
    for j=1:LL
        plot(XV(j),YV(j),'go','MarkerEdgeColor','g',...
            'MarkerFaceColor','g',...
            'MarkerSize',6);
        pause(pause_time);
    end
    return;
end



% ..................... boxes domain plot .................................
if strcmp(plot_mode,'boxes')

    axis off, axis tight, axis equal;
    hold on;

    % plot boxes
    [inside,in_doubts,sp_xnV,sp_xdV,sp_ynV,sp_ydV,boxVx,singular_points_X,...
        singular_points_Y]=inRS([],structure_NURBSPL,Nbox);
    Nboxes=size(boxVx,1);

    for k=1:Nboxes
        a=boxVx(k,1); b=boxVx(k,2); c=boxVx(k,3); d=boxVx(k,4);
        C = [220 220 220]/256;
        fill([a b b a a], [c c d d c],C);
        hold on;
    end


end


% ..................... normal domain plot ................................
% plot domain boundary
axis equal
plot(XV,YV,'k-','MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',1,'LineWidth',1);








