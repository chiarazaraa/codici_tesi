function [xyw,xywc,P,L,subs] = polygcirc3_000_L(n,v,a,b,cc,r,conv,pos)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% Computation of a basic and a compressed positive cubature formula
% on a polygonal element with a circular side, namely the set
% union (convex arc) or difference (concave arc) of a convex
% polygonal element with a circular segment.
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% n: polynomial degree of exactness
% a,b: extrema of the circular arc
% cc,r: circular arc center and radius
% v: polygon vertices (2-column array of coords)
%    note that to "v" are added by default the arc extrema, in
%    counterclockwise order
% conv: conv=1 for a convex arc, conv=0 for a concave arc
%
% WARNING: the figure vertices are a,v(1,:),...,v(end,:),b and MUST BE
% in COUNTERCLOCKWISE order
% (the arc ba is clockwise on the circle if concave and counterclockwise if convex)
%--------------------------------------------------------------------------
% OUTPUT:
%--------------------------------------------------------------------------
% xyw: 3-column array of cubature nodes and positive weights
% xywc: 3-column array of compressed cubature nodes and positive weights
%--------------------------------------------------------------------------
% AUTHORS:
%--------------------------------------------------------------------------
% Authors: E. Artioli, A. Sommariva and M. Vianello
% Written: April 24, 2018
% Revised: December 02, 2021
%--------------------------------------------------------------------------

if nargin < 8, pos=1; end

L=[];
P=[];
subs=[];

if conv==1
    
    P=[a;v;b;a];
    xyw=polygauss_2013(n,P);
    nw=circtrap(n,b,a,b,a,cc,r);
    xyw=[xyw;nw];
    L=[L; size(xyw,1)];
    A=b;
            B=a;
            C=b;
            D=a;
            subs=[subs; A; B; C; D];
    
else %conv==0
    
    k=0;l=0;t=0;
    
    aa=a-cc;bb=b-cc;
    angleab=acos((aa*bb')/(norm(aa)*norm(bb)));
    anglea=angle(aa(1)+sqrt(-1)*aa(2));
    angleb=angle(bb(1)+sqrt(-1)*bb(2));
    
    if angleb<=pi
        clockba=(anglea>=angleb & anglea<angleb+pi);
    else
        clockba=(anglea>angleb-pi & anglea<=angleb);
    end;
    
    if length(v(:,1))==2
        v=[v(1,:);(v(1,:)+v(2,:))/2;v(2,:)];
    end;
    
    vv=v-repmat(cc,length(v(:,1)),1);
    
    for i=1:length(v(:,1))
        angle1=acos((vv(i,:)*aa')/(norm(aa)*norm(vv(i,:))));
        angle2=acos((vv(i,:)*bb')/(norm(bb)*norm(vv(i,:))));
        
        if angle1<angle2 & angle2>angleab & clockba==0
            e(k+1,:)=v(i,:);k=k+1;
        end;
        
        if clockba==1 | (angle1<=angleab & angle2<=angleab)
            z(l+1,:)=v(i,:);
            zz=z(l+1,:)-cc;
            theta=angle(zz(1)+sqrt(-1)*zz(2));
            u(l+1,:)=cc+r*[cos(theta) sin(theta)];
            l=l+1;
        end;
        
        if angle2<angle1 & angle1>angleab & clockba==0
            eta(t+1,:)=v(i,:);t=t+1;
        end;
        
    end;
    
    xyw=[];
    
    if k>1
        P1=[a;e;a];
        xyw=polygauss_2013(n,P1);
        L=[L; size(xyw,1)];
    end;
    
    if k==0 & l>=1
        if norm(u(1,:)-a)>10^(-14)
            nw=circtrap(n,a,u(1,:),z(1,:),z(1,:),cc,r);
            L=[L; size(nw,1)];
            xyw=[xyw;nw];
            A=a;
            B=(u(1,:));
            C=(z(1,:));
            D=(z(1,:));
            subs=[subs; A; B; C; D];
        end;
    end;
    
    if k>=1 & l>=1
        nw=circtrap(n,a,u(1,:),e(k,:),z(1,:),cc,r);
        L=[L; size(nw,1)];
        xyw=[xyw;nw];
        A=a;
            B=(u(1,:));
            C=(e(k,:));
            D=(z(1,:));
            subs=[subs; A; B; C; D];
    end;
    
    for j=1:l-1
        nw=circtrap(n,u(j,:),u(j+1,:),z(j,:),z(j+1,:),cc,r);
        L=[L; size(nw,1)];
        xyw=[xyw;nw];
            A=(u(j,:));
            B=(u(j+1,:));
            C=(z(j,:));
            D=(z(j+1,:));
            subs=[subs; A; B; C; D];
    end;
    
    if t==0 & l>=1
        if norm(u(l,:)-b)>10^(-14)
            nw=circtrap(n,u(l,:),b,z(l,:),z(l,:),cc,r);
            L=[L; size(nw,1)];
            xyw=[xyw;nw];
            A=(u(l,:));
            B=(b);
            C=(z(l,:));
            D=(z(l,:));
            subs=[subs; A; B; C; D];
        end;
    end;
    
    if t>=1 & l>=1
        nw=circtrap(n,u(l,:),b,z(l,:),eta(1,:),cc,r);
        L=[L; size(nw,1)];
        xyw=[xyw;nw];
        A=(u(l,:));
            B=(b);
            C=(z(l,:));
            D=(eta(1,:));
            subs=[subs; A; B; C; D];
    end;
    
    if t>1
        P2=[b;eta;b];
        nw=polygauss_2013(n,P2);
        L=[L; size(nw,1)];
        xyw=[xyw;nw];
    end;
    
    if l==0
        
        if k>=1 & t>=1
            nw=circtrap(n,a,b,e(k,:),eta(1,:),cc,r);
            L=[L; size(nw,1)];
            A=(a);
            B=(b);
            C=(e(k,:));
            D=(eta(1,:));
            subs=[subs; A; B; C; D];
        end;
        
        if k==0 & t>=1
            nw=circtrap(n,a,b,eta(1,:),eta(1,:),cc,r);
            L=[L; size(nw,1)];
            A=(a);
            B=(b);
            C=(eta(1,:));
            D=(eta(1,:));
            subs=[subs; A; B; C; D];
        end;
        
        if k>=1 & t==0
            nw=circtrap(n,a,b,e(k,:),e(k,:),cc,r);
            L=[L; size(nw,1)];
            A=(a);
            B=(b);
            C=(e(k,:));
            D=(e(k,:));
            subs=[subs; A; B; C; D];
        end;
        
        xyw=[xyw;nw];
        
    end;
    
end;

if nargout >= 2
    [pts,w,momerr] = comprexcub(n,[xyw(:,1) xyw(:,2)],xyw(:,3),pos);
    xywc=[pts(:,1) pts(:,2) w];

else
    xywc=[];
end

end

%%%%%%%%%%%%%%%%%%%%

function xyw = circtrap(n,A,B,C,D,CC,r)
% by E. Artioli, A. Sommariva, M. Vianello 
% April 2018

% n: polynomial exactness degree
% A,B: circular arc extrema coords 1 x 2
% C,D: base segment extrema coords 1 x 2, AC and BD are the sides
% CC: circle center coords 1 x 2
% r: circle radius

% xyw: 3-column array xyw(:,1:2) nodes, xyw(:,3) weights



Z=(A(1)-CC(1))+sqrt(-1)*(A(2)-CC(2));
W=(B(1)-CC(1))+sqrt(-1)*(B(2)-CC(2));

az=angle(Z);aw=angle(W);

if aw>=az
    if aw-az<=pi
        alpha=az;beta=aw;U=C;V=D;
    else
        alpha=aw;beta=2*pi+az;U=D;V=C;
    end;
end;

if az>aw
    if az-aw<=pi
        alpha=aw;beta=az;U=D;V=C;
    else
        alpha=az;beta=aw+2*pi;U=C;V=D;
    end;
end;

om=(beta-alpha)/2;g=(beta+alpha)/2;
s=2*sin(om);

A1=[r*cos(g) r*sin(g);0 0];
B1=[-r*sin(g) r*cos(g);(V(1)-U(1))/s (V(2)-U(2))/s];
C1=[CC(1) CC(2);(V(1)+U(1))/2 (V(2)+U(2))/2];

xyw=gqellblend(n,A1,B1,C1,-om,om);

end

%%%%%%%%%%%%%%%%%%%%

function xyw=polygauss_2013(N,polygon_sides,rotation,P,Q)

%--------------------------------------------------------------------------
% REFERENCE PAPER:
% [1] A. SOMMARIVA and M. VIANELLO
% "Gauss-like and triangulation-free cubature over polygons".
%
% INPUT:
%
% N     : DEGREE OF THE 1 DIMENSIONAL GAUSS-LEGENDRE RULE.
%
% polygon_sides: IF THE POLYGON HAS "L" SIDES, "boundary.pts" IS A
%         VARIABLE CONTAINING ITS VERTICES, ORDERED COUNTERCLOCKWISE.
%         AS LAST ROW MUST HAVE THE COMPONENTS OF THE FIRST VERTEX.
%         IN OTHER WORDS, THE FIRST ROW AND LAST ROW ARE EQUAL.
%         "polygon_sides" IS A "L+1 x 2" MATRIX.
%
%            --------- NOT MANDATORY VARIABLES ---------
%
% rotation: 0: NO ROTATION.
%           1: AUTOMATIC.
%           2: PREFERRED DIRECTION ROTATION BY P, Q.
%
% P, Q: DIRECTION THAT FIXES THE ROTATION.
%
% OUTPUT:
%
% xyw     : THE GAUSS LIKE FORMULA PRODUCES THE NODES (xyw(:,1),xyw(:,2))
%           AND THE WEIGHTS xyw(:,3) OF A CUBATURE RULE ON THE POLYGON.
%
%--------------------------------------------------------------------------
% EXAMPLE 1 (NO ROTATION.)
%---------------------------
%
% >> xyw=polygauss_2013(2,[0 0; 1 0; 1 1; 0 1; 0 0],0)
%
% xyw =
%
%     0.2113    0.2113    0.2500
%     0.2113    0.7887    0.2500
%     0.7887    0.2113    0.2500
%     0.7887    0.7887    0.2500
%
% >>
%
%--------------------------------------------------------------------------
% EXAMPLE 2 (AUTO ROTATION.)
%-----------------------------
%
% >> xyw=polygauss_2013(2,[0 0; 1 0; 1 1; 0 1; 0 0])
%
% xyw =
%
%     0.0683    0.0444    0.0078
%     0.3028    0.1972    0.0556
%     0.5374    0.3499    0.0616
%     0.6501    0.4626    0.0616
%     0.8028    0.6972    0.0556
%     0.9556    0.9317    0.0078
%     0.9317    0.9556    0.0078
%     0.6972    0.8028    0.0556
%     0.4626    0.6501    0.0616
%     0.3499    0.5374    0.0616
%     0.1972    0.3028    0.0556
%     0.0444    0.0683    0.0078
%     0.1008    0.0119    0.0078
%     0.4472    0.0528    0.0556
%     0.7935    0.0938    0.0616
%     0.9062    0.2065    0.0616
%     0.9472    0.5528    0.0556
%     0.9881    0.8992    0.0078
%     0.8992    0.9881    0.0078
%     0.5528    0.9472    0.0556
%     0.2065    0.9062    0.0616
%     0.0938    0.7935    0.0616
%     0.0528    0.4472    0.0556
%     0.0119    0.1008    0.0078
%
% >>

%--------------------------------------------------------------------------
%% Copyright (C) 2007-2013 Marco Vianello and Alvise Sommariva
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

%% Authors:
%% Marco Vianello    <marcov@euler.math.unipd.it>
%% Alvise Sommariva  <alvise@euler.math.unipd.it>
%% Date: April 30, 2013.
%--------------------------------------------------------------------------

%----------------------------------------------------------------------
% BOUNDARY PTS.
%----------------------------------------------------------------------
x_bd=polygon_sides(:,1);
y_bd=polygon_sides(:,2);

%----------------------------------------------------------------------
% "MINIMUM" RECTANGLE CONTAINING POLYGON.
%----------------------------------------------------------------------
x_min=min(x_bd); x_max=max(x_bd);
y_min=min(y_bd); y_max=max(y_bd);

%----------------------------------------------------------------------
% SOME AUTOMATIC SETTINGS.
%----------------------------------------------------------------------
if nargin < 3
    rotation=1;
end
cubature_type=4;

%--------------------------------------------------------------------------
% POLYGON ROTATION (IF NECESSARY).
%--------------------------------------------------------------------------
switch rotation
    case 0
        %         fprintf('\n \t [ROTATION]: NO.');
        rot_matrix=eye(2);
        axis_abscissa=[x_min y_max]-[x_min y_min];
    case 1
        %         fprintf('\n \t [ROTATION]: AUTOMATIC');
        [polygon_sides,rot_matrix,rot_angle,axis_abscissa,P,Q]=...
            auto_rotation(polygon_sides,[],[]);
        %         fprintf(' [ANGLE CLOCKWISE (RESPECT Y, IN DEGREES)]: %5.5f',...
        %             rot_angle*180/pi);
    case 2
        %         fprintf('\n \t [ROTATION]: PREFERRED DIRECTION');
        nrm_vect=norm(Q-P);
        if (nrm_vect > 0)
            direction_axis=(Q-P)/nrm_vect;
            [polygon_sides,rot_matrix,rot_angle,axis_abscissa,P,Q]=...
                auto_rotation(polygon_sides,P,Q);
        else
            %             fprintf('\n \t [WARNING]: THE DIRECTION VECTOR IS NULL. ')
            %             fprintf('USING AUTOMATIC ROTATION.');
            [polygon_sides,rot_matrix,rot_angle,axis_abscissa,P,Q]=...
                auto_rotation(polygon_sides,P,Q);
        end
        %         fprintf(' [ANGLE CLOCKWISE (RESPECT Y)]:
        %5.5f',rot_angle*180/pi);
end

%--------------------------------------------------------------------------
% COMPUTE NODES AND WEIGHTS OF 1D GAUSS-LEGENDRE RULE.
% TAKEN FROM TREFETHEN PAPER "Is ... Clenshaw-Curtis?".
%--------------------------------------------------------------------------

% DEGREE "N".
[s_N,w_N]=cubature_rules_1D((N-1),cubature_type);
N_length=length(s_N);

% DEGREE "M".
M=N+1;
[s_M,w_M]=cubature_rules_1D((M-1),cubature_type);

%----------------------------------------------------------------------
% L: NUMBER OF SIDES OF THE POLYGON.
% M: ORDER GAUSS INTEGRATION.
% N: ORDER GAUSS PRIMITIVE.
%----------------------------------------------------------------------
L=length(polygon_sides(:,1))-1;

%a=0.5;
a=axis_abscissa(1);

%----------------------------------------------------------------------
% COMPUTE 2D NODES (nodes_x,nodes_y) AND WEIGHTS "weights".
%----------------------------------------------------------------------

nodes_x=[];
nodes_y=[];
weights=[];

for index_side=1:L
    x1=polygon_sides(index_side,1); x2=polygon_sides(index_side+1,1);
    y1=polygon_sides(index_side,2); y2=polygon_sides(index_side+1,2);
    if ~(x1 == a & x2 == a)
        if (y2-y1) ~=0
            
            if (x2-x1) ~=0
                s_M_loc=s_M;
                w_M_loc=w_M;
            else
                s_M_loc=s_N;
                w_M_loc=w_N;
            end
            
            M_length=length(s_M_loc);
            
            half_pt_x=(x1+x2)/2; half_pt_y=(y1+y2)/2;
            half_length_x=(x2-x1)/2; half_length_y=(y2-y1)/2;
            
            
            % GAUSSIAN POINTS ON THE SIDE.
            x_gauss_side=half_pt_x+half_length_x*s_M_loc; %SIZE: (M_loc,1)
            y_gauss_side=half_pt_y+half_length_y*s_M_loc; %SIZE: (M_loc,1)
            
            scaling_fact_plus=(x_gauss_side+a)/2; %SIZE: (M_loc,1)
            scaling_fact_minus=(x_gauss_side-a)/2;%SIZE: (M_loc,1)
            
            local_weights=...
                (half_length_y*scaling_fact_minus).*w_M_loc;%SIZE:(M_loc,1)
            
            term_1=repmat(scaling_fact_plus,1,N_length); % SIZE: (M_loc,N)
            
            term_2=repmat(scaling_fact_minus,1,N_length); % SIZE: (M_loc,N)
            
            rep_s_N=repmat(s_N',M_length,1);
            
            % x, y ARE STORED IN MATRICES. A COUPLE WITH THE SAME INDEX
            % IS A POINT, i.e. "P_i=(x(k),y(k))" FOR SOME "k".
            x=term_1+term_2.*rep_s_N;
            y=repmat(y_gauss_side,1,N_length);
            
            number_rows=size(x,1);
            number_cols=size(x,2);
            
            x=x(:); x=x';
            y=y(:); y=y';
            
            rot_gauss_pts=rot_matrix'*[x;y]; % THE INVERSE OF A ROTATION
            % MATRIX IS ITS TRANSPOSE.
            
            x_rot=rot_gauss_pts(1,:); % GAUSS POINTS IN THE ORIGINAL SYSTEM.
            y_rot=rot_gauss_pts(2,:);
            
            x_rot=reshape(x_rot',number_rows,number_cols);
            y_rot=reshape(y_rot',number_rows,number_cols);
            
            nodes_x=[nodes_x; x_rot];
            nodes_y=[nodes_y; y_rot];
            weights=[weights; local_weights];
            
            
        end
    end
end

weights=weights*w_N';
weights=weights(:);

nodes_x=nodes_x(:);
nodes_y=nodes_y(:);
xyw=[nodes_x nodes_y weights];

end




%----------------------------------------------------------------------
% FUNCTIONS USED IN THE ALGORITHM.
%----------------------------------------------------------------------


%----------------------------------------------------------------------
% 1. "auto_rotation"
%----------------------------------------------------------------------
function [polygon_bd_rot,rot_matrix,rot_angle,axis_abscissa,vertex_1,vertex_2]=...
    auto_rotation(polygon_bd,vertex_1,vertex_2)


% AUTOMATIC ROTATION OF A CONVEX POLYGON SO THAT "GAUSSIAN POINTS",
% AS IN THE PAPER THEY ARE ALL CONTAINED IN THE CONVEX POLYGON.
% SEE THE PAPER FOR DETAILS.


% FIND DIRECTION AND ROTATION ANGLE.
if length(vertex_1) == 0
    % COMPUTING ALL THE DISTANCES BETWEEN POINTS.A LITTLE TIME CONSUMING
    % AS PROCEDURE.
    distances = points2distances(polygon_bd);
    [max_distances,max_col_comp]=max(distances,[],2);
    [max_distance,max_row_comp]=max(max_distances,[],1);
    vertex_1=polygon_bd(max_col_comp(max_row_comp),:);
    vertex_2=polygon_bd(max_row_comp,:);
    direction_axis=(vertex_2-vertex_1)/max_distance;
else
    direction_axis=(vertex_2-vertex_1)/norm(vertex_2-vertex_1);
end

rot_angle_x=acos(direction_axis(1));
rot_angle_y=acos(direction_axis(2));

if rot_angle_y <= pi/2
    if rot_angle_x <= pi/2
        rot_angle=-rot_angle_y;
    else
        rot_angle=rot_angle_y;
    end
else
    if rot_angle_x <= pi/2
        rot_angle=pi-rot_angle_y;
    else
        rot_angle=rot_angle_y;
    end
end


% CLOCKWISE ROTATION.
rot_matrix=[cos(rot_angle) sin(rot_angle);
    -sin(rot_angle) cos(rot_angle)];

number_sides=size(polygon_bd,1)-1;

polygon_bd_rot=(rot_matrix*polygon_bd')';

axis_abscissa=rot_matrix*vertex_1';

end

%----------------------------------------------------------------------
% 3. "cubature_rules_1D"
%----------------------------------------------------------------------

function [nodes,weights]=cubature_rules_1D(n,cubature_type)

% SEE WALDVOGEL PAPER. ADDED NODES

% Weights of the Fejer2, Clenshaw-Curtis and Fejer1 quadrature by DFTs
% n>1. Nodes: x_k = cos(k*pi/n)

N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';

switch cubature_type
    
    case 1 % FEJER 1.
        v0=[2*exp(i*pi*K/n)./(1-4*K.^2); zeros(l+1,1)];
        v1=v0(1:end-1)+conj(v0(end:-1:2));
        weights=ifft(v1);
        k=(1/2):(n-(1/2)); nodes=(cos(k*pi/n))';
        
    case 2 % FEJER 2.
        v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
        v2=-v0(1:end-1)-v0(end:-1:2);
        wf2=ifft(v2); weights=[wf2;0];
        k=0:n; nodes=(cos(k*pi/n))';
        
    case 3 % CLENSHAW CURTIS.
        g0=-ones(n,1); g0(1+l)=g0(1+l)+n; g0(1+m)=g0(1+m)+n;
        g=g0/(n^2-1+mod(n,2));
        v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
        v2=-v0(1:end-1)-v0(end:-1:2);
        wcc=ifft(v2+g); weights=[wcc;wcc(1,1)];
        k=0:n; nodes=(cos(k*pi/n))';
        
    case 4 % GAUSS LEGENDRE
        beta=0.5./sqrt(1-(2*(1:n)).^(-2));
        T=diag(beta,1)+diag(beta,-1);
        [V,D]=eig(T);
        x=diag(D); [x,index]=sort(x); x=x';
        w=2*V(1,index).^2;
        nodes=x';
        weights=w';
        
end

end

%----------------------------------------------------------------------
% 3. "points2distances"
%----------------------------------------------------------------------

function distances = points2distances(points)

% Create euclidean distance matrix from point matrix.

% Get dimensions.
[numpoints,dim]=size(points);

% All inner products between points.
distances=points*points';

% Vector of squares of norms of points.
lsq=diag(distances);

% Distance matrix.
distances=sqrt(repmat(lsq,1,numpoints)+repmat(lsq,1,numpoints)'-2*distances);

end

%%%%%%%%%%%%%%%%%%%%

function xyw = gqellblend(n,A,B,C,alpha,beta)
% by Gaspare Da Fies, Alvise Sommariva and Marco Vianello
% 
% 2 June 2011

% computes the nodes and weights of a product gaussian formula
% exact on total-degree bivariate polynomials of degree <=n
% on the planar region R obtained by linear blending (convex combination)
% of two trigonometric arcs with parametric equations
% P(theta)=A1*cos(theta)+B1*sin(theta)+C1
% Q(theta)=A2*cos(theta)+B2*sin(theta)+C2
% namely
% R = {(x,y)=t*P(theta)+(1-t)*Q(theta), t in [0,1], theta in [alpha,beta],
% 0<beta-alpha<=2*pi}

% uses the routines:
%
% r_jacobi.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
%
% trigauss.m
% http://www.math.unipd.it/~marcov/mysoft/trigauss.m
% this will be soon substituted by an optimized version

% input:
% n: algebraic degree of exactness
% A,B,C: 2x2 matrices of the parametric arc coefficients:
% A1=A(1,:), B1=B(1,:), C1=C(1,:)
% A2=A(2,:), B2=B(2,:), C2=C(2,:)
% [alpha,beta]: angular interval, 0<beta-alpha<=2*pi

% output:
% xyw: 3 columns array of (xnodes,ynodes,weights)


% computing the algebraic and trigonometric degree increment
S1=abs((A(1,1)-A(2,1))*(B(1,2)-B(2,2))+(A(1,2)-A(2,2))*(B(2,1)-B(1,1)))>10*eps;
S2=abs((C(1,1)-C(2,1))*(B(1,2)-B(2,2))+(C(1,2)-C(2,2))*(B(2,1)-B(1,1)))>10*eps;
S3=abs((A(1,1)-A(2,1))*(C(1,2)-C(2,2))+(A(1,2)-A(2,2))*(C(2,1)-C(1,1)))>10*eps;

if (S1 || S2 || S3)
    h=1;
else
    h=0;
end

S4=abs(A(1,2)*A(2,1)-A(1,1)*A(2,2)-B(1,2)*B(2,1)+B(1,1)*B(2,2))>10*eps;
S5=abs(A(1,2)*B(2,1)-A(1,1)*B(2,2)+B(1,2)*A(2,1)-B(1,1)*A(2,2))>10*eps;
S6=abs(B(2,1)*(C(1,2)-C(2,2))-B(2,2)*(C(1,1)-C(2,1)))>10*eps;
S7=abs(A(2,1)*(C(1,2)-C(2,2))-A(2,2)*(C(1,1)-C(2,1)))>10*eps;
S8=abs((C(1,1)-C(2,1))*(B(1,2)-B(2,2))+(C(1,2)-C(2,2))*(B(2,1)-B(1,1)))>10*eps;
S9=abs((A(1,1)-A(2,1))*(C(1,2)-C(2,2))+(A(1,2)-A(2,2))*(C(2,1)-C(1,1)))>10*eps;

if (S4 || S5)
    k=2;
elseif (S6 || S7 || S8 || S9)
    k=1;
else
    k=0;
end

% trigonometric gaussian formula on the arc
tw=trigauss(n+k,alpha,beta);

% algebraic gaussian formula on [0,1]
ab=r_jacobi(ceil((n+h+1)/2),0,0);
xw=gauss(ceil((n+h+1)/2),ab);
xw(:,1)=xw(:,1)/2+1/2;
xw(:,2)=xw(:,2)/2;

% creating the grid
[t,theta]=meshgrid(xw(:,1),tw(:,1));
[w1,w2]=meshgrid(xw(:,2),tw(:,2));

% nodal cartesian coordinates and weights
s=sin(theta(:));
c=cos(theta(:));
p1=A(1,1)*c+B(1,1)*s+C(1,1);
p2=A(1,2)*c+B(1,2)*s+C(1,2);
q1=A(2,1)*c+B(2,1)*s+C(2,1);
q2=A(2,2)*c+B(2,2)*s+C(2,2);
dp1=-A(1,1)*s+B(1,1)*c;% plot(xyw(:,1)+trasl,xyw(:,2),'.','MarkerSize',6);
dp2=-A(1,2)*s+B(1,2)*c;
dq1=-A(2,1)*s+B(2,1)*c;
dq2=-A(2,2)*s+B(2,2)*c;

xyw(:,1)=p1.*t(:)+q1.*(1-t(:));
xyw(:,2)=p2.*t(:)+q2.*(1-t(:));
xyw(:,3)=abs((p1-q1).*(dp2.*t(:)+dq2.*(1-t(:))) - ...
    (p2-q2).*(dp1.*t(:)+dq1.*(1-t(:)))).* w1(:).*w2(:);


% this part plots the region R and the cubature nodes
% plot(xyw(:,1),xyw(:,2),'.','MarkerSize',6);
% axis equal;
% hold on;
% theta=linspace(alpha,beta,500);
% s=sin(theta);
% c=cos(theta);
% p1=A(1,1)*c+B(1,1)*s+C(1,1);
% p2=A(1,2)*c+B(1,2)*s+C(1,2);
% plot(p1,p2);
% q1=A(2,1)*c+B(2,1)*s+C(2,1);
% q2=A(2,2)*c+B(2,2)*s+C(2,2);
% plot(q1,q2);
% t=linspace(0,1,100);
% seg11=(A(1,1)*cos(alpha)+B(1,1)*sin(alpha)+C(1,1))*t + ...
% (A(2,1)*cos(alpha)+B(2,1)*sin(alpha)+C(2,1))*(1-t);
% seg12=(A(1,2)*cos(alpha)+B(1,2)*sin(alpha)+C(1,2))*t + ...
% (A(2,2)*cos(alpha)+B(2,2)*sin(alpha)+C(2,2))*(1-t);
% plot(seg11,seg12);
% seg21=(A(1,1)*cos(beta)+B(1,1)*sin(beta)+C(1,1))*t + ...
% (A(2,1)*cos(beta)+B(2,1)*sin(beta)+C(2,1))*(1-t);
% seg22=(A(1,2)*cos(beta)+B(1,2)*sin(beta)+C(1,2))*t + ...
% (A(2,2)*cos(beta)+B(2,2)*sin(beta)+C(2,2))*(1-t);
% plot(seg21,seg22);

end

%%%%%%%%%%%%%%%%%%%%

function tw=trigauss(n,alpha,beta)
% by Gaspare Da Fies and Marco Vianello, 
% 8 Nov 2011

% computes the n+1 angles and weights of a trigonometric gaussian
% quadrature formula on [alpha,beta], 0<beta-alpha<=pi

% uses the routines chebyshev.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
% we suggest to put the following statements
% ab = zeros(N,2); sig = zeros(N+1,2*N);
% at the beginning of the body of chebyshev.m to speed-up execution

% input:
% n: trigonometric degree of exactness
% [alpha,beta]: angular interval, 0<beta-alpha<=pi

% output:
% tw: (n+1) x 2 array of (angles,weights)

% the formula integrates the canonical trigonometric basis with accuracy
% from about 10^(-15) (small omega) to about 10^(-13) (omega-->pi)
% up to n=300


% half-length of the angular interval
omega=(beta-alpha)/2;

% modified Chebyshev moments by recurrence
z(1)=2*omega;
z(n+1)=quadgk(@(t)cos(2*n*acos(sin(t/2)/sin(omega/2))),-omega,omega,'MaxIntervalCount',5000);
temp=(2:2:2*n-1);
dl=1/4-1./(4*(temp-1));
dc=1/2-1/sin(omega/2)^2-1./(2*(temp.^2-1));
du=1/4+1./(4*(temp+1));
d=4*cos(omega/2)/sin(omega/2)./(temp.^2-1)';
d(n-1)=d(n-1)-du(n-1)*z(n+1);
z(2:n)=tridisolve(dl(2:n-1),dc(1:n-1),du(1:n-2),d(1:n-1));
mom=zeros(1,2*n+2);
mom(1:2:2*n+1)=z(1:n+1);

% normalization of the moments (monic polynomials)
k=(3:length(mom));
mom(3:end)=exp((2-k)*log(2)).*mom(3:end);

% recurrence coeffs of the monic Chebyshev polynomials
abm(:,1)=zeros(2*n+1,1);
abm(:,2)=0.25*ones(2*n+1,1); abm(1,2)=pi; abm(2,2)=0.5;

% recurrence coeffs for the monic OPS w.r.t. the weight function
% w(x)=2*sin(omega/2)/sqrt(1-sin^2(omega/2)*x^2)
% by the modified Chebyshev algorithm
[ab,normsq]=chebyshev(n+1,mom,abm);

% Gaussian formula for the weight function above
xw=gauss(n+1,ab);

% angles and weights for the trigonometric gaussian formula
tw(:,1)=2*asin(sin(omega/2)*xw(:,1))+(beta+alpha)/2;
tw(:,2)=xw(:,2);

end


function x = tridisolve(a,b,c,d)
%   TRIDISOLVE  Solve tridiagonal system of equations.
% From Cleve Moler's Matlab suite
% http://www.mathworks.it/moler/ncmfilelist.html

%     x = TRIDISOLVE(a,b,c,d) solves the system of linear equations
%     b(1)*x(1) + c(1)*x(2) = d(1),
%     a(j-1)*x(j-1) + b(j)*x(j) + c(j)*x(j+1) = d(j), j = 2:n-1,
%     a(n-1)*x(n-1) + b(n)*x(n) = d(n).
%
%   The algorithm does not use pivoting, so the results might
%   be inaccurate if abs(b) is much smaller than abs(a)+abs(c).
%   More robust, but slower, alternatives with pivoting are:
%     x = T\d where T = diag(a,-1) + diag(b,0) + diag(c,1)
%     x = S\d where S = spdiags([[a; 0] b [0; c]],[-1 0 1],n,n)

x = d;
n = length(x);
for j = 1:n-1
    mu = a(j)/b(j);
    b(j+1) = b(j+1) - mu*c(j);
    x(j+1) = x(j+1) - mu*x(j);
end
x(n) = x(n)/b(n);
for j = n-1:-1:1
    x(j) = (x(j)-c(j)*x(j+1))/b(j);
end
end

%%%%%%%%%%%%%%%%%%%%

% CHEBYSHEV Modified Chebyshev algorithm.
%
%    Given a weight function w encoded by its first 2n modified
%    moments, stored in the (row) vector mom, relative to monic
%    polynomials defined by the (2n-1)x2 array abm of their
%    recurrence coefficients, [ab,normsq]=CHEBYSHEV(n,mom,abm)
%    generates the array ab of the first n recurrence coefficients
%    of the orthogonal polynomials for the weight function w, and
%    the vector normsq of their squared norms. The n alpha-
%    coefficients are stored in the first column, the n beta-
%    coefficients in the second column, of the nx2 array ab. The
%    call [ab,normsq]=CHEBYSHEV(n,mom) does the same, but using the
%    classical Chebyshev algorithm. If n is larger than the sizes
%    of mom and abm warrant, then n is reduced accordingly.
%
function [ab,normsq]=chebyshev(N,mom,abm)
if N<=0, error('N out of range'), end
if N>size(mom,2)/2, N=size(mom,2)/2; end
if nargin<3, abm=zeros(2*N-1,2); end
if N>(size(abm,1)+1)/2, N=(size(abm,1)+1)/2; end
ab(1,1)=abm(1,1)+mom(2)/mom(1); ab(1,2)=mom(1);
if N==1, normsq(1)=mom(1); return, end
sig(1,1:2*N)=0; sig(2,:)=mom(1:2*N);
for n=3:N+1
    for m=n-1:2*N-n+2
        sig(n,m)=sig(n-1,m+1)-(ab(n-2,1)-abm(m,1))*sig(n-1,m) ...
            -ab(n-2,2)*sig(n-2,m)+abm(m,2)*sig(n-1,m-1);
    end
    ab(n-1,1)=abm(n-1,1)+sig(n,n)/sig(n,n-1)-sig(n-1,n-1)/ ...
        sig(n-1,n-2);
    ab(n-1,2)=sig(n,n-1)/sig(n-1,n-2);
end
for n=1:N, normsq(n)=sig(n+1,n); end; normsq=normsq';

end

%%%%%%%%%%%%%%%%%%%%

% R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
%
%    ab=R_JACOBI(n,a,b) generates the first n recurrence
%    coefficients for monic Jacobi polynomials with parameters
%    a and b. These are orthogonal on [-1,1] relative to the
%    weight function w(t)=(1-t)^a(1+t)^b. The n alpha-coefficients
%    are stored in the first column, the n beta-coefficients in
%    the second column, of the nx2 array ab. The call ab=
%    R_JACOBI(n,a) is the same as ab=R_JACOBI(n,a,a) and
%    ab=R_JACOBI(n) the same as ab=R_JACOBI(n,0,0).
%
%    Supplied by Dirk Laurie, 6-22-1998; edited by Walter
%    Gautschi, 4-4-2002.
%
function ab=r_jacobi(N,a,b)
if nargin<2, a=0; end;  if nargin<3, b=a; end
if((N<=0)|(a<=-1)|(b<=-1)) error('parameter(s) out of range'), end
nu=(b-a)/(a+b+2);
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
if N==1, ab=[nu mu]; return, end
N=N-1; n=1:N; nab=2*n+a+b;
A=[nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
n=2:N; nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
ab=[A' [mu; B1; B']];

end

%%%%%%%%%%%%%%%%%%%%

% GAUSS Gauss quadrature rule.
%
%    Given a weight function w encoded by the nx2 array ab of the
%    first n recurrence coefficients for the associated orthogonal
%    polynomials, the first column of ab containing the n alpha-
%    coefficients and the second column the n beta-coefficients,
%    the call xw=GAUSS(n,ab) generates the nodes and weights xw of
%    the n-point Gauss quadrature rule for the weight function w.
%    The nodes, in increasing order, are stored in the first
%    column, the n corresponding weights in the second column, of
%    the nx2 array xw.
%
function xw=gauss(N,ab)
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

end

%%%%%%%%%%%%%%%%%%%%

function [pts,w,momerr] = comprexcub(deg,X,omega,pos)

% compression of bivariate cubature formulas by Tchakaloff points
% or approximate Fekete points
% useful, for example, in node reduction of algebraic cubature formulas
% see the web page: http://www.math.unipd.it/~marcov/Tchakaloff.html

% by Federico Piazzon, Alvise Sommariva and Marco Vianello
% , May 2016


% INPUT:
% deg: polynomial exactness degree
% X: 2-column array of point coordinates
% omega: 1-column array of weights
% pos: NNLS for pos=1, QR with column pivoting for pos=0

% OUTPUT:
% pts: 2-column array of extracted points
% w: 1-column array of corresponding weights (positive for pos=1)
% momerr: moment reconstruction error


% FUNCTION BODY
% total-degree Chebyshev-Vandermonde matrix at X
rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
V=chebvand(deg,X,rect);
% polynomial basis orthogonalization
[Q,R]=qr(V,0);
% tiny imaginary parts could appear
Q=real(Q);
% possible re-orthogonalization
% [Q,R]=qr(Q,0);

% moments of the orthogonal basis
orthmom=Q'*omega;
% weigths computation
switch pos 
    case  1
    % Tchakaloff points (positive weights)
    weights=lsqnonneg(Q',orthmom);
        case  2
    % Tchakaloff points (positive weights)
    weights=LHDM(Q',orthmom);
    otherwise
    % approximate Fekete points (possible negative weights)
    weights=Q'\orthmom;
end;
% indexes of nonvanishing weights and compression
ind=find(abs(weights)>0);
pts=X(ind,:);
w=weights(ind);

% moment reconstruction error
% bivariate Chebyshev basis
% mom=V'*omega;
% momerr=norm(V(ind,:)'*w-mom);
% discrete OP basis
momerr=norm(Q(ind,:)'*w-orthmom);

end

%%%%%%%%%%%%%%%%%%%%

function V = chebvand(deg,X,rect);

% INPUT:
% deg = polynomial degree
% X = 2-column array of point coordinates
% rect = 4-component vector such that the rectangle
% [rect(1),rect(2)] x [rect(3),rect(4)] contains X

% OUTPUT:
% V = Chebyshev-Vandermonde matrix at X, graded lexic. order

% FUNCTION BODY
% rectangle containing the mesh
if isempty(rect)
    rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
end;

% couples with length less or equal to deg
% graded lexicographical order
j=(0:1:deg);
[j1,j2]=meshgrid(j);
dim=(deg+1)*(deg+2)/2;
couples=zeros(dim,2);
for s=0:deg
    good=find(j1(:)+j2(:)==s);
    couples(1+s*(s+1)/2:(s+1)*(s+2)/2,:)=[j1(good) j2(good)];
end;

% mapping the mesh in the square [-1,1]^2
a=rect(1);b=rect(2);c=rect(3);d=rect(4);
map=[(2*X(:,1)-b-a)/(b-a) (2*X(:,2)-d-c)/(d-c)];

% Chebyshev-Vandermonde matrix on the mesh
T1=chebpolys(deg,map(:,1));
T2=chebpolys(deg,map(:,2));
V=T1(:,couples(:,1)+1).*T2(:,couples(:,2)+1);

end

%%%%%%%%%%%%%%%%%%%%

function T=chebpolys(deg,x)

% computes the Chebyshev-Vandermonde matrix on the real line by recurrence

% INPUT:
% deg = maximum polynomial degree
% x = 1-column array of abscissas

% OUTPUT
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg

T=zeros(length(x),deg+1);
t0=ones(length(x),1);
T(:,1)=t0;
t1=x;
T(:,2)=t1;

for j=2:deg
    t2=2*x.*t1-t0;
    T(:,j+1)=t2;
    t0=t1;
    t1=t2;
end;

end






