%% FE code for solving poisson's equation in 2D using linear quadrilateral
% % element

clc
clear 
close all

%% Reading the input file
%Reading the input file

filename = 'Problem_04_Triangular_0.05.inp';

% get node matrix
string1 = '*Node';

% nodes - 2D array with 3 columns - node no., x-coordinate, y-coordinate and rows = no. of nodes 
nodes = readinp(string1,filename);

% get element connectivity matrix for element type CPS8R
string2 = '*Element, type=DC2D6';

% elements - 2D array with columns - element no., nodes in element and connectivity and rows = no. of elements 
elements = readinp(string2,filename);

% get nodes for nset=NodeSetName
string3 = '*Nset, nset=Set-1';

% nsetbc1 - 1D array of node nnumbers in the set  (set-3 for the circle for
% SH Code)
DirBC1 = readinp(string3,filename);

% get elements for elset=ElementSetName           (set 2 for SH Code the
% sides)
string4 = '*Nset, nset=Set-2';

% nsetbc2 - 1D array of node nnumbers in the set 
DirBC2 = readinp(string4,filename); % BC on the left edge

% get elements for elset=ElementSetName
%string5 = '*Nset, nset=Set-5';

% nsetbc2 - 1D array of node nnumbers in the set 
%DirBC3 = readinp(string5,filename); % BC on the right edge

% nsetbc2 - 1D array of element numbers in the set 

NeumannBC1=readmatrix('Problem_04_Element_node_tri_top_0.05.xlsx','Sheet','Sheet1','Range','A1:D14'); % Upper Neumann BC

% nsetbc2 - 1D array of element numbers in the set

NeumannBC2=readmatrix('Problem_04_Element_node_tri_bottom_0.05.xlsx','Sheet','Sheet1','Range','A1:D14'); % Upper Neumann BC

%Conductivity

lambda0=eye(6);

% No. of nodes

n_elem=size(elements,1);  % no. of elements
nNode_elem=6;             % no. of nodes per element
nNodes_tot=size(nodes,1); % Total no. of nodes

%% Quadrature points

% % 2D_Gaussian_Quadrature

qp2d=[0.6,0.2;
      0.2,0.6;
      0.2,0.2;
      1/3,1/3];

qp2d_wi=[25/48;25/48;25/48;-27/48];

%1D Gaussian Quadrature

qp1d=[-1/sqrt(3.0)  1/sqrt(3.0)]';
w1d=[1 1]';

%% Construction of Global stiffness matrix and force vector
%Global stiffness matrix and force vector

kg=zeros(nNodes_tot);
fgs=zeros(nNodes_tot,1);
fgq=zeros(nNodes_tot,1);



% Loop over elements

for k=1:n_elem
    %  Nodal Co-ordinates
    node_no_k=elements(k,2:end);  % Node numbers corresponding to kth  element
    X_coord=nodes(node_no_k,2)';   % X-Coordinates corresponding to nodal numbers of kth element
    Y_coord=nodes(node_no_k,3)';   % Y-Coordinates corresponding to nodal numbers of kth element
    ke=zeros(nNode_elem);         %(6*6 matrix) Elemental stiffness matrix
    fe_s=zeros(nNode_elem,1);
    
    % Loop over integration points
    for qp=1:4
        xi=qp2d(qp,1);
        eta=qp2d(qp,2);
        w2d=qp2d_wi(qp);

     %Evaluating the isoparametric map
     Ns=shapeFn2d(xi,eta); % Ns shape functions in xi and eta
     xqp=Ns*X_coord'; 
     yqp=Ns*Y_coord';

     % Evaluating the Jacobian
     dNsdxideta=dshapeFn2d(xi,eta);
     dxdxi=X_coord*dNsdxideta(1,:)';
     dxdeta=X_coord*dNsdxideta(2,:)';
     dydxi=Y_coord*dNsdxideta(1,:)';
     dydeta=Y_coord*dNsdxideta(2,:)';
     J=[dxdxi,dxdeta;dydxi,dydeta];
     B= (J^-1)'*dNsdxideta;

     % Element stiffness matrix     
     ke=ke+B'*B*lambda0*det(J)*w2d;

     % Elemental source matrix
     s_xy=(0.5/(sqrt(xqp^2+yqp^2)))-4;
     %fe_s(:,1)=fe_s(:,1)+s_xy*Ns'*w2d*det(J);
     fe_s=fe_s+s_xy*Ns'*w2d*det(J);

    end

    % Populating the global stiffness matrix
    for a=1:nNode_elem
        for b=1:nNode_elem
            kg(node_no_k(a),node_no_k(b))=kg(node_no_k(a),node_no_k(b))+ke(a,b);
        end
        fgs(node_no_k(a),1)=fgs(node_no_k(a),1)+fe_s(a,1);
    end

end

%% Applying Neumann BC to the upper edges

for i=1:length(NeumannBC1)
    neumann_node_no_up=NeumannBC1(i,2:end);
    X_coord_neumann_up=nodes(neumann_node_no_up,2)';
    fe_q_up=zeros(length(X_coord_neumann_up),1);

     for qp=1:length(qp1d)
        eta=qp1d(qp);
        w1d=1;
        N1d_up=shapeFn1d(eta);
        X_coord_iso_up=X_coord_neumann_up*N1d_up';
        dN1Ddeta_up=dshapeFn1d(eta);        
        J_1d_up=X_coord_neumann_up*dN1Ddeta_up'; % Here magnitude of Jacobian and dYdeta are the same
        
        % 1D isoparametric shape functions
        q_xy_up=(0.25/sqrt(X_coord_iso_up^2+0.5^2))-1;
        fe_q_up=fe_q_up+q_xy_up*N1d_up'*w1d*abs(J_1d_up);  % check this if the code doesn't work 

     end
   
     for a=1:3
         fgq(neumann_node_no_up(a),1)=fgq(neumann_node_no_up(a),1)+fe_q_up(a);
     end

end

%% Applying Neumann BC to the lower edges

for i=1:length(NeumannBC2)
    neumann_node_no_dn=NeumannBC2(i,2:end);
    X_coord_neumann_dn=nodes(neumann_node_no_dn,2)';
    fe_q_dn=zeros(length(X_coord_neumann_dn),1);

     for qp=1:length(qp1d)
        eta=qp1d(qp);
        w1d=1;
        N1d_dn=shapeFn1d(eta);
        X_coord_iso_dn=X_coord_neumann_dn*N1d_dn';
        dN1Ddeta_dn=dshapeFn1d(eta);        
        J_1d_dn=X_coord_neumann_dn*dN1Ddeta_dn'; % Here Jacobian and dYdeta are the same
        
        % 1D isoparametric shape functions
        q_xy_dn=1-(0.25/sqrt(X_coord_iso_dn^2+0.5^2));
        fe_q_dn=fe_q_dn+q_xy_dn*N1d_dn'*w1d*abs(J_1d_dn);  % check this if the code doesn't work 

     end
   
     for a=1:3
         fgq(neumann_node_no_dn(a),1)=fgq(neumann_node_no_dn(a),1)+fe_q_dn(a);
     end

end

%% Applying Dirichlet BC to the global stiffness matrix

fg=fgq+fgs;

for i=1:length(DirBC1)
    kg(DirBC1(i),DirBC1(i))=kg(DirBC1(i),DirBC1(i))+10^7;
end


for i=1:length(DirBC1)
    fg(DirBC1(i))=0*10^7;
end

% Applying Dirichlet BC to the edges

for i=1:length(DirBC2)
    kg(DirBC2(i),DirBC2(i))=kg(DirBC2(i),DirBC2(i))+10^7;
end


for i=1:length(DirBC2)
    Y_coord_edge=nodes(DirBC2(i),3);
    T_left_edge=0.25^2+0.5^2+Y_coord_edge^2-(0.5*sqrt(0.5^2+Y_coord_edge^2));
    fg(DirBC2(i))=T_left_edge*10^7;
end




%% Solving for the temperatures profile and post processing

tg=kg\fg;


%% Computing the analytical solution
T_an=zeros(length(nodes),1);

for i=1:length(T_an)
    X_an=nodes(nodes(i),2);   % X-Coordinates corresponding to nodal numbers of kth element
    Y_an=nodes(nodes(i),3);   % Y-Coordinates corresponding to nodal numbers of kth element
    T_an(i,1)=(sqrt(X_an^2+Y_an^2)-0.25)^2;
    
end

%% Computing the L2 norm

error=(T_an-tg);
norm_e= norm(error);
norm_tan=norm(T_an);
L2_norm_4=norm_e/norm_tan;

%% postprocessing the solution

% Prepare a grid in zi and eta
[zi, eta]=meshgrid(linspace(0,1,8));

w = 1-(zi+eta);
out = w<0;
zi(out) = nan;
eta(out) = nan;
w(out) = nan;

n_elem=size(elements,1);  % no. of elements

lambda0=@(x,y)(1.0+0.0*x+0.0*y);

Zi= reshape(zi,[numel(zi) 1]);
Eta= reshape(eta,[numel(eta) 1]);

N2d=shapeFn2d(Zi,Eta);

for k=1:n_elem 

    node_no_k=elements(k,2:end);  % Node numbers corresponding to kth  element
    X_coord=nodes(node_no_k,2)';   % X-Coordinates corresponding to nodal numbers of kth element
    Y_coord=nodes(node_no_k,3)';   % Y-Coordinates corresponding to nodal numbers of kth element


    xZi=N2d*X_coord';
    yZi=N2d*Y_coord';


%     %Obtain solution values at those points from (alpha_i*Ni)

    alphas=tg(node_no_k);
    T=N2d*alphas;
 
%     %Surface plot of temperature

    X=reshape(xZi,size(zi));
    Y=reshape(yZi,size(zi));
    Z=reshape(T,size(zi));

    %Plot the solution in this element

    figure(1)
    hold on
    surf(X,Y,Z)
    xlabel('X')
    ylabel('Y')
    title('Temperature field')

    % Gradients of x and y
    lambda=lambda0(xZi,yZi);
    q=zeros(25,2);
    cordXY=[X_coord;Y_coord];
    
    for p=1:length(Zi)
    J=cordXY*dshapeFn2d(Zi(p),Eta(p))';
    dNxy=(J\eye(2))'*dshapeFn2d(Zi(p),Eta(p));
    q(p,:)=-lambda(p)*dNxy*alphas;
%    
    end

    qx=reshape(q(:,1),size(zi));
    qy=reshape(q(:,2),size(zi));

    figure(2)
    hold on
    quiver(X,Y,qx,qy)
    xlabel('X')
    ylabel('Y')
    title('Heat Flux')

end
 
%% Exporting all the graphics
figure(1)
colormap('jet');
colorbar;
exportgraphics(gcf,'Problem_4_temp_field_Sudharshan.png','Resolution',600)

figure(2)
colormap('jet');
exportgraphics(gcf,'Problem_4_heat_flux_Sudharshan.png','Resolution',600)




%% Shape functions and derivatives of the shape functions
% 2D Shape functions
function N=shapeFn2d(xi,eta)
N=[xi.*(2*xi-1),... %N1
   eta.*(2*eta-1),...%N2
   (1-xi-eta).*(1-2*xi-2*eta),...%N3
   4*xi.*eta,...%N4
   4*eta.*(1-xi-eta),...%N5
   4*xi.*(1-xi-eta)];%N6
   end

%1D shape functions
function N=shapeFn1d(xi)
N=[(1/2)*xi*(xi-1),(1-xi)*(1+xi),(1/2)*xi*(1+xi)];
end

% Derivative of 2D shape functions in isoparametric elements

function dN=dshapeFn2d(xi,eta)
% All the shape functions are arranged row-wise
dN=[4*xi-1,0*xi,4*xi+4*eta-3,4*eta,-4*eta,-4*(2*xi+eta-1);
    0*eta,4*eta-1,4*xi+4*eta-3,4*xi,-4*(xi+2*eta-1),-4*xi];
end

%Derivatives of 1D shape functions

function dN=dshapeFn1d(xi)
dN=[(xi-0.5),(-2*xi),(xi+0.5)];
end

end