%% Renan Liupekevicius Carnielli TU/e

% Computational Homogenization of LRAM with viscoelastic constitutive model
clear all; close all;
%% IMPORT MESH
disp('codeblock: SELECT DESIGN')
% select unit cell design
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
design = 'hole';
disp(design);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------



% read mesh text file comsol 5.4
%--------------------------------------------------------------------------
disp('codeblock: READ & CONVERT LINEAR TO QUADRATIC ELEMENT MESH')
% read text file
l   = read_mphtxt_54(design);

% convert to quadratic element
mesh=mesh2d_lin2qua_uc(l);
%--------------------------------------------------------------------------


% read mesh text file comsol 5.6
%--------------------------------------------------------------------------
% disp('codeblock: READ MESH COMSOL 5.6')
% mesh   = read_mphtxt_56(design);
%--------------------------------------------------------------------------

% copy to local variables (please double click on 'mesh' struct in 
% workspace to see the meaning of each cell)
%--------------------------------------------------------------------------
  x    = mesh.x{2,2};    % tensor form
  mx   = mesh.x{2,1};    % matrix form
  conn = mesh.elem{2,1}; % conn matrix of element
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% mesh parameters
  m  = size(  conn, 1); % number of elements in the mesh
  n  = size(  x, 1);    % number of nodes in the mesh 
  fprintf('NUMBER OF MESH ELEMENTS: %d\n', m)
  fprintf('NUMBER OF MESH NODES: %d\n', n)
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit cell size
  a = 2*max(  mx(:,1));
ax=a;
ay=a;

if(max(mx(:,1))+min(mx(:,1))>1e-8); error('RVE is not x centralized');end
if(max(mx(:,1))+min(mx(:,1))>1e-8); error('RVE is not y centralized');end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

%% MATERIAL PROPERTIES

%--------------------------------------------------------------------------
dim   = 2;  % problem dimension
thick = 1;  % thickness (out-of-plane direction), in [m]
%--------------------------------------------------------------------------




    n_solid_phases = 3;
    % define solid phase parameters
    %----------------------------------------------------------------------
    %matrix (part 1)
    matProp{1}.E       =   3.6e9; % Young's modulus
    % matProp{1}.E       =   3.6e6; warning(' test for G.Kuci')
    matProp{1}.nu      =   0.368; % poisson ratio
    matProp{1}.rho     =   1180; % mass density
    % matProp{1}.eta     =  0; %
    matProp{1}.eta     =  0; warning('viscosity of 0 Pa for matrix')
    matProp{1}.G       =   matProp{1}.E/2/(1 +   matProp{1}.nu);
    matProp{1}.kappa   =   matProp{1}.E/3/(1 - 2*matProp{1}.nu);
    matProp{1}.cL      =  sqrt(matProp{1}.E * (1-matProp{1}.nu) ...
                                /(1+matProp{1}.nu)/(1 - 2*matProp{1}.nu)...
                                /matProp{1}.rho);
      
    
    
    % coating (part 2)
    matProp{2}.E       =   0.1175e6; % original
    matProp{2}.nu      =   0.469; %MAKE IT nu=0.3 TO ENSURE NON-COMPRESSIB.
    matProp{2}.rho     =   1300;
    matProp{2}.eta     =   10; warning('viscosity of 10 Pa for rubber')
    matProp{2}.G       =   matProp{2}.E/2/(1 +   matProp{2}.nu);
    matProp{2}.kappa   =   matProp{2}.E/3/(1 - 2*matProp{2}.nu);
    matProp{2}.cL      =   sqrt(matProp{2}.E * (1-matProp{2}.nu) ...
                                /(1+matProp{2}.nu)/(1 - 2*matProp{2}.nu)...
                                /matProp{2}.rho);
    
    % core (part 3)
    matProp{3}.E       =   40.82e9;
    matProp{3}.nu      =   0.37;
    matProp{3}.rho     =   11600;
    matProp{3}.eta     =   0; warning('viscosity of 0 Pa for core')
    matProp{3}.G       =   matProp{3}.E/2/(1 +   matProp{3}.nu);
    matProp{3}.kappa   =   matProp{3}.E/3/(1 - 2*matProp{3}.nu);
    matProp{3}.cL      =  sqrt(matProp{3}.E * (1-matProp{3}.nu) ...
                                /(1+matProp{3}.nu)/(1 - 2*matProp{3}.nu)...
                                /matProp{3}.rho);
    %----------------------------------------------------------------------





% define tensor basis
%--------------------------------------------------------------------------
b  = {'e1'; 'e2'};
ee = cartesianbasis2d(  b{1},   b{2});
e1 = ee(1);
e2 = ee(2);
%--------------------------------------------------------------------------


% calculate fundamental tensors
%--------------------------------------------------------------------------
  I   = identity(2, ee);
  I4  = identity(4, ee);
  I4S = 1/2 * (  I4 + rtranspose( I4));
%--------------------------------------------------------------------------


%%  DESIGNs: LIST OF TAGS 

% MANUALLY import tags from comsol: check acoustic-structure boundary
% list under 'boundary selection'. Same procedure for selecting tag-list of
% the solid phase(s).


switch sscanf(design, '%c')

    case 'LRAM_square_inclusion'
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1:12];
    % solid1 elements tag
    list_solid1_tag = [1 2 7 12];  
    % solid elements tag
    list_solid2_tag = [3 4 8 11];  
    % solid elements tag
    list_solid3_tag = [5 6 9 10];  
    % interface edge elements tag  
    list_itr_edges_tag=[];

    case 'PC'
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1:8];
    % solid1 elements tag
    list_solid1_tag = [1 2 5 8];  
    % solid elements tag
    list_solid2_tag = [];  
    % solid elements tag
    list_solid3_tag = [3 4 6 7];  
    % interface edge elements tag  
    list_itr_edges_tag=[];

    case 'square' 
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1];
    % solid1 elements tag
    list_solid1_tag = [];  
    % solid elements tag
    list_solid2_tag = [1];  
    % solid elements tag
    list_solid3_tag = [];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------

    case 'hole' 
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1];
    % solid1 elements tag
    list_solid1_tag = [];  
    % solid elements tag
    list_solid2_tag = [1];  
    % solid elements tag
    list_solid3_tag = [];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------


    case 'LRAM' 
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1:12];
    % solid1 elements tag
    list_solid1_tag = [1 2 7 12];  
    % solid elements tag
    list_solid2_tag = [3 4 8 11];  
    % solid elements tag
    list_solid3_tag = [5 6 9 10];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------

    case 'LRAM_coarse' 
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1:12];
    % solid1 elements tag
    list_solid1_tag = [1 2 7 12];  
    % solid elements tag
    list_solid2_tag = [3 4 8 11];  
    % solid elements tag
    list_solid3_tag = [5 6 9 10];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------

    case 'LRAM_coarser' 
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1:12];
    % solid1 elements tag
    list_solid1_tag = [1 2 7 12];  
    % solid elements tag
    list_solid2_tag = [3 4 8 11];  
    % solid elements tag
    list_solid3_tag = [5 6 9 10];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------


    case 'LRAM_coarse_shifted' 
    % same as LRAM_coarse shifted (a/2,a/2)
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1:12];
    % solid1 elements tag
    list_solid1_tag = [3 4 7 8 12];  
    % solid elements tag
    list_solid2_tag = [2 5 9 10];  
    % solid elements tag
    list_solid3_tag = [1 6 11];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------
     
    case 'square10x10' 
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1];
    % solid1 elements tag
    list_solid1_tag = [];  
    % solid elements tag
    list_solid2_tag = [1];  
    % solid elements tag
    list_solid3_tag = [];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------

    case 'square20x20' 
    % pure solid phase, heavy core, coating embedded in a matrix
    %----------------------------------------------------------------------
    % solid elements tag
    list_solid_tag = [1];
    % solid1 elements tag
    list_solid1_tag = [];  
    % solid elements tag
    list_solid2_tag = [1];  
    % solid elements tag
    list_solid3_tag = [];  
    % interface edge elements tag  
    list_itr_edges_tag=[];
    %----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add here new case if you designed a new geometry/mesh
% case 'new_mesh_here' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%% SORTING SOLID ELEMENTS

% TOTAL SOLID ELEMENTS DISREGARDING FROM WHICH MATERIAL PHASE
%--------------------------------------------------------------------------    
    % solid elements
      solid_elems =[];
      for i=list_solid_tag
      solid_elems = [solid_elems find( mesh.elem{2,2}==i).'];
      end
    % position vector of the solid nodes
      conns       =   conn(  solid_elems,:);
    % solid nodes
      nodes_s     = unique(reshape(  conns,1,[]));
    % number of solid elements in the mesh
      ms          = size(  conns,1);
    % number of solid nodes
      ns          = length(  nodes_s);
%--------------------------------------------------------------------------



    % SOLID 1
    %--------------------------------------------------------------------------    
        % solid elements
          solid1_elems =[];
          for i=list_solid1_tag
          solid1_elems = [solid1_elems find( mesh.elem{2,2}==i).'];
          end
        % position vector of the solid nodes
          connss{1}       =   conn(  solid1_elems,:);
        % solid nodes
          nodes_ss{1}     = unique(reshape(  connss{1},1,[]));
        % number of solid elements in the mesh
          mss{1}          = size  (  connss{1}, 1);
        % number of solid nodes
          nss{1}          = length(  nodes_ss{1});
    %--------------------------------------------------------------------------
    
    % SOLID 2
    %--------------------------------------------------------------------------    
        % solid elements
          solid2_elems =[];
          for i=list_solid2_tag
          solid2_elems = [solid2_elems find( mesh.elem{2,2}==i).'];
          end
        % position vector of the solid nodes
          connss{2}       =   conn(  solid2_elems,:);
        % solid nodes
          nodes_ss{2}     = unique(reshape(  connss{2},1,[]));
        % number of solid elements in the mesh
          mss{2}          = size  (  connss{2}, 1);
        % number of solid nodes
          nss{2}          = length(  nodes_ss{2});
    %--------------------------------------------------------------------------
    
    % SOLID 3
    %--------------------------------------------------------------------------    
        % solid elements
          solid3_elems =[];
          for i=list_solid3_tag
          solid3_elems = [solid3_elems find( mesh.elem{2,2}==i).'];
          end
        % position vector of the solid nodes
          connss{3}       =   conn(  solid3_elems,:);
        % solid nodes
          nodes_ss{3}     = unique(reshape(  connss{3},1,[]));
        % number of solid elements in the mesh
          mss{3}          = size  (  connss{3},1);
        % number of solid nodes
          nss{3     }     = length(  nodes_ss{3});
    %--------------------------------------------------------------------------

%% PLOT MESH



    %--------------------------------------------------------------------------
    % plot each fluid phase
    figure(2)
    clf
    daspect([1 1 1]);
    hold on;
    plotMesh(  mx,  connss{1},  [0.8 0.8 0.8] );  %gray
    plotMesh(  mx,  connss{2},  [0 0 1] );        %blue
    plotMesh(  mx,  connss{3},  [1 0 0] );       % red
      % plotMesh(  mx,  connss{1},  [1 1 1] );  
    % plotMesh(  mx,  connss{2},  [1 1 1] );        
    % plotMesh(  mx,  connss{3},  [1 1 1] );       
    hold off;
    xlabel('m');
    ylabel('m');
    %--------------------------------------------------------------------------

    % %--------------------------------------------------------------------------
    % % % plot each solid phase
    % figure(2)
    % clf
    % hold on;
    % femplot(  x,  connss{1},'Color' , 'black', 'Nodes', 'off', ...
    % 'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    % femplot(  x,  connss{2},'Color' , 'red', 'Nodes', 'off', ...
    % 'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    % femplot(  x,  connss{3},'Color' , 'blue', 'Nodes', 'off', ...
    % 'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    % % scatter(mx(1:2,1),mx(1:2,2))
    % hold off;
    % set(gcf,'Color','w')
    % set(gca,'XTick',[], 'YTick', [])
    % % exportgraphics(gca,'plot.png','BackgroundColor','none')
    % %--------------------------------------------------------------------------


    % %--------------------------------------------------------------------------
    % % % plot each solid phase
    % figure(2)
    % clf
    % hold on;
    % femplot(  x,  connss{1},'Color' , 'black', 'Nodes', 'off', ...
    % 'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    % femplot(  x,  connss{2},'Color' , 'black', 'Nodes', 'off', ...
    % 'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    % femplot(  x,  connss{3},'Color' , 'black', 'Nodes', 'off', ...
    % 'NodeNumbers' , 'off', 'ElementNumbers', 'off');
    % % scatter(mx(1:2,1),mx(1:2,2))
    % hold off;
    % set(gcf,'Color','w')
    % set(gca,'XTick',[], 'YTick', [])
    % % exportgraphics(gca,'plot.png','BackgroundColor','none')
    % %---------

%% NODE SORTING OF BOUNDARIES
disp('codeblock: NODE SORTING')

% NOTE: this codeblock assumes the top and bottom boundaries are solid
% boundaries.

%--------------------------------------------------------------------------
% mesh precision. Check in COMSOL for correct use of it. Here I set 1e-8 by
% try and error which is not ideal.
precision = 1e-8; 
%--------------------------------------------------------------------------


% sort nodes of the solid boundaries
%--------------------------------------------------------------------------
% nodes on the left boundary of the mesh: line x = -a/2
[  l,~]        = find(abs(  mx(:,1) +   a/2 )<precision ); % 
  left_s       = intersect(  l,  nodes_s)';  % solid boundary

% nodes on the right boundary of the mesh: line x = +a/2                  
[  r,~]        = find(abs(  mx(:,1) -   a/2 )<precision ); % 
  right_s      = intersect(  r,  nodes_s)';  % fluid boundary

% nodes on the bottom boundary of the mesh: line y = -a/2
[  bottom_s,~] = find(abs(  mx(:,2) +   a/2 )<precision ); % 

% nodes on the bottom boundary of the mesh: line y = +a/2
[  top_s,~]    = find(abs(  mx(:,2) -   a/2 )<precision ); % 


% corner nodes
  corners_s = [   intersect(  l,  bottom_s) ...
                  intersect(  r,  bottom_s) ...
                  intersect(  r,  top_s   ) ...
                  intersect(  l,  top_s   )];
% prescribed corners
  pcorners = [  corners_s(1)   corners_s(2)   corners_s(4)];


% exclude corners from solid node lists
  left_s   = setdiff(  left_s  ,   corners_s)';
  right_s  = setdiff(  right_s ,   corners_s)';
  top_s    = setdiff(  top_s   ,   corners_s)';
  bottom_s = setdiff(  bottom_s,   corners_s)';

clear   l; clear   r;
%--------------------------------------------------------------------------


%% INDEX TABLE
% __________________________________________________________________________
% SOLID (S)
%
%        CURRENT_INDEX      
%               |  
%               V
%__________________________________Reordered____________________Reordered__
% CONVERTION |
%____________|__(S)  _|_(S)_ |__(S)______|___(S)_____|__(S)___
%            |    1   |   1  |  dofs_c   |  dofs_pc  |    p 
%            |    2   |   2  |  dofs_un  |  dofs_un  |    f 
%   ORDER    |    .   |   .  |  dofs_in  |  dofs_in  |
%            |    .   |   .  |  dofs_de  |           |
%            |    .   |   .  |           |           | 
%            |    n   |  ns  |           |           |
%____________|___ ____|______|___________|___________|____________


%%  4 INTEGRATION POINTS WITHIN MASTER ELEMENT
%--------------------------------------------------------------------------
% Notation from TensorLab
% quadratic quadrilateral element:
%         4----7----3
%         |         |
%         8         6
%         |         |
%         1----5----2
%      then CONN = [1 2 3 4 5 6 7 8].
%
% The coordinates of integration points x (in the coordinate system of the
% master element). Exact integration up to polynomial of order 5.
%          --------- 
%         | x  x  x |
%         | x  x  x |
%         | x  x  x |
%          ---------
% and the weight factors
%--------------------------------------------------------------------------
xi = [ -sqrt(3/5) *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 + sqrt(3/5) *e2
       -sqrt(3/5) *e1 + sqrt(3/5) *e2
                0 *e1 - sqrt(3/5) *e2
        sqrt(3/5) *e1 +         0 *e2       
                0 *e1 + sqrt(3/5) *e2                
       -sqrt(3/5) *e1 +         0 *e2
                0 *e1 +         0 *e2   ];
w  = [ 25/81
       25/81
       25/81
       25/81
       40/81
       40/81
       40/81
       40/81
       64/81     ];
%--------------------------------------------------------------------------


%% ASSEMBLE STIFFNESS AND MASS MATRICES (SOLID ELEMENTS)
disp('codeblock: ASSEMBLING SOLID MASS AND STIFFNESS MATRICES')
%--------------------------------------------------------------------------
% Knn and Mnn are only temporary matrix and they are not the matrices K and
% M from discretized form in the notes. Note their size is the total number
% of nodes. It's only relevant property is that it preserves the node 
% indexing.
%--------------------------------------------------------------------------


% declare stiffness matrix tensor with initially zero-2nd-order tensors
%--------------------------------------------------------------------------
Knn = zeros(2, ee,   n,   n); 
Mnn = zeros(  n,   n);
Cnn = zeros(2, ee,   n,   n); %Damping matrix
%--------------------------------------------------------------------------



% if multiple solid phases
% SOLID ELEMENTS ( MULTIPLE PHASE MODEL)
%--------------------------------------------------------------------------

%initialize porosity to 1
phi = 1;

%loop
for i =1:n_solid_phases

 % define material parameters from phase i   
 tC    =  matProp{i}.kappa*I*I + 2* matProp{i}.G  *(  I4S - 1/3*  I*  I);
 tEta  =  matProp{i}.eta  *I*I + 2* matProp{i}.eta*(  I4S - 1/3*  I*  I); 
 rhos  =  matProp{i}.rho;

 % initialize area of each solid phase to zero
 V{i}=0;

    for e = 1:mss{i} % loop over all solid elements of phase i
    
    % display computation percentage
    if mod(floor(e/mss{i}*100),10)==0
    fprintf('solid %d: assembled %.2f ', i, e/mss{i}*100); disp('%');
    end
    
    % extract nodes
    iie =   connss{i}(e, :); % nodes of the current element e
    xe  =   x(iie);            % coordinates of these nodes
        
    % compute element matrices in integration point loop
    Ke = zeros(2, ee, 8, 8); % zero matrix of 8x8 2nd order tensors
                             % in basis ee.
    Me = zeros(8,8);         % zero matrix of 8x8 scalars

    Ce = zeros(2, ee, 8, 8); % zero matrix of 8x8 2nd order tensors
                             % in basis ee

    for k = 1:length(w) % loop over 9 integration points 
        
        xi1 = dot(xi(k), e1); % x1 coordinate of the current integr. point
        xi2 = dot(xi(k), e2); % x2 coordinate of the current integr. point

       % column of the shape functions
        Ne = [ -1/4*(1-xi1  )*(1-xi2  )*(1+xi1+xi2)
               -1/4*(1+xi1  )*(1-xi2  )*(1-xi1+xi2)
               -1/4*(1+xi1  )*(1+xi2  )*(1-xi1-xi2)
               -1/4*(1-xi1  )*(1+xi2  )*(1+xi1-xi2)
                1/2*(1-xi1^2)*(1-xi2  )
                1/2*(1+xi1  )*(1-xi2^2)
                1/2*(1-xi1^2)*(1+xi2  )
                1/2*(1-xi1  )*(1-xi2^2)           ];
        
        % column of the gradient of the shape functions
        % with respect to the local coordinates of the mater element
        gradxiNe = [ - 1/4 * (-1 + xi2   )* (2*xi1 +   xi2) *e1 ...
                     - 1/4 * (-1 + xi1   )* (xi1   + 2*xi2) *e2
                     - 1/4 * (2*xi1 - xi2)* (-1    + xi2  ) *e1 ...
                     - 1/4 * (1 + xi1    )* (xi1   - 2*xi2) *e2
                       1/4 * (1 + xi2    )* (2*xi1 + xi2  ) *e1 ...
                     + 1/4 * (1 + xi1    )* (xi1   + 2*xi2) *e2
                       1/4 * (2*xi1 - xi2)* (1     + xi2  ) *e1 ...
                     + 1/4 * (-1 + xi1   )* (xi1 - 2*xi2  ) *e2
                                      xi1 * (-1 + xi2     ) *e1 ...
                                    + 1/2 * (-1 + xi1^2   ) *e2
                                      1/2 * ( 1 - xi2^2   ) *e1 ...
                                    - xi2 * ( 1 + xi1     ) *e2
                                    - xi1 * ( 1 + xi2     ) *e1 ...
                                    + 1/2 * ( 1 - xi1^2   ) *e2
                                      1/2 * (-1 + xi2^2   ) *e1 ...
                                    + xi2 * (-1 + xi1     ) *e2 ];
        % Jacobian
        J = gradxiNe' * xe;
       
        % column of the gradient of the shape functions
        % with respect to the global coordinates of the mesh
        gradNe = dot(inv(J), gradxiNe);
        % element matrix and right hand side
        Me = Me + w(k) *  Ne *  rhos*   Ne'  *   thick * det(J);
        Ke = Ke + w(k) * dot(gradNe,   tC , gradNe') *   thick * det(J);
        Ce = Ce + w(k) * dot(gradNe,   tEta, gradNe') * thick * det(J);
   
    end % end of interation point loop
   
   % assembly
    Knn(iie, iie) = Knn(iie, iie) + Ke;
    Mnn(iie, iie) = Mnn(iie, iie) + Me;
    Cnn(iie,iie)  = Cnn(iie,iie)  + Ce;  


   % get element area
   V{i}=V{i}+get_element_area(xe,ee);

   end % end of element loop

% compute porosity   
phi = phi - V{i}/a^2;

% assignment actual mass, stiffness and damping matrices
K=Knn(  nodes_s,  nodes_s);
M=Mnn(  nodes_s,  nodes_s);
C=Cnn(  nodes_s,  nodes_s);
end

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% elimiante numerical errors (these matrices are already symmetric)
K = (K+K')/2;
M = (M+M')/2;
C = (C+C')/2;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
clear Me Ke Ce iie xe e
clear Mnn Knn Cnn
clear rhos tC teta
clear k xi1 xi2 gradxiNe J w Ne gradNe
%--------------------------------------------------------------------------




%% INDEX OF DISPLACEMENT/PRESSURE DOF UNDER NODES_S/NODES_F ORDERING
%--------------------------------------------------------------------------
% The matrices Mnn, Knn, Ann, Bnn, and Cnn have indices that correspond to
% the node number, i.g., line 1 are the terms related to node 1.
% The matrices M, K, A, B, and C, however, have indices that do not 
% correspond to the node number. Their indices correspond to the
% displacements degree of freedom or pressure degrees of freedom.
% Suppose, for example, you want to find the index in matrix M that
% correspond to the bottom-left corner displacement dof. This corner node
% is stored in the vector [corner], say first entry, corner(1). The index
% we are looking for, sometimes refered as 'new' index, is found using the
% following command:
% 
% find(nodes_s==corner(1)).
%
% get_ind() was implemented for this operation
%--------------------------------------------------------------------------


% warning:  the operation below doesnt do anything if there is no fluid.

%--------------------------------------------------------------------------
% solid 
    % bottom nodes
    dofs_b = get_ind(  nodes_s,   bottom_s ); 
    % top nodes
    dofs_t = get_ind(  nodes_s,   top_s    ); 
    % left nodes
    dofs_l = get_ind(  nodes_s,   left_s   ); 
    % right nodes
    dofs_r = get_ind(  nodes_s,   right_s  ); 
    % corner nodes
    dofs_c = get_ind(  nodes_s,   corners_s);
    % prescribed corners
    dofs_pc = get_ind(  nodes_s,   pcorners);
%--------------------------------------------------------------------------



 
% __________________________________________________________________________
% SOLID (S)
%
%                  CURRENT_INDEX      
%                        |  
%                        V
%______________________________Reordered____________________Reordered__
% CONVERTION |
%____________|__(S)  _|_(S)_ |__(S)______|___(S)_____|__(S)___
%            |    1   |   1  |  dofs_c   |  dofs_pc  |    p 
%            |    2   |   2  |  dofs_un  |  dofs_un  |    f 
%   ORDER    |    .   |   .  |  dofs_in  |  dofs_in  |
%            |    .   |   .  |  dofs_de  |           |
%            |    .   |   .  |           |           | 
%            |    n   |  ns  |           |           |
%____________|___ ____|______|___________|___________|____________



%% ELIMINATING DEPENDENT NODES OF SOLID PHASE


%--------------------------------------------------------------------------
% constrain relation of periodicity 
% u(top)    = u(bottom) +  (u(coners(4))-u(coners(1))) ones()
% u(right)  = u(left)   +  (u(coners(2))-u(coners(1))) ones()
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    % dependent dofs of periodic bc
    dofs_de = [dofs_t dofs_r];    % [top right], excluding corners
    
    % independent nodes of periodic bc
    dofs_in = [dofs_b dofs_l];    % [bottom left], excluding corners
    

    % safety check
    if(length(dofs_in)~= length(dofs_de)) 
        error("The periodic bc must connect same number of nodes");
    end

% unconstrained nodes ( nodes not at the boundary)
dofs_un = setdiff(1:  ns, [dofs_in dofs_de dofs_c]); % exclude corners                           
% include corners as first elements to follow the notes notation, ideally
% technically we should change the variable name here, but we keep the same
% for convenience. The corner nodes are concatenated at the beginning below
dofs_un = [dofs_c dofs_un];           
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% reordering mass, stiffness, and coupling matrix components
%
% order
% c - corner nodes which were grouped with u nodes for simplicity
% u - unconstrained nodes
% i - independent nodes
% d - dependent nodes (on i nodes)

% mass tensor matrix
% M=  Muu Mui Mud   
%     Miu Mii Mid  
%     Mdu Mdi Mdd  
% transform to tensor matrix
M  =   I*M;
% split
Muu=M(dofs_un, dofs_un); Mui=M(dofs_un, dofs_in); Mud=M(dofs_un, dofs_de);  
Miu=M(dofs_in, dofs_un); Mii=M(dofs_in, dofs_in); Mid=M(dofs_in, dofs_de);  
Mdu=M(dofs_de, dofs_un); Mdi=M(dofs_de, dofs_in); Mdd=M(dofs_de, dofs_de);

% stiffness tensor matrix
% K=  Kuu Kui Kud   
%     Kiu Kii Kid  
%     Kdu Kdi Kdd  
Kuu=K(dofs_un, dofs_un); Kui=K(dofs_un, dofs_in); Kud=K(dofs_un, dofs_de);  
Kiu=K(dofs_in, dofs_un); Kii=K(dofs_in, dofs_in); Kid=K(dofs_in, dofs_de);  
Kdu=K(dofs_de, dofs_un); Kdi=K(dofs_de, dofs_in); Kdd=K(dofs_de, dofs_de);

% Damping tensor matrix
% C=  Cuu Cui Cud   
%     Ciu Cii Cid  
%     Cdu Cdi Cdd  
Cuu=C(dofs_un, dofs_un); Cui=C(dofs_un, dofs_in); Cud=C(dofs_un, dofs_de);  
Ciu=C(dofs_in, dofs_un); Cii=C(dofs_in, dofs_in); Cid=C(dofs_in, dofs_de);  
Cdu=C(dofs_de, dofs_un); Cdi=C(dofs_de, dofs_in); Cdd=C(dofs_de, dofs_de); 
%---------------------------------------------------------

%__________________________________________________________________________
% INDEX TABLE
% FLUID (F)
%
%                               CURRENT_INDEX      
%                                      |  
%                                      V
%__________________________________Reordered____________________Reordered__
% CONVERTION | mesh   | nodes_s |   -----      |    dofs_re   |   ---
%    RULE    | mesh   |         |   nodes_f    |              |  nodes_f
%____________|___(S)__|___(S)___|_____(S)______|______(S)_____|____(S)_____
%            |    1   |    1    |   dofs_c     |    dofs_pc   |     p  
%            |    2   |    2    |   dofs_un    |    dofs_un   |     f  
%   ORDER    |    .   |    .    |   dofs_in    |    dofs_in   |
%            |    .   |    .    |   dofs_de    |              |
%            |    .   |    .    |              |              | 
%            |    n   |    ns   |              |              |
%____________|___ ____|_________|______________|______________|____________


%--------------------------------------------------------------------------
% exclude corners again to stick with notes notation, they will be treated
% separately by the projection matrix
dofs_un([1 2 3 4])=[];
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% projection matrix submatrices zeros, ones and identities
O_1u = zeros(1,length(dofs_un));
O_u1 = O_1u';

O_1B = zeros(1,length(dofs_b));
O_B1 = O_1B';

O_1L = zeros(1,length(dofs_l));
O_L1 = O_1L';

O_uB = zeros(length(dofs_un),length(dofs_b));
O_Bu = O_uB';

O_uL = zeros(length(dofs_un),length(dofs_l));
O_Lu = O_uL';

O_BL = zeros(length(dofs_b),length(dofs_l));
O_LB = O_BL';

I_uu = eye(length(dofs_un),length(dofs_un));
I_BB = eye(length(dofs_b),length(dofs_b));
I_LL = eye(length(dofs_l),length(dofs_l));

I_B1 = ones(length(dofs_b),1);
I_L1 = ones(length(dofs_l),1);
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% projection matrix
%     u_1   = T * u_1
%     u_2         u_2
%     u_3         u_4
%     u_4         u_un
%     u_un        u_in
%     u_in
%     u_de

% projection matrix
T  = [  1       0      0      O_1u    O_1B    O_1L
        0       1      0      O_1u    O_1B    O_1L
       -1       1      1      O_1u    O_1B    O_1L
        0       0      1      O_1u    O_1B    O_1L
       O_u1    O_u1   O_u1    I_uu    O_uB    O_uL
       O_B1    O_B1   O_B1    O_Bu    I_BB    O_BL
       O_L1    O_L1   O_L1    O_Lu    O_LB    I_LL
      -I_B1    O_B1   I_B1    O_Bu    I_BB    O_BL
      -I_L1    I_L1   O_L1    O_Lu    O_LB    I_LL];
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% reduced mass tensor matrix (dependent nodes  and corner 3 eliminated)
M_ =  T.'*[Muu Mui Mud;
           Miu Mii Mid;
           Mdu Mdi Mdd]*T;

% reduced stiffness tensor matrix (dependent nodes and corner 3 eliminated)
K_ =  T.'*[Kuu Kui Kud;
           Kiu Kii Kid;
           Kdu Kdi Kdd]*T;

% reduced damping tensor matrix (dependent nodes and corner 3 eliminated)
C_ =  T.'*[Cuu Cui Cud;
           Ciu Cii Cid;
           Cdu Cdi Cdd]*T;
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% set of remaining nodes, this is the vector that contains the rules for
% changing the indices from 'before projection' to 'after projection'.
% reminder: u_3 u_de were eliminated
%
% re - remaining
%
dofs_re = [dofs_c(1)  dofs_c(2)  dofs_c(4)  dofs_un dofs_in];
ns_re   = length(dofs_re);

% consistency check
% length(dofs_re)~=size(M_,1) or size(K_,1) or size(C_,1)
%--------------------------------------------------------------------------

% __________________________________________________________________________
% SOLID (S)
%
%                                       CURRENT_INDEX      
%                                             |  
%                                             V
%______________________________Reordered____________________Reordered__
% CONVERTION |
%____________|__(S)  _|_(S)_ |__(S)______|___(S)_____|__(S)___
%            |    1   |   1  |  dofs_c   |  dofs_pc  |    p 
%            |    2   |   2  |  dofs_un  |  dofs_un  |    f 
%   ORDER    |    .   |   .  |  dofs_in  |  dofs_in  |
%            |    .   |   .  |  dofs_de  |           |
%            |    .   |   .  |           |           | 
%            |    n   |  ns  |           |           |
%____________|___ ____|______|___________|___________|____________


%% PRESCRIBED AND FREE NODES SPLIT
%--------------------------------------------------------------------------
% INDEX OF DISPLACEMENT DOF IN NODES_S (ELIMINATED DEPENDENT DOFS)
% same index change procedure must be executed but only in the solid phase.

%|__NOTES_NOTATION_|____ DESCRIPTION_______________________________________
%|                 |
%|        p        | prescribed nodes in the solid phase
%|_________________|_______________________________________________________
%|                 |
%|        f        | free nodes in the solid phase                                                      
%|                 |
%|_________________|_______________________________________________________
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% solid
p  = get_ind(dofs_re,dofs_pc); % = [1 2 3] by construction
f  = setdiff(1:ns_re,p);

%--------------------------------------------------------------------------

% __________________________________________________________________________
% SOLID (S)
% 
%                                                       CURRENT_INDEX      
%                                                              |  
%                                                              V
%______________________________Reordered____________________Reordered______
% CONVERTION |
%____________|__(S)  _|_(S)_ |__(S)______|___(S)_____|__(S)________________
%            |    1   |   1  |  dofs_c   |  dofs_pc  |    p 
%            |    2   |   2  |  dofs_un  |  dofs_un  |    f 
%   ORDER    |    .   |   .  |  dofs_in  |  dofs_in  |
%            |    .   |   .  |  dofs_de  |           |
%            |    .   |   .  |           |           | 
%            |    n   |  ns  |           |           |
%____________|___ ____|______|___________|___________|_____________________


%% PARTITIONING

%--------------------------------------------------------------------------
% solid
    % mass tensor matrix
    % M =  Mpp Mpf   
    %      Mfp Mff 
    M_p_p   = M_(p , p ); M_p_f   = M_(p , f );
    M_f_p   = M_(f , p ); M_f_f   = M_(f , f );
    % stiffness tensor matrix
    % K =  Kpp Kpf   
    %      Kfp Kff 
    K_p_p   = K_(p , p ); K_p_f   = K_(p , f );
    K_f_p   = K_(f , p ); K_f_f   = K_(f , f );
    % damping tensor matrix
    % C =  Cpp Cpf   
    %      Cfp Cff 
    C_p_p   = C_(p , p ); C_p_f   = C_(p , f );
    C_f_p   = C_(f , p ); C_f_f   = C_(f , f );

%--------------------------------------------------------------------------



%% TRANFORM TENSOR MATRIX TO EQUIVALENT SCALAR MATRIX
disp('codeblock: TRANFORM TENSOR MATRIX TO EQUIVALENT SCALAR MATRIX')

%--------------------------------------------------------------------------

% mass matrix M
sM_p_p = tm2sm( M_p_p, ee , 2 ); sM_p_f = tm2sm( M_p_f, ee , 2 );
sM_f_p = tm2sm( M_f_p, ee , 2 ); sM_f_f = tm2sm( M_f_f, ee , 2 );

% stiffness matrix K
sK_p_p = tm2sm( K_p_p, ee , 2 ); sK_p_f = tm2sm( K_p_f, ee , 2 );
sK_f_p = tm2sm( K_f_p, ee , 2 ); sK_f_f = tm2sm( K_f_f, ee , 2 );

% damping matrix C
sC_p_p = tm2sm( C_p_p, ee , 2 ); sC_p_f = tm2sm( C_p_f, ee , 2 );
sC_f_p = tm2sm( C_f_p, ee , 2 ); sC_f_f = tm2sm( C_f_f, ee , 2 );

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%elimiante numerical errors (these matrices are already symmetric) 
     sM_p_p  = sparse(( sM_p_p +  sM_p_p.' )/2 );
     sM_f_f  = sparse(( sM_f_f +  sM_f_f.' )/2 );
     sK_p_p  = sparse(( sK_p_p +  sK_p_p.' )/2 );
     sK_f_f  = sparse(( sK_f_f +  sK_f_f.' )/2 );
     sM_p_f  = sparse(( sM_p_f +  sM_f_p.' )/2 );
     sM_f_p  =  sM_p_f.';
     sK_f_p  = sparse(( sK_f_p +  sK_p_f.' )/2 );
     sK_p_f  =  sK_f_p.';
%--------------------------------------------------------------------------

%% STATE SPACE
%--------------------------------------------------------------------------
% Eigenvalue problem solving 
O_p_p = sparse(zeros( size(ee,1)*length(p) , size(ee,1)*length(p)));
O_p_f = sparse(zeros( size(ee,1)*length(p) , size(ee,1)*length(f)));
O_f_p = sparse(zeros( size(ee,1)*length(f) , size(ee,1)*length(p)));
O_f_f = sparse(zeros( size(ee,1)*length(f) , size(ee,1)*length(f)));
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
cases = 2; warning( 'state space format [CM;M0]dotw+[K0;0-M]w=0')
% cases = 3; warning( 'state space format [-K0;0M]dotw+[0K;KC]w=0')
switch cases
    case 2
    % lambda submatrices
    lam_B_B = [sC_p_p sM_p_p; sM_p_p O_p_p];
    
    lam_B_I = [sC_p_f sM_p_f; sM_p_f O_p_f];
    
    lam_I_B = [sC_f_p sM_f_p; sM_f_p O_f_p];
    
    lam_I_I = [sC_f_f sM_f_f; sM_f_f O_f_f];
    
    %mu submatrices

    mu_B_B  = [sK_p_p O_p_p; O_p_p -sM_p_p];
    
    mu_B_I  = [sK_p_f O_p_f; O_p_f -sM_p_f];
    
    mu_I_B  = [sK_f_p O_f_p; O_f_p -sM_f_p];
        
    mu_I_I  = [sK_f_f O_f_f; O_f_f -sM_f_f];  

    %Case 3
    case 3  
    % lambda submatrices
    lam_B_B = [-sK_p_p O_p_p; O_p_p sM_p_p];
    
    lam_B_I = [-sK_p_f O_p_f; O_p_f sM_p_f];
    
    lam_I_B = [-sK_f_p O_f_p; O_f_p sM_f_p];
    
    lam_I_I = [-sK_f_f O_f_f; O_f_f sM_f_f];

    %mu submatrices

    mu_B_B  = [O_p_p sK_p_p; sK_p_p sC_p_p];

    mu_B_I  = [O_p_f sK_p_f; sK_p_f sC_p_f];
    
    mu_I_B  = [O_f_p sK_f_p; sK_f_p sC_f_p];
    
    mu_I_I  = [O_f_f sK_f_f; sK_f_f sC_f_f];
end
%--------------------------------------------------------------------------


%% EIGEN: EIGEN STUDY, PRINT EIGENFREQUENCIES AND EIGENMODES
% 
% % 
% %--------------------------------------------------------------------------
% % define desired number of modes to compute/display
% % n_modes =12; %eta=10 
% n_modes =16; %eta=0
% 
% disp('ITERATIVE SOLVER')
% 
% 
%     % n_modes=16
%     if matProp{2}.eta == 10 ; warning('added complex conjugate n_modes/2 to the basis')
%     tic;[phi_I_Q, Dr] = eigs(mu_I_I,lam_I_I ,n_modes/2,-1i*2*pi*(1000 + 1i*10),'SubspaceDimension',2*n_modes,'Display',1);toc; % eta=10
%     Dr       = [Dr zeros(n_modes/2); zeros(n_modes/2) conj(Dr)];
%     phi_I_Q  = [phi_I_Q conj(phi_I_Q)];
% 
%     % mode selection for eta=10
%     % mode_selec = [1 2 3 4 5 8 9 10 11 12 13 16];warning('mode selection')
%     mode_selec = [1 2 4 5 9 10 12 13];warning('mode selection degenerate eta=10')
%     Dr=Dr(mode_selec,mode_selec);
%     phi_I_Q=phi_I_Q(:,mode_selec);
%     n_modes=length(mode_selec);
% 
%     elseif matProp{2}.eta == 0 ; warning('added complex conjugate n_modes/2 to the basis ')
%     tic;[phi_I_Q, Dr] = eigs(mu_I_I,lam_I_I ,n_modes/2,-1i*2*pi*(200 ),'SubspaceDimension',2*n_modes,'Display',1);toc; % eta=0
%     Dr       = [Dr zeros(n_modes/2); zeros(n_modes/2) conj(Dr)];
%     phi_I_Q  = [phi_I_Q conj(phi_I_Q)];
%     mode_selec = [2 3 4 5 7 8 10 11 12 13 15 16];warning('selected only the degenerate pairs eta=0')
%     Dr=Dr(mode_selec,mode_selec);
%     phi_I_Q=phi_I_Q(:,mode_selec);
%     n_modes=length(mode_selec);
%     else
% 
%     % if matProp{2}.eta == 10 ; warning('test few modes')
%     % tic;[phi_I_Q, Dr] = eigs(mu_I_I,lam_I_I ,n_modes/2,-1i*2*pi*(800 + 1i*1000),'SubspaceDimension',2*n_modes,'Display',1);toc; % eta=10
%     % Dr       = [Dr zeros(n_modes/2); zeros(n_modes/2) conj(Dr)];
%     % phi_I_Q  = [phi_I_Q conj(phi_I_Q)];
%     % 
%     % % % mode selection for eta=10
%     % % mode_selec = [1 2 3 4 5 8 9 10 11 12 13 16];warning('mode selection')
%     % % Dr=Dr(mode_selec,mode_selec);
%     % % phi_I_Q=phi_I_Q(:,mode_selec);
%     % % n_modes=length(mode_selec);
%     % else
%         tic;[phi_I_Q, Dr] = eigs(mu_I_I,lam_I_I ,n_modes,-1i*2*pi*(600),'SubspaceDimension',300,'Display',1);toc; warning('eigenvalue lambda= - i omega');
%     % tic;[phi_I_Q, Dr] = eigs(mu_I_I,lam_I_I ,n_modes,-1i*2*pi*(10+1i*1),'SubspaceDimension',5*n_modes,'Display',1);toc; warning('eigenvalue lambda= - i omega');
%     % tic;[phi_I_Q, Dr] = eigs(mu_I_I,lam_I_I ,n_modes,-1i*2*pi*(1000+1i*60),'SubspaceDimension',200,'Display',1);toc; warning('eigenvalue lambda= - i omega');
%     end
% 
% 
% 
%     % define
%     eigfreqs    = diag(Dr/(-1i*2*pi));
%     eigfreqs_el = abs(diag(Dr/(-1i*2*pi)));
% 
% 
%     %display
%     fprintf('eigenfrequencies: \n ');
%     for i=1:n_modes;fprintf('%f + i%f Hz\n ', real(eigfreqs(i)),imag(eigfreqs(i)));end
%     % fprintf('damping ratio: \n ');
%     % fprintf('%f \n ', abs(imag(eigfreqs)./abs(eigfreqs)));
%     % fprintf('Elastic eigenfrequencies: \n ');
%     % fprintf('%f Hz\n ', eigfreqs_el);
% 
    
%% EIGEN: compute & split beteween degenerate or not (all (almost) modes )


%  % tic;[phi_I_Q_FULL, Dr_FULL] = eig(full(mu_I_I),full(lam_I_I));toc;
% % load workspace "eig_FULL_eta1000.mat" instead of runing the line above
% % load eig_FULL_eta00.mat
% 
% % def
% eigfreqs_FULL    = diag(Dr_FULL/(-1i*2*pi));
% 
% % select modes
% % f>1 Hz
% ind              = real(eigfreqs_FULL)>1;
% eigfreqs_FULL    = eigfreqs_FULL(ind);
% phi_I_Q_FULL     = phi_I_Q_FULL(:,ind);
% Dr_FULL          = Dr_FULL(ind,ind);
% 
% 
% % f< 3000 Hz
% ind              = real(eigfreqs_FULL) < 3000;
% eigfreqs_FULL    = eigfreqs_FULL(ind);
% phi_I_Q_FULL     = phi_I_Q_FULL(:,ind);
% Dr_FULL          = Dr_FULL(ind,ind);
% 
% % order
% [~,ind]=sort(real(eigfreqs_FULL), 'ascend');
% eigfreqs_FULL    = eigfreqs_FULL(ind);
% phi_I_Q_FULL     = phi_I_Q_FULL(:,ind);
% Dr_FULL          = Dr_FULL(ind,ind);
% 
% 
% eigfreqs = eigfreqs_FULL;
% phi_I_Q  = phi_I_Q_FULL;
% Dr       = Dr_FULL;
% n_modes  = size(Dr,1);
% 
% Dr       = [Dr zeros(n_modes); zeros(n_modes) conj(Dr)];
% phi_I_Q  = [phi_I_Q conj(phi_I_Q)];
% n_modes  = 2*n_modes;
% eigfreqs    = diag(Dr/(-1i*2*pi));
% 
% fprintf('eigenfrequencies: \n ');
% for i=1:n_modes;fprintf('%f + i%f Hz\n ', real(eigfreqs(i)),imag(eigfreqs(i)));end


%% EIGEN: compute & split beteween degenerate or not (6 modes from paper: 4 modes and 2 rotational)


if matProp{2}.eta==10
n_modes =16; %eta=10 
disp('ITERATIVE SOLVER')



    tic;[phi_I_Q, Dr] = eigs(mu_I_I,lam_I_I ,n_modes/2,-1i*2*pi*(1000 + 1i*10),'SubspaceDimension',2*n_modes,'Display',1);toc; % eta=10
    Dr       = [Dr zeros(n_modes/2); zeros(n_modes/2) conj(Dr)];
    phi_I_Q  = [phi_I_Q conj(phi_I_Q)];

    % mode selection for eta=10
    mode_selec = [1 2 3 4 5 8 9 10 11 12 13 16];warning('mode selection, exclude 1069s')
    Dr=Dr(mode_selec,mode_selec);
    phi_I_Q=phi_I_Q(:,mode_selec);
    n_modes=length(mode_selec);




    % define
    eigfreqs    = diag(Dr/(-1i*2*pi));
    eigfreqs_el = abs(diag(Dr/(-1i*2*pi)));


    %display
    fprintf('eigenfrequencies: \n ');
    for i=1:n_modes;fprintf('%f + i%f Hz\n ', real(eigfreqs(i)),imag(eigfreqs(i)));end
end


if matProp{2}.eta==2
n_modes =12; %eta=2 
disp('ITERATIVE SOLVER')



    tic;[phi_I_Q, Dr] = eigs(mu_I_I,lam_I_I ,n_modes/2,-1i*2*pi*(500 + 1i*10),'SubspaceDimension',2*n_modes,'Display',1);toc; % eta=10
    Dr       = [Dr zeros(n_modes/2); zeros(n_modes/2) conj(Dr)];
    phi_I_Q  = [phi_I_Q conj(phi_I_Q)];

    % mode selection for eta=10
    mode_selec = [1 2 3 4 5 8 9 10 11 12 13 16];warning('mode selection')
    Dr=Dr(mode_selec,mode_selec);
    phi_I_Q=phi_I_Q(:,mode_selec);
    n_modes=length(mode_selec);




    % define
    eigfreqs    = diag(Dr/(-1i*2*pi));
    eigfreqs_el = abs(diag(Dr/(-1i*2*pi)));


    %display
    fprintf('eigenfrequencies: \n ');
    for i=1:n_modes;fprintf('%f + i%f Hz\n ', real(eigfreqs(i)),imag(eigfreqs(i)));end
end


if matProp{2}.eta==0
n_modes =18; %eta=0 
disp('ITERATIVE SOLVER')



    tic;[phi_I_Q, Dr] = eigs(mu_I_I,lam_I_I ,n_modes/2,-1i*2*pi*(600 ),'SubspaceDimension',2*n_modes,'Display',1);toc; 
    Dr       = [Dr zeros(n_modes/2); zeros(n_modes/2) conj(Dr)];
    phi_I_Q  = [phi_I_Q conj(phi_I_Q)];

    % % mode selection for eta=0 (6 modes from paper)
    % mode_selec = [1 2 3 4 5 8 9 10 11 12 13 16];warning('mode selection')
    % Dr=Dr(mode_selec,mode_selec);
    % phi_I_Q=phi_I_Q(:,mode_selec);
    % n_modes=length(mode_selec);
    
    % mode selection for eta=0 (8 modes from paper)
    mode_selec = [1 2 3 4 5 6 8 9 10 11 13 14 15 17 18];warning('mode selection')
    Dr=Dr(mode_selec,mode_selec);
    phi_I_Q=phi_I_Q(:,mode_selec);
    n_modes=length(mode_selec);



    % define
    eigfreqs    = diag(Dr/(-1i*2*pi));
    eigfreqs_el = abs(diag(Dr/(-1i*2*pi)));


    %display
    fprintf('eigenfrequencies: \n ');
    for i=1:n_modes;fprintf('%f + i%f Hz\n ', real(eigfreqs(i)),imag(eigfreqs(i)));end
end




%% EIGEN:  MODE NORMALIZATION WITH RESPESCT TO THE LAMDA MATRIX
%--------------------------------------------------------------------------
vec_norms = sqrt(diag(phi_I_Q.'*lam_I_I*phi_I_Q));
phi_I_Q_n=phi_I_Q./vec_norms.'; %divide each column of phi by element of vec_norms 
%--------------------------------------------------------------------------
% test
% round(abs(100*phi_I_Q_n.'*lam_I_I*phi_I_Q_n))
% phi_I_Q_n.'* mu_I_I*phi_I_Q_n/(-1i*2*pi)



%% separate degenerate
    
if matProp{2}.eta==10
    mode_deg = [1 2 4 5 7 8 10 11];warning('degenerate mode selection degenerate eta=10, exluding 1069s modes');
    % mode_deg = [1 2 4 5 9 10 12 13];warning('degenerate mode selection degenerate eta=10')
end

if matProp{2}.eta==0
    % mode_deg = [1 2 5 6 7 8 11 12];warning('6 mode selection degenerate eta=0')
    mode_deg = [1 2 5 6 9 10 12 13];warning('8 mode selection degenerate eta=0')
end

    Dr_deg=Dr(mode_deg,mode_deg);
    phi_I_Q_n_deg=phi_I_Q_n(:,mode_deg);


%% second eigs call   
%--------------------------------------------------------------------------
warning('second eigenvalue to  othogonalization')
almost_I   = phi_I_Q_n_deg.'*lam_I_I*phi_I_Q_n_deg;
[phi2,D2]=eig(almost_I); warning('eig call changes the order of the eigenvectors, check phi_new_nT*mu_I_I*phi_new_n / (-2*1i*pi)')
vec_norms = sqrt(diag(phi2.'*almost_I*phi2));
phi2_n=phi2./vec_norms.';
phi_new_n = phi_I_Q_n_deg*phi2_n;

% test
% phi_new_n.'*lam_I_I*phi_new_n
% phi_new_n.'*mu_I_I*phi_new_n / (-2*1i*pi)
%--------------------------------------------------------------------------

% override previous degenerate modes
%--------------------------------------------------------------------------
phi_I_Q_n(:,mode_deg)=phi_new_n; % overrride. Note that the order was reshuffled
Dr(mode_deg,mode_deg)=phi_new_n.'*mu_I_I*phi_new_n; % to account for new order of eigenvalues
eigfreqs = diag(phi_I_Q_n.'* mu_I_I*phi_I_Q_n/(-1i*2*pi));
%--------------------------------------------------------------------------


% % % exclude modes that do not create a bandgap (for a test)
% %--------------------------------------------------------------------------
% mode_bg = [1 2 10 11];warning('bangap mode selection eta=10')
% Dr=Dr(mode_bg,mode_bg);
% phi_I_Q_n=phi_I_Q_n(:,mode_bg);
% eigfreqs = diag(phi_I_Q_n.'* mu_I_I*phi_I_Q_n/(-1i*2*pi));
% n_modes=length(eigfreqs);
% %--------------------------------------------------------------------------

eigfreqs

%% override previous modes with only degenerate non-rotational modes
% %--------------------------------------------------------------------------
% warning('only degenerate non-rotational modes')
% phi_I_Q_n=phi_new_n; % overrride. Note that the order was reshuffled
% Dr=phi_new_n.'*mu_I_I*phi_new_n; % to account for new order of eigenvalues
% eigfreqs = diag(phi_I_Q_n.'* mu_I_I*phi_I_Q_n/(-1i*2*pi));
% n_modes=size(Dr,1);
% %--------------------------------------------------------------------------


eigfreqs

%% PLOT EIGEN: ASSIGN DISP. SHAPE MODES TO TENSOR VECTOR & PLOT
%--------------------------------------------------------------------------
% displacement 'ured' is followed is 'u' displacement reduced 'red' thanks
% to the elimitation of the dependend points due to periodic bc
% 'ured' stands for 
%--------------------------------------------------------------------------

% take only the first displacement column (discard velocities)
%--------------------------------------------------------------------------
phi_plot = real(phi_I_Q_n(1:0.5*size(phi_I_Q_n,1),:));
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% prescribed (reduced) displacements are zero
% declare & assign
ured_p = zeros(1, ee, length(p), size(phi_I_Q_n,2) ); 

        
% free (reduced) displacements are the eigenvectors    
ured_f = phi_plot(1:2:2*length(f)-1,:)*e1+ phi_plot(2:2:2*length(f),:)*e2;


% total reduced displacements (indices ordered as [dofs_pc dofs_un dofs_in])
    % declare ured(unconstrained+corners+independent nodes, modes)
    ured = zeros(1, ee, length(p)+length(f), size(phi_I_Q_n,2));
    % ured = zeros(1, ee, length(p)+length(f), size(phiRed,2));
    % assign
    % assign prescribed nodes
    ured(p,:) = ured_p;
    % assign free nodes
    ured(f,:) = ured_f;
    % free unused variables
    clear ured_f; clear uref_p;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% declare displacement tensor vectors (total # of dofs)
vu     = zeros(1, ee,   ns, size(phi_I_Q_n,2) );
% vu     = zeros(1, ee,   ns, size(phiRed,2) );
vu_aux = zeros(1, ee,   ns, size(phi_I_Q_n,2) ); 
% vu_aux = zeros(1, ee,   ns, size(phiRed,2) );
% assign. Note vu_aux is in order [dofs_c dofs_un dofs_in dofs_de]
vu_aux = T*ured;
% assign to vu. Note vu is in nodes_s order
vu([dofs_c dofs_un dofs_in dofs_de],:)  = vu_aux;
% free unused variables
clear vu_aux;

% prepare for plot: length(vu_plot)>length(vu)
    % declare disp. vector for plot -  whole mesh
    vu_plot = zeros(1, ee,   ns, size(phi_I_Q_n,2) );
    % vu_plot = zeros(1, ee,   ns, size(phiRed,2) );
    % assign vu to vector for plot, only solid nodes are non-zero
    vu_plot( nodes_s,:) = vu;
%--------------------------------------------------------------------------




% Plot
% PLOT EIGEN DISP MODE SHAPES (RESCALED DISP.)
%--------------------------------------------------------------------------
disp('codeblock: PLOT')



% % sfreqs=freqs; %sorted frequencies vector
mode=2; 

% for mode=16:1:length(eigfreqs)

% for mode=16:1:17

% plot the solid displacement in case only 3 solid phases
figure
clf
axis equal
vu_plot_n=vu_plot(:,mode)/norm(norm(vu(:,mode)))*a; % rescale
% vu_plot_n=sign*vu_plot(:,mode);
hold on
femplot(  x          ,   conns, 'Color', 'g');
femplot(  x+vu_plot_n,   connss{1}, 'Color', 'k');
femplot(  x+vu_plot_n,   connss{2}, 'Color', 'r');
femplot(  x+vu_plot_n,   connss{3}, 'Color', 'b');
title("Mode shape of eigenfrequency "  + num2str(eigfreqs(mode)) + " Hz" );
hold off
%--------------------------------------------------------------------------
% end



%% LONGWAVE BASIS

% set xR to centroid
xR = 0*e1+0*e2;
ones_p = ones(length(p) ,1);

%--------------------------------------------------------------------------


% prescribed solid points are [corners_s]
%--------------------------------------------------------------------------
    % MODEL INPUT
    % uM     = 0*(e1+e2);
    % graduM = 0.0 *e1*e1 + 0.1   *e1*e2 +...
    %           0.0 *e2*e1 + 0.1   *e2*e2;
    %  uM     = 1e5*(-1*e1-0*e2);
    % graduM = (0*e1*e1 + 1 *e1*e2 +...
    %           1 *e2*e1 + 0  *e2*e2);

    uM     = 0*a*(e1+0.8*e2);
    graduM = (1 *e1*e1 + 1  *e1*e2 +...
              0.0 *e2*e1 + 1   *e2*e2)/a;
    % first-order comp. homogenization, prescribed disp. at the corners
    % note input index of x() is always "mesh" type.
    % delta xs_p
    dxs_p = x(pcorners)-xR*ones_p;
    % u_p    =  uM *ones_p + dot(graduM.', dxs_p);
    u_p    =  uM *ones_p + dot(dxs_p,graduM);
    % transform to equivalent scalar (no need for sparsicity on a vector)
    su_p   = full(tm2sm(u_p, ee, 1));

if length(p)~=length(  pcorners); error('something wrong');end
%--------------------------------------------------------------------------


% Quasistatic basis based on generalized stiffness matrix -inv(mu_ff)*mu_fp
% same approach as in Bram de Kraker
% S_I_B =  [-sK_f_f\sK_f_p                       zeros(dim*length(f),dim*length(p)); 
%           zeros(dim*length(f),dim*length(p))  -sM_f_f\sM_f_p]; warning('generalized stiffness condensation for long wavelength basis')


% % Purely elastic quasistatic basis: -inv(K_ff)*K_fp
S_f_p = -sK_f_f\sK_f_p;
S_I_B = [S_f_p zeros(size(sK_f_f,1),size(sK_f_p,2)); 
        zeros(size(sK_f_f,1),size(sK_f_p,2)) S_f_p]; warning('purely elastic long wavelength basis');

% Quasistatic basis based on generalized stiffness matrix -inv(mu_ff)*mu_fp
% S_I_B = - mu_I_I\mu_I_B ; warning('generalized stiffness condensation for long wavelength basis')

% Quasistatic basis based on generalized mass matrix -inv(lam_ff)*lam_fp
% S_I_B = - lam_I_I\lam_I_B ; warning('generalized mass condensation for long wavelength basis')


%show matrix structure
% showmatrix(S_I_B)

% PLOT LONGWAVE MODE

% right stationary hybrid state vector
sw_f  =  S_I_B(1:0.5*length(S_I_B),1:length(su_p)) * su_p; %Select only upper left submatrix




% %% STATIONARY: ASSIGN DISP. SOLUTION TO TENSOR VECTOR
%--------------------------------------------------------------------------
% 'ured' is  reduced displacement thanks
% to the elimitation of the dependend points due to periodic bc
%
%  length(ured) = length([dofs_pc; dofs_un; dofs_in])
%  length(u   ) = length([dofs_c ; dofs_un; dofs_in; dofs_de])
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% free (reduced) displacements is the right stationary solution
    % declare
    % nodes are between index 1 to 2*length(f)
    ured_f = sw_f(1:2:2*length(f)-1)*e1+ sw_f(2:2:2*length(f))*e2;
% total reduced displacements (indices ordered as [dofs_un dofs_in])
    % declare ured(ordered: [unconstrained independent] nodes)
    ured = zeros(1, ee, length(p)+length(f), 1);
    % assign prescribed nodes
    ured(p,:) = u_p;
    % assign free nodes
    ured(f,:) = ured_f;
    % free unused variables
    clear ured_f; clear u_p;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% declare displacement tensor vectors (total # of dofs)
vu     = zeros(1, ee,   ns, 1);
vu_aux = zeros(1, ee,   ns, 1); 
% assign. Note vu_aux is in order [dofs_c dofs_un dofs_in dofs_de]
vu_aux = T*ured;
% assign to vu. Note vu is in [nodes_s] order
vu([dofs_c dofs_un dofs_in dofs_de])  = vu_aux;
% free unused variables
clear vu_aux; clear ured;

% prepare for plot: length(vu_plot)>length(vu)
    % declare disp. vector for plot -  whole mesh
    vu_plot = zeros(1, ee,   n, 1);
    % assign vu to vector for plot, only solid nodes are non-zero
    vu_plot(  nodes_s) = vu;

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
 figure(5)
    clf
    axis equal
    hold on
    vu_plot_n= vu_plot/norm(norm(vu))*a; % rescale
    % vu_plot_n= vu_plot; % dont't rescale
    femplot(  x,   conns, 'Color', 'k');
    femplot(  x+vu_plot_n,   conns, 'Color', 'r');
    title("Stationary (check if normalized), x_R=(" ...
           +string(dot(xR,e1))+","+string(dot(xR,e2))+")" );
    hold off
%--------------------------------------------------------------------------

%% REDUCED COUPLED DYNAMIC MODEL
% Comments
%--------------------------------------------------------------------------
% 1. check about transposition .' or hermitian transposition '
% regular transposition is being used.
% 2. notation is, for instance, tlam_B_B means matrix with second order
% tensor components of size = (P,P).A
% t - matrix with second order tensor components
% v - matrix with first order   (vector) tensor components 
% s - matrix with scalar components
%--------------------------------------------------------------------------

% compute reduced matrices on B nodes with Q eigenmodes
%--------------------------------------------------------------------------
% quasistatic part
lam_lw  =  full( ...
                 lam_B_B  +  lam_B_I * S_I_B   +   S_I_B.' * lam_I_B  ...
                + S_I_B.' *  lam_I_I * S_I_B ...
                                                                         );
mu_lw   =  full( ...
                 mu_B_B   +   mu_B_I * S_I_B   +   S_I_B.' *  mu_I_B  ...
                + S_I_B.' *   mu_I_I * S_I_B  ...
                                                                         );
lam_Q_B  =   phi_I_Q_n.' * (lam_I_B +  lam_I_I * S_I_B )  ;
 mu_Q_B  =   phi_I_Q_n.' * ( mu_I_B +   mu_I_I * S_I_B )  ;

lam_B_Q  =  lam_Q_B.';
 mu_B_Q  =   mu_Q_B.';


% elimiante numerical errors (these matrices are already symmetric)
lam_lw  = (lam_lw+lam_lw.')/2;
 mu_lw  = ( mu_lw+ mu_lw.')/2;
%--------------------------------------------------------------------------



% show matrices
%--------------------------------------------------------------------------
caxisrange = [-15 5];
figure
subplot(2,2,1);imagesc(log10(abs(lam_lw))); title('lam_{lw}');  caxis(caxisrange);
subplot(2,2,2);imagesc(log10(abs(mu_lw)));  title('mu_{lw}');   caxis(caxisrange);
subplot(2,2,3);imagesc(log10(abs(lam_Q_B)));  title('lam_{QB}');caxis(caxisrange);
subplot(2,2,4);imagesc(log10(abs(mu_Q_B)));  title('mu_{QB}');  caxis(caxisrange);
h = colorbar('Position', [0.93 0.1 0.02 0.8]); % Adjust position and size of the colorbar
 % Set color axis limits common to all plots
% convert equivalent scalar matrix back to tensor matrix form
%--------------------------------------------------------------------------


% extract submatrices in tensor assemblies
%--------------------------------------------------------------------------
tC_lw   =  sm2tm(lam_lw(1:dim*length(p),1:dim*length(p)),ee,2);
tM_lw   =  sm2tm(lam_lw(1:dim*length(p),dim*length(p)+1:end),ee,2);     % tmu_vv equals -tlam_dotvdotu and -tlam_dotudotv
tK_lw   =  sm2tm(mu_lw(1:dim*length(p),1:dim*length(p)),ee,2);
va      =  sm2tm(lam_Q_B(: , 1:dim*length(p)).',ee,1).'; % transpose because contraction happens on the column
vb      =  sm2tm(lam_Q_B(: , dim*length(p)+1:end).',ee,1).'; % transpose because contraction happens on the column
vc      =  sm2tm(mu_Q_B(: , 1:dim*length(p)).',ee,1).';
I_Q_Q   =    phi_I_Q_n.'*lam_I_I*phi_I_Q_n;
LAM_Q_Q =    phi_I_Q_n.'* mu_I_I*phi_I_Q_n;
% -------------------------------------------------------------------

%% HOMGENIZED COEFFICIENTS (uM,vb,eta) old notation


% %--------------------------------------------------------------------------
% % compute homogenized macroscopic momentum (visco-elastic)
% %--------------------------------------------------------------------------
% % long wavelength
% Api        =   (1/a^2) *  ones_p.' * tC_lw * ones_p;
% Bpi        =   (1/a^2) *  ones_p.' * tC_lw * dxs_p;
% Cpi        =   (1/a^2) * (ones_p.' * tM_lw).' ;
% % Dpi doesnt exist keeping v_b internal variables
% Epi        =   (1/a^2) *  ones_p.' * tK_lw * ones_p;
% Fpi        =   (1/a^2) * (ones_p.' * tK_lw * dxs_p).';
% Gpi        =   0*(1/a^2) * (ones_p.' * tK_lw).' ;    % Gpi is identically zero by construction
% % Hpi doesnt exist keeping v_b internal variables
% 
% % local resonance
% a_Q        =    (1/a^2) * ( ones_p.' * va.' ) .' ;
% c_Q        =    (1/a^2) * ( ones_p.' * vc.'    ) .' ;
% 
% %--------------------------------------------------------------------------
% % compute homogenized macroscopic solid stress properties (visco-elastic)
% %--------------------------------------------------------------------------
% % long wavelength
% As        =   (1/a^2) *  dxs_p.' * tC_lw * ones_p;
% Bs        =   (1/a^2) *  dxs_p.' * tC_lw * dxs_p;
% Cs        =   (1/a^2) * (dxs_p.' * tM_lw).' ;
% % Ds doesnt exist keeping v_b internal variables
% Es        =   (1/a^2) *  dxs_p.' * tK_lw * ones_p;
% Fs        =   (1/a^2) * (dxs_p.' * tK_lw * dxs_p).';
% Gs        = 0*(1/a^2) * (dxs_p.' * tK_lw ).' ; % Gs is identically zero by construction
% % Hs doesnt exist keeping v_b internal variables
% 
% % local resonance
% b_Q       =   1/(a^2) * vb * dxs_p;
% d_Q       =   1/(a^2) * vc * dxs_p;
% 
% 
% 
% %--------------------------------------------------------------------------
% % Constraint Equation
% %--------------------------------------------------------------------------
% 
% 
% % tlam_dotvdotv_oV = 1/(a^2) *tlam_dotvdotv; % tlam_dotvdotv is identically zero by construction
% % tmu_vv_oV        = 1/(a^2) *tmu_vv;
% % vm_dotv_Q_oV     = 1/(a^2) *vm_dotv_Q;     
% % vk_v_Q_oV        = 1/(a^2) *vk_v_Q;
% 
% CpiV  =    (ones_p.' * tM_lw).' ;
% CsV   =    (dxs_p.'  * tM_lw).' ;
% GpiV  =    (ones_p.' * tmu_uv).' ;    % GpiV is identically zero by construction
% GsV   =    (dxs_p.'  * tmu_uv).' ;    % GsV is identically zero by construction
% 
% 
% %--------------------------------------------------------------------------
% % Evolution equation
% %--------------------------------------------------------------------------
% I_Q_Q_oV    =   1/(a^2) * phi_I_Q_n.'*lam_I_I*phi_I_Q_n;
% LAM_Q_Q_oV  =   1/(a^2) * phi_I_Q_n.'* mu_I_I*phi_I_Q_n;
% % I_Q_Q_oV    =  1/(a^2) * eye(n_modes);
% % LAM_Q_Q_oV  = 1/(a^2) * Dr;
% 
% 
% alpha = (a^2) * a_Q;
% beta  = (a^2) * b_Q;
% gamma = (a^2) * c_Q;
% delta = (a^2) * d_Q;
% 
% O_Q_2    =    zeros(length(diag(Dr)),2);
% O_vv     =    zeros(dim*length(p));
% %--------------------------------------------------------------------------


%% HOMGENIZED COEFFICIENTS (uM,vb,eta) new notation (paper)


% compute effective material coefficients: contraction due to averaging
%--------------------------------------------------------------------------
% momentum
trho_col_T      =   (1/a^2) * ones_p.' * tM_lw         ; warning('trho_col is a column')
trho_mat        =   (1/a^2) *            tM_lw         ; warning('trho_mat is a matrix')

% stress
t4eta_LT   =   (1/a^2) *  dxs_p.' * tC_lw *  dxs_p;
t4C_LT     =   (1/a^2) *  dxs_p.' * tK_lw *  dxs_p;
% evolution eqs
va1        =     ( va * ones_p ) ;
tax        =     ( va *  dxs_p ) ;
vb1        =     ( vb * ones_p ) ;
% tbx        =     ( vb *  dxs_p ) ;    
% vc1        =     ( vc    * ones_p ) ;
% tcx        =     ( vc    *  dxs_p ) ;

% divide by volume for disperison computation
va1_V          =     (1/a^2) *( va * ones_p ) ;
tax_V          =     (1/a^2) *( va *  dxs_p ) ;
va_V           =     (1/a^2) *( va ) ;
vb_V           =     (1/a^2) *( vb ) ;
I_Q_Q_V        =     (1/a^2) * I_Q_Q;
LAM_Q_Q_V      =     (1/a^2) *LAM_Q_Q;
%--------------------------------------------------------------------------


trho      =   (1/a^2) * ones_p.' * tM_lw * ones_p;

%% HOMGENIZED COEFFICaIENTS (uM,vM,eta) NEW NOTATION

% warning('commented codeblock')

%--------------------------------------------------------------------------
% compute homogenized macroscopic momentum (visco-elastic)
%--------------------------------------------------------------------------
% momentum
teta      =   (1/a^2) * ones_p.' * tC_lw * ones_p;
trho      =   (1/a^2) * ones_p.' * tM_lw * ones_p;
tC        =   (1/a^2) * ones_p.' * tK_lw * ones_p;

% stress
t4eta_LT   =   (1/a^2) *  dxs_p.' * tC_lw *  dxs_p;
t3rho_LT   =   (1/a^2) *  dxs_p.' * tM_lw * ones_p;
t4C_LT     =   (1/a^2) *  dxs_p.' * tK_lw *  dxs_p;

% evolution eqs
va1        =     ( va * ones_p ) ;
tax        =     ( va *  dxs_p ) ;
vb1        =     ( vb * ones_p ) ;
tbx        =     ( vb *  dxs_p ) ;    
vc1        =     ( vc    * ones_p ) ;
tcx        =     ( vc    *  dxs_p ) ;

% divide by volume for disperison computation
va1_V        =     (1/a^2) *( va * ones_p ) ;
tax_V        =     (1/a^2) *( va *  dxs_p ) ;
vb1_V        =     (1/a^2) *( vb * ones_p ) ;
tbx_V        =     (1/a^2) *( vb *  dxs_p ) ;    
vc1_V        =     (1/a^2) *( vc    * ones_p ) ;
tcx_V        =     (1/a^2) *( vc    *  dxs_p ) ;
I_Q_Q_V      =     (1/a^2) * I_Q_Q;
LAM_Q_Q_V    =     (1/a^2) * LAM_Q_Q;

warning('check: teta << t4eta, t3rho << trho and tC<<t4C ')


%% HOMOGENIZED COEFFICIENTS STATE SPACE 3

if cases==3
% check the following matrices are the same (reduced stiffness)
% lam_lw(1:dim*length(p),1:dim*length(p))
% mu_lw(1:dim*length(p),dim*length(p)+1:end)
% mu_lw(dim*length(p)+1:end,1:dim*length(p))

% then the main homogenized properties are
t4C_LT     =   (1/a^2) *  dxs_p.' * tmu_uv *  dxs_p;
t4eta_LT   =   (1/a^2) *  dxs_p.' * tmu_vv *  dxs_p;
trho       =   (1/a^2) * ones_p.' * tlam_dotvdotv * ones_p;

va        =    (1/a^2) * ( va * ones_p ) ;
vb        =    (1/a^2) * ( vb * ones_p ) ;  
% evolution eq.
I_Q_Q_oV   =   (1/a^2) * phi_I_Q_n.'*lam_I_I*phi_I_Q_n;
LAM_Q_Q_oV =   (1/a^2) * phi_I_Q_n.'* mu_I_I*phi_I_Q_n;
end
%% PLOT HOMOGENIZED STIFFNESS AND VISCOSITY

if cases==2
    rhoa = tm2sm((matProp{2}.rho)*I,ee,2);
    Cfour   =  matProp{2}.kappa*I*I + 2* matProp{2}.G*(  I4S - 1/3*  I*  I);
    etafour = matProp{2}.eta*I*I + 2* matProp{2}.eta*(  I4S - 1/3*  I*  I);


    C4_homog   = tm2sm(t4C_LT,ee,4);
    C4_homog_v = t2voigt(C4_homog ,ee,4)
    C4 = tm2sm(Cfour,ee,4);
    C4_v = t2voigt(C4,ee,4)

    eta4_homog = tm2sm(t4eta_LT,ee,4);
    eta4_homog_v = t2voigt(eta4_homog,ee,4)
    D4 = tm2sm(etafour,ee,4);
    D4_v =t2voigt(D4,ee,4)

    full(rhoa)
    full(tm2sm(trho,ee,2))
  
    warning(['microscale quantities compared to homogenized quantities (only works for ' ...
    'homogeneous square unit cell']);
 
end



%% WEIGHTED DENSITIES
%--------------------------------------------------------------------------
% computed with analytical expression
r3 = 5e-3;
r2 = 7.5e-3;
% a  = 20e-3;
V3 = pi*r3^2;
V2 = pi*r2^2-V3;
V1 = a^2-V2-V3;
rho_weighted_ava=(V1*matProp{1}.rho+V2*matProp{2}.rho+V3*matProp{3}.rho)/a^2;
% computed with area of mesh
rho_weighted_av=(V{1}*matProp{1}.rho+V{2}*matProp{2}.rho+V{3}*matProp{3}.rho)/a^2;
%--------------------------------------------------------------------------
 


%% test anisotropy for stiffness

% stiffness
C4  = tm2sm(ltranspose(t4C_LT),ee,4);
C4v = t2voigt(C4 ,ee,4)
% C1111_over_C2121= C4v(1,1)/C4v(3,3)

% % if the following bulk modulus are different then it is not isotropic RVE
BulkModulus1111 = C4v(1,1)-4/3*C4v(3,3)
BulkModulus1122 = C4v(1,2)+2/3*C4v(3,3)

%% test anisotropy for viscosity

% stiffness
eta4  = tm2sm(ltranspose(t4eta_LT),ee,4);
eta4v = t2voigt(eta4 ,ee,4)

% % if the following bulk modulus are different then it is not isotropic RVE
BulkModulus1111 = eta4v(1,1)-4/3*eta4v(3,3)
BulkModulus1122 = eta4v(1,2)+2/3*eta4v(3,3)

%% coupling strength
% % dipolar
% normvb1=abs(norm(vb1))
% [~,ind]=sort(normvb1),'descend')
% [eigfreqs(ind) abs(norm(vb1(ind)))]

% % monopolar
% normtbx=sqrt( abs(dot(e1,tbx,e1)).^2 + abs(dot(e1,tbx,e2)).^2 + abs(dot(e2,tbx,e1)).^2 + abs(dot(e2,tbx,e2)).^2   )
% [~,ind]=sort(normvbx),'descend')
% [eigfreqs(ind) normtbx(ind)]