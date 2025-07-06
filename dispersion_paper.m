% R. Liupekevicius Carnielli 18-09-2024
% paper plots on enriched continuum with 8+N dofs
% of which 1 displacement, 3 corner velocities and N internal variables
% for the modes, i.e. dep. variables uM,v1,v2,v4,xi1,...,xiN

%% DISPERSION BLOCH ANALYSIS COMSOL (REAL)

% figure(23)
hold on
grid on
% scatter(real(dispersion_bloch(:,1)),real(dispersion_bloch(:,2)),20);
% % ylim([0 1000])
scatter(real(dispersion_bloch(:,1)),(real(dispersion_bloch(:,2))), 'Marker', '+', 'MarkerEdgeColor', 'blue', 'LineWidth', 1.5);
ylim([0 1450])
fontsize(gca, 21,'pixels')
% xlim([-1 2])
box on


%% DISPERSION BLOCH ANALYSIS COMSOL damping (IMAG) 


% figure(23)
hold on
grid on
% scatter(real(dispersion_bloch(:,1)),real(dispersion_bloch(:,2)),20);
% ylim([0 1000])
scatter(real(dispersion_bloch(:,1)),(imag(dispersion_bloch(:,2)))./(abs(dispersion_bloch(:,2))), 'Marker', '+', 'MarkerEdgeColor', 'blue', 'LineWidth', 2);
% scatter(real(dispersion_bloch(:,1)),(imag(dispersion_bloch(:,2))), 'Marker', '+', 'MarkerEdgeColor', 'black', 'LineWidth', 2);
ylim([0 1])
% xlim([0 0.15])
fontsize(gca, 21,'pixels')
box on

% ylim([0 .25])
%ylim([0.8 1])

%% DISPERSION CURVE CLASSIC LONGWAVE + DOMINANT LOCAL RESONANCE (REAL)



% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------


% number of depend. variables describing the continuum
% (uM,v1,v2,v4,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 8+n_modes; % propagating 
%-------------------------------------------------------------------------      -

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.001;
range_of_dimless_k = -1:kstep:2;
% range_of_dimless_k = -0.05:kstep:0.1;
% range_of_dimless_k = -0.1:kstep:0.2;
% range_of_dimless_k = 0:kstep:1;
% range_of_dimless_k = 0:kstep:0.2;
% range_of_dimless_k = 0:kstep:0.05;
%--------------------------------------------------------------------------

% % wave number normalized by pi/ax
% %--------------------------------------------------------------------------
% Gamma_index=   round(length(range_of_dimless_k)/3)+1;
% X_index    = 2*round(length(range_of_dimless_k)/3)+1;
% 
% % dimless wavenumber normalized by pi/ax
% k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
%                range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
%                range_of_dimless_k(X_index+1:end    )    *ay/ax ];
% %--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end

    
    
    % matrices -- symmetric --

    Muu = full(tm2sm(dot( k , t4eta_LT , k ),ee,2)) ; 
    Muv = full(tm2sm( trho_col_T ,ee,2))            ;
    Mue = full(tm2sm( va1_V.' ,ee,1))               ;
    
    Mvu = Muv.'                  ; 
    Mvv = zeros(6)               ; % Mvv are indentically zero by construction
    Mve = full(tm2sm( vb_V.' ,ee,1)) ;     

    Meu = Mue.';
    Mev = Mve.';
    Mee = I_Q_Q_V; 

    %stiffness matrix
    Kuu = full(tm2sm(dot( k , t4C_LT ,k),ee,2)) ; 
    Kuv = zeros(2,6); 
    Kue = zeros(2,n_modes);
    
    Kvu =   Kuv.'                     ; 
    Kvv = - full(tm2sm( trho_mat ,ee,2)) ;
    Kve = - full(tm2sm( va_V.' ,ee,1)) ; 
   
    
    Keu = Kue.';
    Kev = Kve.';
    Kee = LAM_Q_Q_V;



    M = full([Muu Muv Mue; Mvu Mvv Mve; Meu Mev Mee]);
    K = full([Kuu Kuv Kue; Kvu Kvv Kve; Keu Kev Kee]);

   % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % % sort d and V with respect to freqs frequencies
    % %--------------------------------------------------------------------------
    % [~,Ind]         = sort(diag(-imag(eigenvalues)));
    % eigenvalues     = eigenvalues(Ind,Ind);
    % eigenvectors    = eigenvectors(:,Ind);
    % %--------------------------------------------------------------------------
    
    % eigen frequencies
    %--------------------------------------------------------------------------
    lambda(:,i) = diag(eigenvalues);
    %--------------------------------------------------------------------------

    i=i+1;
    end
%--------------------------------------------------------------------------
% plot
%--------------------------------------------------------------------------

%define eigenfreqs and discard first two branches (negative freq lw and 0 freq)
%--------------------------------------------------------------------------
eigenfrequencies = lambda/(-1i*2*pi);
%--------------------------------------------------------------------------


%plot one color
%--------------------------------------------------------------------------
% figure(23)
hold on
grid on
for i=1:number_of_waves
% plot(range_of_dimless_k,real(eigenfrequencies(i,:)),"color", 'b', "linewidth", 2)    
scatter(range_of_dimless_k,real(eigenfrequencies(i,:)),'Marker', '.', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
end
hold off
ylim([0 1450])
%--------------------------------------------------------------------------

%% DISPERSION CURVE CLASSIC LONGWAVE + DOMINANT LOCAL RESONANCE purely elastic ionlcuding A (REAL) test 1500 Hz mode




% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------


% number of depend. variables describing the continuum
% (uM,v1,v2,v4,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 8+n_modes; % propagating 
%-------------------------------------------------------------------------      -

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.001;
range_of_dimless_k = -1:kstep:2;
% range_of_dimless_k = -0.05:kstep:0.1;
% range_of_dimless_k = -0.1:kstep:0.2;
% range_of_dimless_k = 0:kstep:1;
% range_of_dimless_k = 0:kstep:0.2;
% range_of_dimless_k = 0:kstep:0.05;
%--------------------------------------------------------------------------

% % wave number normalized by pi/ax
% %--------------------------------------------------------------------------
% Gamma_index=   round(length(range_of_dimless_k)/3)+1;
% X_index    = 2*round(length(range_of_dimless_k)/3)+1;
% 
% % dimless wavenumber normalized by pi/ax
% k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
%                range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
%                range_of_dimless_k(X_index+1:end    )    *ay/ax ];
% %--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end

    
    
    % matrices -- symmetric --

    Muu = full(tm2sm(dot( k , t4eta_LT , k ),ee,2)) ; 
    Muv = full(tm2sm( trho_col_T ,ee,2))            ;
    Mue = full(tm2sm( va1_V.'  ,ee,1)) + 1i*tm2sm(dot(k,tax_V.'),ee,1)   ;
    
    Mvu = Muv.'                  ; 
    Mvv = zeros(6)               ; % Mvv are indentically zero by construction
    Mve = full(tm2sm( vb_V.' ,ee,1)) ;     

    Meu = Mue.';
    Mev = Mve.';
    Mee = I_Q_Q_V; 

    %stiffness matrix
    Kuu = full(tm2sm(dot( k , t4C_LT ,k),ee,2)) ; 
    Kuv = zeros(2,6); 
    % Kue = zeros(2,n_modes);
    Kue = full(tm2sm( vc1_V.'  ,ee,1)) + 1i*tm2sm(dot(k,tcx_V.'),ee,1) ;

    Kvu =   Kuv.'                     ; 
    Kvv = - full(tm2sm( trho_mat ,ee,2)) ;
    Kve = - full(tm2sm( va_V.' ,ee,1)) ; 
   
    
    Keu = Kue.';
    Kev = Kve.';
    Kee = LAM_Q_Q_V;



    M = full([Muu Muv Mue; Mvu Mvv Mve; Meu Mev Mee]);
    K = full([Kuu Kuv Kue; Kvu Kvv Kve; Keu Kev Kee]);

   % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % % sort d and V with respect to freqs frequencies
    % %--------------------------------------------------------------------------
    % [~,Ind]         = sort(diag(-imag(eigenvalues)));
    % eigenvalues     = eigenvalues(Ind,Ind);
    % eigenvectors    = eigenvectors(:,Ind);
    % %--------------------------------------------------------------------------
    
    % eigen frequencies
    %--------------------------------------------------------------------------
    lambda(:,i) = diag(eigenvalues);
    %--------------------------------------------------------------------------

    i=i+1;
    end
%--------------------------------------------------------------------------
% plot
%--------------------------------------------------------------------------

%define eigenfreqs and discard first two branches (negative freq lw and 0 freq)
%--------------------------------------------------------------------------
eigenfrequencies = lambda/(-1i*2*pi);
%--------------------------------------------------------------------------


%plot one color
%--------------------------------------------------------------------------
% figure(23)
hold on
grid on
for i=1:number_of_waves
% plot(range_of_dimless_k,real(eigenfrequencies(i,:)),"color", 'b', "linewidth", 2)    
scatter(range_of_dimless_k,real(eigenfrequencies(i,:)),'Marker', '.', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
end
hold off
ylim([0 1450])
%--------------------------------------------------------------------------

%% DISPERSION CURVE CLASSIC LONGWAVE + DOMINANT LOCAL RESONANCE (IMAG)



% RVE geometry
%--------------------------------------------------------------------------
ax=a;
ay=a;
%--------------------------------------------------------------------------


% number of depend. variables describing the continuum
% (uM,v1,v2,v4,eta1,...,etaQ)
%--------------------------------------------------------------------------
number_of_waves = 8+n_modes; % propagating 
%-------------------------------------------------------------------------      -

% dimensionless wave number definition (for loop)
%--------------------------------------------------------------------------
kstep=0.001;
% range_of_dimless_k = -1:kstep:2;
range_of_dimless_k = -0.1:kstep:0.2;
% range_of_dimless_k = 0:kstep:1;
% range_of_dimless_k = 0:kstep:0.2;
% range_of_dimless_k = 0:kstep:0.05;
%--------------------------------------------------------------------------

% % wave number normalized by pi/ax
% %--------------------------------------------------------------------------
% Gamma_index=   round(length(range_of_dimless_k)/3)+1;
% X_index    = 2*round(length(range_of_dimless_k)/3)+1;
% 
% % dimless wavenumber normalized by pi/ax
% k_plot     = [ range_of_dimless_k(1        :Gamma_index)*sqrt(ax^2+ay^2)/ax ...
%                range_of_dimless_k(Gamma_index+1:X_index)*ax/ax ...
%                range_of_dimless_k(X_index+1:end    )    *ay/ax ];
% %--------------------------------------------------------------------------

% initialize zero-matrix for storaging the eigenvalues of given k
%--------------------------------------------------------------------------
lambda = zeros(number_of_waves,length(range_of_dimless_k)); % 
%--------------------------------------------------------------------------


% loop in k
%--------------------------------------------------------------------------
    i=1;
    for dimless_k = range_of_dimless_k
    
    % wavenumber
    if(dimless_k<0)
    k =  - dimless_k    *pi/ax *e1 -  dimless_k    *pi/ay *e2;
    elseif(dimless_k<1)
    k =    dimless_k    *pi/ax *e1 +                 0    *e2;
    else
    k =          1      *pi/ax *e1 + (dimless_k-1) *pi/ay *e2;
    end

    
    
    % matrices -- symmetric --

    Muu = full(tm2sm(dot( k , t4eta_LT , k ),ee,2)) ; 
    Muv = full(tm2sm( trho_col_T ,ee,2))            ;
    Mue = full(tm2sm( va1_V.' ,ee,1))               ;
    
    Mvu = Muv.'                  ; 
    Mvv = zeros(6)               ; % Mvv are indentically zero by construction
    Mve = full(tm2sm( vb_V.' ,ee,1)) ;     

    Meu = Mue.';
    Mev = Mve.';
    Mee = I_Q_Q_V; 

    %stiffness matrix
    Kuu = full(tm2sm(dot( k , t4C_LT ,k),ee,2)) ; 
    Kuv = zeros(2,6); 
    Kue = zeros(2,n_modes);
    
    Kvu =   Kuv.'                     ; 
    Kvv = - full(tm2sm( trho_mat ,ee,2)) ;
    Kve = - full(tm2sm( va_V.' ,ee,1)) ; 
   
    
    Keu = Kue.';
    Kev = Kve.';
    Kee = LAM_Q_Q_V;



    M = full([Muu Muv Mue; Mvu Mvv Mve; Meu Mev Mee]);
    K = full([Kuu Kuv Kue; Kvu Kvv Kve; Keu Kev Kee]);

   % compute eigevalues and eigenvectors
    %--------------------------------------------------------------------------
    [eigenvectors,eigenvalues] = eig(K,M);
    %--------------------------------------------------------------------------

    % % sort d and V with respect to freqs frequencies
    % %--------------------------------------------------------------------------
    % [~,Ind]         = sort(diag(-imag(eigenvalues)));
    % eigenvalues     = eigenvalues(Ind,Ind);
    % eigenvectors    = eigenvectors(:,Ind);
    % %--------------------------------------------------------------------------
    
    % eigen frequencies
    %--------------------------------------------------------------------------
    lambda(:,i) = diag(eigenvalues);
    %--------------------------------------------------------------------------

    i=i+1;
    end
%--------------------------------------------------------------------------
% plot
%--------------------------------------------------------------------------

%define eigenfreqs and discard first two branches (negative freq lw and 0 freq)
%--------------------------------------------------------------------------
eigenfrequencies = lambda/(-1i*2*pi);
%--------------------------------------------------------------------------


%plot one color
%--------------------------------------------------------------------------
% figure(23)
hold on
grid on
for i=1:number_of_waves
scatter(range_of_dimless_k,imag(eigenfrequencies(i,:))./abs(eigenfrequencies(i,:)),'Marker', '.', 'MarkerEdgeColor', [0.4 0.4 0.4], 'LineWidth', 2);
end
hold off
%--------------------------------------------------------------------------