% DESCRIPTION
% Renan Liupekevicius Carnielli TU/e
% start in 27-08-2023

% Compute frequency-dependent material properties.
V = a^2; warning('overrriding the struct V with phase volumes')



% %% EFFECTIVE MACRO DENSITIES rho1, rho2, rho3 and rho4 (12 MODES)
% 
% % extract the local resonant modes
% omega_lr = diag(phi_I_Q_n.'*mu_I_I*phi_I_Q_n /(-1i));
% 
% % def frequency range
% f     = (0:1:1450)'; %column
% % f     = (0:1:10084)'; %column
% omega = 2*pi*f;
% 
% 
% % using for loop is less efficient but for variable number of modes. Note
% % that the implementation is vectorized for frequency axis and the for
% % loop is only for the modes.
% 
% % initialize densities with frequency independent terms
% rho1 = zeros(2, ee, length(omega), 1); 
% rho2 = trho*ones(length(omega),1); 
% rho4 = trho*ones(length(omega),1);
% 
% for j=1:length(omega_lr) %number o modes n_modes
% 
% rho1 = rho1 + 1/V./( 1i*omega-1i*omega_lr(j) )   * (-1)       *  va1(j)*va1(j) ;
% 
% rho2 = rho2 + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)*  va1(j)*vb1(j) ...       
%             + 1/V./( 1i*omega-1i*omega_lr(j) )                *  va1(j)*va1(j);
% 
% rho4 = rho4 + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-omega.^2)*  vb1(j)*vb1(j) ...
%             + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)* (va1(j)*vb1(j)+vb1(j)*va1(j)) ...
%             + 1/V./( 1i*omega-1i*omega_lr(j) )                *  va1(j)*va1(j);   
% end
% 
% 
% 
% 
% %% EFFECTIVE RHO using vM (SIMPLIFIED EQS)
% 
% rho_eff_vM = zeros(2, ee, length(omega), 1); 
% 
% for i=1:length(omega) 
%    rho_eff_vM(i) = rho1(i) + dot(rho2(i), inv(rho4(i)) , rho2(i).' ) ;
% end


%% EFFECTIVE rho1(2ND ORDER TENSOR), rho2, rho3 (column of 2ND ORDER TENSORS) and rho4 (matrix assembly of 2ND ORDER TENSORS)

% local resonant modes
omega_lr = diag(phi_I_Q_n.'*mu_I_I*phi_I_Q_n /(-1i));

% frequency axis
f     = (0:1:1450); %line
% f     = (0:1:4000); %line
% f     = (1500:0.1:1600); %line

omega = 2*pi*f;
% f     = f.'    ;
% omega = omega.';

% using for loop is less efficient but for variable number of modes. Note
% that the implementation is vectorized for frequency axis and the for
% loop is only for the modes.

% initialize tensor rho1, tensor column rho2/3 n tensor matrix assembly rho4
% initialize densities with frequency independent terms

% rho1 is a 2nd order tensor for each omega
rho_  =    zeros(2, ee, 1,length(omega) );

% % rho2 is 1x3 column for each omega
% rho_1  = trho*ones( 1,length(omega) );
% rho_2  = trho*ones( 1,length(omega) );
% rho_3  = trho*ones( 1,length(omega) );
% 
% % rho4 is 3x3 symmetric matrix for each omega
% rho_11 = trho*ones( 1,length(omega) );
% rho_12 = trho*ones( 1,length(omega) );
% rho_13 = trho*ones( 1,length(omega) );
% rho_22 = trho*ones( 1,length(omega) );
% rho_23 = trho*ones( 1,length(omega) );
% rho_33 = trho*ones( 1,length(omega) );

% rho2 is 1x3 column for each omega
rho_1  = trho_col_T(1)*ones( 1,length(omega) );
rho_2  = trho_col_T(2)*ones( 1,length(omega) );
rho_3  = trho_col_T(3)*ones( 1,length(omega) );

% rho4 is 3x3 symmetric matrix for each omega
rho_11 = trho_mat(1,1)*ones( 1,length(omega) );
rho_12 = trho_mat(1,2)*ones( 1,length(omega) );
rho_13 = trho_mat(1,3)*ones( 1,length(omega) );
rho_22 = trho_mat(2,2)*ones( 1,length(omega) );
rho_23 = trho_mat(2,3)*ones( 1,length(omega) );
rho_33 = trho_mat(3,3)*ones( 1,length(omega) );




for j=1:length(omega_lr) %loop on number o modes n_modes

% rho1 <-paper notation
rho_      = rho_ ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )   * (-1)*         va1(j)*va1(j) ;

% undetilde{rho2} components <-paper notation
rho_1  = rho_1 ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)*  va1(j)*vb(j,1) ...       
           + 1/V./( 1i*omega-1i*omega_lr(j) )   *               va1(j)*va(j,1);

rho_2  = rho_2 ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)*  va1(j)*vb(j,2) ...       
           + 1/V./( 1i*omega-1i*omega_lr(j) )   *               va1(j)*va(j,2);

rho_3  = rho_3 ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)*  va1(j)*vb(j,3) ...       
           + 1/V./( 1i*omega-1i*omega_lr(j) )   *               va1(j)*va(j,3);

% undetline{rho4} components <-paper notation
rho_11 = rho_11 ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-omega.^2)*  vb(j,1)*vb(j,1) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)*( va(j,1)*vb(j,1) + vb(j,1)*va(j,1) ) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )   *               va(j,1)*va(j,1);

rho_12 = rho_12 ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-omega.^2)*  vb(j,1)*vb(j,2) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)*( va(j,1)*vb(j,2) + vb(j,1)*va(j,2) ) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )   *               va(j,1)*va(j,2);   

rho_13 = rho_13 ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-omega.^2)*  vb(j,1)*vb(j,3) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)*( va(j,1)*vb(j,3) + vb(j,1)*va(j,3) ) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )   *               va(j,1)*va(j,3);   

rho_22 = rho_22 ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-omega.^2)*  vb(j,2)*vb(j,2) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)*( va(j,2)*vb(j,2) + vb(j,2)*va(j,2) ) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )   *               va(j,2)*va(j,2);  

rho_23 = rho_23 ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-omega.^2)*  vb(j,2)*vb(j,3) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)*( va(j,2)*vb(j,3) + vb(j,2)*va(j,3) ) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )   *               va(j,2)*va(j,3);   

rho_33 = rho_33 ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-omega.^2)*  vb(j,3)*vb(j,3) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )  .* (-1i*omega)*( va(j,3)*vb(j,3) + vb(j,3)*va(j,3) ) ...
           + 1/V./( 1i*omega-1i*omega_lr(j) )   *               va(j,3)*va(j,3);   

end




%% EFFECTIVE RHO


% inv rho4

rho_eff_vb = zeros(2, ee, length(omega), 1); 

for i=1:length(omega)
   
   rho_line_i      = [rho_1(i) rho_2(i) rho_3(i)];
   rho_matr_i      = [rho_11(i)   rho_12(i)   rho_13(i)
                      rho_12(i).' rho_22(i)   rho_23(i)
                      rho_13(i).' rho_23(i).' rho_33(i)];

   inv_rho_matr_i  = sm2tm( inv(tm2sm( rho_matr_i, ee , 2 )) , ee ,2);
              

   rho_eff_vb(i) = rho_(i) + dot(rho_line_i, inv_rho_matr_i , rho_line_i.' ) ;
   % rho_eff(i) = rho1(i) + dot(rho2(i), inv(rho4(i)) , rho2(i).' ) -rho1(1) ;

end



%% Plot vM and vb

% 
% figure
% hold on
% plot(f,real(dot(e1,rho_eff_vb,e1)),'b',LineWidth=1.5)
% plot(f,imag(dot(e1,rho_eff_vb,e1)),'b--',LineWidth=1.5)
% plot(f,real(dot(e1,rho_eff_vM,e1)),'r')
% plot(f,imag(dot(e1,rho_eff_vM,e1)),'r--')
% % ylim([-100 100])
% 
% hold off
% grid on
% legend('real vb','imaginary vb','real vM','imaginary vM') 
% title('effective macroscopic density 11')
% 
% 
% figure
% hold on
% plot(f,real(dot(e1,rho_eff_vb,e2)),'b',LineWidth=1.5)
% plot(f,imag(dot(e1,rho_eff_vb,e2)),'b--',LineWidth=1.5)
% plot(f,real(dot(e1,rho_eff_vM,e2)),'r')
% plot(f,imag(dot(e1,rho_eff_vM,e2)),'r--')
% 
% % ylim([-100 100])
% hold off
% grid on
% legend('real vb','imaginary vb','real vM','imaginary vM') 
% title('effective macroscopic density 12')
% 
% figure
% hold on
% plot(f,real(dot(e2,rho_eff_vb,e1)),'b',LineWidth=1.5)
% plot(f,imag(dot(e2,rho_eff_vb,e1)),'b--',LineWidth=1.5)
% plot(f,real(dot(e2,rho_eff_vM,e1)),'r')
% plot(f,imag(dot(e2,rho_eff_vM,e1)),'r--')
% 
% % ylim([-100 100])
% hold off
% grid on
% legend('real vb','imaginary vb','real vM','imaginary vM') 
% title('effective macroscopic density 21')
% 
% figure
% hold on
% plot(f,real(dot(e2,rho_eff_vb,e2)),'b',LineWidth=1.5)
% plot(f,imag(dot(e2,rho_eff_vb,e2)),'b--',LineWidth=1.5)
% plot(f,real(dot(e2,rho_eff_vM,e2)),'r')
% plot(f,imag(dot(e2,rho_eff_vM,e2)),'r--')
% hold off
% grid on
% legend('real vb','imaginary vb','real vM','imaginary vM') 
% title('effective macroscopic density 22')
% 
% 
% 
% 
% %% Plot vM 
% 
% 
% figure
% plot(f,real(dot(e1,rho_eff_vM,e1)),'r')
% hold on
% plot(f,imag(dot(e1,rho_eff_vM,e1)),'r--')
% % ylim([-100 100])
% hold off
% grid on
% legend('real vM','imaginary vM') 
% title('effective macroscopic density 11')
% 
% 
% figure
% plot(f,real(dot(e1,rho_eff_vM,e2)),'r')
% hold on
% plot(f,imag(dot(e1,rho_eff_vM,e2)),'r--')
% % ylim([-100 100])
% hold off
% grid on
% legend('real vM','imaginary vM') 
% title('effective macroscopic density 12')
% 
% figure
% plot(f,real(dot(e2,rho_eff_vM,e1)),'r')
% hold on
% plot(f,imag(dot(e2,rho_eff_vM,e1)),'r--')
% % ylim([-100 100])
% hold off
% grid on
% legend('real vM','imaginary vM') 
% title('effective macroscopic density 21')
% 
% figure
% plot(f,real(dot(e2,rho_eff_vM,e2)),'r')
% hold on
% plot(f,imag(dot(e2,rho_eff_vM,e2)),'r--')
% hold off
% grid on
% legend('real vM','imaginary vM')  
% title('effective macroscopic density 22')

%% Plot vb


figure
hold on
% ylim([-100 100])
plot(f,real(dot(e1,rho_eff_vb,e1)),'r')
plot(f,imag(dot(e1,rho_eff_vb,e1)),'r--')
hold off
grid on
% legend('real vb','imaginary vb') 
title('effective macroscopic density 11')


figure
hold on
plot(f,real(dot(e1,rho_eff_vb,e2)),'r')
plot(f,imag(dot(e1,rho_eff_vb,e2)),'r--')
% ylim([-100 100])
hold off
grid on
legend('real vb','imaginary vb') 
title('effective macroscopic density 12')

figure
hold on
plot(f,real(dot(e2,rho_eff_vb,e1)),'r')
plot(f,imag(dot(e2,rho_eff_vb,e1)),'r--')
% ylim([-100 100])
hold off
grid on
legend('real vb','imaginary vb') 
title('effective macroscopic density 21')

figure
hold on
plot(f,real(dot(e2,rho_eff_vb,e2)),'r')
plot(f,imag(dot(e2,rho_eff_vb,e2)),'r--')
hold off
grid on
legend('real vb','imaginary vb') 
title('effective macroscopic density 22')


%% Plot vb paper


figure
hold on
ylim([0 1450])
plot(real(dot(e1,rho_eff_vb,e1)),f,'r',"LineWidth",2)
plot(imag(dot(e1,rho_eff_vb,e1)),f,'r--',"LineWidth",2)
hold off
grid on
% legend('real vb','imaginary vb') 
title('effective macroscopic density 11')
xlim([-1.5e4 1e4])
fontsize(gca, 21,'pixels')
box on


%% write

writematrix([f' real(dot(e1,rho_eff_vb,e1))],'f_rhoR.txt', 'Delimiter','space')
writematrix([f' imag(dot(e1,rho_eff_vb,e1))],'f_rhoI.txt', 'Delimiter','space')
writematrix([f' real(dot(e1,rho_eff_vb,e2))],'f_rho12R.txt', 'Delimiter','space')
writematrix([f' imag(dot(e1,rho_eff_vb,e2))],'f_rho12I.txt', 'Delimiter','space')

% writematrix([f' real(dot(e1,rho_eff_vb,e1))],'f_rhoR_67.txt', 'Delimiter','space')
% writematrix([f' imag(dot(e1,rho_eff_vb,e1))],'f_rhoI_67.txt', 'Delimiter','space')

% writematrix([f' real(dot(e1,rho_eff_vb,e1))],'f_rhoR.txt', 'Delimiter','space')



%% EFFECTIVE MACRO DENSITIES (2ND ORDER TENSORS) rho1, rho2, rho3 and rho4 (4 MODES)

% bottle neck is still that for loop on omega so improvements below are
% uselless


% omega_lr = diag(phi_I_Q_n.'*mu_I_I*phi_I_Q_n /(-1i));
% 
% f     = (0:0.1:1084)'; %column
% % f     = (0:1:10084)'; %column
% % 
% omega = 2*pi*f;
% 
% rho1 = +1/V./( 1i*omega-1i*omega_lr(1) )    *  (-1)* va1(1)*va1(1) ...
%        +1/V./( 1i*omega-1i*omega_lr(2) )    *  (-1)* va1(2)*va1(2) ...
%        +1/V./( 1i*omega-1i*omega_lr(3) )    *  (-1)* va1(3)*va1(3) ...
%        +1/V./( 1i*omega-1i*omega_lr(4) )    *  (-1)* va1(4)*va1(4)  ;
% 
% 
% rho2 =  trho*ones(length(omega),1)  ...
%        +1/V./( 1i*omega-1i*omega_lr(1) )   .*  (-1i*omega)* va1(1)*vb1(1)       +    1/V./( 1i*omega-1i*omega_lr(1) )     *  va1(1)*va1(1)...
%        +1/V./( 1i*omega-1i*omega_lr(2) )   .*  (-1i*omega)* va1(2)*vb1(2)       +    1/V./( 1i*omega-1i*omega_lr(2) )     *  va1(2)*va1(2)...
%        +1/V./( 1i*omega-1i*omega_lr(3) )   .*  (-1i*omega)* va1(3)*vb1(3)       +    1/V./( 1i*omega-1i*omega_lr(3) )     *  va1(3)*va1(3)...
%        +1/V./( 1i*omega-1i*omega_lr(4) )   .*  (-1i*omega)* va1(4)*vb1(4)       +    1/V./( 1i*omega-1i*omega_lr(4) )     *  va1(4)*va1(4) ;
% 
% 
% rho4 =  trho*ones(length(omega),1)  ...
%        +1/V./( 1i*omega-1i*omega_lr(1) )   .*(-omega.^2)  *  vb1(1)*vb1(1)                +             1/V./( 1i*omega-1i*omega_lr(1) )  .*(-1i*omega)*(va1(1)*vb1(1)+vb1(1)*va1(1))                +         1/V./( 1i*omega-1i*omega_lr(1) )   *  va1(1)*va1(1) ...
%        +1/V./( 1i*omega-1i*omega_lr(2) )   .*(-omega.^2)  *  vb1(2)*vb1(2)                +             1/V./( 1i*omega-1i*omega_lr(2) )  .*(-1i*omega)*(va1(2)*vb1(2)+vb1(2)*va1(2))                +         1/V./( 1i*omega-1i*omega_lr(2) )   *  va1(2)*va1(2) ...
%        +1/V./( 1i*omega-1i*omega_lr(3) )   .*(-omega.^2)  *  vb1(3)*vb1(3)                +             1/V./( 1i*omega-1i*omega_lr(3) )  .*(-1i*omega)*(va1(3)*vb1(3)+vb1(3)*va1(3))                +         1/V./( 1i*omega-1i*omega_lr(3) )   *  va1(3)*va1(3) ...
%        +1/V./( 1i*omega-1i*omega_lr(4) )   .*(-omega.^2)  *  vb1(4)*vb1(4)                +             1/V./( 1i*omega-1i*omega_lr(4) )  .*(-1i*omega)*(va1(4)*vb1(4)+vb1(4)*va1(4))                +         1/V./( 1i*omega-1i*omega_lr(4) )   *  va1(4)*va1(4) ;
% 

% Optimized implementation for 12 modes (very quick)

% rho1 = +1/V./( 1i*omega-1i*omega_lr(1) )    *  (-1)* va1(1)*va1(1) ...
%        +1/V./( 1i*omega-1i*omega_lr(2) )    *  (-1)* va1(2)*va1(2) ...
%        +1/V./( 1i*omega-1i*omega_lr(3) )    *  (-1)* va1(3)*va1(3) ...
%        +1/V./( 1i*omega-1i*omega_lr(4) )    *  (-1)* va1(4)*va1(4) ...
%        +1/V./( 1i*omega-1i*omega_lr(5) )    *  (-1)* va1(5)*va1(5) ...
%        +1/V./( 1i*omega-1i*omega_lr(6) )    *  (-1)* va1(6)*va1(6) ...
%        +1/V./( 1i*omega-1i*omega_lr(7) )    *  (-1)* va1(7)*va1(7) ...
%        +1/V./( 1i*omega-1i*omega_lr(8) )    *  (-1)* va1(8)*va1(8) ...
%        +1/V./( 1i*omega-1i*omega_lr(9) )    *  (-1)* va1(9)*va1(9) ...
%        +1/V./( 1i*omega-1i*omega_lr(10) )    *  (-1)* va1(10)*va1(10) ...
%        +1/V./( 1i*omega-1i*omega_lr(11) )    *  (-1)* va1(11)*va1(11) ...
%        +1/V./( 1i*omega-1i*omega_lr(12) )    *  (-1)* va1(12)*va1(12) ;
% 
% 
% rho2 =  trho*ones(length(omega),1)  ...
%        +1/V./( 1i*omega-1i*omega_lr(1) )   .*  (-1i*omega)* va1(1)*vb1(1)       +1/V./( 1i*omega-1i*omega_lr(1) )     *  va1(1)*va1(1)...
%        +1/V./( 1i*omega-1i*omega_lr(2) )   .*  (-1i*omega)* va1(2)*vb1(2)       +1/V./( 1i*omega-1i*omega_lr(2) )     *  va1(2)*va1(2)...
%        +1/V./( 1i*omega-1i*omega_lr(3) )   .*  (-1i*omega)* va1(3)*vb1(3)       +1/V./( 1i*omega-1i*omega_lr(3) )     *  va1(3)*va1(3)...
%        +1/V./( 1i*omega-1i*omega_lr(4) )   .*  (-1i*omega)* va1(4)*vb1(4)       +1/V./( 1i*omega-1i*omega_lr(4) )     *  va1(4)*va1(4) ...
%        +1/V./( 1i*omega-1i*omega_lr(5) )   .*  (-1i*omega)* va1(5)*vb1(5)       +1/V./( 1i*omega-1i*omega_lr(5) )     *  va1(5)*va1(5)...
%        +1/V./( 1i*omega-1i*omega_lr(6) )   .*  (-1i*omega)* va1(6)*vb1(6)       +1/V./( 1i*omega-1i*omega_lr(6) )     *  va1(6)*va1(6) ...
%        +1/V./( 1i*omega-1i*omega_lr(7) )   .*  (-1i*omega)* va1(7)*vb1(7)        +1/V./( 1i*omega-1i*omega_lr(7) )     *  va1(7)*va1(7) ...
%        +1/V./( 1i*omega-1i*omega_lr(8) )   .*  (-1i*omega)* va1(8)*vb1(8)       +1/V./( 1i*omega-1i*omega_lr(8) )     *  va1(8)*va1(8) ...
%        +1/V./( 1i*omega-1i*omega_lr(9) )   .*  (-1i*omega)* va1(9)*vb1(9)       +1/V./( 1i*omega-1i*omega_lr(9) )     *  va1(9)*va1(9) ...
%        +1/V./( 1i*omega-1i*omega_lr(10) )  .*  (-1i*omega)* va1(10)*vb1(10)        +1/V./( 1i*omega-1i*omega_lr(10) )     *  va1(10)*va1(10) ...
%        +1/V./( 1i*omega-1i*omega_lr(11) )  .*  (-1i*omega)* va1(11)*vb1(11)        +1/V./( 1i*omega-1i*omega_lr(11) )     *  va1(11)*va1(11)...
%        +1/V./( 1i*omega-1i*omega_lr(12) )  .*  (-1i*omega)* va1(12)*vb1(12)        +1/V./( 1i*omega-1i*omega_lr(12) )     *  va1(12)*va1(12);
% 
% 
% rho4 =  trho*ones(length(omega),1)  ...
%        +1/V./( 1i*omega-1i*omega_lr(1) )   .*(-omega.^2)  *  vb1(1)*vb1(1)                +             1/V./( 1i*omega-1i*omega_lr(1) )  .*(-1i*omega)*(va1(1)*vb1(1)+vb1(1)*va1(1))                +         1/V./( 1i*omega-1i*omega_lr(1) )   *  va1(1)*va1(1) ...
%        +1/V./( 1i*omega-1i*omega_lr(2) )   .*(-omega.^2)  *  vb1(2)*vb1(2)                +             1/V./( 1i*omega-1i*omega_lr(2) )  .*(-1i*omega)*(va1(2)*vb1(2)+vb1(2)*va1(2))                +         1/V./( 1i*omega-1i*omega_lr(2) )   *  va1(2)*va1(2) ...
%        +1/V./( 1i*omega-1i*omega_lr(3) )   .*(-omega.^2)  *  vb1(3)*vb1(3)                +             1/V./( 1i*omega-1i*omega_lr(3) )  .*(-1i*omega)*(va1(3)*vb1(3)+vb1(3)*va1(3))                +         1/V./( 1i*omega-1i*omega_lr(3) )   *  va1(3)*va1(3) ...
%        +1/V./( 1i*omega-1i*omega_lr(4) )   .*(-omega.^2)  *  vb1(4)*vb1(4)                +             1/V./( 1i*omega-1i*omega_lr(4) )  .*(-1i*omega)*(va1(4)*vb1(4)+vb1(4)*va1(4))                +         1/V./( 1i*omega-1i*omega_lr(4) )   *  va1(4)*va1(4) ...
%        +1/V./( 1i*omega-1i*omega_lr(5) )   .*(-omega.^2)  *  vb1(5)*vb1(5)                +             1/V./( 1i*omega-1i*omega_lr(5) )  .*(-1i*omega)*(va1(5)*vb1(5)+vb1(5)*va1(5))                +         1/V./( 1i*omega-1i*omega_lr(5) )   *  va1(5)*va1(5) ...
%        +1/V./( 1i*omega-1i*omega_lr(6) )   .*(-omega.^2)  *  vb1(6)*vb1(6)                +             1/V./( 1i*omega-1i*omega_lr(6) )  .*(-1i*omega)*(va1(6)*vb1(6)+vb1(6)*va1(6))                +         1/V./( 1i*omega-1i*omega_lr(6) )   *  va1(6)*va1(6) ...
%        +1/V./( 1i*omega-1i*omega_lr(7) )   .*(-omega.^2)  *  vb1(7)*vb1(7)                +             1/V./( 1i*omega-1i*omega_lr(7) )  .*(-1i*omega)*(va1(7)*vb1(7)+vb1(7)*va1(7))                +         1/V./( 1i*omega-1i*omega_lr(7) )   *  va1(7)*va1(7) ...
%        +1/V./( 1i*omega-1i*omega_lr(8) )   .*(-omega.^2)  *  vb1(8)*vb1(8)                +             1/V./( 1i*omega-1i*omega_lr(8) )  .*(-1i*omega)*(va1(8)*vb1(8)+vb1(8)*va1(8))                +         1/V./( 1i*omega-1i*omega_lr(8) )   *  va1(8)*va1(8) ...
%        +1/V./( 1i*omega-1i*omega_lr(9) )   .*(-omega.^2)  *  vb1(9)*vb1(9)                +             1/V./( 1i*omega-1i*omega_lr(9) )  .*(-1i*omega)*(va1(9)*vb1(9)+vb1(9)*va1(9))                +         1/V./( 1i*omega-1i*omega_lr(9) )   *  va1(9)*va1(9) ...
%        +1/V./( 1i*omega-1i*omega_lr(10) )   .*(-omega.^2)  *  vb1(10)*vb1(10)                +             1/V./( 1i*omega-1i*omega_lr(10) )  .*(-1i*omega)*(va1(10)*vb1(10)+vb1(10)*va1(10))                +         1/V./( 1i*omega-1i*omega_lr(10) )   *  va1(10)*va1(10) ...
%        +1/V./( 1i*omega-1i*omega_lr(11) )   .*(-omega.^2)  *  vb1(11)*vb1(11)                +             1/V./( 1i*omega-1i*omega_lr(11) )  .*(-1i*omega)*(va1(11)*vb1(11)+vb1(11)*va1(11))                +         1/V./( 1i*omega-1i*omega_lr(11) )   *  va1(11)*va1(11) ...
%        +1/V./( 1i*omega-1i*omega_lr(12) )   .*(-omega.^2)  *  vb1(12)*vb1(12)                +             1/V./( 1i*omega-1i*omega_lr(12) )  .*(-1i*omega)*(va1(12)*vb1(12)+vb1(12)*va1(12))                +         1/V./( 1i*omega-1i*omega_lr(12) )   *  va1(12)*va1(12) ;

