% Renan Liupekevicius Carnielli TU/e
% 18/10/2024



%% eta =10
f10h    = FRF_eta10_HOM(:,1);
f10d    = FRF_eta10_DNS(:,1);
% homogenized
FRFHOM10 = FRF_eta10_HOM(:,3);
% direct numerical simulation
FRFDNS10 = FRF_eta10_DNS(:,3);


figure
semilogx(f10d,FRFDNS10, 'r.', f10h,FRFHOM10 , 'r');
% ylim([0 1])
grid on


%% eta =0
f00h    = FRF_eta00_HOM(:,1);
f00d    = FRF_eta00_DNS(:,1);
% homogenized
FRFHOM00 = FRF_eta00_HOM(:,3);
% direct numerical simulation
FRFDNS00 = FRF_eta00_DNS(:,3);


figure
semilogx(f00d,FRFDNS00, 'b.', f00h,FRFHOM00 , 'b');
% ylim([0 1])
grid on


%% compare

figure
semilogx(f00d,FRFDNS00, 'b.', f00h,FRFHOM00 , 'b',f10d,FRFDNS10, 'r.', f10h,FRFHOM10 , 'r');
ylim([-14 -5])
grid on
box on
fontsize(gca, 21,'pixels')