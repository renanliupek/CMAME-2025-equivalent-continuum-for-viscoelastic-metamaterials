% plot transmission curves
% 11/10/2024


%% eta =10
f10h    = fRT_eta10_HOM(:,1);
f10d    = fRT_eta10_DNS(:,1);
% homogenized
RHOM10 = fRT_eta10_HOM(:,3);
THOM10 = fRT_eta10_HOM(:,4);
AHOM10  = 1-abs(RHOM10).^2-abs(THOM10).^2;
TLHOM10 = 20*log10(1./abs(THOM10));
% direct numerical simulation
RDNS10=fRT_eta10_DNS(:,3);
TDNS10=fRT_eta10_DNS(:,4);
ADNS10  = 1-abs(RDNS10).^2-abs(TDNS10).^2;
TLDNS10 = 20*log10(1./abs(TDNS10));


figure
semilogx(f10d,abs(RDNS10), 'r.', f10h,abs(RHOM10), 'r',f10d,abs(TDNS10), 'b.', f10h,abs(THOM10), 'b',f10d,(ADNS10), 'g.', f10h,(AHOM10), 'g')
% semilogx(f10,abs(ADNS10), 'r.', f10,abs(AHOM10), 'r')
ylim([0 1])
grid on

figure
semilogx(f10d, TLDNS10, 'b.', f10h, TLHOM10, 'b')
grid on
 % xlim([0 1450])

% figure
% hold on
% semilogx(f10,angle(RDNS10), 'r.')
% semilogx(f10,angle(RHOM10), 'r')

%% eta = 0
f00h    = fRT_eta00_HOM(:,1);
f00d    = fRT_eta00_DNS(:,1);
% homogenized
RHOM00 = fRT_eta00_HOM(:,3);
THOM00 = fRT_eta00_HOM(:,4);
AHOM00  = 1-abs(RHOM00).^2-abs(THOM00).^2;
TLHOM00 = 20*log10(1./abs(THOM00));
% direct numerical simulation
RDNS00=fRT_eta00_DNS(:,3);
TDNS00=fRT_eta00_DNS(:,4);
ADNS00  = 1-abs(RDNS00).^2-abs(TDNS00).^2;
TLDNS00 = 20*log10(1./abs(TDNS00));


figure
semilogx(f00d,abs(RDNS00), 'r.', f00h,abs(RHOM00), 'r',f00d,abs(TDNS00), 'b.', f00h,abs(THOM00), 'b',f00d,(ADNS00), 'g.', f00h,(AHOM00), 'g')
% semilogx(f00,abs(ADNS00), 'r.', f00,abs(AHOM00), 'r')
ylim([0 1])
grid on

figure
semilogx(f00d, TLDNS00, 'b.', f00h, TLHOM00, 'b')
grid on

% figure
% hold on
% semilogx(f,angle(RDNS), 'r.')
% semilogx(f,angle(RHOM), 'r')


%% compare


figure
semilogx(f00d, TLDNS00, 'b.', ...
         f00h, TLHOM00, 'b', ...
         f10d, TLDNS10, 'r.', ...
         f10h, TLHOM10, 'r');
% slg(1).LineWidth = 1.5;

ylim([0 50])
% xlim([100 1450])
xlim([100 4000])
grid on
box on
fontsize(gca, 21,'pixels')


figure
% slg=semilogx(f00d(1:3:end), ADNS00(1:3:end), 'b.', ...
%              f00h, AHOM00, 'b', ...
%              f10d(1:3:end), ADNS10(1:3:end), 'r.', ...
%              f10h, AHOM10, 'r');
slg=semilogx(f00d, ADNS00, 'b.', ...
             f00h, AHOM00, 'b', ...
             f10d, ADNS10, 'r.', ...
             f10h, AHOM10, 'r');
slg(1).LineWidth = 4;
% slg(1).Marker = 'o';
% slg(1).MarkerEdgeColor = 'black';
% slg(1).MarkerFaceColor = 'blue';
ylim([0 1])
% xlim([100 1450])
xlim([100 4000])
grid on
box on
fontsize(gca, 21,'pixels')
