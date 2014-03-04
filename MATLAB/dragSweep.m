%-------- dragSweep --------%
%                           %
%       Todd Reichert       %
%        Feb 1, 2011        %
%                           %
%---------------------------%

% Performs a sweep over Reynolds number, thickness and fraction of laminar
% flow, computing the drag of an airfoil. Used for validation and as a
% reference.
clear all
close all

NRe = 20;
Re = linspace(50000,500000,NRe);
xtcU = [0.15 0.5 0.6];
xtcL = [0.15 1 1];
tc = 0.15;
for i = 1:NRe
    for j = 1:length(xtcU)
        Cd(i,j) = dragCoefficient(Re(i),tc,xtcU(j),xtcL(j));
    end
end
for i = 1:NRe
    for j = 1:length(xtcU)
        Cd2(i,j) = dragCoefficientFit(Re(i),tc,xtcU(j),xtcL(j));
    end
end

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)*0.45 scrsz(4)*0.5 scrsz(3)*0.5 scrsz(4)*0.4])
plot(Re,Cd,Re,Cd2,'--')
xlabel('Re')
ylabel('C_d')
legend('x_t/c = 0.15','x_t/c = 0.5','x_t/c = 0.6')
hold on
Re = [100000 200000 300000 500000];
DAE31 = [0.0281 0.0153 0.0117 0.0088];
BE10056 = [0.0212 0.0135 0.0109 0.0082];
BE10759 = [0.0248 0.0148 0.0114 0.009];
DAE31trip = [0.0281 0.0212 0.0190 0.0171];
BE10056trip = [0.0263 0.0213 0.0189 0.0169];
BE10759trip = [0.0263 0.0212 0.0201 0.0174];
plot(Re,DAE31,'or',Re,BE10056,'sr',Re,BE10759,'dr')
plot(Re,DAE31trip,'ob',Re,BE10056trip,'sb',Re,BE10759trip,'db')
legend('x_t/c = 0.15','x_t/c = 0.5','x_t/c = 0.6','DAE31','BE10056','BE10759')
axis([0 500000 0.005 0.03])


% Ntc = 20;
% Nxtc = 7;
% tc = linspace(0.06,0.20,Ntc);
% xtc = linspace(0,0.6,Nxtc);
% Re = 250000;
% for i = 1:Ntc
%     for j = 1:Nxtc
%         Cd(i,j) = dragCoefficient(Re,tc(i),xtc(j),xtc(j));
%     end
% end
% subplot(1,2,2);
% plot(tc,Cd)
% xlabel('t/c')
% ylabel('C_d')
% legend('x_t/c = 0','x_t/c = 0.1','x_t/c = 0.2','x_t/c = 0.3','x_t/c = 0.4','x_t/c = 0.5','x_t/c = 0.6','location','SouthEast')
% hold on
% xx = [0.06 0.09 0.12 0.15];
% yy = [0.01289 0.01383 0.01493 0.01631];
% plot(xx,yy,'^')
% xx = [0.06 0.09 0.12 0.15];
% yy = [0.00838 0.00888 0.00947 0.0102];
% plot(xx,yy,'k^')

set(gcf,'color',[1 1 1])