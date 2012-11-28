%%%%%%%%%%%%%%%%%%%
%%%   SYSTEMS   %%%
%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%% systems definition
G1 = tf(1,[1 2 1]);
G2 = tf(1,[1 1 1]);
G3 = tf(1,[1 0.5 1]);

figure(1);
    step(G1);
    hold on;
    step(G2);
    step(G3);
    grid on;
legend('G_1','G_2','G_3',4);
FigChar(figure(1),'t [s]','y [-]',12,'',14,3);
print(1, '-depsc2', '../Figures/Steps');