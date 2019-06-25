%% Setup 
clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked') 

% Read in
density=dlmread('density.txt');

% Split up roads
road01=density(:,1:25);
road02=density(:,26:50);
road03=density(:,51:75);
road04=density(:,76:100);
road05=density(:,101:125);
road06=density(:,126:150);
road07=density(:,151:175);
road08=density(:,176:200);
road09=density(:,201:225);
road10=density(:,226:250);
road11=density(:,251:275);
road12=density(:,276:300);
road13=density(:,301:325);
road14=density(:,326:350);
road15=density(:,351:375);
road16=density(:,376:400);

clear density

% Plot density profile

close all
clc

 %%
 
 siz=50;
 
for i=1:338
    subplot(12,6,[13 19])
        plot(1:25,road01(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 1')
    subplot(12,6,[25 31])
        plot(1:25,road02(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 2')
    subplot(12,6,[3 9])
        plot(1:25,road03(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 3')
    subplot(12,6,[14 20])
        plot(1:25,road04(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 4')
    subplot(12,6,[26 32])
        plot(1:25,road05(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 5')    
    subplot(12,6,[16 22])
        plot(1:25,road06(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 6') 
    subplot(12,6,[17 23])
        plot(1:25,road07(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 7') 
    subplot(12,6,[34 40])
        plot(1:25, road08(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 8')
    subplot(12,6,[58 64])
        plot(1:25, road09(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 9')
    subplot(12,6,[29 35])
        plot(1:25, road10(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 10')
    subplot(12,6,[41 47])
        plot(1:25, road11(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 11')
    subplot(12,6,[53 59])
        plot(1:25, road12(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 12')
    subplot(12,6,[65 71])
        plot(1:25, road13(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 13')
    subplot(12,6,[24 30])
        plot(1:25, road14(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 14')
    subplot(12,6,[48 54])
        plot(1:25, road15(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 15')
    subplot(12,6,[66 72])
        plot(1:25, road16(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 16')
    
    set(gcf,'PaperSize',[siz*0.8,siz*1.2]);
    filename=sprintf('frame%04i',i);
    print(filename,'-dpng')
    %pause(0.001)
end

%% Plot Quasi microscopic emulation
%  Overlayed density profile

figure
for i=1:338
    spacing=cumsum(1./density(i,:));
    spacing=spacing+(1.25-spacing(25))*ones(1,length(spacing));
    plot(spacing,25*ones(1,length(spacing)),'ok')
    hold on
    plot(spacing,density(i,:))
    hold off
    xlim([-0.7 1.3])
    ylim([10 105])
    pause(0.0001)
end

clear i spacing 


%% Save figure as PDF
% figure(2);



