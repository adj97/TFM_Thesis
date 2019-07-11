%% Setup 
clear density frames i length
close all
clc
set(0,'DefaultFigureWindowStyle','docked') 

% Read in
density=dlmread('density.txt');

% split
wo5_road1 = density(:,1:50);
wo5_road2 = density(:, 51:100);
wo5_road3 = density(:,101:150);
wo5_road4 = density(:, 151:200);
wo5_road5 = density(:,201:250);
wo5_road6 = density(:, 251:300);
clear density

%%

% Simple Plot density profile
[frames, length]=size(fo_road1);
maximum = max(max(fo_road1));

siz = 30;

for i=1:400
    subplot(2,4,[1 5])
    plot(0:length-1,fo_road1(i,:),'b-')
    hold on    
    plot(0:length-1,so_road1(i,:),'g-')
    plot(0:length-1,M_road1(i,:),'r-')
    plot(0:length-1,wo5_road1(i,:),'k-')
    hold off
    legend('1st', '2nd','MUSCL3','WENO5')
    title(sprintf('Frame %d', i))
    xlim([0 length])
    ylim([0 maximum])
    
    subplot(2,4,2)
    plot(0:length-1,fo_road2(i,:),'b-')
    hold on
    plot(0:length-1,so_road2(i,:),'g-')
    plot(0:length-1,M_road2(i,:),'r-')
    plot(0:length-1,wo5_road2(i,:),'k-')
    hold off
    xlim([0 length])
    ylim([0 maximum])
    
    subplot(2,4,6)
    plot(0:length-1,fo_road3(i,:),'b-')
    hold on
    plot(0:length-1,so_road3(i,:),'g-')
    plot(0:length-1,M_road3(i,:),'r-')
    plot(0:length-1,wo5_road3(i,:),'k-')
    hold off
    xlim([0 length])
    ylim([0 maximum])
    
    subplot(2,4,3)
    plot(0:length-1,fo_road4(i,:),'b-')
    hold on
    plot(0:length-1,so_road4(i,:),'g-')
    plot(0:length-1,M_road4(i,:),'r-')
    plot(0:length-1,wo5_road4(i,:),'k-')
    hold off
    xlim([0 length])
    ylim([0 maximum])
    
    subplot(2,4,7)
    plot(0:length-1,fo_road5(i,:),'b-')
    hold on
    plot(0:length-1,so_road5(i,:),'g-')
    plot(0:length-1,M_road5(i,:),'r-')
    plot(0:length-1,wo5_road5(i,:),'k-')
    hold off
    xlim([0 length])
    ylim([0 maximum])
    
    subplot(2,4,[4 8])
    plot(0:length-1,fo_road6(i,:),'b-')
    hold on
    plot(0:length-1,so_road6(i,:),'g-')
    plot(0:length-1,M_road6(i,:),'r-')
    plot(0:length-1,wo5_road6(i,:),'k-')
    hold off
    xlim([0 length])
    ylim([0 maximum])
    
    set(gcf,'PaperSize',[siz*0.8,siz*1.2]);
    filename=sprintf('frame%04i',i);
    print(filename,'-dpng')
end

 %% network plot
 
 %siz=50;
 
for i=1:338
    subplot(5, 5, 2)
        plot(1:25,road01(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 1')
        camroll(-90)
    subplot(5, 5, 4)
        plot(1:25,road02(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 2')
        camroll(90)
    subplot(5, 5, 6)
        plot(1:25,road03(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 3')
    subplot(5, 5, 8)
        plot(1:25,road04(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 4')
    subplot(5, 5, 10)
        plot(1:25,road05(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 5')    
    subplot(5, 5, 12)
        plot(1:25,road06(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 6')
        camroll(-90)
    subplot(5, 5, 14)
        plot(1:25,road07(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 7') 
        camroll(90)
    subplot(5, 5, 16)
        plot(1:25, fliplr(road08(i,:)),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 8')
    subplot(5, 5, 18)
        plot(1:25, fliplr(road09(i,:)),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 9')
    subplot(5, 5, 20)
        plot(1:25, fliplr(road10(i,:)),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 10')
    subplot(5, 5, 22)
        plot(1:25, road11(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 11')
        camroll(-90)
    subplot(5, 5, 24)
        plot(1:25, road12(i,:),'o')
        ylim([0 15])
        xlim([0 25])
        title('Road 12')
        camroll(90)
    

    pause()
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



