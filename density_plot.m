%% Setup 
clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked') 

% Read in
density=dlmread('density.txt');

% Simple Plot density profile
[frames, length]=size(density);
maximum=max(max(density));

for i=1:frames
    plot(1:length,density(i,:),'o')
    xlim([0 length])
    ylim([0 maximum*1.1])
    pause()
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
    
    %set(gcf,'PaperSize',[siz*0.8,siz*1.2]);
    %filename=sprintf('frame%04i',i);
    %print(filename,'-dpng')
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



