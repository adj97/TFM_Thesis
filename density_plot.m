%% Setup 
clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked') 

% Read in
density=dlmread('density.txt');

% Split up roads
road1=density(:,1:25);
road2=density(:,26:50);
road3=density(:,51:75);
road4=density(:,76:100);
road5=density(:,101:125);
road6=density(:,126:150);

clear density

%% Plot density profile

close all
clc

dx=0.2;
siz=18;

for i=1:338
    subplot(2,4,[1 5])
        plot(1:25,road1(i,:),'o')
        ylim([0 113])
        title('Road 1')
        text(5,100,sprintf('Frame %i',i))
    subplot(2,4,2)
        plot(1:25,road2(i,:),'o')
        ylim([0 113])
        title('Road 2')
    subplot(2,4,6)
        plot(1:25,road3(i,:),'o')
        ylim([0 113])
        title('Road 3')
    subplot(2,4,3)
        plot(1:25,road4(i,:),'o')
        ylim([0 113])
        title('Road 4')
    subplot(2,4,7)
        plot(1:25,road5(i,:),'o')
        ylim([0 113])
        title('Road 5')
    subplot(2,4,[4 8])
        plot(1:25,road6(i,:),'o')
        ylim([0 113])
        title('Road 6')
    %set(gcf,'PaperSize',[siz,siz]);
    %filename=sprintf('frame%04i',i);
    %print(filename,'-dpng')
    pause()
end

clear dx ans dimension

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



