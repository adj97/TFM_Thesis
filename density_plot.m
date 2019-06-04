%% Setup 
clear all
close all
clc

% Read in
density=dlmread('density.txt');
%struct=load('rho.mat');
%reference_density=struct.rho;
%clear struct

% Split up roads
road1=density(:,1:25);
road2=density(:,26:50);
road3=density(:,51:75);

clear density

%% Plot density profile

close all

dx=0.2;

for i=1:338
    subplot(3,2,1)
        plot(1,1);
        text(1,1,sprintf('Frame %i',i))
        box off
        axis off
    subplot(3,2,3)
        plot(1:25,road1(i,:),'o')
        ylim([0 113])
        title('Road 1')
    subplot(3,2,2)
         plot(1:25,road2(i,:),'o')
         ylim([0 113])
         title('Road 2')
    subplot(3,2,6)
         plot(1:25,road3(i,:),'o')
         ylim([0 113])
         title('Road 3')
    pause()
end

clear dx i ans dimension

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

