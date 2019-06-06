%% Setup 
clear all
close all
clc

% Read in
density=dlmread('density.txt');
struct=load('ref_density.mat');
ref_density=struct.ref_density;
clear struct

% Split up reference roads
ref_road1=ref_density(:,1:25);
ref_road2=ref_density(:,26:50);
ref_road3=ref_density(:,51:75);

% Split up roads
road1=density(:,1:25);
road2=density(:,26:50);
road3=density(:,51:75);

if sum([sum(road1~=ref_road1),sum(road2~=ref_road2),sum(road3~=ref_road3)])==0
    'Correct profiles'
else
    'Incorrect profiles'
end

clear ref_density density ans

%% Plot density profile

close all
clc

dx=0.2;

for i=1:338
    subplot(7,2,[1,3])
        plot(1,2,'ro');
        hold on
        plot(3,4,'bx')
        text(1.5,2.5,sprintf('Frame %i',i))
        legend('Updated', 'Original Reference')
        box off
        axis off
        hold off
    subplot(7,2,[5,7,9,11])
        plot(1:25,road1(i,:),'ro')
        hold on
        plot(1:25,ref_road1(i,:),'bx')
        hold off
        ylim([0 113])
        title('Road 1')
    subplot(7,2,[2,4,6])
        plot(1:25,road2(i,:),'ro')
        hold on
        plot(1:25,ref_road2(i,:),'bx')
        hold off
        ylim([0 113])
        title('Road 2')
    subplot(7,2,[10,12,14])
        plot(1:25,road3(i,:),'ro')
        hold on
        plot(1:25,ref_road3(i,:),'bx')
        hold off
        ylim([0 113])
        title('Road 3')
    pause(0.000000001)
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

