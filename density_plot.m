%% Setup 
clear all
close all
clc

% Read in
density=dlmread('density.txt');
struct=load('rho.mat');
reference_density=struct.rho;
clear struct

%% Plot density profile

close all

dx=0.2;
dimension=size(density);

for i=1:length(density)
    subplot(2,1,1)
        plot(1:dimension(2),density(i,:),'o')
        ylim([0 113])
        title('My Density')
    subplot(2,1,2)
        plot(1:dimension(2),reference_density(i,:),'o')
        ylim([0 113])
        title('Reference Density')
    pause(0.1)
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

