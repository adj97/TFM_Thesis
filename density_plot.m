%% Setup 
clear density
close all
clc
set(0,'DefaultFigureWindowStyle','docked') 

% Read in
density=dlmread('density.txt');

dx=0.05;

% split
r1_fine = density(:,1:1*5/dx);
r2_fine = density(:, 1*5/dx+1:2*5/dx);
r3_fine = density(:,2*5/dx+1:3*5/dx);
r4_fine = density(:, 3*5/dx+1:4*5/dx);
r5_fine = density(:,4*5/dx+1:5*5/dx);
r6_fine = density(:, 5*5/dx+1:6*5/dx);
r7_fine = density(:,6*5/dx+1:7*5/dx);
r8_fine = density(:, 7*5/dx+1:8*5/dx);

maximum = max(max(max(density)), maximum)*1.1;
clear density dx


%% plot

close all
figure

% Simple Plot density profile
[nframes, clength]=size(r1_coarse);
[nframes, mlength]=size(r1_med);
[nframes, flength]=size(r1_fine);

length=1;

fsize = 12;
siz = 20;

for i=12
    subplot(8,8,[25 26 33 34])
    plot((0:clength-1)/(clength-1),r1_coarse(i,:))
    hold on
    plot((0:mlength-1)/(mlength-1),r1_med(2*i-1,:))
    plot((0:flength-1)/(flength-1),r1_fine(4*i-2,:))
    hold off
    title('q_1 = 0.5')
    xlim([0 length])
    ylim([0 maximum])
    pbaspect([1 1 1])
    
    subplot(8,8,[31 32 39 40])
    plot((0:clength-1)/(clength-1),r2_coarse(i,:))
    hold on
    plot((0:mlength-1)/(mlength-1),r2_med(2*i-1,:))
    plot((0:flength-1)/(flength-1),r2_fine(4*i-2,:))
    hold off
    xlim([0 length])
    ylim([0 maximum])
    pbaspect([1 1 1])
    
    subplot(8,8,[52 53 60 61])
    plot((0:clength-1)/(clength-1),r3_coarse(i,:))
    hold on
    plot((0:mlength-1)/(mlength-1),r3_med(2*i-1,:))
    plot((0:flength-1)/(flength-1),r3_fine(4*i-2,:))
    hold off
    xlim([0 length])
    ylim([0 maximum])
    pbaspect([1 1 1])
    
    subplot(8,8,[4 5 12 13])
    plot((0:clength-1)/(clength-1),r4_coarse(i,:))
    hold on
    plot((0:mlength-1)/(mlength-1),r4_med(2*i-1,:))
    plot((0:flength-1)/(flength-1),r4_fine(4*i-2,:))
    hold off
    xlim([0 length])
    ylim([0 maximum])
    pbaspect([1 1 1])
    
    subplot(8,8,[35 36 43 44])
    plot((0:clength-1)/(clength-1),r5_coarse(i,:))
    hold on
    plot((0:mlength-1)/(mlength-1),r5_med(2*i-1,:))
    plot((0:flength-1)/(flength-1),r5_fine(4*i-2,:))
    hold off
    xlim([0 length])
    ylim([0 maximum])
    pbaspect([1 1 1])
 
    subplot(8,8,[37 38 45 46])
    plot((0:clength-1)/(clength-1),r6_coarse(i,:))
    hold on
    plot((0:mlength-1)/(mlength-1),r6_med(2*i-1,:))
    plot((0:flength-1)/(flength-1),r6_fine(4*i-2,:))
    hold off
    xlim([0 length])
    ylim([0 maximum])
    pbaspect([1 1 1])
    
    subplot(8,8,[21 22 29 30])
    plot((0:clength-1)/(clength-1),r7_coarse(i,:))
    hold on
    plot((0:mlength-1)/(mlength-1),r7_med(2*i-1,:))
    plot((0:flength-1)/(flength-1),r7_fine(4*i-2,:))
    hold off
    xlim([0 length])
    ylim([0 maximum])
    pbaspect([1 1 1])
 
    subplot(8,8,[19 20 27 28])
    plot((0:clength-1)/(clength-1),r8_coarse(i,:))
    hold on
    plot((0:mlength-1)/(mlength-1),r8_med(2*i-1,:))
    plot((0:flength-1)/(flength-1),r8_fine(4*i-2,:))
    hold off
    xlim([0 length])
    ylim([0 maximum])
    pbaspect([1 1 1])
    
    subplot(8,8,[49 50 57 58])
    text(0, 0, sprintf('Time t=0.024752'))%, (i-1)*0.5/(nframes-1)))
    text(0, 0.3, sprintf('Frame %d of 405', i))
    axis off;
    
    subplot(8,8,[6 7 8 14 15 16])
    plot([0 1],[0 1])
    hold on
    plot([1 0],[0 1])
    plot([1 0],[0 1])
    hold off
    xlim([-1 0])
    ylim([-1 0])
    legend('Coarse', 'Medium', 'Fine', ...
           'location', 'southwest')
    axis off;
    %pause()
   %
    
    set(gcf,'PaperSize',[siz*1.2,siz]);
    %filename=sprintf('frame%04i',i);
    filename = 'trafficCircles3';
    print(filename,'-dpdf')
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



