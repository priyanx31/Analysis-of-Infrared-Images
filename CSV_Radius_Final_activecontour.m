clear all
close all
clc
tic

%% Directory

dir_main = 'C:\Users\Lenovo PC\Downloads\Example1\Example1';
cd(dir_main)
dir_data = dir_main;

load rainbowpalette.mat
load modifiedJet.mat

%%

water = 0;
ethanol = abs(water-100);
temperature = 50;
concentration = 0.3;
casenum = 1;

Tmin = 46;
Tmax = 57;

sf = 1;  %%%% scaling factor pixel to mm

nstart = 551;
nend = 660;
nlife = 840;

fnum = sprintf('%d-%d_%d_%0.1f_%d_%d.csv', water, ethanol, temperature, concentration, casenum, nlife);
dir_path=fullfile(dir_data, fnum);
BG = importdata(dir_path);

 figure(1), imagesc(BG), colormap(palette), caxis([Tmin Tmax]), colorbar, axis image

[~,ncnt] = size(nstart:5:nend);
RZ = zeros(ncnt,4);
nimg = 1;

for ii = nstart+10
    
    fnum = sprintf('%d-%d_%d_%0.1f_%d_%d.csv', water, ethanol, temperature, concentration, casenum, ii);
    dir_path=fullfile(dir_data, fnum);
    M = importdata(dir_path);
    
    figure(2), imagesc(M), colormap(palette), caxis([Tmin Tmax]), colorbar, axis image
    
    I = (M - BG)+temperature;
    
    I = medfilt2(I, [7 7]);
  %  I = medfilt2(I, [4 4]);
    I = imsharpen(I);
    J = imresize(I, 0.25);
    
    [ncols, nrows] = size(J);
    [X,Y] = meshgrid(1:nrows,1:ncols);
    N = imgaussfilt(J, 5.0);
   
   % figure(3), surf(X,Y,N), colormap(palette), colorbar
    
    figure(3), imagesc(J), colormap(palette), colorbar, axis image
    
    figure(4), imagesc(N), colormap(palette), caxis([Tmin Tmax]), colorbar, axis image
    
    BW = edge(N, 'Prewitt');   %%%%% Prewitt Canny
    I1 = imclearborder(BW);
    
    figure(5), imagesc(I1), colormap(gray), colorbar, axis image
    
    dsize = 10;
    se = strel('disk',dsize,0);
    I2 = imclose(I1, se);
    I3 = imfill(I2, 'holes');
    
     figure(6), imagesc(I3), colormap(gray), colorbar, axis image
        
        I3_new=imresize(I3,4);
        mask = zeros(size(I3_new));
        mask(25:end-25,25:end-25) = 1;
        bw = activecontour(I3_new,mask,'edge');
  
    
figure(7), imagesc(bw), colormap(gray), colorbar, axis image
    
    
    
    [B, ~] = bwboundaries(bw);
    [ns1, ~] = size(B);
    nval = zeros(ns1,1);
    for jj=1:ns1
        [ns2, ns3] =size(B{jj});
        nval(jj) = ns2;
    end
    [~, ptm] = max(nval);
    Bxy = B{ptm};
    
        figure(8);
        hold on
        imagesc(M), colormap(palette), caxis([Tmin Tmax]), colorbar, axis image
        plot(Bxy(:,2), Bxy(:,1), '-r');
        hold off
        xlim([1 640]), ylim([1 512]), ax = gca; ax.YDir = 'reverse';
        set(gca, 'box', 'on');
    
    clear x
    
    x(:,1) = Bxy(:,2);
    x(:,2) = Bxy(:,1);
    
    [z, r] = fitcircle(x);
    
    RZ(nimg, 1) = ii;
    RZ(nimg, 2) = r;
    RZ(nimg, 3) = z(1);
    RZ(nimg, 4) = z(2);
    
    t = linspace(0, 2*pi, 100);
    
    X = (Bxy(:,2)-z(1))/sf;
    Y = (Bxy(:,1)-z(2))/sf;
    
    Xf = ((z(1)  + r  * cos(t)) - z(1))/sf;
    Yf = ((z(2)  + r * sin(t)) - z(2))/sf;
    
    hfig=figure(9);
    hold on
    plot(X, Y, 'ob','markersize',4.0,'markerfacecolor','b');
    plot(Xf, Yf, 'r-','linewidth',3);
    set(gca,'FontSize',28);
    hold off
    xlim([-150 150]);
    ylim([-150 150]);
    axis square
    
    fnum = sprintf('C_%d-%d_%d_%0.1f_%d_%d.png', water, ethanol, temperature, concentration, casenum, ii);
    dir_path=fullfile(dir_data, fnum);
    print(hfig,'-dpng',dir_path);
    
    fnum = sprintf('B_%d-%d_%d_%0.1f_%d_%d.mat', water, ethanol, temperature, concentration, casenum, ii);
    dir_path=fullfile(dir_data, fnum);
    save(dir_path, 'Bxy');
    
    %close all
    
    nimg = nimg+1;
end

fnum = sprintf('RZ_%d-%d_%d_%0.1f_%d.mat', water, ethanol, temperature, concentration, casenum);
dir_path=fullfile(dir_data, fnum);
save(dir_path, 'RZ');

%%
toc