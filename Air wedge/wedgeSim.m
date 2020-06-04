% The numerical simulation for the air wedge: 
% We use the parameter extracted from the experimenats (Cf. wedgePreProc.m)
% to conduct a simulation for comparison
% Authors: Yuye Ling and Mengyuan Wang

clear all;
% close all;
load('wedgePreProc.mat');
linData = h5read('rawSpectrumAirWedgeThorlab.h5','/rawData');

% numSpec defines the sampling number of the spectral interferogram
% numObj defines the grid size of the original function
numSpec = 2048;
numObj = 70000;
factor = 8;

% The starting and ending wavelength are obtained directly from Thorlabs
% support department
lambdaSt = 791.6e-9;
lambdaEnd = 994e-9;
kSt = 2 * pi / lambdaEnd;
kEnd = 2 * pi / lambdaSt;
k = linspace(kEnd, kSt, numSpec)';
% Sk.mat contains the emission spectrum of the light source after the
% lambda-to-k interpolation
load('Sk.mat');


% We use the data (sep) up to 2449th column to fit the angle between the
% top and bottom interface
% The lateral pixel size is 2 um (3000 pixels for 6 mm FOV)
% The axial pixel size is 1.9438 um
dx = 2e-6;

% dz_fft gives the axial pixel size (digital resolution) after conventional
% IDFT processing (w/o zero padding)
dzFFT = 0.5 * 1 / (1 / lambdaSt - 1 / lambdaEnd);
% dzObj gives the discretization grid size in the object domain
% This value was set to be extremely small so that the small tilting angle
% of the top layer could be smoothly discretized
dzObj = 0.1e-7;
zObj = linspace(0, (numObj - 1) * dzObj, numObj)';
% Construct transform matrix
[zGridObj, kGrid] = meshgrid(zObj, k);
matFourObj = real(exp(2j * zGridObj .* kGrid));
specObj = repmat(Sk, 1 , numObj);

% Obtain fitted separation, top layer location and bottom layer location
lateralGrid = [2351: 1: 2500];
% separation = (thickFit.p1 * lateralGrid + thickFit.p2) * dzFFT;
topTrue = topFit.p1 * lateralGrid + topFit.p2;
botTrue = botFit.p1 * lateralGrid + botFit.p2;
wedgeTrue = zeros(numObj, length(topTrue));


for iCol = 1: length(topTrue)
    wedgeTrue(round((topTrue(iCol) - 1) * dzFFT / 8 / dzObj + 1), iCol) = 1;
    wedgeTrue(round((botTrue(iCol) - 1) * dzFFT / 8 / dzObj + 1), iCol) = mean(peakBot(2351: 2450)./peakTop(2351: 2450));
    
%     wedgeTrue(round(322 * dzFFT / dz0), iCol) = 1;
%     wedgeTrue(round(322 * dzFFT / dz0) + round(separation(iCol) / dz0), iCol) = 0.7;

    matTranObj = specObj .* matFourObj;

    fringe(:, iCol) = matTranObj * wedgeTrue(:, iCol);

    % The simulation is zero padded by "factor"
    simImg (:, iCol) = abs(ifft(fringe(:, iCol), factor * numSpec));
end 
img = abs(ifft(linData(:, 2351: 2500), factor * numSpec));
figure('Name','The experimental measured image (IFFT reconstruction w/ 8x padding)');
set(gcf,'outerposition',get(0,'screensize'));
imagesc(lateralGrid * dx * 1e6, [2551: 2650] * dzFFT * 1e6, 10 .* log10(img(2551: 2650, 91: end))); 
colormap hot; colorbar; caxis([15, 21])
text(4700,5032,'1', 'Color','white','FontSize',20)
text(4725,5032,'2', 'Color','white','FontSize',20)
text(4745,5032,'3', 'Color','white','FontSize',20)
text(4765,5032,'4', 'Color','white','FontSize',20)
text(4790,5032,'5', 'Color','blue','FontSize',20)
text(4810,5032,'6', 'Color','blue','FontSize',20)
text(4830,5032,'7', 'Color','blue','FontSize',20)
text(4853,5032,'8', 'Color','blue','FontSize',20)
text(4875,5032,'9', 'Color','blue','FontSize',20)
text(4896,5032,'10', 'Color','blue','FontSize',20)
text(4920,5032,'11', 'Color','blue','FontSize',20)

figure('Name','The simulated image (IFFT reconstruction w/ 8x padding)');
set(gcf,'outerposition',get(0,'screensize'));
imagesc(lateralGrid * dx * 1e6, [2551: 2650] * dzFFT * 1e6, 10 .* log10(simImg(2551: 2650, 91: end))); 
colormap hot; colorbar; caxis([26.5, 33])
text(4700,5032,'1', 'Color','white','FontSize',20)
text(4725,5032,'2', 'Color','white','FontSize',20)
text(4745,5032,'3', 'Color','white','FontSize',20)
text(4765,5032,'4', 'Color','white','FontSize',20)
text(4790,5032,'5', 'Color','blue','FontSize',20)
text(4810,5032,'6', 'Color','blue','FontSize',20)
text(4830,5032,'7', 'Color','blue','FontSize',20)
text(4853,5032,'8', 'Color','blue','FontSize',20)
text(4875,5032,'9', 'Color','blue','FontSize',20)
text(4896,5032,'10', 'Color','blue','FontSize',20)
text(4920,5032,'11', 'Color','blue','FontSize',20)

save('wedgeSim.mat', 'wedgeTrue', 'simImg', 'fringe', 'topTrue', 'botTrue', 'wedgeTrue', 'dx');