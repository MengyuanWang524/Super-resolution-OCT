% Recontructing the air wedge: 
% We recontruct the air wedge from both the simulated and experimental
% data by using the proposed algorithm. For comparison, we also used
% Lucy-Richardson deconvolution algorithm which is available in MATLAB.
% Authors: Yuye Ling and Mengyuan Wang

load('wedgePreProc.mat');
load('wedgeSim.mat');
% numSpec defines the sampling number of the spectral interferogram
% numRecn defines the grid size of the original function
numSpec = 2048;
numRecn = 70000;

% The starting and ending wavelength are obtained directly from Thorlabs
% support department
lambdaSt = 791.6e-9;
lambdaEnd = 994e-9;
kSt = 2 * pi / lambdaEnd;
kEnd = 2 * pi / lambdaSt;
k = linspace(kEnd, kSt, numSpec)';
load('Sk.mat');
Sk = Sk ./ (sum(Sk)) * numRecn;
fringe = linData(:, 2351: 2500, 1);

% The lateral pixel size is 2 um (3000 pixels for 6 mm FOV)
% The axial pixel size is 1.9438 um
dx = 2e-6;

% dz_fft gives the axial pixel size (digital resolution) after conventional
% IDFT processing (w/o zero padding)
dzFFT = 0.5 * 1 / (1 / lambdaSt - 1 / lambdaEnd);
% T defines the grid size of the reconstructed function
numRecn = 2048 * factor / 2;

lambda = 1000;
% dz = 0.1e-6;
dzRecn = dzFFT / factor;
zRecn = linspace(0, (numRecn - 1) * dzRecn, numRecn)';
%z = linspace(0, 4095 * dz, T);

% We first use Lucy-Richardson algorithm to deconvolve the zero padded
% image
PSF = abs(ifftshift(ifft(Sk, factor * numSpec)));
lateralGrid = [2351: 1: 2500];
img = img(:, lateralGrid);
for iCol = 1: length(topTrue)
    simDeconvImg(:, iCol) = deconvlucy(simImg(:, iCol), PSF);
    deconvImg(:, iCol) = deconvlucy(img(:, iCol), PSF);
end

figure('Name','The deconvolved image (Lucy-Richardson algorithm w/ 8x padding)');
set(gcf,'outerposition',get(0,'screensize'));
imagesc(lateralGrid * dx * 1e6, [2551: 2650] * dzFFT * 1e6, 10 .* log10(deconvImg(2551: 2650, 91: end))); 
colormap hot; colorbar; caxis([-5, 6])

% We use the same technique to fit the deconvolved image (experimental) to
% obtain the location of top and bottom interfaces
axialGrid = [2551: 2650]';
thickTrue = thickFit.p1 .* lateralGrid + thickFit.p2; 
for iCol = 1: 150
    if (iCol > 1)
        options.StartPoint = [gauss2Fit{iCol - 1}.a1 gauss2Fit{iCol - 1}.b1...
        gauss2Fit{iCol - 1}.c1 gauss2Fit{iCol - 1}.a2 ...
        gauss2Fit{iCol - 1}.b2 gauss2Fit{iCol - 1}.c2];
    end
    gauss2Fit{iCol} = fit(axialGrid, deconvImg(2551: 2650, iCol), 'gauss2');
    locTopDeconv(iCol) = gauss2Fit{iCol}.b1;
    peakTopDeconv(iCol) = gauss2Fit{iCol}.a1;
    locBotDeconv(iCol) = gauss2Fit{iCol}.b2;
    peakBotDeconv(iCol) = gauss2Fit{iCol}.a2;
    thickAirWedgeDeconv(iCol) = abs(locTopDeconv(iCol) - locBotDeconv(iCol));
    errDeconv(iCol) = thickAirWedgeDeconv(iCol) - thickTrue(iCol);
end
% Fit the deconvolved image (simulated)
for iCol = 1: 150
    if (iCol > 1)
        options.StartPoint = [gauss2Fit{iCol - 1}.a1 gauss2Fit{iCol - 1}.b1...
        gauss2Fit{iCol - 1}.c1 gauss2Fit{iCol - 1}.a2 ...
        gauss2Fit{iCol - 1}.b2 gauss2Fit{iCol - 1}.c2];
    end
    gauss2Fit{iCol} = fit(axialGrid, simDeconvImg(2551: 2650, iCol), 'gauss2');
    locTopSimDeconv(iCol) = gauss2Fit{iCol}.b1;
    peakTopSimDeconv(iCol) = gauss2Fit{iCol}.a1;
    locBotSimDeconv(iCol) = gauss2Fit{iCol}.b2;
    peakBotSimDeconv(iCol) = gauss2Fit{iCol}.a2;
    thickAirWedgeSimDeconv(iCol) = abs(locTopSimDeconv(iCol) - locBotSimDeconv(iCol));
    errSimDeconv(iCol) = thickAirWedgeSimDeconv(iCol) - thickTrue(iCol);
end

% The following codes reconstruct the image by using the proposed
% optimization technique
% WARNING: depending on the chosen parameters, the program might take hours
% to days to complete. It is suggested to skip this part and directly load
% the results
% [zGridRecn, kGrid] = meshgrid(zRecn, k);
% matFourRecn = exp(2j * zGridRecn .* kGrid);
% specRecn = repmat(Sk, 1 , numRecn);
% matTranRec = specRecn .* matFourRecn;
% D = eye(numRecn);
% for iCol = 1: size(fringe, 2)
%     [recImg(:, iCol), history] = lasso(matTranRec, fringe(:, iCol), D, lambda, 10, 1.2);
% end
% str = sprintf('lambda_%d_factor_%d.mat', lambda, factor);
% save(str, 'recImg','-v6');

recImg = h5read('recExpImg_factor_8.h5','/rawData');

figure('Name','The reconstructed image (Proposed method)');
set(gcf,'outerposition',get(0,'screensize'));
imagesc(lateralGrid * dx * 1e6, [2551: 2650] * dzFFT * 1e6, 10 .* log10(recImg(2551: 2650, 91: end))); 
colormap hot; colorbar; caxis([12, 22])

axialGrid = [2551: 2650]';
for iCol = 1: 150
    if (iCol > 1)
        options.StartPoint = [gauss2Fit{iCol - 1}.a1 gauss2Fit{iCol - 1}.b1...
        gauss2Fit{iCol - 1}.c1 gauss2Fit{iCol - 1}.a2 ...
        gauss2Fit{iCol - 1}.b2 gauss2Fit{iCol - 1}.c2];
    end
    gauss2Fit{iCol} = fit(axialGrid, recImg(2551: 2650, iCol), 'gauss2');
    locTopRecn(iCol) = gauss2Fit{iCol}.b1;
    peakTopRecn(iCol) = gauss2Fit{iCol}.a1;
    locBotRecn(iCol) = gauss2Fit{iCol}.b2;
    peakBotRecn(iCol) = gauss2Fit{iCol}.a2;
    thickAirWedgeRecn(iCol) = abs(locTopRecn(iCol) - locBotRecn(iCol));
    errRecn(iCol) = thickAirWedgeRecn(iCol) - thickTrue(iCol);
end