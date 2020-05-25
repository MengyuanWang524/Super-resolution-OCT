% The preprocessing of the air wedge data: 
% Extract the location of the top and bottom interfaces of the air wedge
% by using a two-term Gaussian model
% Author: Yuye Ling

linData = h5read('rawSpectrumAirWedgeThorlab.h5','/rawData');
img = abs(ifft(linData(:, :), size(linData, 1)));
figure('Name','Original air wedge image'); 
imagesc(10.*log10(img(1: 1024, :)));colormap hot; colorbar
caxis([0, 20])

cropImg = img(281: 480, :);
figure('Name','Cropped image'); 
imagesc(10.*log10(cropImg));colormap hot; colorbar
caxis([0, 20])
index = [281: 1: 480]';

options = fitoptions('gauss2', 'Robust', 'LAR', 'MaxFunEvals', 6000, ...
    'MaxIter', 4000);
% The start point is estimated by inspection
options.StartPoint = [550 307 0.8 50 440 0.9];
options.Lower = [0 0 0 0 0 0];
for iCol = 1: size(cropImg, 2)
    if (iCol > 1)
        % The starting point for the new fitting will be set to the fitting
        % parameters from the last one
        options.StartPoint = [gauss2Fit{iCol - 1}.a1 gauss2Fit{iCol - 1}.b1...
            gauss2Fit{iCol - 1}.c1 gauss2Fit{iCol - 1}.a2 ...
            gauss2Fit{iCol - 1}.b2 gauss2Fit{iCol - 1}.c2];
    end
    gauss2Fit{iCol} = fit(index, cropImg(:, iCol), 'gauss2', options);
    locTop(iCol) = gauss2Fit{iCol}.b1;
    peakTop(iCol) = gauss2Fit{iCol}.a1;
    locBot(iCol) = gauss2Fit{iCol}.b2;
    peakBot(iCol) = gauss2Fit{iCol}.a2;
    thickAirWedge(iCol) = abs(locTop(iCol) - locBot(iCol));
end

figure('Name','Exemplary fitting results'); 
plot(index, cropImg(:, 2000));
hold on; plot(gauss2Fit{2000})

% We use the data (sep) up to 2475th column to fit the angle between the
% top and bottom interface
% The lateral pixel size is 2 um (3000 pixels for 6 mm FOV)
% The axial pixel size is 1.9438 um
dx = 2e-6;

lateralGrid = [2000: 1: 2449];
thick = thickAirWedge(2001:2450);
[thickFit, paraFit] = fit(lateralGrid', thick', 'poly1');
angleWedge = atand(thickFit.p1);

figure('Name','Fitting the angle of the wedge'); 
plot(lateralGrid, thick);
hold on; plot(thickFit)
text(4.72e-3,5,sprintf('Rsq = %g\nAdj = %g',paraFit.rsquare,paraFit.adjrsquare))

% The top and bottom interfaces are fitted as well
% We only use the last 100 columns before the tipping point of the wedge
lateralGrid = [2100 * dx: dx: dx * 2249];
top = locTop(2101:2250);
bot = locBot(2101:2250);
[topFit, paraFit] = fit(lateralGrid', top', 'poly1');
[botFit, paraFit] = fit(lateralGrid', bot', 'poly1');

save('wedgePreProc.mat');