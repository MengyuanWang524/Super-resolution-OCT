% The preprocessing of the air wedge data: 
% Extract the location of the top and bottom interfaces of the air wedge
% by using a two-term Gaussian model
% Author: Yuye Ling

factor = 8;
linData = h5read('rawSpectrumAirWedgeThorlab.h5','/rawData');
img = abs(ifft(linData(:, :), size(linData, 1) * factor));
figure('Name','Original air wedge image'); 
imagesc(10 .* log10(img(1: 1024 * factor, :))); colormap hot; colorbar
caxis([0, 20])

axialGrid = [280 * factor + 1: 480 * factor]';
cropImg = img(axialGrid, :);
figure('Name','Cropped image'); 
imagesc(10.*log10(cropImg)); colormap hot; colorbar
caxis([0, 20])


options = fitoptions('gauss2', 'Robust', 'LAR', 'MaxFunEvals', 6000, ...
    'MaxIter', 4000);
% The start point is estimated by inspection
options.StartPoint = [848 / factor 317 * factor 0.8 500 / factor 346 * factor 0.9];
options.Lower = [0 0 0 0 0 0];
for iCol = 1: size(cropImg, 2)
    if (iCol > 1)
        % The starting point for the new fitting will be set to the fitting
        % parameters from the last one
        options.StartPoint = [gauss2Fit{iCol - 1}.a1 gauss2Fit{iCol - 1}.b1...
            gauss2Fit{iCol - 1}.c1 gauss2Fit{iCol - 1}.a2 ...
            gauss2Fit{iCol - 1}.b2 gauss2Fit{iCol - 1}.c2];
    end
    gauss2Fit{iCol} = fit(axialGrid, cropImg(:, iCol), 'gauss2', options);
    locTop(iCol) = gauss2Fit{iCol}.b1;
    peakTop(iCol) = gauss2Fit{iCol}.a1;
    locBot(iCol) = gauss2Fit{iCol}.b2;
    peakBot(iCol) = gauss2Fit{iCol}.a2;
    thickAirWedge(iCol) = abs(locTop(iCol) - locBot(iCol));
end

figure('Name','Exemplary fitting results'); 
plot(axialGrid, cropImg(:, 2000));
hold on; plot(gauss2Fit{2000})

lateralGrid = [2001: 1: 2450];
thick = thickAirWedge(lateralGrid);
[thickFit, paraFit] = fit(lateralGrid', thick', 'poly1', 'Exclude', [398: 400]);
angleWedge = atand(thickFit.p1);
% lateralGrid = [2351: 2500]';


figure('Name','Fitting the angle of the wedge');
plot(lateralGrid, thick);
hold on; plot(thickFit)
text(2020,15,sprintf('Rsq = %g\nAdj = %g',paraFit.rsquare,paraFit.adjrsquare))

% The top and bottom interfaces are fitted as well
% We only use the last column 2321 to 2470 before the tipping point of the
% wedge, becuase the interfaces by themselves have aberrations
lateralGrid = [2351: 1: 2450];
top = locTop(lateralGrid);
bot = locBot(lateralGrid);
[topFit, paraFitTop] = fit(lateralGrid', top', 'poly1');
[botFit, paraFitBot] = fit(lateralGrid', bot', 'poly1', 'Exclude', [47: 49]);

figure('Name','Fitting the top interface of the wedge');
plot(lateralGrid, top);
hold on; plot(topFit)
text(2375,2569,sprintf('Rsq = %g\nAdj = %g',paraFitTop.rsquare,paraFitTop.adjrsquare))

figure('Name','Fitting the bottom interface of the wedge');
plot(lateralGrid, bot);
hold on; plot(botFit)
text(2375,2585,sprintf('Rsq = %g\nAdj = %g',paraFitBot.rsquare,paraFitBot.adjrsquare))

save('wedgePreProc.mat');