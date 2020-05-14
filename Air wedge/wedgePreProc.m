% The preprocessing of the air wedge data: 
% Extract the location of the top and bottom interfaces of the air wedge
% by using a two-term Gaussian model
% Authors: Yuye Ling

RawDataNFrame = h5read('raw_spectrum_air_wedge_Thorlab.h5','/rawData');
img = abs(fft(RawDataNFrame(:, :)));
figure('Name','Original air wedge image'); 
imagesc(10.*log10(img(1: 1024, :)));colormap hot; colorbar
caxis([40, 60])

cropImg = img(281: 480, :);
figure('Name','Cropped image'); 
imagesc(10.*log10(cropImg));colormap hot; colorbar
caxis([40, 60])
index = [281: 1: 480]';

options = fitoptions('gauss2', 'Robust', 'LAR', 'MaxFunEvals', 6000, ...
    'MaxIter', 4000);
% The start point is estimated by inspection
options.StartPoint = [7e5 307 0.8 8e4 440 0.9];
options.Lower = [0 0 0 0 0 0];
for iCol = 1: size(cropImg, 2)
    if (iCol > 1)
        % The starting point for the new fitting will be set to the fitting
        % parameters from the previous one
        options.StartPoint = [gauss2Fit{iCol - 1}.a1 gauss2Fit{iCol - 1}.b1...
            gauss2Fit{iCol - 1}.c1 gauss2Fit{iCol - 1}.a2 ...
            gauss2Fit{iCol - 1}.b2 gauss2Fit{iCol - 1}.c2];
    end
    gauss2Fit{iCol} = fit(index, cropImg(:, iCol), 'gauss2', options);
    locTop(iCol) = gauss2Fit{iCol}.b1;
    locBot(iCol) = gauss2Fit{iCol}.b2;
    thickAirWedge(iCol) = abs(locTop(iCol) - locBot(iCol));
end

figure('Name','Exemplary fitting results'); 
plot(index, cropImg(:, 2000));
hold on; plot(gauss2Fit{2000})

% We use the data (sep) up to 2475th column to fit the angle between the
% top and bottom interface
% The lateral pixel size is 2 um (3000 pixels for 6 mm FOV)
% The axial pixel size is 1.9438 um

lateralGrid = [0: 2e-6: 2e-6 * 2474];
thick = thickAirWedge(1: 2475) * 1.9438e-6;
linearFit = fit(lateralGrid', thick', 'poly1');
angleWedge = atand(linearFit.p1);

figure('Name','Fitting the angle of the wedge'); 
plot(lateralGrid, thick);
hold on; plot(linearFit)