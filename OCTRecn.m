%% OCTRecn
% INPUT: inputFringe -- the measured spectral interferogram
%        lambda -- the Lagrange multiplier
%        factor -- super-resolving factor
%        options -- struct that saves device/simulation specific paramters
% OUTPUT: recImg -- reconstructed image
% AUTHORS: Yuye Ling and Mengyuan Wang
% HISTORY: Created 2020/06/12

%%
function [recImg] = OCTRecn(inputFringe, lambda, factor, options)
    % numSpec defines the sampling number of the spectral interferogram
    % numRecn defines the grid size of the original function
    numRecn = options.numSpec;
    dzRecn = options.dzFFT / factor;
    zRecn = linspace(0, (numRecn - 1) * dzRecn, numRecn)';

%     fringe = h5read('rawSpectrumOnionThorlab.h5','/rawData');
    [zGridRecn, kGrid] = meshgrid(zRecn, options.k);
    matFourRecn = exp(2j * zGridRecn .* kGrid);
    specRecn = repmat(options.Sk, 1 , numRecn);
    matTranRec = specRecn .* matFourRecn;
    D = eye(numRecn);
    recImg = zeros(numRecn, size(inputFringe, 2));
    for iCol = 1: size(inputFringe, 2)
        iCol
        [recImg(:, iCol), history] = lasso(matTranRec, inputFringe(:, iCol), D, lambda, 10, 1.2);
    end
    str = sprintf('onion_lambda_%d_factor_%d.mat', lambda, factor);
    save(str, 'recImg','-v6');
end