% Numerical simulation for superresolution
% Authors: Mengyuan Wang and Yuye Ling
function[] = superresolutionFunction(SNR , separation ,lambda,dz ,numIter)
% Function superresolution requires input of parameters below:
% SNR: for additive noise in dB
% separation: defined actual separation for xtrue, measured in ¦Ìm
% lambda: Lagrange mutiplier in l1-optimization.
% dz:self-defined grid space, measured in meter. e.g. dz = 1e-07
% numIter: number of iteration in one trial.

% M defines the sampling number of the spectral interferogram
% N defines the grid size of the original function
% T defines the grid size of the reconstructed function
    M = 400;
    N = 1000;
    T = 500;
% Simulation set-up
    lambda0 = 1310e-9;
    FWHM_lambda = 40e-9;
    lambda_st = lambda0 - 50e-9;
    lambda_end = lambda0 + 50e-9;
%     k_st = 2 * pi / lambda_end;
%     k_end = 2 * pi / lambda_st;
    k0 = 2 * pi / lambda0;
    delta_k = (pi / sqrt(log(2))) * (FWHM_lambda / lambda0^2);
    k = linspace(2 * pi / lambda_st, 2 * pi / lambda_end, M)';
    Sk = exp(-((k-k0)/delta_k).^2);
    dz0 = 0.1e-6;
% Define two grids: grids in object domain and reconstruction domain
    gridObj = linspace(0, (N - 1) * dz0, N )';
    gridRec = linspace(0, (T - 1) * dz, T )';
    [X0, Y0] = meshgrid(gridObj, k);
    [X, Y] = meshgrid(gridRec, k);
    matFourObj =(exp(2j * X0 .* Y0));
    matFourRec = (exp(2j .* X .* Y));
    specObj = repmat( Sk, 1 , N);
    specRec = repmat( Sk, 1 , T);
% Initialize a three dimension variable for storing result   
    result = zeros(2, N, numIter);

%%
    for i = 1:numIter
    noise =  randn(M,1);
    rtrue = zeros(N, 1);
    % xtrue(61:90) = 100;
    startpoint = randi([100 800],1,1); 
    rtrue( startpoint ) = 100;
    rtrue( startpoint + separation) = 100;

    %temp = xtrue;
    % noise = noise .* sqrt((sum(xtrue.^2) / 10.^(20/10)) ./ sum(noise.^2));
    % xtrue = xtrue + noise;
    matTranObj = specObj .* matFourObj;
    matTranRec = specRec .* matFourRec;
    %  
    b = matTranObj * rtrue;
    noise = noise .* sqrt((sum(abs(b).^2) / 10^(SNR/10)) ./ sum(abs(noise).^2));
    b = b + noise;
    % b = hilbert(b);
    D = eye(T);
   
    [r, history] = lasso((matTranRec), b, D, lambda, 10, 1);

%%
% Store the result
    result(1,:, i) = [r', zeros( 1, N-T) ];
    result(2,:, i) = rtrue;

end
%%
% Save the result into '.mat', use 'Plot.m' to plot
filename = sprintf('twospike_%g_%u_%g_%g',SNR,lambda,dz,separation);
save ( [filename, '.mat'] ,'result','-v6')
end
