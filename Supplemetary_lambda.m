% Simulation data for supplemetary Figure3 & 4: 
% Authors: Mengyuan Wang and Yuye Ling
% 
% M defines the sampling number of the spectral interferogram
% N defines the grid size of the original function
% T defines the grid size of the reconstructed function
close all;
    SNR =  15; 
    separation= 80; 
    lambda= [10,100,500,1000]; 
    dz = 1e-6
%     dz = [1,4,8]* 1e-6;
    M = 500;
    N = 1000;
    T = 400;
 % Simulation set-up
    lambda0 = 1310e-9;
    FWHM_lambda = 30e-9;
    lambda_st = lambda0 - 50e-9;
    lambda_end = lambda0 + 50e-9;
    k_st = 2 * pi / lambda_end;
    k_end = 2 * pi / lambda_st;
    k0 = 2 * pi / lambda0;
    delta_k = (pi / sqrt(log(2))) * (FWHM_lambda / lambda0^2);
    k = linspace(2 * pi / lambda_st, 2 * pi / lambda_end, M)';
    Sk = exp(-((k-k0)/delta_k).^2);
    dz_fft = 0.5 * 1 / (1 / lambda_st - 1 / lambda_end);
%     dz = dz_fft
    dz0 = 0.1e-6;
     for p= 1: length(dz)
%define two grids:grids in object domain and reconstruction domain
    gridObj = linspace(0, (N - 1) * dz0, N )';
    gridRec = linspace(0, (T - 1) * dz(p), T )';
    [X0, Y0] = meshgrid(gridObj, k);
    [X, Y] = meshgrid(gridRec, k);
    matFourObj =(exp(2j * X0 .* Y0));
    matFourRec = (exp(2j .* X .* Y));
    specObj = repmat( Sk, 1 , N);
    specRec = repmat( Sk, 1 , T);

    result = zeros(2, N);

 
% Random noise: 
    noise =  randn(M,1);
% Define the actual axial function as xtrue
    rtrue = zeros(N, 1);
%   xtrue(61:90) = 100;
    startpoint = randi([100 800],1,1); 
    rtrue( startpoint ) = 100;
    rtrue( startpoint + separation) = 100;
    z_range = (N - 1) * 1 * dz0;
    Original_z = linspace(0, z_range, N) * 1e6;  
    
    figure
    plot(Original_z,rtrue );
     hold on
    matTranObj = specObj .* matFourObj;
    matTranRec = specRec .* matFourRec; 
    b = matTranObj * rtrue;
    i = 1;
    noise = noise .* sqrt((sum(abs(b).^2) / 10^(SNR/10)) ./ sum(abs(noise).^2));
    bn = b + noise;
    D = eye(T);
    
 for s =1: length(lambda) 
    [x, history] = lasso((matTranRec), bn, D, lambda(s), 10, 1);
    Reconstructed_z = linspace(0, dz(p) * (T-1), T) * 1e6;
   plot(Reconstructed_z, abs(x'));
   str{s}=['lambda=',num2str(lambda(s))];
%    ans = SNR(s)
end
 legend('Object',str)
end