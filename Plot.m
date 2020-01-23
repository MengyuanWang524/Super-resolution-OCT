clear all
% file name 
    load ('twospike_30_500_1e-06_100.mat')
% I equals to the iteration number
    I = 1 ;
    T =500;
    N = 1000;
% dz0 equals to that in supperresolutionFunction.m
    dz0 = 0.1e-6;
 % self-defined grid. 8.568e-6 is grid defined by Fourier relationship
    dz =1e-6;
    z_range = (N - 1) * 1 * dz0;
    Original_z = linspace(0, z_range, N) * 1e6;
    Reconstructed_z = linspace(0, dz * (T-1), T) * 1e6;
    for i = 1:I
    figure; plot( Reconstructed_z, abs(result(1, 1:T, i)))
    hold on; plot( Original_z, abs(result(2, :, i)) / 2)
    end
