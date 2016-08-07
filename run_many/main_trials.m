N = 64;
L = 12;
rho = N/L^2;
T = 1;
Nsteps = 100000;
maxdr = 1;
rCutoff = 2.5;
initialConfig = 'random';

close all;
[finalU,finalConfiguration] = ...
    MonteCarlo2DLJ(N,T,rho,Nsteps,maxdr,initialConfig,rCutoff);

plotParticles(finalConfiguration,L,2^(1/6)/2);
title('end');
