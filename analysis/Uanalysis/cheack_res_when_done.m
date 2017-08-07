%% check results when run is done
%s = '_hcr_';

N = M.simulationParam.N;
T = M.simulationParam.T;
rho = M.simulationParam.rho;
m = M.simulationParam.m;

% U,P vs. step
[M, h] = M.plotPropVsStep('U');
saveas(h,['UvsStep' s M.fileName(startind:end-4) '.fig']);
saveas(h,['UvsStepN' s M.fileName(startind:end-4) 'angdep' '.jpg']);
close all;
    
[M, h] = M.plotPropVsStep('P');
saveas(h,['PvsStepN' s M.fileName(startind:end-4) '.fig']);
saveas(h,['PvsStepN' s M.fileName(startind:end-4) '.jpg']);
close all;

disp('ploted U,P vs step');

% var U,P vs. step
[M, varU, varP, steps, varVarU, varVarP] = M.varOfvar('plotVarVsStep',false,...
    'plotVarVarVsStep',false,'saveFig',false);

figure;
hold on;
plot(steps,varU);
plot(steps,varP,'r');
legend({'Energy variance','Pressure variance'});
xlabel('steps');
ylabel('Energy or Pressure variance');
title(['variance for Energy and Pressure Vs. steps, T = '...
   num2str(M.simulationParam.T) ' N = '...
   num2str(M.simulationParam.N) ' \rho = '... 
   num2str(M.simulationParam.rho) s]);

name = ['varUPvsSteps' s M.fileName(startind:end-4)];
saveas(gcf,[name '.fig']);
saveas(gcf,[name '.jpg']);
close all;

disp('ploted var U,P vs step');

% var var U,P vs. step
figure;
hold on;
plot(steps,varVarU);
plot(steps,varVarP,'r');
legend({'Energy variance of variance','Pressure variance of variance'});
xlabel('steps');
ylabel('Energy or Pressure variance of variance');
title(['variance of variance for Energy and Pressure Vs. steps, T = '...
   num2str(M.simulationParam.T) ' N = '...
   num2str(M.simulationParam.N) ' \rho = '... 
   num2str(M.simulationParam.rho) s]);

name = ['varVarUPvsSteps' s M.fileName(startind:end-4) ];
saveas(gcf,[name '.fig']);
saveas(gcf,[name '.jpg']);
close all;

disp('ploted var var U,P vs step');

% RDF
plot(M.data.RDFbins,mean(M.data.RDFhisto,3));
name = ['RDF' s M.fileName(startind:end-4)];
saveas(gcf,[name '.fig']);
saveas(gcf,[name '.jpg']);
close all;

disp('ploted RDF');


% snapshots
name = ['snapshot' s M.fileName(startind:end-4)];
M.showStep('first');
saveas(gcf,[name '1.fig']);
saveas(gcf,[name '1.jpg']);
close all;
M.showStep('mid');
saveas(gcf,[name '2.fig']);
saveas(gcf,[name '1.jpg']);
close all;
M.showStep('last');
saveas(gcf,[name '3.fig']);
saveas(gcf,[name '3.jpg']);
close all;


disp('snapshots');

% inefficiancy
%[M,tauP,tauU] = M.inefficiency(10);


%disp('ploted inefficiancy');
