%% check results when run is done

N = M.simulationParam.N;
T = M.simulationParam.T;
rho = M.simulationParam.rho;
m = M.simulationParam.m;

% U,P vs. step
[M, h] = M.plotPropVsStep('U');
saveas(h,['UvsStepN' num2str(N) 'T' my_num2str(T) 'm' num2str(m) 'angdep' '.fig']);
saveas(h,['UvsStepN' num2str(N) 'T' my_num2str(T) 'm' num2str(m) 'angdep' '.jpg']);
close all;
    
[M, h] = M.plotPropVsStep('P');
saveas(h,['PvsStepN' num2str(N) 'T' my_num2str(T) 'm' num2str(m) 'angdep' '.fig']);
saveas(h,['PvsStepN' num2str(N) 'T' my_num2str(T) 'm' num2str(m) 'angdep' '.jpg']);
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
   num2str(M.simulationParam.rho) ' angle dep']);

name = ['varUPvsSteps_T' my_num2str(M.simulationParam.T)...
   'N' my_num2str(M.simulationParam.N) 'rho'...
   my_num2str(M.simulationParam.rho) 'angdep'];
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
   num2str(M.simulationParam.rho) ' angle dep']);

name = ['varVarUPvsSteps_T' my_num2str(M.simulationParam.T)...
   'N' my_num2str(M.simulationParam.N) 'rho'...
   my_num2str(M.simulationParam.rho) 'angdep'];
saveas(gcf,[name '.fig']);
saveas(gcf,[name '.jpg']);
close all;

disp('ploted var var U,P vs step');

% RDF
plot(M.data.RDFbins,mean(M.data.RDFhisto,3));
name = ['RDF_T' my_num2str(M.simulationParam.T)...
   'N' my_num2str(M.simulationParam.N) 'rho'...
   my_num2str(M.simulationParam.rho) 'angdep'];
saveas(gcf,[name '.fig']);
saveas(gcf,[name '.jpg']);
close all;

disp('ploted RDF');


% snapshots
name = ['snapshot_T' my_num2str(M.simulationParam.T)...
   'N' my_num2str(M.simulationParam.N) 'rho'...
   my_num2str(M.simulationParam.rho) 'angdep'];
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
