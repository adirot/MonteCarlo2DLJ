
list = dir('hcr/N*T0_0*mat');
list = {list.name};

for i = 1:length(list)
    M = MC2DLJoutput(['hcr/' list{1,i}]);
    M.showStep('last');
   saveas(gcf,['hcr/snapshotN' num2str(M.simulationParam.N) 'T' my_num2str(M.simulationParam.T) 'rho' my_num2str(M.simulationParam.rho) 'm' num2str(M.simulationParam.m) 'hcr_last.fig']);
   disp(i);
end
