folderName = '';
list = dir([folderName 'N625*T0_*rho0_1*m6' '*mat']);
list = {list.name};

for i = 1:length(list)
    M = MC2DLJoutput([folderName list{1,i}]);
    M.showStep('last');
    title([folderName 'snapshot N = ' num2str(M.simulationParam.N) ' T = ' my_num2str(M.simulationParam.T) ' \rho = ' my_num2str(M.simulationParam.rho) ' m = ' num2str(M.simulationParam.m) ' sweeps ' num2str(M.indIndata)]);
   saveas(gcf,[folderName 'snapshotN' num2str(M.simulationParam.N) 'T' my_num2str(M.simulationParam.T) 'rho' my_num2str(M.simulationParam.rho) 'm' num2str(M.simulationParam.m) '_last.fig']);
   disp(i);
   close all;
end

for i = 1:length(list)
    M(i) = MC2DLJoutput([folderName '/' list{1,i}]);
   indIndata(i) = M.indIndata;
   disp(indIndata(i));
end

