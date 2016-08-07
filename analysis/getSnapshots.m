
list = dir('N*T*rho0_1*m3*mat');
list = {list.name};

for i = 1:length(list)
    M = MC2DLJoutput(list{1,i});
    M.showStep('last');
   save(gcf,['snapshotN' num2srt(M.simulationParam.N) 'T' my_num2str(M.simulationParam.T) 'rho' my_num2str(M.simulationParam.rho) 'm' num2str(M.simulationParam.m) 'last.fig']);
end
