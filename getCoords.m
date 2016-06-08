
list = dir('N*T1rho*m3*mat');
list = {list.name}'

for i = 1:length(list)
    M = MC2DLJoutput(list{i,1});
    coords = M.currentCoords;
    T = M.simulationParam.T;
    rho = M.simulationParam.rho;
    m = M.simulationParam.m;
   save(['coordsN625T' my_num2str(T) 'rho' my_num2str(rho) 'm' num2str(m) 'last.mat'],'coords');
end
