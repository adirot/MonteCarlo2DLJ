% find stdstd for all files in folder

T = [0.45,0.6,0.8,1,1.5,2];
rho = [0.005,0.01,0.05,0.1];
m = [3,4,5,6];

list = dir('N*mat');
list = {list.name}';
stepRDFstdstd0_1 = zeros(6,4,4);
stepRDFstdstd0_01 = zeros(6,4,4);

for i = 1:length(list)
        varVarPeaks = [];
        M = MC2DLJoutput(list{i,1});
        t = M.simulationParam.T;
        r = M.simulationParam.rho;
        mp = M.simulationParam.m;
        
        RDFp = M.data.RDFhisto(:,30);
        for j = 2:length(RDFp)
            stdRDFp(j-1) = std(RDFp(1:j))/mean(RDFp(1:j));
        end
        for j = 2:length(stdRDFp)
            varVarPeaks(j-1) = std(stdRDFp(1:j))/mean(stdRDFp(1:j));
        end

        [~,ind] = min(abs(varVarPeaks - 0.1));
        [~,indt] = find(T==t);
        [~,indr] = find(rho==r);
        [~,indm] = find(m==mp);
       stepRDFstdstd0_1(indt,indr,indm) = ind;
       [~,ind] = min(abs(varVarPeaks - 0.01));
       stepRDFstdstd0_01(indt,indr,indm) = ind;
        
        disp(i);
end
save('all_stdstdRDFp.mat');