%% minimal save mode: rhoDistrib %%
r = 0.1; m = 3; sq = 121; sqind = sqrt(sq);
list = dir(['N*T' my_num2str(t) 'rho' my_num2str(r) '*_m' num2str(m) '*mat']);
list = {list.name};
if length(list) > 1
    error('more than one file found!');
end

M = MC2DLJoutput(list{1,1});

load('rhoDisribDiffnumOfSquaresOptionb.mat','meanhistnumOfCellsInSquare');
sumhistnumOfPartInSquare =...
    meanhistnumOfCellsInSquare{tind,1,1,sqind}*sq*M.indIndata;
meanhistnumOfPartInSquare_inStep{tind,1} = meanhistnumOfCellsInSquare{tind,1,1,sqind};
steps(1,1) = M.indIndata;

totSec = 0;
for i = 1:100000
    tic;
    steps(1,1+i) = steps(1,1) + i; 
    M = M.MonteCarlo(10,0,'save2data',false);
    [M, histxnumOfPartInSquare, histnumOfPartInSquare,...
               numOfPartInSquare] =...
               M.calcRhoDistrib(sq,'coordsFromObj',true);
    
    sumhistnumOfPartInSquare = sumhistnumOfPartInSquare +...
        histnumOfPartInSquare;
    meanhistnumOfPartInSquare_inStep{tind,i+1} = sumhistnumOfPartInSquare/steps(1,1+i);

    if mod(i,500)
        save(['minimal_save_rhoDist' my_num2str(t) '.mat'],'meanhistnumOfPartInSquare_inStep','steps','-v7.3');
    end
    
    % Create log file
    if i > 1
        delete(lastFileName);
    end
    
    totSec = totSec + toc;
    lastFileName = ['T' my_num2str(t) 'rho' my_num2str(r)...
        'm' num2str(m) 'stepsDone' num2str(i) 'secPassed'...
        my_num2str(totSec) '.txt'];
    fileID = fopen(lastFileName, 'w');
    fclose(fileID);
    
    
end


save(['minimal_save_rhoDist' my_num2str(t) '.mat'],'meanhistnumOfPartInSquare_inStep','steps','-v7.3');