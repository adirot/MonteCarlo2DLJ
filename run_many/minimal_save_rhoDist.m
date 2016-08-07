%% minimal save mode: rhoDistrib %%
r = 0.1; m = 3; sq = 256; sqind = sqrt(sq); t = 2; simind = 9;
%list = dir(['N*T' my_num2str(t) 'rho' my_num2str(r) '*_m' num2str(m) '*mat']);
%list = {list.name};
%if length(list) > 1
%    error('more than one file found!');
%end

%M = MC2DLJoutput(list{1,1});

%load('rhoDisribDiffnumOfSquaresOptionb.mat','meanhistnumOfCellsInSquare');
%sumhistnumOfPartInSquare =...
%    meanhistnumOfCellsInSquare{tind,1,1,sqind}*sq*M.indIndata;
%meanhistnumOfPartInSquare_inStep{tind,1} = meanhistnumOfCellsInSquare{tind,1,1,sqind};

steps(1,1) = 1;
[M, histXnumOfPartInSquare, histnumOfPartInSquare,...
               numOfPartInSquare] =...
               M.calcRhoDistrib(sq,'coordsFromObj',true);
sumhistnumOfPartInSquare = histXnumOfPartInSquare;
           
totSec = 0;
for i = 1:1000
    tic;
    steps(1,1+i) = steps(1,1) + i; 
    M = M.MonteCarlo(10,0,'save2data',false);
    [M, histXnumOfPartInSquare, histnumOfPartInSquare,...
               numOfPartInSquare] =...
               M.calcRhoDistrib(sq,'coordsFromObj',true);
    
    sumhistnumOfPartInSquare = sumhistnumOfPartInSquare +...
        histnumOfPartInSquare;
    meanhistnumOfPartInSquare_inStep{simind,i+1} =...
        sumhistnumOfPartInSquare/steps(1,1+i);

    if mod(i,100)
        save(['minimal_save_rhoDist_t' my_num2str(t) 'rho'...
            my_num2str(r) 'simind' my_num2str(simind) '.mat'],...
            'meanhistnumOfPartInSquare_inStep','steps','-v7.3');
    end
    
    % Create log file
    if i > 1
        delete(lastFileName);
    end
    
    totSec = totSec + toc;
    lastFileName = ['T' my_num2str(t) 'rho' my_num2str(r)...
        'm' num2str(m) 'simind' my_num2str(simind)...
        'stepsDone' num2str(i) 'secPassed'...
        my_num2str(totSec) '.txt'];
    fileID = fopen(lastFileName, 'w');
    fclose(fileID);
    
    
end


save(['minimal_save_rhoDist_t' my_num2str(t) 'rho'...
            my_num2str(r) 'simind' my_num2str(simind) '.mat'],...
            'meanhistnumOfPartInSquare_inStep','steps','-v7.3');