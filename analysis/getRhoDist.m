%% get all rho distributions, mean over a few simulations
%  list = dir('N625*T0_0*mat'); list = {list.name}; for ii = 1:length(list) M = MC2DLJoutput([list{1,ii}]); M.data.allDists = []; disp(ii); end;
% getRhoDist(18^2,'',t,0.1,6)
function [PL, PLN ,PLavg, PLNavg , cellsInSq] =...
    getRhoDist(numOfSquares,folderName,inputT,inputRho,inputm)

totSec = 0;
tic; 
tind = 0;
fileListind = 1;

for t = inputT
    tind = tind + 1;
    rind = 0;
    for r = inputRho
        rind = rind + 1;
        mind = 0;
        for m = inputm
            disp([num2str(tind) num2str(rind) num2str(mind)]);
            mind = mind + 1;
            
            list = dir([folderName 'N625*T' my_num2str(t) 'rho' my_num2str(r) '*_m'...
                my_num2str(m) '*.mat']);
            list = {list.name}';      
            numOfRuns = length(list);

            for ii = 1:numOfRuns
                try
                    M = MC2DLJoutput([folderName list{ii,1}]);
                    filelist{1,fileListind} = [folderName list{ii,1}];
                    fileListind = fileListind + 1;
                    
                    [M, rho{tind,rind,mind,ii}, PL{tind,rind,mind,ii},...
                    cellsInSq{tind,rind,mind,ii}] =...
                    M.calcRhoDistrib(numOfSquares,'startFrom',4000);
                    PLN{tind,rind,mind,ii} = ...
                        PL{tind,rind,mind,ii}.*(0:M.simulationParam.N);

                    %avg
                    if ii == 1
                        PLavg{tind,rind,mind} = PL{tind,rind,mind,ii};
                        PLNavg{tind,rind,mind} = PLN{tind,rind,mind,ii};
                    else
                        PLavg{tind,rind,mind} =...
                            PLavg{tind,rind,mind} + PL{tind,rind,mind,ii}; 
                        PLNavg{tind,rind,mind} = ...
                            PLNavg{tind,rind,mind} + PLN{tind,rind,mind,ii};
                    end

                    %create log file
                    logFileInit = [folderName 'rhoDistrib_T' my_num2str(t)...
                        'rho' my_num2str(r) 'm' num2str(m)...
                        'numOfSquares' num2str(numOfSquares) list{ii,1}];
                    if and(tind == 1,and(rind == 1,and(mind == 1,ii == 1)))
                        oldLogFileName = [];
                    else
                        oldLogFileName = newLogFileName;
                    end
                    totSec = totSec + toc;
                    newLogFileName =...
                        createLogFile(logFileInit,oldLogFileName,totSec);

                catch
                    if isempty(list)
                        disp(['Missing file: T = ' num2str(t) ' rho = '...
                            num2str(r) ' m = ' num2str(m)]);
                        
                        %create log file
                        logFileInit = [folderName 'rhoDistrib_T' my_num2str(t)...
                            'rho' my_num2str(r) 'm' num2str(m)...
                            'numOfSquares' num2str(numOfSquares) 'missingFile' list{ii,1}];
                        if and(tind == 1,and(rind == 1,and(mind == 1,ii == 1)))
                            oldLogFileName = [];
                        else
                            oldLogFileName = newLogFileName;
                        end
                        totSec = totSec + toc;
                        newLogFileName =...
                            createLogFile(logFileInit,oldLogFileName,totSec);
                    else
                        disp(['Unknown error in : T = ' num2str(t) ' rho = '...
                            num2str(r) ' m = ' num2str(m)]);
                        
                        %create log file
                        logFileInit = [folderName 'rhoDistrib_T' my_num2str(t)...
                            'rho' my_num2str(r) 'm' num2str(m)...
                            'numOfSquares' num2str(numOfSquares) 'unknownError' list{ii,1}];
                        if and(tind == 1,and(rind == 1,and(mind == 1,ii == 1)))
                            oldLogFileName = [];
                        else
                            oldLogFileName = newLogFileName;
                        end
                        totSec = totSec + toc;
                        newLogFileName =...
                            createLogFile(logFileInit,oldLogFileName,totSec);
                    end
                end
            end
            % avg
            PLavg{tind,rind,mind} = PLavg{tind,rind,mind}/numOfRuns;
            PLNavg{tind,rind,mind} = PLNavg{tind,rind,mind}/numOfRuns;
            
            % plot
            plotArea =...
                sum(rho{tind,rind,mind,1}(2)*PLavg{tind,rind,mind});
            
            plot(rho{tind,rind,mind,1},PLavg{tind,rind,mind}/plotArea);
            
            name = [folderName 'rhoDistrib_avg_Psubsys_T' my_num2str(t)...
                'rho' my_num2str(r) 'm' my_num2str(m)...
                'numOfSquares' num2str(numOfSquares)];
            title(name);
            saveas(gcf, [name '.fig']);
            saveas(gcf, [name '.jpg']);
            close all;
            
            plotArea =...
                sum(rho{tind,rind,mind,1}(2)*PLNavg{tind,rind,mind});
            
            plot(rho{tind,rind,mind,1},PLNavg{tind,rind,mind}/plotArea);
            
            name = [folderName 'rhoDistrib_avg_Pcells_T' my_num2str(t)...
                'rho' my_num2str(r) 'm' my_num2str(m)];
            title(name);
            saveas(gcf, [name '.fig']);
            saveas(gcf, [name '.jpg']);
            close all;
            
        end
    end
end

save([folderName 'rhoDistrib_date' nowdatetimestr() '.mat'],...
    'inputT','inputRho','inputm','rho','numOfSquares',...
    'PL', 'PLN' ,'PLavg', 'PLNavg' , 'cellsInSq','filelist','-v7.3');