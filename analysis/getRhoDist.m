%% get all rho distributions, mean over a few simulations
numOfSquares = 16^2;
folderName = 'rhoDist/';

totSec = 0;
tic; 
tind = 0;
for t = [1 1.5 2]
    tind = tind + 1;
    rind = 0;
    for r = 0.1
        rind = rind + 1;
        mind = 0;
        for m = 3
            disp([num2str(tind) num2str(rind) num2str(mind)]);
            mind = mind + 1;
            
            list = dir([folderName 'N*T' my_num2str(t) 'rho' my_num2str(r) '*_m'...
                my_num2str(m) '*.mat']);
            list = {list.name}';      
            numOfRuns = length(list);

            for ii = 1:numOfRuns
                try
                    M = MC2DLJoutput([folderName list{ii,1}]);
                catch
                    if isempty(list)
                        disp(['Missing file: T = ' num2str(t) ' rho = '...
                            num2str(r) ' m = ' num2str(m)]);
                    else
                        disp(['Unknown error in : T = ' num2str(t) ' rho = '...
                            num2str(r) ' m = ' num2str(m)]);
                    end
                end
                [M, rho{tind,rind,mind,ii}, PL{tind,rind,mind,ii},...
                cellsInSq{tind,rind,mind,ii}] =...
                M.calcRhoDistrib(numOfSquares);
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
                    'rho' my_num2str(r) 'm' num2str(m)];
                if and(tind == 1,and(rind == 1,and(mind == 1,ii == 1)))
                    oldLogFileName = [];
                else
                    oldLogFileName = newLogFileName;
                end
                totSec = totSec + toc;
                newLogFileName =...
                    createLogFile(logFileInit,oldLogFileName,totSec);
            end
            % avg
            PLavg{tind,rind,mind} = PLavg{tind,rind,mind}/numOfRuns;
            PLNavg{tind,rind,mind} = PLNavg{tind,rind,mind}/numOfRuns;
            
            % plot
            plotArea =...
                sum(rho{tind,rind,mind,1}(2)*PLavg{tind,rind,mind});
            
            plot(rho{tind,rind,mind,1},PLavg{tind,rind,mind}/plotArea);
            
            name = [folderName 'rhoDistrib_avg_Psubsys_T' my_num2str(t)...
                'rho' my_num2str(r) 'm' my_num2str(m)];
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

save([folderName 'rhoDistrib.mat'],...
    'rho', 'PL', 'PLN' ,'PLavg', 'PLNavg' , 'cellsInSq');