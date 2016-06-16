%% get all rho distributions
numOfSquares = 900;

tind = 0;
for t = [0.45 0.6 0.8 1 1.5 2]
    tind = tind + 1;
    rind = 0;
    for r = [0.005 0.01 0.05 0.1]
        rind = rind + 1;
        mind = 0;
        for m = 3:6
            disp([num2str(tind) num2str(rind) num2str(mind)]);
            mind = mind + 1;
            
            list = dir(['N*T' my_num2str(t) 'rho' my_num2str(r) '*_m'...
                my_num2str(m) '*.mat']);
            list = {list.name}';      
            if length(list) > 1
                error('more than one file found');
            end

            try
                M = MC2DLJoutput(list{1,1});
            catch
                if isempty(list)
                    disp(['Missing file: T = ' num2str(t) ' rho = '...
                        num2str(r) ' m = ' num2str(m)]);
                else
                    disp(['Unknown error in : T = ' num2str(t) ' rho = '...
                        num2str(r) ' m = ' num2str(m)]);
                end
            end
            [M, rho{tind,rind,mind}, PL{tind,rind,mind},...
                cellsInSq{tind,rind,mind}] =....
                M.calcRhoDistrib(numOfSquares);
            cellArea{tind,rind,mind} = (M.simulationParam.L)^2/numOfSquares;
            
            plot(rho{tind,rind,mind},...
                cellsInSq{tind,rind,mind}.*PL{tind,rind,mind}/625);
            
            name = ['rhoDistrib_T' my_num2str(t) 'rho' my_num2str(r)...
                'm' my_num2str(m)];
            save(gcf, [name '.fig']);
            save(gcf, [name '.jpg']);
            close all;
            
        end
    end
end

save('rhoDistrib.mat', 'rho', 'PL', 'cellsInSq');