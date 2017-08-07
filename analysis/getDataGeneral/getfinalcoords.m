%% get final coords

tind = 0;
for t = [0.45 0.6 0.8 1 1.5 2]
    tind = tind + 1;
    rind = 0;
    for r = [0.005 0.01 0.05 0.1]
        rind = rind + 1;
        mind = 0;
        for m = 3:6
            mind = mind + 1;
            
            list = dir(['N*T' my_num2str(t) 'rho' my_num2str(r) '*_m'...
                my_num2str(m) '*.mat']);
            list = {list.name}';      
            if length(list) > 1
                error('more than one file found');
            end
            
            try
                data = matfile(list{1,1},'Writable',true);
                lastind = data.indIndata;
                coords{tind,rind,mind} = data.allCoords(1:2,1:625,lastind);
            catch
                if isempty(list)
                    disp(['Missing file: T = ' num2str(t) ' rho = '...
                        num2str(r) ' m = ' num2str(m)]);
                else
                    disp(['Unknown error in : T = ' num2str(t) ' rho = '...
                        num2str(r) ' m = ' num2str(m)]);
                end
            end
        end
    end
end

save('finalCoords.mat', 'coords');