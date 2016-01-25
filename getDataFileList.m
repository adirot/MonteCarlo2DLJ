    function fileListOrgbyT = getDataFileList(inputFileList)
        fileListOrgbyT = {};
        fileList = inputFileList;
        rho = struct;
        i = 1;
        
        while ~isempty(fileList);
            rho(i).rho = [];
            data = matfile(fileList{1,1});
            sim = data.simulationParam;
            T(i) = sim.T;
            rho(i).rho = [rho(i).rho sim.rho];
            fileListOrgbyT{i,1} = fileList{1,1};
            
            ind = 2;
            indnewlist = 1;
            newfileList = fileList;
            for j = 2:length(fileList)
		
                data = matfile(fileList{1,j});
                sim = data.simulationParam;
                thisT = sim.T;
                
                if thisT == T(i)
                    fileListOrgbyT{i,ind} = fileList{1,j};
                    newfileList = ...
                        {newfileList{1,1:(indnewlist-1)}...
                        newfileList{1,(indnewlist+1):end}};
                    rho(i).rho = [rho(i).rho sim.rho];
                    indnewlist = indnewlist - 1;
                    ind = ind + 1;

                end
                indnewlist = indnewlist + 1;
            end
            fileList = newfileList;
            fileList = fileList(1,2:end);
            i = i + 1;
        end
        clear i;
        
        % sort by Temperature
        [~, indsorted] = sort(T); 
        fileListOrgbyT(indsorted,:) = fileListOrgbyT(:,:);
        
        % sort by density
        for i = 1:length(T)
            [~, indsorted] = sort(unique(rho(i).rho));
            fileListOrgbyT(i,indsorted) =...
                fileListOrgbyT(i,1:length(unique(rho(i).rho)));
        end
        
    end
