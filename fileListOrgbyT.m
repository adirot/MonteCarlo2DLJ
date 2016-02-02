function fileListOrgbyT = fileListOrgbyT(varargin)
    

    % check input, build isotherm object from data files if necessary
    if nargin == 0
        isothermsOrFileList = dir('*mat');
        isothermsOrFileList = {isothermsOrFileList.name};
    else
        isothermsOrFileList = varargin{:};
    end
    
    if isobject(isothermsOrFileList)
        isotherms = isothermsOrFileList;
    else
        fileListOrgbyT = {};
        fileList = isothermsOrFileList;
        i = 1;
        while ~isempty(fileList);
            data = matfile(fileList{1,1});
            sim = data.simulationParam;
            T(i) = sim.T;
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
        
        % sort by Temprature
        [~, indsorted] = sort(T); 
        fileListOrgbyT(indsorted,:) = fileListOrgbyT(:,:);
        
    end
    %colorPlot(1./rho,pressure,'addLegend',leg);
end