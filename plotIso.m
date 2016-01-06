function isotherms = plotIso(varargin)
    

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
        
        
        [Niso, ~] = size(fileListOrgbyT);
        
        for i = 1:Niso
            if i == 1
                isotherms = isotherm(fileListOrgbyT{i,:});
            else
                isotherms(i) = isotherm(fileListOrgbyT{i,:});
            end
        end
    end
                    
                

    Niso = length(isotherms);
    Nrho = length(isotherms(1).rho);
    pressure = zeros(Niso,Nrho);
    rho = zeros(Niso,Nrho);
    
    for i = 1:Niso
        rho(i,:) = isotherms(i).rho;
        pressure(i,:) = isotherms(i).pressure;
        leg{1,i} = ['T = ' num2str(isotherms(i).T)];
        
    end
    
    colorPlot(rho,pressure,'addLegend',leg);
    colorPlot(1./rho,pressure,'addLegend',leg);
end