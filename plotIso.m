function [isotherms,fit,canGetUfromgRind] = plotIso(varargin)
    
    % check input, build isotherm object from data files if necessary

        p = inputParser();
        addOptional(p, 'N', []); 
        addOptional(p, 'Visible', 'on');
        addOptional(p, 'saveFig', 'off');
        addOptional(p, 'fitprop', []);
        addOptional(p, 'residuals', []);
        addOptional(p, 'fileNameEnd', '');
        addOptional(p, 'isotherms', []);
        addOptional(p, 'fit', []);
        addOptional(p, 'canGetUfromgRind', []);
        addOptional(p, 'fileListOrgbyT', []);
        parse(p, varargin{:});
        Results = p.Results;
        N = Results.N;
        Visible = Results.Visible;
        saveFig = Results.saveFig;
        fitprop = Results.fitprop;
        residuals = Results.residuals;
        fileNameEnd = Results.fileNameEnd;
        isotherms = Results.isotherms;
        fit = results.fit;
        canGetUfromgRind = Results.canGetUfromgRind;
        fileListOrgbyT = Results.fileListOrgbyT;
        
        if isempty(N) && isempty(isotherms) && isempty(fileListOrgbyT)
            error(['provide number of particles,'...
                'for example: plotIso(''N'',625), or isotherm object, or fileListOrgbyT']);
        else
            isothermsOrFileList = dir(['N' num2str(N) '*mat']);
            isothermsOrFileList = {isothermsOrFileList.name};
        end
            
    
    if isempty(isotherms)
        if isempty(fileListOrgbyT)
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
        
        [Niso, ~] = size(fileListOrgbyT);
        
        for i = 1:Niso
            if i == 1
                isotherms = isotherm(fileListOrgbyT{i,:});
            else
                isotherms(i) = isotherm(fileListOrgbyT{i,:});
            end
        end
    end
                    
                
    if isempty(Niso)
        Niso = length(isotherms);
    end
    Nrho = length(isotherms(1).rho);
    pressure = zeros(Niso,Nrho);
    rho = zeros(Niso,Nrho);
    
    for i = 1:Niso
        rho(i,:) = isotherms(i).rho;
        pressure(i,:) = isotherms(i).pressure;
        leg{1,i} = ['T = ' num2str(isotherms(i).T)];
        
    end
    
    h1 = figure('Visible',Visible);
    colorPlot(rho,pressure,'addLegend',leg,'figHandle',h1);
    title(['isotherms for N = ' num2str(N)]);
    xlabel('density, reduced units');
    ylabel('pressure, reduced units');
    
    if ~isempty(fitprop) && isempty(fit)
        for i = 1:length(fitprop)
            fit{1,i} = fitIso(isotherms,fitprop{i},'figureHandle',h1);
        end
    end
    
    if ~isempty(residuals)
        for i = 1:length(fit)
            if residuals{1,i}
                for iso = 1:length(isotherms)
                    f = fit{1,i}(1,iso);
                    res = my_residuals(pressure(iso,:),f.P);
                    for j = 1:length(res)
                        text(rho(iso,j),f.P(j),num2str(res(j)));
                    end
                end
            end
        end
    end
    
    if isempty(canGetUfromgRind)
        canGetUfromgRind = canGetUfromgR(isotherms,fit,0.1);
    end
    
    for i = 1:length(isotherms)
        plot(isotherms(i).rho(1,canGetUfromgRind{1,i}),...
            isotherms(i).pressure(1,canGetUfromgRind{1,i}),'m*');
    end
    
    if saveFig
        saveas(h1,['isotherms_N' num2str(N) 'P_rho' fileNameEnd '.fig']);
        saveas(h1,['isotherms_N' num2str(N) 'P_rho' fileNameEnd '.jpg']);
    end

    
    h2 = figure('Visible',Visible);
    colorPlot(1./rho,pressure,'addLegend',leg,'figHandle',h2);
    title(['isotherms for N = ' num2str(N)]);
    xlabel('volume, reduced units');
    ylabel('pressure, reduced units');
    
   
    if saveFig
        saveas(h2,['isotherms_N' num2str(N) 'P_V' fileNameEnd '.fig']);
        saveas(h2,['isotherms_N' num2str(N) 'P_V' fileNameEnd '.jpg']);
    end
    
end