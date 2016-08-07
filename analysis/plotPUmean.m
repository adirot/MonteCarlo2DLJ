function plotPUmean(varargin)
    
        % check input

        p = inputParser();
        addOptional(p, 'N', []); 
        addOptional(p, 'Visible', 'on');
        addOptional(p, 'saveFig', false);
        addOptional(p, 'fileNameEnd', '');
        addOptional(p, 'fileList', []);
        parse(p, varargin{:});
        Results = p.Results;
        N = Results.N;
        Visible = Results.Visible;
        saveFig = Results.saveFig;
        fileNameEnd = Results.fileNameEnd;
        fileList = Results.fileList;
        
        if isempty(fileList)
            fileList = dir(['N' N '*.mat']);
            fileList = {fileList.name};
        end
        
        for i = 1:length(fileList)
            data = matfile(fileList{1,i});
            sim = data.simulationParam;
            T = num2str(sim.T);
            if strcmp(N,'')
                N = sim.N;
                N = num2str(N);
            end
            rho = num2str(sim.rho);
            if ~existInMatfile(fileList{1,i},'meanP') 
                M = MC2DLJoutput(fileList{1,i});
                M.meanProp();
            end
            
            h1 = figure;
            h1 = plot(data.meanP,'Visible',Visible);
            title(['mean pressure for T = ' T ' N = ' N ' \rho = ' rho]);
            xlabel('step in simulation');
            ylabel('mean pressure (reduced units)');
            
            h2 = figure;
            h2 = plot(data.meanU,'Visible',Visible);
            title(['mean energy for T = ' T ' N = ' N ' \rho = ' rho]);
            xlabel('step in simulation');
            ylabel('mean energy (reduced units)');
            
            if saveFig 
                saveas(h1,['meanP_T'...
                    T 'N' N 'rho' my_num2str(sim.rho) fileNameEnd '.fig']);
                saveas(h1,['meanP_T'...
                    T 'N' N 'rho' my_num2str(sim.rho) fileNameEnd '.jpg']);
                saveas(h2,['meanU_T'...
                    T 'N' N 'rho' my_num2str(sim.rho) fileNameEnd '.fig']);
                saveas(h2,['meanU_T'...
                    T 'N' N 'rho' my_num2str(sim.rho) fileNameEnd '.jpg']);
                close all;
            end
        end

end