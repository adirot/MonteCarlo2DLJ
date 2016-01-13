function [RDF,fileListOrgbyT,MC2DLJs,length2plot,legT,legrho] = plotRDF(varargin)
    

    % check input, build MC2DLJoutput object from data files if necessary

        p = inputParser();
        addOptional(p, 'N', []); 
        addOptional(p, 'Visible', 'on');
        addOptional(p, 'saveFig', 'off');
        parse(p, varargin{:});
        Results = p.Results;
        N = Results.N;
        Visible = Results.Visible;
        saveFig = Results.saveFig;
        
        if isempty(N)
            error(['provide number of particles,'...
                'for example: plotRDF(''N'',625)']);
        else
            MC2DLJOrFileList = dir(['N' num2str(N) '*mat']);
            MC2DLJOrFileList = {MC2DLJOrFileList.name};
        end
    
        fileListOrgbyT = {};
        fileList = MC2DLJOrFileList;
        i = 1;
        rhoInd = 1;
        while ~isempty(fileList);
            data = matfile(fileList{1,1});
            sim = data.simulationParam;
            T(i) = sim.T;
            rho(rhoInd) = sim.rho;
            rhoInd = rhoInd + 1;
            fileListOrgbyT{i,1} = fileList{1,1};
            
            ind = 2;
            indnewlist = 1;
            newfileList = fileList;
            for j = 2:length(fileList)
                data = matfile(fileList{1,j});
                sim = data.simulationParam;
                thisT = sim.T;
                rho(rhoInd) = sim.rho;
                rhoInd = rhoInd + 1;
                
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
        
        % sort by Temperature
        [~, indsorted] = sort(T); 
        fileListOrgbyT(indsorted,:) = fileListOrgbyT(:,:);
        
        % sort by density
        [~, indsorted] = sort(unique(rho));
        fileListOrgbyT(:,indsorted) = fileListOrgbyT(:,:);
        
        [Niso, Nrho] = size(fileListOrgbyT);
        
        for i = 1:Niso
            for j = 1:Nrho
                if i == 1 && j == 1
                    MC2DLJs = MC2DLJoutput(fileListOrgbyT{1,1});
                else
                    MC2DLJs(i,j) = MC2DLJoutput(fileListOrgbyT{i,j}); 
                end
            end
        end
    
                    
                
    MC2DLJs(1,1) = MC2DLJs(1,1).calcRDF(10,300);
    numOfBins = length(MC2DLJs(1,1).bins);
    histo = zeros(Niso,Nrho,numOfBins);
    bins = zeros(Niso,Nrho,numOfBins);
    length2plot = zeros(Niso,Nrho);
    
    for i = 1:Niso
        for j = 1:Nrho
            j
            i
            MC2DLJs(i,j) = MC2DLJs(i,j).calcRDF(10,300);
            length2plot(i,j) = length(MC2DLJs(i,j).bins); 
            bins(i,j,1:length2plot(i,j)) = MC2DLJs(i,j).bins;
            histo(i,j,1:length2plot(i,j)) = MC2DLJs(i,j).histo;
            legrho{1,j} = ['\rho = '...
                num2str(MC2DLJs(i,j).simulationParam.rho)];
        end
        legT{1,i} = ['T = ' num2str(MC2DLJs(i,1).simulationParam.T)];
    end
    
    RDF.histo = histo;
    RDF.bins = bins;
    
    % plot for one temprature in the same figure,different densities
%    [~, b, c] = size(RDF.bins(1,:,:));
%     x = zeros(b,c);
%     y = zeros(b,c);
%     for i = 1:Niso
%        x(:,:) = RDF.bins(i,:,:);
%        y(:,:) = RDF.histo(i,:,:);
%        colorPlot(x,y,'addLegend',legrho,...
%            'lineStyle','-');
%        title(['T = ' num2str(MC2DLJs(i,1).simulationParam.T)]);
%        xlabel('distance, reduced units');
%        ylabel('g(r)');
%        saveas(gcf,['RDF_T' my_num2str(MC2DLJs(i,1).simulationParam.T)...
%             '.fig']);
%        saveas(gcf,['RDF_T' my_num2str(MC2DLJs(i,1).simulationParam.T)...
%             '.jpg']);
%        close all;
% 
%     end
%     
   % plot for one density in the same figure, different Temperatures
%     for j = 1:Nrho
%          colorPlot(RDF.bins(:,j,:),RDF.histo(:,j,:),'addLegend',legT,...
%              'lineStyle','-','length2plot',length2plot(:,j));
%         title(['\rho = ' num2str(MC2DLJs(1,j).simulationParam.rho)]);
%         xlabel('distance, reduced units');
%         ylabel('g(r)');
%         saveas(gcf,['RDF_rho' my_num2str(MC2DLJs(1,j).simulationParam.rho)...
%             '.fig']);
%         saveas(gcf,['RDF_rho' my_num2str(MC2DLJs(1,j).simulationParam.rho)...
%             '.jpg']);
%         close all;
%     end
    
    
end