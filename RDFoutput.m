classdef RDFoutput

    properties
        data, dataFileList, maxDist,numOfBins, MC2DLJs;
        length2plot, legrho, legT;
    end
    
    methods
        
        %constructor
        function obj = RDFoutput(varargin)
              p = inputParser();
              addOptional(p, 'dataFileList', []);
              addOptional(p, 'RDFdata', []);
              parse(p, varargin{:});
              Results = p.Results;
              dataFileList = Results.dataFileList;
              RDFdata = Results.RDFdata;
              
              if isempty(RDFdata) && isempty(dataFileList)
                  
                  % get all data from files in folder
                  fileList = dir('N*mat');
                  fileList = {fileList.name};
                  
                  dataFileList = getDataFileList(fileList);
                  obj.dataFileList = dataFileList;
                  
                  % build MC2DLJoutput objects
                  MC2DLJs = getMC2DLJs(obj.dataFileList);
                  obj.MC2DLJs = MC2DLJs;
                  
                  % create RDF data matfile
                  save('RDF.mat','MC2DLJs','dataFileList','-v7.3');
                  obj.data = matfile('RDF.mat');
                  obj.data = matfile('RDF.mat','Writable',true);
                  clear MC2DLJs dataFileList;
                  
              else
                 % add constructor - if data file exists 
                 
              end
        end
        
        function obj = calcAllRDF(obj,maxDist,numOfBins,varargin)
                
                p = inputParser();
                addOptional(p, 'talk', false);
                parse(p, varargin{:});
                Results = p.Results;
                talk = Results.talk;

                
                [Niso, Nrho] = size(obj.MC2DLJs);
                obj.MC2DLJs(1,1) =...
                    obj.MC2DLJs(1,1).calcRDF(maxDist,numOfBins);
                obj.data.histo = zeros(Niso,Nrho,numOfBins);
                obj.data.bins = zeros(Niso,Nrho,numOfBins);
                obj.length2plot = zeros(Niso,Nrho);

                for i = 1:Niso
                    for j = 1:Nrho
                        if talk
                            j
                            i
                        end
                        obj.MC2DLJs(i,j) =...
                            obj.MC2DLJs(i,j).calcRDF(maxDist,numOfBins);
                        obj.length2plot(i,j) = length(obj.MC2DLJs(i,j).RDFbins); 
                        obj.data.bins(i,j,1:obj.length2plot(i,j)) =...
                            obj.MC2DLJs(i,j).RDFbins;
                        obj.data.histo(i,j,1:obj.length2plot(i,j)) =...
                            obj.MC2DLJs(i,j).RDFhisto;
                        obj.legrho{1,j} = ['\rho = '...
                            num2str(obj.MC2DLJs(i,j).simulationParam.rho)];
                    end
                    obj.legT{1,i} = ['T = '...
                        num2str(obj.MC2DLJs(i,1).simulationParam.T)];
                end

        end
        
        function obj = plotRDFT(obj,varargin)
        % plot for one temprature in the same figure,different densities
           
           p = inputParser();
           addOptional(p, 'saveFig', true);
           addOptional(p, 'keepFigOpen', true);
           parse(p, varargin{:});
           Results = p.Results;
           saveFig = Results.saveFig;
           keepFigOpen = Results.keepFigOpen;

           [Niso, ~] = size(obj.MC2DLJs);
           
           [~, b, c] = size(obj.data.bins(1,:,:));
            x = zeros(b,c);
            y = zeros(b,c);
            for i = 1:Niso
                   x(:,:) = obj.data.bins(i,:,:);
                   y(:,:) = obj.data.histo(i,:,:);
                   colorPlot(x,y,'addLegend',obj.legrho,...
                       'lineStyle','-');
                   title(['T = ' num2str(obj.MC2DLJs(i,1).simulationParam.T)]);
                   xlabel('distance, reduced units');
                   ylabel('g(r)');

                   if saveFig
                       saveas(gcf,['RDF_T'...
                           my_num2str(obj.MC2DLJs(i,1).simulationParam.T)...
                            '.fig']);
                       saveas(gcf,['RDF_T'...
                           my_num2str(obj.MC2DLJs(i,1).simulationParam.T)...
                            '.jpg']);
                   end

                   if ~keepFigOpen
                        close all;
                   end

            end


        end
        
        function obj = plotRDFrho(obj,varargin)
        % plot for one density in the same figure, different Temperatures

               p = inputParser();
               addOptional(p, 'saveFig', true);
               addOptional(p, 'keepFigOpen', true);
               parse(p, varargin{:});
               Results = p.Results;
               saveFig = Results.saveFig;
               keepFigOpen = Results.keepFigOpen;

               [~, Nrho] = size(obj.MC2DLJs);
               
                for j = 1:Nrho
                     colorPlot(obj.data.bins(:,j,:)...
                         ,obj.data.histo(:,j,:),'addLegend',obj.legT,...
                         'lineStyle','-',...
                         'length2plot',obj.length2plot(:,j));
                    title(['\rho = '...
                        num2str(obj.MC2DLJs(1,j).simulationParam.rho)]);
                    xlabel('distance, reduced units');
                    ylabel('g(r)');
                    
                    if saveFig
                        saveas(gcf,['RDF_rho'...
                            my_num2str(obj.MC2DLJs(1,j).simulationParam.rho)...
                            '.fig']);
                        saveas(gcf,['RDF_rho'...
                            my_num2str(obj.MC2DLJs(1,j).simulationParam.rho)...
                            '.jpg']);
                    end
                    
                    if ~keepFigOpen
                        close all;
                    end
                end



        end
            
    end
    

end

    function fileListOrgbyT = getDataFileList(inputFileList)
        fileListOrgbyT = {};
        fileList = inputFileList;
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
        
    end
            
    function MC2DLJs = getMC2DLJs(fileListOrgbyT)
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
    end