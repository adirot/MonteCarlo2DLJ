classdef RDFoutput

    properties
        data, dataFileList, maxDist,numOfBins, MC2DLJs, N;
        length2plot, legrho, legT;
        
    end
    
    methods
        
        %constructor
        function obj = RDFoutput(varargin)
              p = inputParser();
              addOptional(p, 'dataFileList', []);
              addOptional(p, 'RDFdata', []);
              addOptional(p, 'N',[]);
              parse(p, varargin{:});
              Results = p.Results;
              dataFileList = Results.dataFileList;
              RDFdata = Results.RDFdata;
              N = Results.N;
              
              if isempty(N)
                    error(['provide number of particles,'...
                            'for example: RDFoutput(''N'',625)']);
              end
              
              obj.N = N;
              
              if isempty(RDFdata) 
                  if isempty(dataFileList)
                  
                      % get all data from files in folder
                      dataFileList = dir(['N' num2str(N) '*mat']);
                      dataFileList = {dataFileList.name};
                  end
                  
                  dataFileList = getDataFileList(dataFileList);
                  obj.dataFileList = dataFileList;
                  
                  % build MC2DLJoutput objects
                  MC2DLJs = getMC2DLJs(obj.dataFileList);
                  obj.MC2DLJs = MC2DLJs;
                  
                  % create RDF data matfile
                  save(['RDF_N' num2str(N) '.mat']...
                      ,'MC2DLJs','dataFileList','-v7.3');
                  obj.data = matfile(['RDF_N' num2str(N) '.mat']);
                  obj.data = matfile(['RDF_N' num2str(N) '.mat'],...
                      'Writable',true);
                  clear MC2DLJs dataFileList;
                  
                  % get old RDF data if it exists
                  
                        [Niso, Nrho] = size(obj.MC2DLJs);
                        [~,numOfBins] = size(obj.MC2DLJs(1,1).data.RDFbins);
                        obj.numOfBins = numOfBins;
                        obj.data.histo = zeros(Niso,Nrho,numOfBins);
                        obj.data.bins = zeros(Niso,Nrho,numOfBins);
                        for iso = 1:Niso
                            for rho = 1:Nrho
                                if ~isempty(whos(obj.MC2DLJs(iso,rho).data,'RDFhisto'))
                                    obj.data.histo(iso,rho,:) = reshape(...
                                        mean(obj.MC2DLJs(iso,rho).data.RDFhisto),...
                                        [1,1,numOfBins]);
                                    obj.data.bins(iso,rho,:) = reshape(...
                                        obj.MC2DLJs(iso,rho).data.RDFbins,...
                                        [1,1,numOfBins]);
                                end
                            end
                        end
                  

                  
              else
                 % add constructor - if data file exists 
                 
              end
        end
        
        function obj = calcAllRDF(obj,maxDist,numOfBins,varargin)
                
                p = inputParser();
                addOptional(p, 'talk', false);
                addOptional(p, 'skipExisting', false);
                parse(p, varargin{:});
                Results = p.Results;
                talk = Results.talk;
                skipExisting = Results.skipExisting;
                
                [Niso, Nrho] = size(obj.MC2DLJs);
                obj.data.histo = zeros(Niso,Nrho,numOfBins);
                obj.data.bins = zeros(Niso,Nrho,numOfBins);
                obj.length2plot = zeros(Niso,Nrho);
                obj.numOfBins = numOfBins;
                
                        

                for i = 1:Niso
                    for j = 1:Nrho
                        if talk
                            j
                            i
                        end
                        
                        if skipExisting
                            if isempty(whos(obj.MC2DLJs(i,j).data,'RDFhisto'))
                                obj.MC2DLJs(i,j) =...
                                obj.MC2DLJs(i,j).calcRDF(maxDist,numOfBins);
                                
                            end
                        else
                           obj.MC2DLJs(i,j) =...
                                obj.MC2DLJs(i,j).calcRDF(maxDist,numOfBins); 
                        end
                        
                        obj.length2plot(i,j) =...
                            length(obj.MC2DLJs(i,j).data.RDFbins); 
 
                        a = obj.MC2DLJs(i,j).data.RDFbins(1,...
                            1:obj.length2plot(i,j));
                        obj.data.bins(i,j,1:obj.length2plot(i,j)) =...
                            reshape(a,[1,1,obj.length2plot(i,j)]);
                        
                        a = obj.MC2DLJs(i,j).data.RDFhisto(1,...
                            1:obj.length2plot(i,j));
                        obj.data.histo(i,j,1:obj.length2plot(i,j)) =...
                            reshape(a,[1,1,obj.length2plot(i,j)]);
                        
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
           addOptional(p, 'Visible', 'on');
           addOptional(p, 'plotLog', false);
           parse(p, varargin{:});
           Results = p.Results;
           saveFig = Results.saveFig;
           keepFigOpen = Results.keepFigOpen;
           Visible = Results.Visible;
           plotLog = Results.plotLog;
           

           [Niso, Nrho] = size(obj.MC2DLJs);
           
           [~, b, numOfBins] = size(obj.data.bins(1,:,:));
            x = zeros(b,numOfBins);
            y = zeros(b,numOfBins);
            
            if plotLog
                   obj.data.logRDF = zeros(Niso,Nrho,numOfBins);
                   obj.data.length2plotlog = zeros(Niso,Nrho);
            end
            
            for i = 1:Niso
                    if plotLog
                        
                        % don't log zeros
                         length2plotlog = zeros(Nrho);
                         for j = 1:Nrho
                            nonZeroInd = find(obj.data.histo(i,j,:));
                            length2plotlog(j) =...
                                length(obj.data.bins(i,j,nonZeroInd));
                            y(j,1:length2plotlog(j)) =...
                                -log(obj.data.histo(i,j,nonZeroInd));
                            x(j,1:length2plotlog(j)) =...
                                obj.data.bins(i,j,nonZeroInd);
                            length2plotlog(j) =...
                                length(obj.data.bins(i,j,nonZeroInd));
                            
                                                        
                            % save log data
                            obj.data.logRDF(i,j,1:numOfBins)...
                             = reshape(y(j,:),[1,1,numOfBins]);
                            obj.data.length2plotlog(i,j) =...
                                length2plotlog(j); 

                         end
                    else
                        x(:,:) = obj.data.bins(i,:,:);
                        y(:,:) = obj.data.histo(i,:,:);
                    end
                   
                   
                   
                   h = figure('Visible', Visible);
                   colorPlot(x,y,'addLegend',obj.legrho,...
                       'lineStyle','-','figHandle',h,...
                       'length2plot',length2plotlog);
                   if plotLog
                       ylabel('-log(g(r))');
                   else
                       ylabel('g(r)');
                   end
                   
                   title(['T = ' num2str(obj.MC2DLJs(i,1).simulationParam.T)]);
                   xlabel('distance, reduced units');
                   

                   if saveFig
                       if plotLog
                           saveas(gcf,['logRDF_T'...
                           my_num2str(obj.MC2DLJs(i,1).simulationParam.T)...
                            'N' num2str(obj.N)  '.fig']);
                           saveas(gcf,['logRDF_T'...
                               my_num2str(obj.MC2DLJs(i,1).simulationParam.T)...
                                'N' num2str(obj.N)  '.jpg']);
                       else
                           saveas(gcf,['RDF_T'...
                               my_num2str(obj.MC2DLJs(i,1).simulationParam.T)...
                                'N' num2str(obj.N)  '.fig']);
                           saveas(gcf,['RDF_T'...
                               my_num2str(obj.MC2DLJs(i,1).simulationParam.T)...
                                'N' num2str(obj.N)  '.jpg']);
                       end
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
               addOptional(p, 'Visible' , 'on');
               addOptional(p, 'plotLog' , false);
               parse(p, varargin{:});
               Results = p.Results;
               saveFig = Results.saveFig;
               keepFigOpen = Results.keepFigOpen;
               Visible = Results.Visible; 
               plotLog = Results.plotLog;
               
               [Niso, Nrho] = size(obj.MC2DLJs);
               if plotLog
                   obj.data.logRDF = zeros(Niso,Nrho,obj.numOfBins);
                   obj.data.length2plotlog = zeros(Niso,Nrho);
               end
               
                for j = 1:Nrho
                     
                     h = figure('Visible',Visible);
                     if plotLog
                         % don't log zeros
                         x = zeros(Niso,obj.numOfBins);
                         y = zeros(Niso,obj.numOfBins);
                         length2plotlog = zeros(Niso);
                         for i = 1:Niso 
                            nonZeroInd = find(obj.data.histo(i,j,:));
                            length2plotlog(i) =...
                                length(obj.data.bins(i,j,nonZeroInd));
                            y(i,1:length2plotlog(i)) =...
                                -log(obj.data.histo(i,j,nonZeroInd));
                            x(i,1:length2plotlog(i)) =...
                                obj.data.bins(i,j,nonZeroInd);
                            length2plotlog(i) =...
                                length(obj.data.bins(i,j,nonZeroInd));
                            
                            % save log data
                            obj.data.logRDF(i,j,1:obj.numOfBins)...
                             = reshape(y(i,:),[1,1,obj.numOfBins]);
                            obj.data.length2plotlog(i,j) = length2plotlog(i); 
                         end
                         
                         colorPlot(x,y,'addLegend',obj.legT,...
                         'lineStyle','-',...
                         'length2plot',length2plotlog...
                         ,'figHandle',h);
                         ylabel('-log(g(r))');
                         
                         
                     else
                         colorPlot(obj.data.bins(:,j,:)...
                             ,obj.data.histo(:,j,:),'addLegend',obj.legT,...
                             'lineStyle','-',...
                             'length2plot',obj.length2plot(:,j)...
                             ,'figHandle',h);
                         ylabel('g(r)');
                     end
                    
                    title(['\rho = '...
                        num2str(obj.MC2DLJs(1,j).simulationParam.rho)]);
                    xlabel('distance, reduced units');
                    ylabel('g(r)');
                    
                    if saveFig
                        if log
                            logstr = 'log';
                        else
                            logstr = '';
                        end
                        
                        saveas(gcf,[logstr 'RDF_rho'...
                            my_num2str(obj.MC2DLJs(1,j).simulationParam.rho)...
                             'N' num2str(obj.N) '.fig']);
                        saveas(gcf,[logstr 'RDF_rho'...
                            my_num2str(obj.MC2DLJs(1,j).simulationParam.rho)...
                             'N' num2str(obj.N) '.jpg']);
                    end
                    
                    if ~keepFigOpen
                        close all;
                    end
                end



        end
        
        function obj = fitlogRDF(obj,varargin)
            
            [Niso, Nrho] = size(obj.MC2DLJs);
            ftnm = fittype('4*(x^-n - x^m)');
            ftm = fittype('4*(x^-12 - x^m)');
            
            for iso = 1:Niso
                for rho = 1:Nrho
                    fitsnm{iso,rho} =...
                        fit(obj.data.bins(iso,rho,...
                            1:obj.data.length2plotlog(iso,rho)),...
                        obj.data.logRDF(iso,rho,...
                            1:obj.data.length2plotlog(iso,rho)),ftnm);
                    fitsm{iso,rho} =...
                        fit(obj.data.bins(iso,rho,...
                            1:obj.data.length2plotlog(iso,rho)),...
                        obj.data.logRDF(iso,rho,...
                            1:obj.data.length2plotlog(iso,rho)),ftm);
                end
                obj.data.fitsnm = fitsnm;
                obj.data.fitsm = fitsm;
            end
        end
        
            
        
    end
    

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