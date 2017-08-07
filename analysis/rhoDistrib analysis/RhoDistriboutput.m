classdef RhoDistriboutput

    properties
        data, dataFileList, rhoDistribParam,numOfBins, MC2DLJs, length2plot,...
            legrho, legT;
    end
    
    methods
        
        %constructor
        function obj = RhoDistriboutput(varargin)
            
              p = inputParser();
              addOptional(p, 'dataFileList', []);
              addOptional(p, 'fileListOrg', []);
              addOptional(p, 'N', '');              
              addOptional(p, 'rhoDistribdata', []);
              addOptional(p, 'dataFileNameEnd', '');
              parse(p, varargin{:});
              Results = p.Results;
              dataFileList = Results.dataFileList;
              fileListOrg = Results.fileListOrg;
              N = Results.N;
              rhoDistribdata = Results.rhoDistribdata;
              dataFileNameEnd = Results.dataFileNameEnd;
              
              if isempty(rhoDistribdata)
                    
                  if isempty(dataFileList) && isempty(fileListOrg)
                      % get all data from files in folder
                      dataFileList = dir(['N' num2str(N) 'T*.mat']);
                      dataFileList = {dataFileList.name};
                      dataFileList = getDataFileList(dataFileList);
                      
                  else
                      if ~isempty(fileListOrg)
                          dataFileList = fileListOrg;
                      else
                          dataFileList = getDataFileList(dataFileList); 
                      end
                  end
                  
                  
                  obj.dataFileList = dataFileList;
                  
                  % build MC2DLJoutput objects
                  MC2DLJs = getMC2DLJs(obj.dataFileList);
                  obj.MC2DLJs = MC2DLJs;
                  
                  % create rhoDistrib data matfile
                  save(['rhoDistribN' num2str(N) dataFileNameEnd '.mat']...
                      ,'MC2DLJs','dataFileList','-v7.3');
                  obj.data = matfile(['rhoDistribN' num2str(N) dataFileNameEnd '.mat']);
                  obj.data = matfile(['rhoDistribN' num2str(N) dataFileNameEnd...
                      '.mat'],'Writable',true);
                  clear MC2DLJs dataFileList;
                  
              else
                 if ~isempty(rhoDistribdata)
                        obj.data = matfile(rhoDistribdata);
                        obj.data = matfile(rhoDistribdata,...
                            'Writable',true);
                        obj.MC2DLJs = obj.data.MC2DLJs;
                        obj.numOfBins = size(obj.data.bins);
                        MC2DLJs1 = obj.MC2DLJs(1,1);
                        obj.rhoDistribParam = MC2DLJs1.data.rhoDistribParam;
                        clear MC2DLJs1;
                        obj.dataFileList = obj.data.dataFileList;
                        obj.length2plot = obj.data.length2plot;
                        obj.legT = obj.data.legT;
                        obj.legrho = obj.data.legrho;
                 end
                 
              end
        end
        
        function obj = calcAllrhoDistrib(obj,squares,numOfBins,varargin)
                
                p = inputParser();
                addOptional(p, 'talk', false);
                addOptional(p, 'startFrom', 1);
                parse(p, varargin{:});
                Results = p.Results;
                talk = Results.talk;
                startFrom = Results.startFrom;
                
                [Niso, Nrho] = size(obj.MC2DLJs);
                [obj.MC2DLJs(1,1), rhoNorm, PL] =...
                    obj.MC2DLJs(1,1).calcRhoDistrib(squares,...
                    numOfBins,'startFrom',startFrom);
                obj.data.histo = zeros(Niso,Nrho,numOfBins);
                obj.data.bins = zeros(Niso,Nrho,numOfBins);
                obj.data.histo(1,1,1:numOfBins) =...
                    reshape(PL,[1,1,numOfBins]);
                obj.data.bins(1,1,1:numOfBins) =...
                    reshape(rhoNorm,[1,1,numOfBins]);
                obj.data.squares = squares;
                obj.data.numOfBins = numOfBins;
                
                for i = 1:Niso
                    for j = 1:Nrho
                        if talk
                            j
                            i
                        end
                        if ~and(i == 1,j==1)
                            [obj.MC2DLJs(i,j), rhoNorm, PL] =...
                                obj.MC2DLJs(i,j).calcRhoDistrib(squares,...
                                    numOfBins,'startFrom',startFrom); 
                            obj.data.bins(i,j,1:numOfBins) =...
                                reshape(rhoNorm, [1,1,numOfBins]);
                            obj.data.histo(i,j,1:numOfBins) =...
                                reshape(PL, [1,1,numOfBins]);
                        end
                        obj.legrho{1,j} = ['\rho = '...
                            num2str(obj.MC2DLJs(i,j).simulationParam.rho)];
                    end
                    obj.legT{1,i} = ['T = '...
                        num2str(obj.MC2DLJs(i,1).simulationParam.T)];
                end
                
                obj.data.legT = obj.legT;
                obj.data.legrho = obj.legrho;

        end
        
        function obj = plotrhoDistribT(obj,varargin)
        % plot for one temprature in the same figure,different densities
           
           p = inputParser();
           addOptional(p, 'saveFig', true);
           addOptional(p, 'keepFigOpen', true);
           addOptional(p, 'Visible', 'on');
           addOptional(p, 'addStr2title', '');
           addOptional(p, 'addFileNameEnd', '');
           parse(p, varargin{:});
           Results = p.Results;
           saveFig = Results.saveFig;
           keepFigOpen = Results.keepFigOpen;
           Visible = Results.Visible;
           addStr2title = Results.addStr2title;
           addFileNameEnd = Results.addFileNameEnd; 
           
           [Niso, ~] = size(obj.MC2DLJs);
           
           [~, b, c] = size(obj.data.bins(1,:,:));
            x = zeros(b,c);
            y = zeros(b,c);
            for i = 1:Niso
                   x(:,:) = obj.data.bins(i,:,:);
                   y(:,:) = obj.data.histo(i,:,:);
                   h = figure('Visible',Visible); 
                   colorPlot(x,y,'addLegend',obj.legrho,...
                       'lineStyle','-','figHandle',h);
                   title(['T = '...
                       num2str(obj.MC2DLJs(i,1).simulationParam.T)...
                        addStr2title]);
                   xlabel('distance, reduced units');
                   ylabel('PL');

                   if saveFig
                       saveas(gcf,['rhoDistrib_N'...
                           num2str(obj.MC2DLJs(i,1).simulationParam.N) 'T' ...
                           my_num2str(obj.MC2DLJs(i,1).simulationParam.T)...
                             addFileNameEnd '.fig']);
                       saveas(gcf,['rhoDistrib_N'...
                           num2str(obj.MC2DLJs(i,1).simulationParam.N) 'T' ...
                           my_num2str(obj.MC2DLJs(i,1).simulationParam.T)...
                             addFileNameEnd '.jpg']);
                   end

                   if ~keepFigOpen
                        close all;
                   end

            end


        end
        
        function obj = plotrhoDistribrho(obj,varargin)
        % plot for one density in the same figure, different Temperatures

               p = inputParser();
               addOptional(p, 'saveFig', true);
               addOptional(p, 'keepFigOpen', true);
               addOptional(p, 'Visible','on');
               addOptional(p, 'addStr2title', '');
               addOptional(p, 'addFileNameEnd', '');
               parse(p, varargin{:});
               Results = p.Results;
               saveFig = Results.saveFig;
               keepFigOpen = Results.keepFigOpen;
               Visible = Results.Visible; 
               addStr2title = Results.addStr2title;
               addFileNameEnd = Results.addFileNameEnd;
               
               [~, Nrho] = size(obj.MC2DLJs);
               
                for j = 1:Nrho
                     h = figure('Visible',Visible);
                     colorPlot(obj.data.bins(:,j,:)...
                         ,obj.data.histo(:,j,:),'addLegend',obj.legT,...
                         'lineStyle','-',...
                         'figHandle',h);
                    title(['\rho = '...
                        num2str(obj.MC2DLJs(1,j).simulationParam.rho)...
                         addStr2title]);
                    xlabel('distance, reduced units');
                    ylabel('g(r)');
                    
                    if saveFig
                        saveas(gcf,['RDF_N'...
                             num2str(obj.MC2DLJs(1,j).simulationParam.N) '_rho'...
                            my_num2str(obj.MC2DLJs(1,j).simulationParam.rho)...
                             addFileNameEnd '.fig']);
                        saveas(gcf,['RDF_N'...
                            num2str(obj.MC2DLJs(1,j).simulationParam.N) 'rho'...
                            my_num2str(obj.MC2DLJs(1,j).simulationParam.rho)...
                             addFileNameEnd '.jpg']);
                    end
                    
                    if ~keepFigOpen
                        close all;
                    end
                end



        end
            
    end
    

end
            
