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
              addOptional(p, 'MC2DLJs', []);
              addOptional(p, 'RDFdata', []);
              addOptional(p, 'numOfBins', 300);
              addOptional(p, 'N',[]);
              addOptional(p, 'dataFileListOrgbyT',[]);
              parse(p, varargin{:});
              Results = p.Results;
              dataFileList = Results.dataFileList;
              MC2DLJs = Results.MC2DLJs;
              RDFdata = Results.RDFdata;
              numOfBins = Results.numOfBins;
              N = Results.N;
              dataFileListOrgbyT = Results.dataFileListOrgbyT;
              
              if isempty(N)
                    error(['provide number of particles,'...
                            'for example: RDFoutput(''N'',625)']);
              end
              
              obj.N = N;

              
              if isempty(RDFdata) 
                    if isempty(dataFileListOrgbyT)
                        if isempty(dataFileList) 

                                % get all data from files in folder
                                dataFileList = dir(['N' num2str(N) '*mat']);
                                dataFileList = {dataFileList.name};
                        end

                        dataFileList = getDataFileList(dataFileList);
                    else
                        dataFileList = dataFileListOrgbyT;
                    end
                    obj.dataFileList = dataFileList;
                  
                    % build MC2DLJoutput objects
                    if isempty(MC2DLJs)
                        MC2DLJs = getMC2DLJs(obj.dataFileList);
                        obj.MC2DLJs = MC2DLJs;
                    else
                        obj.MC2DLJs = MC2DLJs;
                    end
                  
                    % create RDF data matfile
                    save(['RDF_N' num2str(N) '.mat']...
                      ,'MC2DLJs','dataFileList','-v7.3');
                    obj.data = matfile(['RDF_N' num2str(N) '.mat']);
                    obj.data = matfile(['RDF_N' num2str(N) '.mat'],...
                      'Writable',true);
                    clear MC2DLJs dataFileList;
                  
              else
                    obj.data = matfile(RDFdata,'Writable',true);
                    obj.MC2DLJs = obj.data.MC2DLJs;
                    obj.dataFileList = obj.data.dataFileList;
                    
                    if ~isempty(whos(obj.data,'length2plot'))
                        obj.data.length2plot = length2plot;
                    end
                    
                 
              end
              
              % get old RDF data if it exists
                  
              [Niso, Nrho] = size(obj.MC2DLJs);
              obj.numOfBins = numOfBins;
              
              if isempty(whos(obj.data,'histo'))
                    histo = zeros(Niso,Nrho,numOfBins);
                    bins = zeros(Niso,Nrho,numOfBins);
              else
                    histo = obj.data.histo;
                    bins = obj.data.bins;
              end
              
              for iso = 1:Niso
                    for rho = 1:Nrho
                        if sum(histo(iso,rho,:)) == 0
                            if ~isempty(obj.MC2DLJs(iso,rho).fileName)
                            		if ~isempty(whos(obj.MC2DLJs(iso,rho).data,'RDFhisto'))
                                		histo(iso,rho,:) = reshape(mean(obj.MC2DLJs(iso,rho).data.RDFhisto),...
                                            		[1,1,numOfBins]);
                               			bins(iso,rho,:) = reshape(...
                                            		obj.MC2DLJs(iso,rho).data.RDFbins,...
                                            		[1,1,numOfBins]);
                                    end
                            end
                        end
                    end
              end
              
              obj.data.histo = histo;
              obj.data.bins = bins;
              
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
			if ~isempty(obj.MC2DLJs(i,j).indIndata)	                        
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
				if ~isempty(obj.dataFileList{i,j})
                            nonZeroInd = find(obj.data.histo(i,j,:));
                            bins = obj.data.bins(i,j,1:numOfBins);
                            bins = bins(1,1,nonZeroInd);
                            histo = obj.data.histo(i,j,1:numOfBins);
                            histo = histo(1,1,nonZeroInd);
                            length2plotlog(j) = length(bins);
                            y(j,1:length2plotlog(j)) = -log(histo);
                            x(j,1:length2plotlog(j)) = bins;
                            length2plotlog(j) = length(bins);
                                    
                            % save log data
                            obj.data.logRDF(i,j,1:numOfBins)...
                             = reshape(y(j,:),[1,1,numOfBins]);
                            obj.data.length2plotlog(i,j) =...
                                length2plotlog(j); 
				end
                         end
                    else
                        x(:,:) = obj.data.bins(i,:,:);
                        y(:,:) = obj.data.histo(i,:,:);
                    end
                   
                   
                   
                   if plotLog
                        h = figure('Visible', Visible);
                            colorPlot(x,y,'addLegend',obj.legrho,...
                            'lineStyle','-','figHandle',h,...
                            'length2plot',length2plotlog);

                       ylabel('-log(g(r))');
                   else
                        h = figure('Visible', Visible);
                            colorPlot(x,y,'addLegend',obj.legrho,...
                            'lineStyle','-','figHandle',h);

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
                     if ~isempty(obj.data.bins(:,j,:))
                         h = figure('Visible',Visible);
                         if plotLog
                             % don't log zeros
                             x = zeros(Niso,obj.numOfBins);
                             y = zeros(Niso,obj.numOfBins);
                             length2plotlog = zeros(Niso);
                             for i = 1:Niso 

                                histo = obj.data.histo(i,j,1:obj.numOfBins);
                                bins = obj.data.bins(i,j,1:obj.numOfBins);
                                nonZeroInd = find(histo);
                                histo = histo(1,1,nonZeroInd);
                                bins = bins(1,1,nonZeroInd);
                                length2plotlog(i) = length(bins);
                                y(i,1:length2plotlog(i)) = -log(histo);
                                x(i,1:length2plotlog(i)) = bins;
                                length2plotlog(i) = length(bins);

                                % save log data
                                obj.data.logRDF(i,j,1:obj.numOfBins)...
                                 = reshape(y(i,:),[1,1,obj.numOfBins]);
                                obj.data.length2plotlog(i,j) =...
                                    length2plotlog(i);

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
                            if plotLog
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



        end
        
%         function obj = fitlogRDF(obj,varargin)
%             
%             [Niso, Nrho] = size(obj.MC2DLJs);
%             ftnm = fittype('4*(x^-n - x^m)');
%             ftm = fittype('4*(x^-12 - x^m)');
%             start = 1;
%             
%             for iso = 1:Niso
%                 for rho = 1:Nrho
%                     len2plot = obj.data.length2plotlog(iso,rho);
%                     
%                     x = reshape...
%                         (obj.data.bins(iso,rho,1:len2plot),[len2plot,1]);
%                     y = reshape...
%                         (obj.data.histo(iso,rho,1:len2plot),[len2plot,1]);
%                     
%                     
%                     success = false;
%                     
%                     while and(start < len2plot,~success)
%                         try
%                             fitsm{iso,rho} = fit(x,y,ftm);
%                             if confint(
%                             success = true;
%                         catch
%                             start = start + 10;
%                         end
%                     end
%                     
%                     if ~success
%                         fitsm{iso,rho} = [];
%                     end
%                     
%                 end
%                 obj.data.fitsnm = fitsnm;
%                 obj.data.fitsm = fitsm;
%             end
%         end
%         
% %        function obj = plotGoodFits(obj,
%         
%             
%         
     end
     
 
 end
            
    function MC2DLJs = getMC2DLJs(fileListOrgbyT)
        [Niso, Nrho] = size(fileListOrgbyT);
        
        for i = 1:Niso
            for j = 1:Nrho
                if i == 1 && j == 1
                    MC2DLJs = MC2DLJoutput(fileListOrgbyT{1,1});
                else
			if exist(fileListOrgbyT{i,j},'file') == 2
                    		MC2DLJs(i,j) = MC2DLJoutput(fileListOrgbyT{i,j});
			end
                end
            end
        end
    end
    
%     function start = find_best_start_for_fit(x,y,ft,minconf)
%             
%             % check input: 
%             if length(minconf) ~= length(coeffnames(ft))
%                 error(['minconf must be a vector with the values'...
%                     ' of minconf for each coefficient in ft']);
%             end
%             
%             success = false;
%             
%             while and(start < len2plot, ~success)
%                     try
%                         fitsnm{iso,rho} = fit(x,y,ft);
%                         if 
%                         success = true;
%                     catch
%                         start = start + 10;
%                     end
%             end
%                     
%             if ~success
%                     fitsnm{iso,rho} = [];
%             end
%  
%     end
% 
    
