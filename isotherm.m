classdef isotherm
    
    properties
        pressure, rho, T;
        MC2DLJ, datafileList;
        simulationParam;
        cutEquilirization, firstSteps2ignore;
    end
    
    methods
        
        % constructor
        function obj = isotherm(varargin)
                    
                    p = inputParser();
                    addOptional(p, 'verelet', []);
                    addOptional(p, 'N', []);
                    addOptional(p, 'T', []);
                    addOptional(p, 'rho', []);
                    addOptional(p, 'initialmaxdr', []);
                    addOptional(p, 'initialConfig', []);
                    addOptional(p, 'rCutoff', []);
                    addOptional(p, 'r', []);
                    addOptional(p, 'fileNameInit', '');
                    addOptional(p, 'cutEquilirization', true);
                    addOptional(p, 'firstSteps2ignore', 50);
                    addOptional(p, 'datafileList', []);
                    parse(p, varargin{:});
                    Results = p.Results;
                    rl = Results.verelet;
                    N = Results.N;
                    T = Results.T;
                    rho = Results.rho;
                    initialmaxdr = Results.initialmaxdr;
                    initialConfig = Results.initialConfig;
                    rCutoff = Results.rCutoff;
                    r = Results.r;
                    rl = Results.verelet;
                    fileNameInit = Results.fileNameInit;
                    cutEquilirization = Results.cutEquilirization;
                    datafileList = Results.datafileList;
                    firstSteps2ignore  = Results.firstSteps2ignore;
                    
                    obj.cutEquilirization = cutEquilirization;

            
            if ~isempty(datafileList)
                
                    % create object from a given file list
                    obj.datafileList = datafileList;
                    data = matfile(obj.datafileList{1,1});
                    sim = data.simulationParam;
                    obj.T = sim.T;
                    obj.simulationParam = data.simulationParam;

                    for i = 1:length(obj.datafileList)
                            if i == 1
                                obj.MC2DLJ = ...
                                    MC2DLJoutput(obj.datafileList{i,1});
                            else
                                obj.MC2DLJ(i) =...
                                    MC2DLJoutput(obj.datafileList{i,1});
                            end
                            
                            if cutEquilirization
                                
                               % check if meanPlrcEq needs to be calculated
                               % calculate it if it was nit calculated yet
                               if isempty(whos(obj.MC2DLJ(i).data,...
                                       'meanPlrcEq'))
                                   obj.MC2DLJ(i) =...
                                       obj.MC2DLJ(i).calcMeanWithoutFirstSteps(...
                                       firstSteps2ignore);
                               end
                                   
                               obj.pressure(i) =...
                                   obj.MC2DLJ(i).data.meanPlrcEq;
                            else
                               [~,s] = size(obj.MC2DLJ(i).data.meanPlrc);
                               obj.pressure(i) =...
                                   obj.MC2DLJ(i).data.meanPlrc(1,s);
                            end
                            obj.rho(i) = obj.MC2DLJ(i).simulationParam.rho;
                    end
                
                
            else
          
                    
                    obj.simulationParam.N = N;
                    obj.T = T;
                    obj.rho = rho;
                    obj.simulationParam.L = sqrt(N./rho);
                    obj.simulationParam.initialmaxdr = initialmaxdr;
                    obj.simulationParam.initialConfig = initialConfig;
                    obj.simulationParam.rCutoff = rCutoff;
                    obj.simulationParam.r = r;
                    obj.simulationParam.rl = rl;
                    
                    for i = 1:length(rho)
                        if i == 1
                            obj.MC2DLJ = MC2DLJoutput(N,T,rho(i),initialmaxdr,...
                                initialConfig,rCutoff,r,...
                                'verelet',rl,'pressure',true,...
                                'fileNameInit',fileNameInit);
                        else
                            obj.MC2DLJ(i) = MC2DLJoutput(N,T,rho(i),initialmaxdr,...
                                initialConfig,rCutoff,r,...
                                'verelet',rl,'pressure',true,...
                                'fileNameInit',fileNameInit);
                        end
                        obj.datafileList{1,i} = obj.MC2DLJ(i).fileName;
                    end
            end    
        end
        
        function obj = calcIso(obj,Nsteps,saveEvery)
                 
                 M = obj.MC2DLJ;
                 parfor i = 1:length(obj.rho)
                           M(i) = M(i).MonteCarlo(Nsteps,saveEvery);
                           M(i) = M(i).meanProp();
                 end
                    
                 obj.MC2DLJ = M;
                 clear M;
                 
                 for i = 1:length(obj.rho)
                        
                        if obj.cutEquilirization
                                % check if meanPlrcEq needs to be calculated
                                % calculate it if it was nit calculated yet
                            if isempty...
                                    (whos(obj.MC2DLJ(i).data,'meanPlrcEq'))
                                obj.MC2DLJ(i) =...
                                obj.MC2DLJ(i).calcMeanWithoutFirstSteps(...
                                               obj.firstSteps2ignore);
                            end

                            obj.pressure(i) =...
                                       obj.MC2DLJ(i).data.meanPlrcEq;
                        else
                            [~,s] = size(obj.MC2DLJ(i).data.meanPlrc);
                                obj.pressure(i) =...
                                       obj.MC2DLJ(i).data.meanPlrc(1,s);
                        end
                            obj.rho(i) = obj.MC2DLJ(i).simulationParam.rho;
                            [~,s] = size(obj.MC2DLJ(i).data.meanPlrc);
                            obj.pressure(i) =...
                                obj.MC2DLJ(i).data.meanPlrc(1,s);
                        
                 end
             
        end
        
        function [obj, figHandle] = plotPropVsStep(obj,prop,varargin)
            
            p = inputParser();
            addOptional(p, 'figHandle', []);
            addOptional(p, 'normalizeByMean', false);
            parse(p, varargin{:});
            Results = p.Results;
            figHandle = Results.figHandle;
            normalizeByMean = Resuslts.normalizeByMean;
                    
            if isempty(figHandle)
                figHandle = figure();
            end
                    
            indIndata = obj.MC2DLJ(1).indIndata;
            rhoN = length(obj.rho);
            x = zeros(rhoN,indIndata);
            y = x;
            switch prop
                case 'P'
                    for i = 1:rhoN
                        y(i,:) = obj.MC2DLJ(i).data.allPlrc;
                        
                        if normalizeByMean
                            y(i,:) = y(i,:)/mean(y(i,:));
                        end
                        
                        x(i,:) = 1:indIndata;
                        leg{1,i} = ['rho = ' num2str(obj.rho(i))];
                    end
                    
                case 'U'
                    for i = 1:rhoN
                        y(i,:) = obj.MC2DLJ(i).data.allUlrc;
                        
                        if normalizeByMean
                            y(i,:) = y(i,:)/mean(y(i,:));
                        end
                        
                        x(i,:) = 1:indIndata;
                        leg{1,i} = ['rho = ' num2str(obj.rho(i))];
                    end
            end
            
            figHandle = colorPlot(x,y,'lineStyle','-','addLegend',leg,...
                'figHandle',figHandle);
            title([prop ' Vs. steps for T = ' num2str(obj.T)]);
            xlabel('step');
            ylabel(prop);
        end
    end
end