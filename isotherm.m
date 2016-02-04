classdef isotherm
    
    properties
        pressure, rho, T;
        MC2DLJ, datafileList;
        cv;
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
                    addOptional(p, 'm', 6);
                    addOptional(p, 'fileNameInit', '');
                    addOptional(p, 'cutEquilirization', false);
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
                    m = Results.m;
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
                               % calculate it if it was not calculated yet
                               if isempty(whos(obj.MC2DLJ(i).data,...
                                       'meanPlrcEq'))
                                   
                                   obj.MC2DLJ(i) =...
                                       obj.MC2DLJ(i).calcMeanWithoutFirstSteps(...
                                       firstSteps2ignore);
                               end
                                   
                               obj.pressure(i) =...
                                   obj.MC2DLJ(i).data.meanPlrcEq;
                            else
                                if existInMatfile(obj.MC2DLJ(i).data,...
                                        'allPlrc')
                                    if ~existInMatfile(obj.MC2DLJ(i).data,...
                                            'meanPlrc')
                              
                                        obj.MC2DLJ(i).meanProp('P');
                                    end
                                   
                                else
                                    obj.MC2DLJ(i).calcPressureAfterRun();
                                    obj.MC2DLJ(i).meanProp('P');
                                end
                                
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
                    obj.simulationParam.m = m;
                    
                    for i = 1:length(rho)
                        if i == 1
                            if strcmp(initialConfig,'auto')
                                if rho(i) > 0.4
                                    init = 'hex';
                                else
                                    init = 'random';
                                end
                            else
                                init = initialConfig;
                            end
                            obj.MC2DLJ = MC2DLJoutput(N,T,rho(i),initialmaxdr,...
                                init,rCutoff,r,...
                                'verelet',rl,'pressure',true,...
                                'fileNameInit',fileNameInit,'m',m);
                        else
                            if strcmp(initialConfig,'auto')
                                if rho(i) > 0.4
                                    init = 'hex';
                                else
                                    init = 'random';
                                end
                            else
                                init = initialConfig;
                            end
                            obj.MC2DLJ(i) = MC2DLJoutput(N,T,rho(i),initialmaxdr,...
                                init,rCutoff,r,...
                                'verelet',rl,'pressure',true,...
                                'fileNameInit',fileNameInit,'m',m);
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
            normalizeByMean = Results.normalizeByMean;
                    
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
        
        function obj = getSnapShots(obj,step,varargin)
           p = inputParser();
           addOptional(p, 'saveFig', true);
           addOptional(p, 'keepFigOpen', true);
           addOptional(p, 'rhos', 'all');
           addOptional(p, 'rhosInd', []);           
           parse(p, varargin{:});
           Results = p.Results;
           saveFig = Results.saveFig;
           keepFigOpen = Results.keepFigOpen;
           rhos = Results.rhos;
           rhosInd = Results.rhosInd;
           
           if ischar(rhos)
               if strcmp(rhos,'all')
                   rhos = obj.rho;
               end
           end
           
           if isempty(rhosInd)
               for i = 1:length(rhos)
                    obj.MC2DLJ(obj.rho == rhos(i)).showStep(step,...
                        'saveFig',saveFig,'keepFigOpen',keepFigOpen);
               end
           else
               for i = 1:length(rhosInd)
                    obj.MC2DLJ(rhosInd(i)).showStep(step,...
                        'saveFig',saveFig,'keepFigOpen',keepFigOpen);
               end
           end
        end
        
        function obj = calcMeanWithoutFirstSteps(obj, firstSteps2ignore)
            for i = 1:length(obj.rho)
                obj.MC2DLJ(i).calcMeanWithoutFirstSteps(firstSteps2ignore);
            end
        end
        
        function obj = calcCv(obj)
            for i = 1:length(obj.rho)
                obj.MC2DLJ(i) = obj.MC2DLJ(i).calcCv();
                obj.cv(i) = obj.MC2DLJ(i).data.cv;
            end
        end
    end
end