classdef isotherm
    
    properties
        pressure, rho, T;
        MC2DLJ, datafileList;
        simulationParam;
    end
    
    methods
        
        % constructor
        function obj = isotherm(varargin)
            
            if iscellstr(varargin) && ~isempty(varargin)
                if exist(varargin{1},'file')
                    % create object from a given file list
                    obj.datafileList = varargin(:);
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
                           [~,s] = size(obj.MC2DLJ(i).data.meanPlrc);
                           obj.pressure(i) = obj.MC2DLJ(i).data.meanPlrc(1,s);
                           obj.rho(i) = obj.MC2DLJ(i).simulationParam.rho;
                    end
                end
                
            else
          
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
                           [~,s] = size(obj.MC2DLJ(i).data.meanPlrc);
                           obj.pressure(i) = obj.MC2DLJ(i).data.meanPlrc(1,s);
                 end
             
        end
        
        function [obj, figHandle] = plotPropVsStep(obj,prop)
            indIndata = obj.MC2DLJ(1).indIndata;
            rhoN = length(obj.rho);
            x = zeros(rhoN,indIndata);
            y = x;
            switch prop
                case 'P'
                    for i = 1:rhoN
                        y(i,:) = obj.MC2DLJ(i).data.allPlrc;
                        x(i,:) = 1:indIndata;
                        leg{1,i} = ['rho = ' num2str(obj.rho(i))];
                    end
                    
                case 'U'
                    for i = 1:rhoN
                        y(i,:) = obj.MC2DLJ(i).data.allUlrc;
                        x(i,:) = 1:indIndata;
                        leg{1,i} = ['rho = ' num2str(obj.rho(i))];
                    end
            end
            
            figHandle = colorPlot(x,y,'lineStyle','-','addLegend',leg);
            title([prop ' Vs. steps for T = ' num2str(obj.T)]);
            xlabel('step');
            ylabel(prop);
        end
    end
end