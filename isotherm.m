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
                    parse(p, varargin{:});
                    Results = p.Results;
                    rl = Results.verelet;
                    N = Results.N;
                    T = Results.T;
                    rho = Results.rho;
                    initialmaxdr = Results.initialmaxdr;
                    initialconfig = Results.initialconfig;
                    rCutoff = Results.rCutoff;
                    r = Results.r;
                    rl = Results.verelet;
                    
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
                                'verelet',rl,'pressure',true);
                        else
                            obj.MC2DLJ(i) = MC2DLJoutput(N,T,rho(i),initialmaxdr,...
                                initialConfig,rCutoff,r,...
                                'verelet',rl,'pressure',true);
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
        
    end
end