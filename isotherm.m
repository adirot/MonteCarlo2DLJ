classdef isotherm
    
    properties
        pressure, rho, T;
        MC2DLJ, datafileList;
        simulationParam;
    end
    
    methods
        
        % constructor
        function obj = isotherm(N,T,rho,initialmaxdr,initialConfig,...
                rCutoff,r,varargin)
            
            if nargin == 1
                % create object from a given file list
                obj.datafileList = varargin{1};
                data = matfile(obj.datafileList{1,1});
                obj.T = data.T;
                obj.simulationParam = data.simulationParam;
                
                for i = 1:length(obj.datafileList)
                       obj.MC2DLJ(i) = MC2DLJoutput(obj.datafileList{1,i});
                       [~,s] = size(obj.MC2DLJ(i).data.meanPlrc);
                       obj.pressure(i) = obj.MC2DLJ(i).data.meanPlrc(1,s);
                       obj.rho(i) = obj.MC2DLJ(i).simulationParam.rho;
                end
                
                
            else
            
                    p = inputParser();
                    addOptional(p, 'verelet', []);
                    parse(p, varargin{:});
                    Results = p.Results;
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