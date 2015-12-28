function MC2DLJ = MonteCarlo2DLJ(varargin)

%% Monte-Carlo in NVT ensemble for Lennard-Jonse potantioal in 2D %%

% usage options:
% 1. new Monte Carlo simulation: 
%    MC2DLJ = MonteCarlo2DLJ(N,T,rho,Nsteps,maxdr,...
%               initialConfig,rCutoff,saveEvery)
% 2. continue Nstep in an old Monte Carlo simulation
%    MC2DLJ = MonteCarlo2DLJ(fileName,Nstep)
%       where fileName is the matlab file contaning the old simulation data
%       constructed with the class "MC2DLJoutput"

% get input

switch nargin
    case 2
        fileName = varargin{1};
        Nsteps = varargin{2};
        data = matfile(fileName);
        data = matfile(fileName,'Writable',true);
        
    case 8
        N = varargin{1};
        T = varargin{2};
        rho = varargin{3};
        Nsteps = varargin{4};
        maxdr = varargin{5};
        initialConfig = varargin{6};
        rCutoff = varargin{7};
        saveEvery = varargin{8};
end

r = 2^(1/6)/2; % particle radius in reduced units 

%% create initial configuration, find all pair distances %%
MC2DLJ = MC2DLJoutput(N,T,rho,maxdr,initialConfig,rCutoff,r);
 
%MC2DLJ.showStep(0);

%% monte carlo %%

% perform Nsteps steps of monte carlo, save every N*saveEvery steps
MC2DLJ = MC2DLJ.MonteCarlo(Nsteps,saveEvery);

end
