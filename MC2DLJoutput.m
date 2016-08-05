classdef MC2DLJoutput
% an object that holds simulation data and parameter.    

% properties:

% simulationParam - simulation parmeter. includes:
%       simulationParam.N - number of particles
%       simulationParam.rho - reduced density
%       simulationParam.L - reduced board size
%       simulationParam.T - reduced temperature
%       simulationParam.r - reduced radius of a particle
%       simulationParam.initialmaxdr - initial maximum particle displacement
%       simulationParam.initialConfig - initial configuration of particles: 'random' or
%       'hcp' (6 nieghbors for each patrticle)
%       simulationParam.hcr - hard core repulsion: true or false 
%       simulationParam.rCutoff - the cutoff distance for the energy
%       simulationParam.rl - rl for verelet algorithm. if empty - verelet
%                               algorithm will not be used.
%       simulationParam.pressure - true or false, calculate the pressure or
%                                   not.
%       simulationParam.angleDependent - true or false (if the simulation
%               particles have an angle that defines them)
%       simulationParam.angleDependence - an anonymous function of two angles,
%               describing the angle dependence of the potantial. 
%       simulationParam.maxdAng - maximum change in angle in an MC step 


% currentU - the energy in the last step calculated
% Ulrc - the energy in the last step calculated, with long range
%           corrections
% currentVir - the virial in the last step calculated
% Plrc - the pressure in the last step calculated, with long range
%           corrections
% currentCoords - the coordinates of all particles in the last step 
%           calculated (2 by N matrix)
% currentDists - the pair distances of all particles in the last step 
%           calculated (N by N matrix)
% currentAngs - list of current cell angles
% currentAlphas - list of current alpha angles (alpha(i,j) is the angle
%       between the oriantations of cells i,j)
% currentThetas - list of current Theta angles (theta(i,j) is the angle
%       between the oriantation of cell i and the connecting vector of
%       cells i,j.  
% moveCount - counts accepted moves
% currentPressure - the pressure in the last step calculated
% currentmaxdr - the current maximum displacement of the Monte Carlo step.
% currentSweep - how many steps calculated so far
% indIndata - the number of steps saved so far (not every step is saved,
%               so this may be lower than currentSweep)
% fileName - the name of the matfile where all the data is stored. the file
%           name tells you the simulation parameters like N,T, etc. the last 
%           numbers in the file name are creation date and time
% histo - histogram of the radial distribution function (mean)
% bins - bins of the radial distribution function

% to save time, only the current step is saved in the object. the rest of
% the data is stored in a mat file. you can excess the data useing obj.data.

% data includes:
%       data.allCoords - all the coordinates in each step saved (size: 2 by N by
%           indIndata)
%       data.allDists - all the pair distances in each step saved (size: N by N by
%           indIndata)
%       data.allU - the energy in each step saved (size: 1 by indIndata)
%       data.allUlrc - the energy in each step saved, with long range correction (size: 1 by indIndata)
%       data.allV - the virial in each step saved (size: 1 by indIndata)
%       data.allPlrc - pressure in each step saved, with long range correction (size: 1 by indIndata)
%       data.sweepInd - the number of step clulated up to each saved point
%       data.moveCount - counts accepted moves
%       data.currentmaxdr - the current maximum displacement of the Monte
%               Carlo step.
%       data.simulationParam - all prameters of the simulation.
%       data.RDFhisto - histogram of the radial distribution function (for every
%                   step saved)
%       data.RDFbins - bins of the radial distribution function
%       data.allAng - all cell oriantations in each step
%       data.allAlphas - all alpha angles (allAlphas(i,j,k) is the angle between
%                   the oriantations of cells i,j in the k'th step)
%       data.allThetas - all Theta angles (allThetas(i,j,k) is the angle between
%                   the oriantation of cell i and the vector connecting
%                   cells i,j in the k'th step)


% usage example: 

% create new simulation output object:
%    N = 100; T = 0.1; rho = 0.1; initialmaxdr = 1; initialConfig = 'random';
%    rCutoff = 2.5; r = 2^(1/6)/2; rl = 3;
%    MC2DLJ = MC2DLJoutput(N,T,rho,initialmaxdr,initialConfig,rCutoff...
%               ,r,'verelet',rl,'pressure',true);

% monte carlo for 2000 steps, save every 5*N steps
%    MC2DLJ = MC2DLJ.MonteCarlo(1000,5)

% show a snapshot of some steps:
%   MC2DLJ.showStep('first');
%   MC2DLJ.showStep('mid');
%   MC2DLJ.showStep('last');
%   MC2DLJ.showStep(500); % if step 500 was not saved, this will show the
%                               closest step to 500.

% monte carlo for another 1000 steps, save every 5*N steps (this will
% continue the previus simulation)
%    MC2DLJ = MC2DLJ.MonteCarlo(1000,5)

% get all the energies for your data:
%   allU = MC2DLJ.data.allU
% get the coordinates of the 6'th step
%   coords = MC2DLJ.data.allCoords(:,:,6);

% delete MC2DLJ object, and get it back from your saved data
%   fileName = MC2DLJ.fileName;
%   clear MC2DLJ;
%   MC2DLJ = MC2DLJoutput(fileName);


   properties
      
      simulationParam = struct;
      currentmaxdr,...
          moveCount,currentCoords,currentDists,currentU,currentPressure,...
          currentVir,Ulrc,Plrc,currentSweep,currentStep,...
          currentThetas,currentAlphas,currentAngs,...
          currentBettas;
      fileName,data,indIndata;
      RDFhisto, RDFbins;
      runNum;
      RDFpeakLocs;
      
      
   end
    methods
        
        % constructor 
       function obj = MC2DLJoutput(varargin)
            
            % constructors for MC2DLJoutput. usage options:
            % 1. create new simulation output object: 
            %    obj = MC2DLJoutput(N,T,rho,initialmaxdr,...
            %           initialConfig,rCutoff,r);
            % 2. create a simulation output object linked to an old
            %   simulation, using the old simulations data file name:
            %   obj = MC2DLJoutput(fileName);
            
            
 
                if nargin == 1
                    % create simulation output object from old data
                    obj.fileName = varargin{1};
                    obj.data = matfile(obj.fileName);
                    obj.data = matfile(obj.fileName,'Writable',true);
                    obj.simulationParam = obj.data.simulationParam;
                    obj.currentmaxdr = obj.data.currentmaxdr;
                    obj.moveCount = obj.data.moveCount;
                    obj.indIndata = obj.data.indIndata;
                    
                    obj.currentCoords = obj.data.allCoords(:,:,obj.indIndata);
                    obj.currentDists =...
                        calcDists(obj.currentCoords,obj.simulationParam.L);
                    obj.currentU = obj.data.allU(1,obj.indIndata);
                    
                    if ~isfield(obj.simulationParam,'angleDependent')
                        obj.simulationParam.angleDependent = false;
                        obj.simulationParam.angleDependence = [];
                        obj.simulationParam.maxdAng = [];
                    end
                    
                    if ~isfield(obj.simulationParam,'ufunc')
                        obj.simulationParam.ufunc = [];
                    end
                    
                    if ~isfield(obj.simulationParam,'hcr')
                        obj.simulationParam.hcr = false;
                    end
                    
                    if obj.simulationParam.angleDependent
                   
                        obj.currentAngs = obj.data.allAngs(:,:,obj.indIndata);
%                        obj.currentAlphas = obj.data.allAlphas(:,:,obj.indIndata);
%                        obj.currentThetas = obj.data.allThetas(:,:,obj.indIndata);
                        obj.currentBettas = obj.data.allBettas(:,:,obj.indIndata);
                    end
                    
                    if ~isempty(obj.data.allV)
                        obj.currentVir = obj.data.allV(1,obj.indIndata);
                        obj.currentPressure = obj.data.allP(1,obj.indIndata);
                    end
                    
                    if existInMatfile(obj.fileName,'runNum')
                        obj.runNum = obj.data.runNum;
                    end
                    
                    if existInMatfile(obj.fileName,'RDFhisto')
                        obj.RDFhisto = mean(obj.data.RDFhisto);
                        obj.RDFbins = mean(obj.data.RDFbins);
                    end
                    
                    if existInMatfile(obj.fileName,'stepInd')
                        obj.currentStep = obj.data.stepInd(1,obj.indIndata);
                    else
                        obj.currentSweep = obj.data.sweepInd(1,obj.indIndata);
                    end

                end
                    
                if nargin >= 7
                    
                    
                    p = inputParser();
                    addOptional(p, 'verelet', []);
                    addOptional(p, 'pressure', false);
                    addOptional(p, 'fileNameInit', '');
                    addOptional(p, 'm', 6);
                    addOptional(p, 'runNum', []);
                    addOptional(p, 'angleDependent', false);
                    addOptional(p, 'angleDependence', []);
                    addOptional(p, 'maxdAng', []);
                    addOptional(p, 'ufunc', []);
                    addOptional(p, 'hcr', false);
                    addOptional(p, 'dontSaveDists', false);
                    parse(p, varargin{8:end});
                    Results = p.Results;
                    rl = Results.verelet;
                    pressure = Results.pressure;
                    fileNameInit = Results.fileNameInit;
                    m = Results.m;
                    runNum = Results.runNum;
                    angleDependent = Results.angleDependent;
                    angleDependence = Results.angleDependence;
                    maxdAng = Results.maxdAng;
                    ufunc = Results.ufunc;
                    hcr = Results.hcr;
                    dontSaveDists = Results.dontSaveDists;
                    
                    if dontSaveDists
                        pressure = true;
                    end
                    
                    if isempty(ufunc)
                        if hcr
                            ufunc = @(r) (-((1./r).^m));
                            ufuncstr = '_hcr_';
                        else
                            ufunc = @(r) (((1./r).^12)-((1./r).^m));
                            ufuncstr = '';
                        end
                    else
                        ufuncstr ='_ufunccustum_';
                    end
                    
                    % create a simulation output object for a new
                    % simulation
                    N = varargin{1};
                    T = varargin{2};
                    rho = varargin{3};
                    initialmaxdr = varargin{4};
                    initialConfig = varargin{5};
                    rCutoff = varargin{6};
                    r = varargin{7};
                    
                    obj.simulationParam.N = N;
                    obj.simulationParam.T = T;
                    obj.simulationParam.rho = rho;
                    obj.simulationParam.L = sqrt(N/rho);
                    obj.simulationParam.initialmaxdr = initialmaxdr;
                    obj.simulationParam.initialConfig = initialConfig;
                    obj.simulationParam.rCutoff = rCutoff;
                    obj.simulationParam.r = r;
                    obj.simulationParam.rl = rl;
                    obj.simulationParam.pressure = pressure;
                    obj.simulationParam.m = m;
                    obj.simulationParam.angleDependent = angleDependent;
                    obj.simulationParam.angleDependence = angleDependence;
                    obj.simulationParam.maxdAng = maxdAng;
                    obj.simulationParam.ufunc = ufunc;
                    obj.simulationParam.hcr = hcr;
                    obj.simulationParam.dontSaveDists = dontSaveDists;
                    
                    obj.currentmaxdr = obj.simulationParam.initialmaxdr;
                    obj.moveCount = 0;
                    obj.indIndata = 0;
                    obj.runNum = runNum;
                    
                    if ~isempty(rl)
                        vereletstr = ['verelet' my_num2str(rl)];
                    else
                        vereletstr = '';
                    end
                    
                    if ~isempty(pressure)
                        pressurestr = 'pressure';
                    else
                        pressurestr = '';
                    end
                    
                    if angleDependent
                        angleDependencestr = 'angleDependent_';
                    else
                        angleDependencestr = '';
                    end
                    
                    if angleDependent
                        maxdAngstr = ['maxdAng' my_num2str(maxdAng) '_'];
                    else
                        maxdAngstr = '';
                    end
                    
                    % make sure that L/2 is larger than rCutoff:
                    if obj.simulationParam.L/2 < rCutoff
                        error(['choose a larger number of particles, so that'...
                             'L/2 will be larger than rCutoff']);
                    end
                    
                    obj.fileName = [fileNameInit 'N' num2str(N)...
                        'T' my_num2str(T)...
                        'rho' my_num2str(rho) 'initialmaxdr'...
                        my_num2str(initialmaxdr) 'initialconfig_'...
                        initialConfig 'rCutoff'...
                         my_num2str(rCutoff) vereletstr '_' pressurestr...
                          '_m' num2str(m) ...
                          angleDependencestr ...
                          maxdAngstr ...
                          ufuncstr ...
                          'runNum' num2str(runNum) ...
                          'date'...
                        nowdatetimestr()];
                        
                    % create initial configuration
                    allCoords = zeros(2,N,1);
                    allDists = zeros(N,N,1);
                    [allDists(:,:,1),allCoords(:,:,1)] = ...
                            createInitialConfig(obj.simulationParam.L,N,r...
                            ,initialConfig);
                        
                    if angleDependent
                        allAngs = zeros(1,N,1);
%                         allAlphas = zeros(N,N,1);
%                         allThetas = zeros(N,N,1);
                        
                        allAngs(:,:,1) = rand(1,N,1)*pi;
                        %allAngs(:,:,1) = zeros(1,N,1);
                        allAngs(:,:,1) = rand(1,N,1)*pi;
%                         allAlphas(:,:,1) =...
%                             tril(bsxfun(@minus,allAngs(:,:,1),allAngs(:,:,1)'),-1);

                         xdist = bsxfun(@minus,allCoords(1,:,1),allCoords(1,:,1)');
                         ydist = bsxfun(@minus,allCoords(2,:,1),allCoords(2,:,1)');
%                         allThetas(:,:,1) = ...
%                             tril(allAlphas(:,:,1)-atan(ydist./xdist),-1);
                         allBettas = zeros(N,N,1);
                         allBettas(:,:,1) = tril(atan(ydist./xdist),-1);
                         
                        allAlphas = [];
                        allThetas = [];
                        
                    else
                        allAngs = [];
                        allAlphas = [];
                        allThetas = [];
                        allBettas = [];
                        
                    end

                        
                    % calculate initial energy
                    d = reshapeDist(allDists);
%                     ang = [reshapeDist(allAlphas),...
%                             reshapeDist(allThetas)];
%                     allU = pairU(d,rCutoff,m,...
%                         'angleDependence',angleDependence,...
%                        'relativeCellAngles',ang,'ufunc',ufunc);
                    [ang1, ang2] = reshapeAng(allAngs);
                    ang{1,1} = reshapeDist(allBettas);
                    ang{1,2} = ang1;
                    ang{1,3} = ang2;

                      allU = pairU(d,rCutoff,m,...
                        'angleDependence',angleDependence,...
                        'relativeCellAngles',ang,...
                        'numOfrelativeCellAngles',3,'ufunc',ufunc);
%                     [ang1, ang2] = reshapeAng(allAngs);
%                     ang{1,1} = reshapeDist(allBettas);
%                     ang{1,2} = ang1;
%                     ang{1,3} = ang2;
% 
%                       allU = pairU(d,rCutoff,m,...
%                         'angleDependence',angleDependence,...
%                         'relativeCellAngles',ang,...
%                         'numOfrelativeCellAngles',3,'ufunc',ufunc);
                      allU = pairU(d,rCutoff,m,...
                        'angleDependence',angleDependence,'ufunc',ufunc);
                    %%% this is the long range correction - should be fixed
                    %%% for the case of angle dependence
                    %allUlrc = allU - pi*rho*N/rCutoff^4;
                    
                    % calculate initial virial and pressere
                    if pressure
                        allV = ...
                            calcVirial(d,rho,12,m,N,rCutoff);
                        allP = T*rho + allV;
                        allPlrc = allP - 3*pi*rho^2/rCutoff^4;
                    else
                        allV = [];
                        allP = [];
                        allPlrc = [];
                    end
                    clear d; clear ang;
                    
                    % save to data file
                    sweepInd = 0;
                    moveCount = 0;
                    indIndata = 1;
                    currentmaxdr = initialmaxdr;
                    simulationParam = obj.simulationParam;
%                     save(obj.fileName, 'allDists','allCoords',...
%                             'allU','allUlrc','allV','allP','allPlrc','sweepInd',...
%                             'moveCount','indIndata'...
%                             ,'currentmaxdr','simulationParam','runNum',...
%                             'allAlphas','allThetas',...
%                             'allAngs','allBettas','-v7.3');
                    save(obj.fileName, 'allDists','allCoords',...
                            'allU','allV','allP','allPlrc','sweepInd',...
                            'moveCount','indIndata'...
                            ,'currentmaxdr','simulationParam','runNum',...
                            'allAlphas','allThetas',...
                            'allAngs','allBettas','-v7.3');
                    obj.data = matfile(obj.fileName);
                    obj.data = matfile(obj.fileName,'Writable',true);
                        
                    obj.currentCoords = allCoords;
                    obj.currentDists = allDists;
                    obj.currentAngs = allAngs;
                    obj.currentAlphas = allAlphas;
                    obj.currentThetas = allThetas;
                    obj.currentBettas = allBettas;
                    obj.currentU = allU;
                    %obj.Ulrc = allUlrc;
                    obj.currentVir = allV;
                    obj.currentPressure = allP;
                    obj.Plrc = allPlrc;
                    obj.currentSweep = sweepInd;
                    obj.indIndata = 1;
                end
            
        end
        
       function obj = MonteCarlo(obj,Nsweeps,saveEvery,varargin)
           p = inputParser();
           addOptional(p, 'TalkEvery', []);
           addOptional(p, 'save2data', true);
           addOptional(p, 'logFile', true);
           addOptional(p, 'logFileInit', '');
           parse(p, varargin{:});
           R = p.Results;
           TalkEvery = R.TalkEvery;
           save2data = R.save2data;
           logFile = R.logFile;
           logFileInit = R.logFileInit;
           
            N = obj.simulationParam.N; 
            rho = obj.simulationParam.rho;
            rCutoff = obj.simulationParam.rCutoff;        
            m = obj.simulationParam.m;
            
            angleDependent = obj.simulationParam.angleDependent;
            angleDependence = obj.simulationParam.angleDependence;
            maxdAng = obj.simulationParam.maxdAng;
            ufunc = obj.simulationParam.ufunc;
            hcr = obj.simulationParam.hcr;
            if hcr
                hardCoreRepRad = 1;
                %if obj.simulationParam.angleDependent
                    %hardCoreRepRad = (m/12)^(1/(m-12))/2;
                    
                %end
            else 
                hardCoreRepRad = 0;
            end
            
            if saveEvery == 0
                % Only save last step
                numOfruns2save = N*Nsweeps;
            else
                numOfruns2save = N*saveEvery;
            end
            
            totSec = 0;
            tic;
            sweepCount = 0;
            totSec = 0;
            
            while(sweepCount < Nsweeps)
                tic;
                [finalU,finalV,finalPressure,...
                    finalConfiguration,finalDistances,...
                    currentmoveCount,...
                    finalAngs,...
                    finalBettas] = ... 
                    MonteCarlo2DLJHeart(...
                    N,...
                    obj.simulationParam.T,...
                    rho,...
                    numOfruns2save,...
                    obj.currentmaxdr,...
                    obj.currentCoords,...
                    rCutoff,...
                    obj.currentDists,...
                    obj.currentU,...
                    'verelet',obj.simulationParam.rl,...
                    'virial',obj.currentVir,...
                    'm',m,...
                    'angleDependent',angleDependent,...
                    'angleDependence',angleDependence,...
                    'initialAng',obj.currentAngs,...
                    'initialBettas', obj.currentBettas,...
                    'maxdAng',maxdAng,...
                    'ufunc',ufunc,...
                    'hardCoreRepRad',hardCoreRepRad,...
                    'TalkEvery',TalkEvery);
                
                sweepCount = sweepCount + numOfruns2save/N;
                obj.currentSweep = obj.currentSweep + numOfruns2save/N;
                obj.moveCount = obj.moveCount + currentmoveCount;
                obj.currentU = finalU;
                obj.currentVir = finalV;
                obj.currentPressure = finalPressure;
                obj.currentCoords = finalConfiguration;
                obj.currentDists = finalDistances;
                obj.currentAngs = finalAngs;
%                obj.currentAlphas = finalAlphas;
%                obj.currentThetas = finalThetas;
                obj.currentBettas = finalBettas;
            
                obj.currentAlphas = finalAlphas;
                obj.currentThetas = finalThetas;
                
                % long range correction
                obj.Ulrc = obj.currentU + 4*pi*N*rho/((2-m)*rCutoff^(m-2));
                if obj.simulationParam.pressure
                    obj.Plrc = obj.currentPressure -...
                        2*m*pi*rho^2/((2-m)*rCutoff^(m-2));
                end
            
                if ~isempty(TalkEvery)
                    talk = true;
                else
                    talk = false;
                end
                
                
                
                
                if save2data
%                     obj = obj.addStep2data(obj.currentSweep,finalConfiguration,...
%                         finalDistances,finalU,finalV,finalPressure,...
%                         obj.moveCount,obj.currentmaxdr,obj.Ulrc,obj.Plrc,...
%                         angleDependent,obj.currentAngs,...
%                         obj.currentAlphas,obj.currentThetas,...
%                         obj.simulationParam.dontSaveDists,talk);
                    obj = obj.addStep2data(obj.currentSweep,finalConfiguration,...
                        finalDistances,finalU,finalV,finalPressure,...
                        obj.moveCount,obj.currentmaxdr,obj.Ulrc,obj.Plrc,...
                        angleDependent,obj.currentAngs,...
                        obj.currentBettas,...
                        obj.simulationParam.dontSaveDists,talk);

                end
                
                

                
                % Create log file
                if logFile
                    
                    if sweepCount > numOfruns2save/N
                        delete(lastFileName);
                    end

                    totSec = totSec + toc;
                    lastFileName = [logFileInit ...
                        'T' my_num2str(obj.simulationParam.T)...
                        'rho' my_num2str(obj.simulationParam.rho)...
                        'm' num2str(obj.simulationParam.m)...
                        'sweepsDone' num2str(sweepCount) 'secPassed'...
                        my_num2str(totSec) 'date' nowdatetimestr() '.txt'];
                    fileID = fopen(lastFileName, 'w');
                    fclose(fileID);
                end
                
                clear finalU finalV finalConfiguration finalDistances...
                    currentmoveCount finalAngs finalAlphas finalThetas...
                    finalBettas
                    
            end
            
        end
        
%        function obj = addStep2data(obj,newInd,newCoords,newDists,newU,newV...
%                 ,newP,moveCount,currentmaxdr,newUlrc,newpressurelrc,...
%                 angleDependent,newAngs,newAlphas,newThetas,dontSaveDists,...
%                 varargin)
       function obj = addStep2data(obj,newInd,newCoords,newDists,newU,newV...
                ,newP,moveCount,currentmaxdr,newUlrc,newpressurelrc,...
                angleDependent,newAngs,newBettas,dontSaveDists,...
                varargin)
          
            p = inputParser();
            addOptional(p, 'talk', false);
            parse(p, varargin{:});
            R = p.Results;
            talk = R.talk;
            
            if talk
                disp(['adding to data. size of allU: '...
                    num2str(size(obj.data.allU))]);
            end
            
            obj.data.sweepInd(1,obj.indIndata+1) = newInd;
            if obj.indIndata == 1
                s = zeros(2,obj.simulationParam.N,2);
                s(:,:,2) = newCoords(:,:,1);
                s(:,:,1) = obj.data.allCoords(:,:);
                obj.data.allCoords = s;
                s = zeros(obj.simulationParam.N,obj.simulationParam.N,2);
                s(:,:,2) = newDists(:,:,1);
                s(:,:,1) = obj.data.allDists(:,:);
                obj.data.allDists = s;
                clear s;
                
                if dontSaveDists
                    [bins, histo] =...
                        calculateRDF(newDists,obj.simulationParam.N,...
                        obj.simulationParam.rho,300,10);
                    obj.data.RDFbins = bins;
                    obj.data.RDFhisto = zeros(1,300,2);
                    obj.data.RDFhisto(1,1:300,1) = histo;
                end
            else
                obj.data.allCoords(:,:,obj.indIndata+1) = newCoords(:,:,1);
                if dontSaveDists
                    [~, histo] =...
                        calculateRDF(newDists,obj.simulationParam.N,...
                        obj.simulationParam.rho,300,10);
                    obj.data.RDFhisto(1,1:300,obj.indIndata+1) = histo;
                else
                    obj.data.allDists(:,:,obj.indIndata+1) = newDists(:,:,1);
                end
            end
            
            obj.data.allU(1,obj.indIndata+1) = newU;
            obj.data.allUlrc(1,obj.indIndata+1) = newUlrc;
            
            if angleDependent
                if obj.indIndata == 1
                    s = zeros(1,obj.simulationParam.N,2);
                    s(:,:,2) = newAngs(:,:,1);
                    s(:,:,1) = obj.data.allAngs(:,:);
                    obj.data.allAngs = s;
                    s = zeros(obj.simulationParam.N,obj.simulationParam.N,2);
 
 %                   s(:,:,2) = newAlphas(:,:,1);
 %                   s(:,:,1) = obj.data.allAlphas(:,:);
 %                   obj.data.allAlphas = s;
 %                   s = zeros(obj.simulationParam.N,obj.simulationParam.N,2);
 %                   s(:,:,2) = newThetas(:,:,1);
 %                   s(:,:,1) = obj.data.allThetas(:,:);
 %                   obj.data.allThetas = s;
                   s(:,:,2) = newBettas(:,:,1);
                   s(:,:,1) = obj.data.allBettas(:,:);
                   obj.data.allBettas = s;
 
                   s(:,:,2) = newAlphas(:,:,1);
                   s(:,:,1) = obj.data.allAlphas(:,:);
                   obj.data.allAlphas = s;
                   s = zeros(obj.simulationParam.N,obj.simulationParam.N,2);
                   s(:,:,2) = newThetas(:,:,1);
                   s(:,:,1) = obj.data.allThetas(:,:);
                   obj.data.allThetas = s;
 
 
                    clear s;
                else
                    obj.data.allAngs(:,:,obj.indIndata+1) = newAngs(:,:,1);
 
%                    obj.data.allAlphas(:,:,obj.indIndata+1) = newAlphas(:,:,1);
%                    obj.data.allThetas(:,:,obj.indIndata+1) = newThetas(:,:,1);
                    obj.data.allBettas(:,:,obj.indIndata+1) = newBettas(:,:,1);
 
                   obj.data.allAlphas(:,:,obj.indIndata+1) = newAlphas(:,:,1);
                   obj.data.allThetas(:,:,obj.indIndata+1) = newThetas(:,:,1);
 

                end
            end
            
            if obj.simulationParam.pressure
                obj.data.allV(1,obj.indIndata+1) = newV;
                obj.data.allP(1,obj.indIndata+1) = newP;
                obj.data.allPlrc(1,obj.indIndata+1) = newpressurelrc;
            end
            obj.indIndata = obj.indIndata + 1;
            obj.data.indIndata = obj.indIndata;
            obj.data.moveCount = moveCount;
            obj.data.currentmaxdr = currentmaxdr;
            
            
            if talk
                disp(['done adding to data. size of allU: '...
                    num2str(size(obj.data.allU))]);
            end
            
       end
       
       function [obj, bins, RDFhisto] =...
               calcRDF(obj,maxDist,numOfBins,varargin)
            %% calculate 2D radial distribution function, with pediodic boundary
            %% conditions.
            
            %         input:
            %         ~~~~~~
            
            %         NumOfBins - number of bins in the histogram
            %
            %         output:
            %         ~~~~~~~
            %         obj.data.bins - x axis of the RDF histogram
            %         obj.data.histo - each colmun is the RDF axis of a diffrent step in the
            %         monte carlo simulation. (so to plot the RDF for the 10'th step we
            %         need to write: plot(bins,histo(:,10));

            %       method:
            %       ~~~~~~~
            %       the RDF is calculated by binnig all pair partical distances into 
            %       a histogram, and normalizing each bin with it's Ideal gas number of
            %       particals. 
            %       when binning the pair distances, we take into acount Periodic
            %       Boudary Conditions (PBC)
            %       finaly, to ansure that when r->infinity : RDF ->1 , we
            %       normalize by multiplying RDF*2/(N-1) where N is the number of
            %       particals. 
            %       for more information
            %   http://www2.msm.ctw.utwente.nl/sluding/TEACHING/APiE_Script_v2011.pdf
            %       page 48 - "Radial distribution function"
            
            p = inputParser();
            addOptional(p, 'save2data', true);
            addOptional(p, 'startFrom', 1);
            addOptional(p, 'talk', false);
            parse(p, varargin{:});
            Results = p.Results;
            save2data = Results.save2data;
            startFrom = Results.startFrom;
            talk = Results.talk;
            
              N = obj.simulationParam.N;
              rho = obj.simulationParam.rho;
              
              bins = linspace(0,maxDist,numOfBins); 
              if save2data
                    obj.data.RDFhisto = ...
                        zeros(obj.indIndata-startFrom+1,numOfBins);      
                    obj.data.RDFbins = bins;
                    obj.RDFhisto = zeros(1,numOfBins);
              else 
                    RDFhisto = zeros(obj.indIndata-startFrom+1,numOfBins);
              end
              
              for step = startFrom:obj.indIndata
                    if talk
                        disp(['calcRDF step: ' num2str(step)]);
                    end
                    dist = obj.data.allDists(:,:,step);
                    [~, histo] = calculateRDF(dist,N,rho,numOfBins,maxDist);
                    if save2data
                            obj.data.RDFhisto(step,:) = histo;
                            obj.RDFhisto = obj.RDFhisto + histo;
                    else
                            RDFhisto(step,:) = histo;
                    end
                            
              end
              
              if save2data
                    obj.RDFbins = bins;
                    obj.RDFhisto = obj.RDFhisto/(obj.indIndata-startFrom+1);
              end
       end
       
       function [obj, RDFlocs] = findRDFpeaks(obj,varargin)
           
           p = inputParser();
           addOptional(p, 'minPeakHeight', 2);
           parse(p, varargin{:});
           Results = p.Results;
           minPeakHeight = Results.minPeakHeight;
           
           % find peaks
           RDFlocs = [];RDFlocsmin = [];
           %  find peaks and valies larger than
           % minPeakHight
           [~, RDFlocs] =...
               findpeaks(obj.RDFhisto,'MINPEAKHEIGHT',minPeakHeight);
           [~, RDFlocsmin] =...
               findpeaks(2-obj.RDFhisto,...
               'MINPEAKHEIGHT',minPeakHeight);
           [~, RDFind] = min(abs(2^(1/6) - obj.RDFbins)); % find location
                                                       % of the first peak
           RDFlocs = unique([RDFlocs RDFind RDFlocsmin]);
           obj.RDFpeakLocs = RDFlocs;
                     
       end
       
       function [obj, varVarPeaks] = RDFvarVar(obj,varargin)
           % calculates the variance and variance of the variance for
           % spesific points in the RDF:  the first peak. 
           
           p = inputParser();
           addOptional(p, 'plotFig', false);
           addOptional(p, 'saveFig', false);
           addOptional(p, 'add2savedFigName', '');
           addOptional(p, 'keepFigOpen', false);
           addOptional(p, 'firstSteps2ignore', 0);
           addOptional(p, 'peackLoc', 39);
           parse(p, varargin{:});
           Results = p.Results;
           plotFig = Results.plotFig;
           saveFig = Results.saveFig;
           add2savedFigName = Results.add2savedFigName;
           keepFigOpen = Results.keepFigOpen;
           firstSteps2ignore = Results.firstSteps2ignore;
           peackLoc = Results.peackLoc;
           
           % check if RDF was calculated allready
           if ~existInMatfile(obj.data,'RDFhisto')
               error('calculate RDF first');
           end
           
           [~, ~, Nsteps] = size(obj.data.RDFhisto);
           
           if firstSteps2ignore == 0
               minInd = 1;
           else
               % find the relevent initial index
               [~, minInd] = min(abs(firstSteps2ignore - obj.data.sweepInd));
           end
           
           % find peaks
%            [obj, locs] = obj.findRDFpeaks();
            locs = peackLoc;
            peaks = zeros(length(locs),Nsteps-minInd+1);
            steps = zeros(length(locs),Nsteps-minInd+1);
            

           for i = minInd:Nsteps
               for k = 1:length(locs)
                   peaks(k,i-minInd+1) =...
                       obj.data.RDFhisto(1,locs(k),i);
               end
               
               % find variance and variance of variance in the peaks
               for j = 1:length(locs)
                    varPeaks(j,i-minInd+1) = var(peaks(j,1:(i-minInd+1)));
                    varVarPeaks(j,i-minInd+1) =...
                        var(varPeaks(j,1:(i-minInd+1)));
               end
           end
           
           for j = 1:length(locs)
               steps(j,:) = obj.data.sweepInd(1,minInd:Nsteps);
           end
           
           if plotFig
               
               for j = 1:length(locs)
                   leg{1,j} = ['peak distance '...
                       num2str(obj.data.RDFbins(1, locs(j)))]; 
               end
               
               colorPlot(steps,varPeaks,'addLegend',leg);
               xlabel('steps');
               ylabel('RDF peaks variance');
               
               if saveFig
                    name = ['varRDFpeaksVsSteps_T'...
                    my_num2str(obj.simulationParam.T)...
                    'N' my_num2str(obj.simulationParam.N) 'rho'...
                    my_num2str(obj.simulationParam.rho)];
                    saveas(gcf,[name add2savedFigName '.fig']);
                    saveas(gcf,[name add2savedFigName '.jpg']);
               end
             
               colorPlot(steps,varVarPeaks,'addLegend',leg);
               xlabel('steps');
               ylabel('RDF peaks variance of variance');
               
               if saveFig
                    name = ['varVarRDFpeaksVsSteps_T'...
                    my_num2str(obj.simulationParam.T)...
                    'N' my_num2str(obj.simulationParam.N) 'rho'...
                    my_num2str(obj.simulationParam.rho)];
                    saveas(gcf,[name add2savedFigName '.fig']);
                    saveas(gcf,[name add2savedFigName '.jpg']);
               end
               
               if ~keepFigOpen
                   close all;
               end
           end
               
       end

       function [obj, histxnumOfPartInSquare, histnumOfPartInSquare,...
               numOfPartInSquare] =...
               calcRhoDistrib(obj, numOfSquares, varargin)
       % We devide every step to numOfSquares subsystem. We count how many
       % particles are in each subsystem, in every step.
       % histnumOfCellsInSquare is the histogram of densities in 
       % the system, considering all steps.
       % To get the probability of a spesific density in a subsystem:
       % PL(rho) = histnumOfPartInSquare(rho)/(indIndata*numOfSquares);
       % To get the probability that a particle will be found in a
       % subsystem with a spesific density:
       % PLN(rho) = numOfPartInSquare(rho)*histnumOfPartInSquare(rho)/N
       
       p = inputParser();
       addOptional(p, 'plotHist', false); 
       addOptional(p, 'plotHist_times_partNum', false);
       addOptional(p, 'startFrom', 1);
       addOptional(p, 'endAt', []);
       addOptional(p, 'coordsFromObj', false);
       parse(p, varargin{:});
       Results = p.Results;
       plotHist = Results.plotHist;
       plotHist_times_partNum = Results.plotHist_times_partNum;
       startFrom = Results.startFrom;
       endAt = Results.endAt;
       coordsFromObj = Results.coordsFromObj;
       
       L = obj.simulationParam.L;
       if isempty(endAt)
            endAt = obj.indIndata;
       end
       N = obj.simulationParam.N;    
           
       if coordsFromObj
           numOfPartInSquare(1,1:numOfSquares) =...
                numOfCellsDistribution(obj.currentCoords,...
                L,numOfSquares);
       else
           numOfPartInSquare = zeros(1,numOfSquares*(endAt-startFrom+1));

           for i = 1:(endAt-startFrom+1)
               numOfPartInSquare(1,(numOfSquares*(i-1)+1):(numOfSquares*i)) =...
                    numOfCellsDistribution(obj.data.allCoords(1:2,1:N,i+startFrom-1),...
                    L,numOfSquares);
           end
       end
       
       histxnumOfPartInSquare = 0:N;
       histnumOfPartInSquare =...
           hist(numOfPartInSquare,histxnumOfPartInSquare);
       % Normalize
       histxnumOfPartInSquare = histxnumOfPartInSquare / ((L^2/numOfSquares)*(N/L^2));       
       
       if plotHist
           figure;
           bin = histxnumOfPartInSquare(2);
           areaUnderPlot = sum(bin*histnumOfPartInSquare.*histxnumOfPartInSquare);
           plot(histxnumOfPartInSquare,histnumOfPartInSquare.*histxnumOfPartInSquare/areaUnderPlot);
           xlabel('Density');
           ylabel('Probability of a particle to be in a square with a certain density');
           title('Density histogram');
       end
       
       if plotHist_times_partNum
           figure;
           bin = histxnumOfPartInSquare(2);
           areaUnderPlot = sum(bin*histnumOfPartInSquare);
           plot(histxnumOfPartInSquare,histnumOfPartInSquare/areaUnderPlot);
           xlabel('Density');
           ylabel('Probability of a subsystem to be in a square with a certain density');
           title('Density histogram');
       end
       
       end

       
       function showStep(obj,step)
           
          if isnumeric(step)
                  try
                      ind = find(obj.data.sweepInd == step);
                  catch
                      ind = find(obj.data.stepInd == step);
                  end
            
                  if isempty(ind) %find the closest index saved
                      try
                          [~,ind] = min(abs(obj.data.sweepInd - step));
                          step = obj.data.sweepInd(1,ind);
                      catch
                          [~,ind] = min(abs(obj.data.stepInd - step));
                          step = obj.data.stepInd(1,ind);
                      end
                  end
          else
              
              switch step
                  case 'first'
                      ind = 1;
                      step = 0;
                  case 'mid'
                     
                      ind = ceil(obj.indIndata/2);
                      try
                          step = obj.data.sweepInd(1,ind);
                      catch
                          step = obj.data.stepInd(1,ind);
                      end
                  case 'last'
                      ind  = obj.indIndata;
                      try
                          step = obj.data.sweepInd(1,ind);
                      catch
                          step = obj.data.stepInd(1,ind);
                      end
              end
              
          end 
          
          plotParticles(obj.data.allCoords(:,:,ind),obj.simulationParam.L...
              ,obj.simulationParam.r,obj.data.allAngs(1,:,ind));
          title(['snapshot of step: ' num2str(step)]);
       end
       
       function obj = meanProp(obj)
            % makes mean vectors for U,P,RDFpeaks.
            % meanProp(i) is the mean of all the values of property on
            % steps 1 to i.
            
            obj.data.meanU = my_mean(obj.data.allU);
            obj.data.meanUlrc = my_mean(obj.data.allUlrc);
            obj.data.meanP = my_mean(obj.data.allP);
            obj.data.meanPlrc = my_mean(obj.data.allPlrc);
            [obj, locs] = obj.findRDFpeaks;
            if isempty(obj.RDFhisto)
                obj.RDFhisto = mean(obj.data.RDFhisto);
                obj.RDFbins = obj.data.RDFbins;
            end
            obj.data.meanRDFpeaks = obj.RDFhisto(1,locs);
            
       end
       
       function [obj, figHandle] = plotPropVsStep(obj,prop,varargin)
           
            p = inputParser();
            addOptional(p, 'figHandle', []);
            addOptional(p, 'lineStyle', []);
            parse(p, varargin{:});
            Results = p.Results;
            figHandle = Results.figHandle;
            lineStyle = Results.lineStyle;
            
            if isempty(figHandle)
                figHandle = figure;
                hold on;
            else
                figure(figHandle);
                hold on;
            end
            
            switch prop
                case 'P'
                    if isempty(lineStyle)
                        plot(obj.data.allPlrc);
                    else
                        plot(obj.data.allPlrc,'lineStyle',lineStyle);
                    end
                case 'U'
                    if isempty(lineStyle)
                        plot(obj.data.allUlrc);
                    else
                        plot(obj.data.allUlrc,'lineStyle',lineStyle);
                    end
            end
            
            figHandle = gcf;
       end
       
       function obj = calcMeanWithoutFirstSteps(obj, firstSteps2ignore)
           obj.data.firstSteps2ignore = firstSteps2ignore;
           
           if ~existInMatfile(obj.data,'allPlrc')
                obj = obj.calcPressureAfterRun();
           else
               [~,s] = size(obj.data.allPlrc);
               if s == 0
                     obj = obj.calcPressureAfterRun();
               end
           end 
           
           obj.data.meanPlrcEq =...
               mean(obj.data.allPlrc(1,firstSteps2ignore:obj.indIndata));
           obj.data.meanUlrcEq =...
               mean(obj.data.allUlrc(1,firstSteps2ignore:obj.indIndata));
           
           if isempty(obj.RDFpeakLocs)
               [obj, ~] = obj.findRDFpeaks();
           end
           
           obj.data.meanRDFpeaksEq = ...
               mean(obj.data.RDFhisto(firstSteps2ignore:obj.indIndata,...
               obj.RDFpeakLocs));
       end
        
       function [obj,tauP,tauU] = inefficiency(obj,n,varargin)
           % for error calculations. see computer simulation of liquids
           % page 192
           
           p = inputParser();
           addOptional(p, 'plotSVsSqrtTau', true);
           addOptional(p, 'saveFigSVsSqrtTau', true);           
           addOptional(p, 'firstSteps2ignore', []);
           addOptional(p, 'keepFigOpen', false);
           parse(p, varargin{:});
           Results = p.Results;
           plotSVsSqrtTau = Results.plotSVsSqrtTau;
           firstSteps2ignore = Results.firstSteps2ignore;
           saveFigSVsSqrtTau = Results.saveFigSVsSqrtTau;
           keepFigOpen = Results.keepFigOpen;
           
           
           % calculate <A>
           if isempty(firstSteps2ignore)
                if existInMatfile(obj.fileName,'meanPlrcEq')
                    meanP = obj.data.meanPlrcEq;
                    meanU = obj.data.meanUlrcEq;
                    meanRDF = obj.data.meanRDFpeaksEq;
                    firstSteps2ignore = obj.data.firstSteps2ignore;
                else
                    error(['you must provide a number of first steps'... 
                        ' to ignore or calculate the mean in advence'... 
                        ' with calcMeanWithoutFirstSteps']);
                end
           else
               obj = calcMeanWithoutFirstSteps(obj, firstSteps2ignore);
               meanP = obj.data.meanPlrcEq;
               meanU = obj.data.meanUlrcEq;
               meanRDFpeaks = obj.data.meanRDFpeaksEq;
           end
           
           P = obj.data.allPlrc(1,firstSteps2ignore:obj.data.indIndata);
           U = obj.data.allUlrc(1,firstSteps2ignore:obj.data.indIndata);
           RDFpeaks = obj.RDFhisto(firstSteps2ignore:obj.data.indIndata,...
               obj.RDFlocs);
           [~, Npeaks] = size(RDFpeaks);

           nt = (obj.indIndata - firstSteps2ignore);
           
           for j = 1:length(n)
                % calculate <A>b
                ind = 1;
                for i = 1:n(j):(nt-n(j)+1)
                   
                    Pmeanb(ind) = mean(P(i:(i+n(j)-1)));
                    Umeanb(ind) = mean(U(i:(i+n(j)-1)));
                    RDFpeaksMeanb(ind,1:Npeaks) =...
                        mean(RDFpeaks(i:(i+n(j)-1),1:Npeaks));
                    ind = ind + 1;
                end
           
                % calculate sigma^2(<A>b) for all the different tau values
                
                nb(j) = length(Pmeanb);
                varMeanP(j) = mean((Pmeanb - mean(Pmeanb)).^2);
                varMeanU(j) = mean((Umeanb - mean(Umeanb)).^2);
                varMeanRDFpeaks(j,1:Npeaks) =...
                    mean((RDFpeaksMeanb - mean(RDFpeaksMeanb)).^2);
                Pmeanb = [];
                Umeanb = [];
                RDFpeaksMeanb = [];
                
           end
            
           %calculate sigma^2(A)
           varP = mean((P(:) - meanP).^2);
           varU = mean((U(:) - meanU).^2);
           varRDFpeaks = mean((RDFpeaks - meanP).^2);

           % calculate tau
           tauP = n.*varMeanP/varP;
           tauU = n.*varMeanU/varU;
           
           for i = 1:Npeaks
                tauRDFpeaks(i,:) = n.*varMeanRDFpeaks(:,i)/varRDFpeaks(:,i);
                tauRDFpeaksX(i,:) = sqrt(n);
                RDFpeaksLeg{1,i} = ['peak num: ' num2str(i)];
           end
           
           if saveFigSVsSqrtTau
               plotSVsSqrtTau = true;
           end
           
           if plotSVsSqrtTau
               figure;
               plot(sqrt(n),tauP);
               hold on;
               title('$$t_A^c$$ for the pressure and energy',...
                   24,'Interpreter','latex');
               xlabel('$$\sqrt{t_b}$$','FontSize',24,'Interpreter','latex');
               ylabel('$$t_A^c=\frac{t_b \sigma^2 (<A>_b)}{\sigma^2(A)}$$',...
                   'FontSize',24,'Interpreter','latex');
               plot(sqrt(n),tauU,'r');
               legend('pressure','energy');
               
               if saveFigSVsSqrtTau
                   fileName = ['SVsSqrtTauN' num2str(obj.simulationParam.N)... 
                       'T' my_num2str(obj.simulationParam.T)...
                       'rho' my_num2str(obj.simulationParam.rho)];
                   saveas(gcf,[fileName '.fig']);
                   saveas(gcf,[fileName '.jpg']);
               end

               
               colorPlot(tauRDFpeaksX,tauRDFpeaks,'addLegend',RDFpeaksLeg);
               title('$$t_A^c$$ for RDF peaks',...
                   24,'Interpreter','latex');
               xlabel('$$\sqrt{t_b}$$','FontSize',24,'Interpreter','latex');
               ylabel('$$t_A^c=\frac{t_b \sigma^2 (<A>_b)}{\sigma^2(A)}$$',...
                   'FontSize',24,'Interpreter','latex');
               
               if saveFigSVsSqrtTau
                   fileName = ['SVsSqrtTauNRDF' num2str(obj.simulationParam.N)... 
                       'T' my_num2str(obj.simulationParam.T)...
                       'rho' my_num2str(obj.simulationParam.rho)];
                   saveas(gcf,[fileName '.fig']);
                   saveas(gcf,[fileName '.jpg']);
               end
               
           end
           
           
           if ~keepFigOpen
               close gcf;
           end
       end
       
       function [obj, varU, varP, steps, varVarU, varVarP] =...
               varOfvar(obj,varargin)
           % calculate the variance of the variance of P,U as a function of
           % steps in the simulation.
           
           p = inputParser();
           addOptional(p,'startFromStep',0);
           addOptional(p, 'plotVarVsStep', true);
           addOptional(p, 'plotVarVarVsStep', true);
           addOptional(p, 'saveFig', true);
           addOptional(p, 'add2savedFigName', '');
           addOptional(p, 'keepFigOpen', false);
           parse(p, varargin{:});
           Results = p.Results;
           plotVarVsStep = Results.plotVarVsStep;
           plotVarVarVsStep = Results.plotVarVarVsStep;
           saveFig = Results.saveFig;
           add2savedFigName = Results.add2savedFigName;
           keepFigOpen = Results.keepFigOpen;
           startFromStep = Results.startFromStep;
           
           if startFromStep == 0
               minInd = 1;
           else
               % find the relevent initial index
               [~, minInd] = min(abs(startFromStep - obj.data.sweepInd));
           end
           
           
           for i = minInd:obj.data.indIndata
               
               varU(i-minInd+1) = var(obj.data.allUlrc(1,minInd:i));
               varVarU(i-minInd+1) = var(varU);
               varP(i-minInd+1) = var(obj.data.allPlrc(1,minInd:i));
               varVarP(i-minInd+1) = var(varP);
           end
           
           obj.data.varU = varU;
           obj.data.varP = varP;
           obj.data.varVarU = varVarU;
           obj.data.varVarP = varVarP;
           
           steps = obj.data.sweepInd(1,minInd:obj.indIndata);
           
           if plotVarVsStep
               figure;
               hold on;
               plot(steps,varU);
               plot(steps,varP,'r');
               legend({'Energy variance','Pressure variance'});
               xlabel('steps');
               ylabel('Energy or Pressure variance');
               title(['variance for Energy and Pressure Vs. steps, T = '...
                   num2str(obj.simulationParam.T) ' N = '...
                   num2str(obj.simulationParam.N) ' \rho = '... 
                   num2str(obj.simulationParam.rho)]);
           end
           
           if saveFig
               name = ['varUPvsSteps_T' my_num2str(obj.simulationParam.T)...
                   'N' my_num2str(obj.simulationParam.N) 'rho'...
                   my_num2str(obj.simulationParam.rho)];
               saveas(gcf,[name add2savedFigName '.fig']);
               saveas(gcf,[name add2savedFigName '.jpg']);
           end
           
           if plotVarVarVsStep
               figure;
               hold on;
               plot(steps,varVarU);
               plot(steps,varVarP,'r');
               legend({'Energy variance of variance',...
                   'Pressure variance of variance'});
               xlabel('steps');
               ylabel('Energy or Pressure variance of variance');
               title(['variance of variance '...
                   'for Energy and Pressure Vs. steps, T = '...
                   num2str(obj.simulationParam.T) ' N = '...
                   num2str(obj.simulationParam.N) ' \rho = '... 
                   num2str(obj.simulationParam.rho)]);
           end
           
           if saveFig
               name = ['varVarUPvsSteps_T' my_num2str(obj.simulationParam.T)...
                   'N' my_num2str(obj.simulationParam.N) 'rho'...
                   my_num2str(obj.simulationParam.rho)];
               saveas(gcf,[name add2savedFigName '.fig']);
               saveas(gcf,[name add2savedFigName '.jpg']);
           end
           
           if ~keepFigOpen
               close all;
           end
       end
        
       function obj = calcPressureAfterRun(obj)
           N = obj.simulationParam.N;
           rho = obj.simulationParam.rho;
           rCutoff = obj.simulationParam.rCutoff;
           T = obj.simulationParam.T;
           Pfix = 3*pi*rho^2/rCutoff^4;
           
           for i = 1:obj.indIndata
               dists = obj.data.allDists(1:N,1:N,i);
               virial = calcVirial(dists,rho,12,6,N,rCutoff);
               obj.data.allV(1,i) = virial;
               P = T*rho + virial;
               obj.data.allP(1,i) = P;
               obj.data.allPlrc(1,i) = P - Pfix;
           end
           
           obj.currentPressure = P;
           obj.currentVir = virial;
           obj.Plrc = mean(obj.data.allPlrc(1,1:obj.indIndata));
       end
    
       function obj = calcCv(obj)
           % calculate the energy flucuations:
           firstSteps2ignore = obj.data.firstSteps2ignore;
           N = obj.simulationParam.N;
           T = obj.simulationParam.T;
           endAt = obj.indIndata;
           
           meanU = obj.data.meanUlrcEq;
           U = obj.data.allUlrc(1,firstSteps2ignore:endAt);
           sigU = mean((U - meanU).^2);

           obj.data.sigU = sigU;
           obj.data.cv = N + (sigU/T^2);
       end
       
       function obj = deleteFirstRuns(obj,firstRuns2Delete)
           % delete from data some of the first runs
           if obj.simulationParam.angleDependent
%                part1 = floor(obj.data.indIndata/2);
%                obj.data.a1 = obj.data.allAlphas(:,:,(firstRuns2Delete+1):part1);
%                obj.data.a2 = obj.data.allAlphas(:,:,(part1+1):obj.data.indIndata);
%                obj.data.allAlphas = obj.data.a1;
%                obj.data.allAlpha(:,:,(part1+1):(obj.data.indIndata - firstRuns2Delete));
%                obj.data.a1 = [];
%                obj.data.a2 = [];
               
               a = obj.data.allAngs(:,:,(firstRuns2Delete+1):obj.data.indIndata);
               obj.data.allAngs = a;
               clear a;
               
%                part1 = floor(obj.data.indIndata/2);
%                obj.data.a1 = obj.data.allThetas(:,:,(firstRuns2Delete+1):part1);
%                obj.data.a2 = obj.data.allThetas(:,:,(part1+1):obj.data.indIndata);
%                obj.data.allThetas = obj.data.a1;
%                obj.data.allThetas(:,:,(part1+1):(obj.data.indIndata - firstRuns2Delete));
%                obj.data.a1 = [];
%                obj.data.a2 = [];
               
           end
           
           a = obj.data.allCoords(:,:,(firstRuns2Delete+1):obj.data.indIndata);
           obj.data.allCoords = a;
           clear a;
           
%            part1 = floor(obj.data.indIndata/2);
%            obj.data.a1 = obj.data.allDists(:,:,(firstRuns2Delete+1):part1);
%            obj.data.a2 = obj.data.allDists(:,:,(part1+1):obj.data.indIndata);
%            obj.data.allDists = obj.data.a1;
%            obj.data.allDists(:,:,(part1+1):(obj.data.indIndata - firstRuns2Delete));
%            obj.data.a1 = [];
%            obj.data.a2 = [];
           
           a = obj.data.allU(:,(firstRuns2Delete+1):obj.data.indIndata);
           obj.data.allU = a;
           clear a;
           a = obj.data.sweepInd(:,(firstRuns2Delete+1):obj.data.indIndata);
           obj.data.sweepInd = a;
           clear a;
           
           
           if ~isempty(obj.data.allP)
                a = obj.data.allP(:,(firstRuns2Delete+1):obj.data.indIndata);
                obj.data.allP = a;
                clear a;
           end
           
           if ~isempty(obj.data.allUlrc)
                a = obj.data.allUlrc(:,(firstRuns2Delete+1):obj.data.indIndata);
                obj.data.allUlrc = a;
                clear a;
           end
           
           
           if ~isempty(obj.data.allPlrc)
                a = obj.data.allPlrc(:,(firstRuns2Delete+1):obj.data.indIndata);
                obj.data.allPlrc = a;
                clear a;
           end
           
           obj.data.indIndata = obj.data.indIndata - firstRuns2Delete;
           obj.indIndata = obj.data.indIndata;
           obj.data.firstRuns2Delete = firstRuns2Delete;
           
       end
       
       function objCopy = copyOutputFile(obj, copyNum)
           % copy MC output flie. adds the string 'copy' and the number
           % copyNum to the file name of the copied object.
                      
           N = obj.simulationParam.N;
           endAt = obj.indIndata;
           
           % create a new data file
           simulationParam = obj.simulationParam;
           dataCopyStr = [obj.fileName 'copy' copyNum '.mat'];
           save(dataCopyStr, 'simulationParam','-v7.3');
           newdata = matfile(dataCopyStr,'Writable',true);
           
           %add the last step of obj to the new data file
           newdata.allCoords = zeros(2,N,2);
           newdata.allCoords(1:2,1:N,1) = obj.data.allCoords(1:2,1:N,endAt);
           [~, ~, s] = size(obj.data.allDists);
           newdata.allDists = zeros(N,N,2);
           newdata.allDists(1:N,1:N,1) = obj.data.allDists(1:N,1:N,s);
           newdata.allU = zeros(1,2);
           newdata.allU(1,1) = obj.data.allU(1,endAt);
           newdata.currentmaxdr = obj.data.currentmaxdr;
           newdata.moveCount = 0;
           newdata.sweepInd = zeros(1,2);
           newdata.sweepInd(1,1) = obj.data.sweepInd(1,endAt);
           newdata.indIndata = 1;
           
           if existInMatfile(obj.data,'allP')
               newdata.allP = zeros(1,2);
               newdata.allP(1,1) = obj.data.allP(1,endAt);
           end
           
           if existInMatfile(obj.data,'allPlrc')
               newdata.allPlrc = zeros(1,2);
               newdata.allPlrc(1,1) = obj.data.allPlrc(1,endAt);
           end
           
           if existInMatfile(obj.data,'allUlrc')
               newdata.allUlrc = zeros(1,2);
               newdata.allUlrc(1,1) = obj.data.allUlrc(1,endAt);
           end
           
           if existInMatfile(obj.data,'allV')
               newdata.allV = zeros(1,2);
               newdata.allV(1,1) = obj.data.allV(1,endAt);
           end
           
           
           if existInMatfile(obj.data,'RDFbins')
               newdata.RDFbins = obj.data.RDFbins;
               [~, s, ~] = size(obj.data.RDFhisto);
               newdata.RDFhisto = zeros(1,s,2);
               newdata.RDFhisto(1,1:s,1) = obj.data.RDFhisto(1,1:s,endAt);
           end
           
           if existInMatfile(obj.data,'allAlphas')
               newdata.allAlphas = zeros(N,N,2);
               newdata.allAlphas(1:N,1:N,1) = obj.data.allAlphas(1:N,1:N,endAt);
               newdata.allAngs = zeros(1,N,2);
               newdata.allAngs(1,1:N,1) = obj.data.allAngs(1,1:N,endAt);
               newdata.allThetas = zeros(N,N,2);
               newdata.allThetas(1:N,1:N,1) = obj.data.allThetas(1:N,1:N,endAt);
           end
           
           % make a copy of the object
           objCopy = MC2DLJoutput(dataCopyStr);
           
           
               
       end
       
       function [obj, varUfromRDF] =...
               getVarUfromRDF(obj, varargin)
       % Integrate over the RDF to get the energy in each step, than get
       % the std of the energy
       
            p = inputParser();
            addOptional(p, 'plotFig', false);
            addOptional(p, 'freen', false);
            addOptional(p, 'freeTandn', false);
            addOptional(p, 'freeTnbound', false);
            addOptional(p, 'freeTnset', false);
            addOptional(p, 'nset', []);
            addOptional(p, 'use_m_n_T', []);
            addOptional(p, 'talk', false);
            parse(p, varargin{:});
            Results = p.Results;
            plotFig = Results.plotFig;
            freen = Results.freen;
            freeTandn = Results.freeTandn;
            freeTnbound = Results.freeTnbound;
            freeTnset = Results.freeTnset;
            nset = Results.nset; 
            use_m_n_T = Results.use_m_n_T;
            talk = Results.talk;
            
            [obj, UfromRDF] = getUfromRDF(obj, 'plotFig'...
                ,plotFig,'freen',freen,'freeTandn',freeTandn,...
                'freeTnbound',freeTnbound,'freeTnset',freeTnset,...
                'nset',nset,'use_m_n_T',use_m_n_T);
            
            % Get U for every step
            [~ , numOfRDFbins] = size(obj.data.RDFhisto);
            inc = obj.data.RDFbins(1,2) - obj.data.RDFbins(1,1); 
            obj.data.UfromRDFinStep = zeros(1,obj.indIndata);
            x = obj.data.RDFbins;
            rho = obj.simulationParam.rho;
            for ii = 1:obj.indIndata
                if talk
                    tic;
                    disp(['calc' num2str(ii) ' of ' num2str(obj.indIndata)]);
                end
                obj.data.UfromRDFinStep(1,ii) =...
                    rho*pi*sum(UfromRDF(2:numOfRDFbins).*obj.data.RDFhisto(ii,2:numOfRDFbins).*x(2:numOfRDFbins)*inc);
                if talk
                    toc;
                end
            end
            
            varUfromRDF = var(obj.data.UfromRDFinStep);
       end
       
       function [obj, UfromRDF] = getUfromRDF(obj, varargin)
       % Get the pair potantial from the RDF
       
        p = inputParser();
        addOptional(p, 'plotFig', false);
        addOptional(p, 'freen', false);
        addOptional(p, 'freeTandn', false);
        addOptional(p, 'freeTnbound', false);
        addOptional(p, 'freeTnset', false);
        addOptional(p, 'nset', []);
        addOptional(p, 'use_m_n_T', []);
        parse(p, varargin{:});
        Results = p.Results;
        plotFig = Results.plotFig;
        freen = Results.freen;
        freeTandn = Results.freeTandn;
        freeTnbound = Results.freeTnbound;
        freeTnset = Results.freeTnset;
        nset = Results.nset; 
        use_m_n_T = Results.use_m_n_T;
        
        % Check if we need to calculate a fit for the RDF
        fitDone = false;
        if use_m_n_T
            minput = use_m_n_T(1);
            ninput = use_m_n_T(2);
            Tinput = use_m_n_T(3);
        else
            if freen
                if existInMatfile(obj.fileName,'fitObjlogRDFfreen')
                    fitDone = true;
                end
            else
                if freeTandn
                    if existInMatfile(obj.fileName,'fitObjlogRDFfreeTandn')
                        fitDone = true;
                    end
                else
                    if freeTnbound
                        if existInMatfile(obj.fileName,...
                                'fitObjlogRDFfreeTnbound')
                            fitDone = true;
                        end
                    else
                        if freeTnset
                            if existInMatfile(obj.fileName,...
                                    'fitObjlogRDFfreeTnset')
                                fitDone = true;
                            end
                        end
                    end
                end
            end
        end
        
        if ~fitDone
            [obj, fitresult, mfit, mError, nfit, nError,...
                Tfit, TError, gof] = pairUfromRDF(obj, 'plotFig'...
                ,plotFig,'freen',freen,'freeTandn',freeTandn,...
                'freeTnbound',freeTnbound,'freeTnset',freeTnset,...
                'nset',nset);
        end
        
        % Find U from the fit results
        UfromRDF = zeros(1, obj.indIndata);
        x = obj.data.RDFbins;
        T = obj.simulationParam.T;
        if use_m_n_T
            UfromRDF = 4*(1/Tinput)*((x.^-ninput) - (x.^-minput));
            obj.data.UfromRDFnmTset = UfromRDF;            
        else
            if freen
                fitObjlogRDFfreen = obj.data.fitObjlogRDFfreen;
                nfit = fitObjlogRDFfreen.nfit;
                mfit = fitObjlogRDFfreen.mfit;
                UfromRDF = 4*(1/T)*((x.^-nfit) - (x.^-mfit));
                obj.data.UfromRDFnfree = UfromRDF;
            else
                if freeTandn
                    fitObjlogRDFfreeTandn = obj.data.fitObjlogRDFfreeTandn;
                    nfit = fitObjlogRDFfreeTandn.nfit;
                    mfit = fitObjlogRDFfreeTandn.mfit;
                    Tfit = fitObjlogRDFfreeTandn.Tfit;
                    UfromRDF = 4*(1/Tfit)*((x.^-nfit) - (x.^-mfit));
                    obj.data.UfromRDFfreeTandn = UfromRDF;
                else
                    if freeTnbound
                        fitObjlogRDFfreeTnbound =...
                            obj.data.fitObjlogRDFfreeTnbound;
                        nfit = fitObjlogRDFfreeTnbound.nfit;
                        mfit = fitObjlogRDFfreeTnbound.mfit;
                        Tfit = fitObjlogRDFfreeTnbound.Tfit;
                        UfromRDF = 4*(1/Tfit)*((x.^-nfit) - (x.^-mfit));
                        obj.data.UfromRDFfreeTnbound = UfromRDF;
                    else
                        if freeTnset
                            fitObjlogRDFfreeTnset =...
                                obj.data.fitObjlogRDFfreeTnset;

                            mfit = fitObjlogRDFfreeTnset.mfit;
                            Tfit = fitObjlogRDFfreeTnset.Tfit;
                            UfromRDF = 4*(1/Tfit)*((x.^-nset) - (x.^-mfit));
                            obj.data.UfromRDFfreeTnset = UfromRDF;
                        end
                    end
                end
            end
        end
        
           
       end
       
       function [obj, fitresult, mfit, mError, nfit, nError,...
               Tfit, TError, gof] = pairUfromRDF(obj, varargin)
       % Fit -log(RDF) to the function funForFit
       
        p = inputParser();
        addOptional(p, 'plotFig', false);
        addOptional(p, 'freen', false);
        addOptional(p, 'freeTandn', false);
        addOptional(p, 'freeTnbound', false);
        addOptional(p, 'freeTnset', false);
        addOptional(p, 'nset', []);
        addOptional(p, 'numOfBinsRDF', 300);
        addOptional(p, 'rCutoffRDF', 10);
        parse(p, varargin{:});
        Results = p.Results;
        plotFig = Results.plotFig;
        freen = Results.freen;
        freeTandn = Results.freeTandn;
        freeTnbound = Results.freeTnbound;
        freeTnset = Results.freeTnset;
        nset = Results.nset; 
        numOfBinsRDF = Results.numOfBinsRDF;
        rCutoffRDF = Results.rCutoffRDF;
        
        
        try
            x = obj.data.RDFbins;
        catch
            obj = obj.calcRDF(rCutoffRDF, numOfBinsRDF);
            x = obj.data.RDFbins;
        end
        
        steps{1,1} = ['steps: ' num2str(obj.indIndata)]; 
        
       [fitresult, mfit, mError, nfit, nError, Tfit, TError, gof] =...
            createFitRDF(x, -log(mean(obj.data.RDFhisto)),...
            obj.simulationParam.T, obj.simulationParam.rho,...
            obj.simulationParam.m , steps, 'plotFig',plotFig,...
            'freen',freen,'freeTandn',freeTandn,...
            'freeTnbound',freeTnbound,'freeTnset',freeTnset,'nset',nset);
        
        fitObj.fitresult = fitresult;
        fitObj.mfit = mfit;
        fitObj.mError = mError;
        fitObj.nfit = nfit;
        fitObj.nError = nError;
        fitObj.Tfit = Tfit;
        fitObj.TError = TError;
        fitObj.gof = gof;
        
        % Save fit to data
        if freen
            obj.data.fitObjlogRDFfreen = fitObj;
        else
            if freeTandn
                obj.data.fitObjlogRDFfreeTandn = fitObj;
            else
                if freeTnbound
                    obj.data.fitObjlogRDFfreeTnbound = fitObj;
                else
                    if freeTnset
                        obj.data.fitObjlogRDFfreeTnset = fitObj;
                    end
                end
            end
        end
        
        
       end
       
        
    end    

end


function [dist,particlesPosition] = ...
    createInitialConfig(L,N,r,initialConfig)

    possibleInitialConfigs = {'random','hex','auto','2par_close'};
    initialConfigInd = strcmp(initialConfig,possibleInitialConfigs);
    % check if input is valid:
    if sum(initialConfigInd) ~= 1
        error(['choose one of the initial configurations: '...
            my_cell2str(possibleInitialConfigs)]);
    else
        switch find(initialConfigInd)
            case 1 % random initial configuration
                [dist,particlesPosition] = randomStart(L,N,r);
            case 2 % hexagonal initial configuration
                [dist,particlesPosition] = hcp(L,N,r);
            case 3 % choose automatically
                if N/L^2 > 0.4
                    [dist,particlesPosition] = hcp(L,N,r);
                else
                    [dist,particlesPosition] = randomStart(L,N,r);
                end
            case 4 % two particle close to each other
                dist = [0,0;2,0];
                particlesPosition = [0,0;2,0];
        end
    end
end

function [dist,particlesPosition] = randomStart(L,N,r)
        % randomize first particle possition in the box 
        % [-L/2,L/2] x [-L/2,L/2]
        particlesPosition(1,1) = L*rand - (L/2);
        particlesPosition(2,1) = L*rand - (L/2);
        dist = zeros(N);

        for j = 2:N

              % choose random possition
              particlesPosition(1,j) = L*rand - (L/2);
              particlesPosition(2,j) = L*rand - (L/2);

              % calculate PBC distances
              xj = particlesPosition(1,j);
              yj = particlesPosition(2,j);
              dist(j,1:j) = distPBC(xj,yj,particlesPosition,L);

              % check for piriodic boundary condition overlaps,
              % randomize new possition if overlaps are found.
              overlapPBC = sum(dist(j,1:(j-1)) < 2*r) > 0;
              countTry = 0;
              while overlapPBC
                      particlesPosition(1,j) = L*rand - (L/2);
                      particlesPosition(2,j) = L*rand - (L/2);

                      % calculate PBC distances
                      xj = particlesPosition(1,j);
                      yj = particlesPosition(2,j);
                      dist(j,1:j) = distPBC(xj,yj,particlesPosition,L);

                      overlapPBC = sum(dist(j,1:(j-1)) < 2*r) > 0;
                      countTry = countTry + 1;
                      if countTry > 1000
                          error(['it is difficult to generate a'...
                              'random distribution with a density '...
                              'of ' num2str(N/L^2) '. try a lower '...
                              'density, or a hexagonal initial'...
                              ' configuration.']);
                      end
              end

        end

end

function [dist,particlesPosition] = hcp(L,N,r)

        % Make sure N is a perfect square
        intRoot = floor(sqrt(N));
        if (sqrt(N) - intRoot) > 1e-7
            % Display an error message
            disp('Number of particles should be a perfect square');
            particlesPosition = [];
            return
        end

        % Calculate the seperation length between particles centers
        sepDist = L/sqrt(N);

        % Make sure the density is not too high
%         if sepDist < 2*r
%             % Display an error message
%             disp('density is too high');
%             particlesPosition = [];
%             return
%         end

        % Find the box size
        Lx = sepDist * sqrt(N);

        % Create a vector of linearly spaced points along the
        % x-direction
        xPos = linspace(sepDist/2, Lx-sepDist/2, sqrt(N));
        % And find the corresponsing y-direction increments
        yPos = xPos;

        % Create a matrix with all combinations of x and y
        [X,Y] = meshgrid(xPos,yPos);
        % Shift coordinates to the be at the center of each
        % particle
        X(1:2:end,:) = X(1:2:end,:) + sepDist/2;

        % Reshape the matrix to be 1D in X and 1D in Y
        % (numel returns the number of elements in a given array)
        particlesPosition =...
            [reshape(X,1,numel(X));reshape(Y,1,numel(Y))];

        % make the board in: [-L/2 L/2]x[-L/2 L/2]
        particlesPosition = particlesPosition - L/2;

                % calculate all pair distances
        dist = zeros(N);
        for par = 1:N

            x = particlesPosition(1,par);
            y = particlesPosition(2,par);
            dist((par+1):N,par) = ...
                distPBC(x,y,particlesPosition(:,(par+1):N),L);
        end


end

function d = reshapeDist(allDists)
    [N,~] = size(allDists);
    d = [];
    for i=1:N
        d = [d ; allDists((i+1):N,i)];
    end
end



function [ang1, ang2] = reshapeAng(allAngs)
    N = length(allAngs);
    ang1 = [];
    ang2 = [];
    for i = 1:(N-1)
        ang1 = [ang1 ones(1,N-i)*allAngs(i)]; 
        ang2 = [ang2 allAngs(1,(i+1):N)];
    end
    ang1 = ang1';
    ang2 = ang2';
end

function [ang1, ang2] = reshapeAng(allAngs)
    N = length(allAngs);
    ang1 = [];
    ang2 = [];
    for i = 1:(N-1)
        ang1 = [ang1 ones(1,N-i)*allAngs(i)]; 
        ang2 = [ang2 allAngs(1,(i+1):N)];
    end
end

function meanProp = my_mean(prop)
    % input: some property of the simulation in all steps.
    % output: meanProp(i) is the mean of all the values of property on
    % steps 1 to i.

        lenProp = length(prop);
        meanProp = zeros(1,lenProp); 
        meanProp(1) = prop(1);
        for i = 2:lenProp
            meanProp(i) = meanProp(i-1) + prop(i);
        end
        one2len = 1:lenProp;
        meanProp = meanProp./one2len;
end  

function dist = calcDists(particlesPosition,L)
    % calculate distance matrix from coordinates
        [~,N] = size(particlesPosition);
        dist = zeros(N);

        for j = 2:N

              % calculate PBC distances
              xj = particlesPosition(1,j);
              yj = particlesPosition(2,j);
              dist(j,1:(j-1)) =...
                  distPBC(xj,yj,particlesPosition(:,1:(j-1)),L);

        end

end

function [bins, histo] =...
    calculateRDF(dists,N,rho,numOfBins,maxDist)
%% calculate 2D radial distribution function, with pediodic boundary
%% conditions.

%         input:
%         ~~~~~~

%         NumOfBins - number of bins in the histogram
%
%         output:
%         ~~~~~~~
%         obj.data.bins - x axis of the RDF histogram
%         obj.data.histo - each colmun is the RDF axis of a diffrent step in the
%         monte carlo simulation. (so to plot the RDF for the 10'th step we
%         need to write: plot(bins,histo(:,10));

%       method:
%       ~~~~~~~
%       the RDF is calculated by binnig all pair partical distances into 
%       a histogram, and normalizing each bin with it's Ideal gas number of
%       particals. 
%       when binning the pair distances, we take into acount Periodic
%       Boudary Conditions (PBC)
%       finaly, to ansure that when r->infinity : RDF ->1 , we
%       normalize by multiplying RDF*2/(N-1) where N is the number of
%       particals. 
%       for more information
%   http://www2.msm.ctw.utwente.nl/sluding/TEACHING/APiE_Script_v2011.pdf
%       page 48 - "Radial distribution function"


bins = linspace(0,maxDist,numOfBins);  

d = reshape(dists,1,[]);
d = nonzeros(d);
d = d(d < maxDist);

histo = hist(d,bins);
increment = bins(2) - bins(1);

% each bin should be normalized according to its volume
for bin = 1:numOfBins

        % histo(bin) is the number of particles in some layer of area 
        % 2pi*rVal*dr, a distance rVal from the central cell
        rVal = bins(bin);
        next_rVal = increment + rVal;

        % Calculate the area of the bin (a ring of radii r,
        % r+dr)
        ereaBin = pi*next_rVal^2 - pi*rVal^2;

        % Calculate the number of particles expected in this bin in
        % the ideal case
        nIdeal = ereaBin*rho;

        % Normalize the bin
        histo(bin) =...
            histo(bin) / nIdeal;

end
histo = 2*histo/(N-1);


end


function Ni = numOfCellsDistribution(coords,L,numOfSquares)
%% Find densities in a montecarlo step

% Given coordinates of N particles in 2D box, we devide the box to numOfSquares
% squares and build a histogram from the densities in each square.

% Input:
% ~~~~~
% coords       : Coordinates in 2D, 2 by N matrix. the center of the box
%                must be at (0,0);
% L            : Size of box side.
% numOfSquares : Number of squares we devide the board into. the root of
%                this number must be real
% numOfBins    : Number of bins in the histogram

% Output:
% ~~~~~
% Ni    : all number of cells in subsystems 

% Parse input:
[~, N] = size(coords); 
m = sqrt(numOfSquares);
if isinteger(m)
    error('The square root of numOfSquares must be an intiger');
end


% Get the density in each square  
squareArea = L^2/numOfSquares;
squareSide = sqrt(squareArea);
Ni = zeros(m,m);
indi = 0;

for i = (-L/2):squareSide:(L/2 - squareSide)
    indi = indi + 1;
    indj = 0;
    for j = (-L/2):squareSide:(L/2 - squareSide)
        indj = indj + 1;
        Ni(indi,indj) =...
            countParticlesInSquare(i,j,coords(1,:),coords(2,:),squareSide);
    end
end

%densities = densities/squareArea;
Ni = reshape(Ni,[1,m^2]);

%numOfCellsInSquare = densities;

% Bin the results to a histogram
%rho = linspace(min(densities),max(densities),numOfBins);
%[PL, rho] = hist(densities,numOfBins);

% Normalize with the hompgeneus density
%rho = rho / (N/L^2);


function N = countParticlesInSquare(i,j,x,y,squareSide)
% Count the number of particles in square i,j
    N = sum(and(and(y > j,y < j+squareSide), and(x > i,x < i+squareSide)));
end
end  
        
