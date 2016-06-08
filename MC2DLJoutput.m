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
% currentStep - how many steps calculated so far
% indIndata - the number of steps saved so far (not every step is saved,
%               so this may be lower than currentStep)
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
%       data.stepInd - the number of step clulated up to each saved point
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
          currentVir,Ulrc,Plrc,currentStep,...
          currentThetas,currentAlphas,currentAngs;
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
                    obj.currentDists = obj.data.allDists(:,:,obj.indIndata);
                    obj.currentU = obj.data.allU(1,obj.indIndata);
                    
                    if ~isfield(obj.simulationParam,'angleDependent')
                        obj.simulationParam.angleDependent = false;
                    end
                    
                    if obj.simulationParam.angleDependent
                   
                        obj.currentAngs = obj.data.allAngs(:,:,obj.indIndata);
                        obj.currentAlphas = obj.data.allAlphas(:,:,obj.indIndata);
                        obj.currentThetas = obj.data.allThetas(:,:,obj.indIndata);
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
                    
                    obj.currentStep = obj.data.stepInd(1,obj.indIndata);

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
                        allAlphas = zeros(N,N,1);
                        allThetas = zeros(N,N,1);
                        allAngs(:,:,1) = rand(1,N,1)*pi;
                        allAlphas(:,:,1) =...
                            tril(bsxfun(@minus,allAngs(:,:,1),allAngs(:,:,1)'),-1);

                         xdist = bsxfun(@minus,allCoords(1,:,1),allCoords(1,:,1)');
                         ydist = bsxfun(@minus,allCoords(2,:,1),allCoords(2,:,1)');
                        allThetas(:,:,1) = ...
                            tril(allAlphas(:,:,1)-atan(ydist./xdist),-1);
                    else
                        allAngs = [];
                        allAlphas = [];
                        allThetas = [];
                        
                    end

                        
                    % calculate initial energy
                    d = reshapeDist(allDists);
                    ang = [reshapeDist(allAlphas),...
                            reshapeDist(allThetas)];
                    allU = pairU(d,rCutoff,m,...
                        'angleDependence',angleDependence,...
                        'relativeCellAngles',ang,'ufunc',ufunc);
                    
                    %%% this is the long range correction - should be fixed
                    %%% for the case of angle dependence
                    allUlrc = allU - pi*rho*N/rCutoff^4;
                    
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
                    stepInd = 0;
                    moveCount = 0;
                    indIndata = 1;
                    currentmaxdr = initialmaxdr;
                    simulationParam = obj.simulationParam;
                    save(obj.fileName, 'allDists','allCoords',...
                            'allU','allUlrc','allV','allP','allPlrc','stepInd',...
                            'moveCount','indIndata'...
                            ,'currentmaxdr','simulationParam','runNum',...
                            'allAlphas','allThetas',...
                            'allAngs','-v7.3');
                    obj.data = matfile(obj.fileName);
                    obj.data = matfile(obj.fileName,'Writable',true);
                        
                    obj.currentCoords = allCoords;
                    obj.currentDists = allDists;
                    obj.currentAngs = allAngs;
                    obj.currentAlphas = allAlphas;
                    obj.currentThetas = allThetas;
                    obj.currentU = allU;
                    obj.Ulrc = allUlrc;
                    obj.currentVir = allV;
                    obj.currentPressure = allP;
                    obj.Plrc = allPlrc;
                    obj.currentStep = stepInd;
                    obj.indIndata = 1;
                end
            
        end
        
       function obj = MonteCarlo(obj,Nsteps,saveEvery)
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
            else 
                hardCoreRepRad = 0;
            end
            
            stepCount = 0;
            while(stepCount < Nsteps)
                [finalU,finalV,finalPressure,...
                    finalConfiguration,finalDistances,...
                    currentmoveCount,...
                    finalAngs,finalAlphas,finalThetas] = ...
                    MonteCarlo2DLJHeart(...
                    N,...
                    obj.simulationParam.T,...
                    rho,...
                    obj.simulationParam.N*saveEvery,...
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
                    'initialAlphas',obj.currentAlphas,...
                    'initialThetas',obj.currentThetas,...
                    'maxdAng',maxdAng,...
                    'ufunc',ufunc,...
                    'hardCoreRepRad',hardCoreRepRad);
                
                stepCount = stepCount + obj.simulationParam.N*saveEvery;
                obj.currentStep = obj.currentStep...
                    + obj.simulationParam.N*saveEvery;
                obj.moveCount = obj.moveCount + currentmoveCount;
                obj.currentU = finalU;
                obj.currentVir = finalV;
                obj.currentPressure = finalPressure;
                obj.currentCoords = finalConfiguration;
                obj.currentDists = finalDistances;
                obj.currentAngs = finalAngs;
                obj.currentAlphas = finalAlphas;
                obj.currentThetas = finalThetas;
            
                % long range correction
                obj.Ulrc = obj.currentU + 4*pi*N*rho/((2-m)*rCutoff^(m-2));
                if obj.simulationParam.pressure
                    obj.Plrc = obj.currentPressure -...
                        2*m*pi*rho^2/((2-m)*rCutoff^(m-2));
                end
            
                obj = obj.addStep2data(obj.currentStep,finalConfiguration,...
                    finalDistances,finalU,finalV,finalPressure,...
                    obj.moveCount,obj.currentmaxdr,obj.Ulrc,obj.Plrc,...
                    angleDependent,obj.currentAngs,...
                    obj.currentAlphas,obj.currentThetas);
                
                clear finalU finalV finalConfiguration finalDistances...
                    currentmoveCount finalAngs finalAlphas finalThetas
            end
            
        end
        
       function obj = addStep2data(obj,newInd,newCoords,newDists,newU,newV...
                ,newP,moveCount,currentmaxdr,newUlrc,newpressurelrc,...
                angleDependent,newAngs,newAlphas,newThetas)
            
            obj.data.stepInd(1,obj.indIndata+1) = newInd;
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
            else
                obj.data.allCoords(:,:,obj.indIndata+1) = newCoords(:,:,1);
                obj.data.allDists(:,:,obj.indIndata+1) = newDists(:,:,1);
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
                    d = reshape(dist,1,[]);
                    d = nonzeros(d);
                    d = d(d < maxDist);
                    
                    histo = hist(d,bins);
                    increment = bins(2) - bins(1);
                   

                    % each bin should be normalized according to its volume
                    for bin = 1:numOfBins

                            % rVal is the number of particles in some layer of area 
                            % da(r)=2pi*r*dr, a distance r from the central cell
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
           addOptional(p, 'keepFigOpen', false);
           addOptional(p, 'firstSteps2ignore', 0);
           parse(p, varargin{:});
           Results = p.Results;
           plotFig = Results.plotFig;
           saveFig = Results.saveFig;
           keepFigOpen = Results.keepFigOpen;
           firstSteps2ignore = Results.firstSteps2ignore;
           
           % check if RDF was calculated allready
           if ~existInMatfile(obj.data,'RDFhisto')
               error('calculate RDF first');
           end
           
           [Nsteps, ~] = size(obj.data.RDFhisto);
           
           if firstSteps2ignore == 0
               minInd = 1;
           else
               % find the relevent initial index
               [~, minInd] = min(abs(firstSteps2ignore - obj.data.stepInd));
           end
           
           % find peaks
%            [obj, locs] = obj.findRDFpeaks();
            locs = 30;
            peaks = zeros(length(locs),Nsteps-minInd+1);
            steps = zeros(length(locs),Nsteps-minInd+1);
            

           for i = minInd:Nsteps
               for k = 1:length(locs)
                   peaks(k,i-minInd+1) =...
                       obj.data.RDFhisto(i,locs(k));
               end
               
               % find variance and variance of variance in the peaks
               for j = 1:length(locs)
                    varPeaks(j,i-minInd+1) = var(peaks(j,1:(i-minInd+1)));
                    varVarPeaks(j,i-minInd+1) =...
                        var(varPeaks(j,1:(i-minInd+1)));
               end
           end
           
           for j = 1:length(locs)
               steps(j,:) = obj.data.stepInd(1,minInd:Nsteps);
           end
           
           if plotFig
               
               for j = 1:length(locs)
                   leg{1,j} = ['peak distance '...
                       num2str(obj.RDFbins(locs(j)))]; 
               end
               
               colorPlot(steps,varPeaks,'addLegend',leg);
               xlabel('steps');
               ylabel('RDF peaks variance');
               
               if saveFig
                    name = ['varRDFpeaksVsSteps_T'...
                    my_num2str(obj.simulationParam.T)...
                    'N' my_num2str(obj.simulationParam.N) 'rho'...
                    my_num2str(obj.simulationParam.rho)];
                    saveas(gcf,[name '.fig']);
                    saveas(gcf,[name '.jpg']);
               end
             
               colorPlot(steps,varVarPeaks,'addLegend',leg);
               xlabel('steps');
               ylabel('RDF peaks variance of variance');
               
               if saveFig
                    name = ['varVarRDFpeaksVsSteps_T'...
                    my_num2str(obj.simulationParam.T)...
                    'N' my_num2str(obj.simulationParam.N) 'rho'...
                    my_num2str(obj.simulationParam.rho)];
                    saveas(gcf,[name '.fig']);
                    saveas(gcf,[name '.jpg']);
               end
               
               if ~keepFigOpen
                   close all;
               end
           end
               
       end

       function [obj, rhoNorm, PL] =...
               calcRhoDistrib(obj,squares,numOfbins,varargin)
       % try:PL is the prob to get a distribution in one squre, considering all
       % steps. we avg on all subsys


       % breaks the board to "squares" pieces, and finds the density in each piece.

            % input:
            % allCoords is a matrix of all coordinates in each step of a monte carlo
            % simulation (a '2' x 'nuber of particles' x 'number of steps' matrix).
            % N must have a natural root.
            % L is the size of the board

            % output:
            % the results is a histogram of densities. rho is a 1 x bins matrix.
            % rhoNorm in the x axis (the bins) of the histogram. Its the densities,
            % normalized with the total desity on the board.
            % PL is the probability of finding a block with some density (the y axis of
            % the histogram.
            % hist(j,i) is the number of particles in block i, in step j of the
            % simulation (the histogram for the coordinates: allCoords(:,:,j)

            % this function uses the function 'my_hist'. get it in this repository: 
            % https://github.com/adirot/general-matlab-functions.git

            p = inputParser();
            addOptional(p, 'calcPL', true);
            addOptional(p, 'startFrom', 1);
            addOptional(p, 'save2data', false);
            parse(p, varargin{:});
            Results = p.Results;
            calcPL = Results.calcPL;
            startFrom = Results.startFrom;
            save2data = Results.save2data;
            
            if save2data
                rhoDistribParam.squares = squares;
                rhoDistribParam.startFrom = startFrom;
                obj.data.rhoDistribParam = rhoDistribParam;  
                obj.simulationParam.rhoDistribParam = rhoDistribParam;
            end
            
            [~,N,numberOfSteps] = size(obj.data.allCoords);
            row = sqrt(squares);
            L = obj.simulationParam.L;
            inc = L/row;
            boxSize = inc^2;
            %hist = zeros(numberOfSteps,bins);
            hist = zeros(squares,numOfbins);
            allRho = zeros(numberOfSteps,squares);

            for step = 1:numberOfSteps
                n = 1;
                for x = 1:row
                    for y = 1:row
            %                 plot(allCoords(1,:,step),allCoords(2,:,step),'+');
            %              hold on; xlim([-25 25]);ylim([-25 25]);
                        % find the number of particles in the box
                        inBoxX = and((obj.data.allCoords(1,:,step) < -L/2 ...
                        + inc*x), ...
                            (obj.data.allCoords(1,:,step) > -L/2 ...
                        + inc*(x-1)));
            %             plot([-L/2 + inc*x ,-L/2 + inc*x] ,[-25 25],'r');
            %             plot([-L/2 + inc*(x-1), -L/2 + inc*(x-1)] ,[-25 25],'r');
                        inBoxY = and((obj.data.allCoords(2,:,step) < -L/2 ...
                        + inc*y) ,...
                            (obj.data.allCoords(2,:,step) > -L/2 ...
                            + inc*(y-1)));
            %             plot( [-25 25], [-L/2 + inc*y ,-L/2 + inc*y],'r');
            %             plot( [-25 25],[-L/2 + inc*(y-1), -L/2 + inc*(y-1)],'r');
                        particlesInBox = sum(inBoxX.*inBoxY);

                        allRho(step,n) = particlesInBox/boxSize;
                        n = n + 1;
                        %hold off;
                    end
                end
            end

                % create a histogram from all distributions.
                rho = linspace(0,max(max(allRho)),numOfbins);

            %     for step = 1:numberOfSteps
            %         hist(step,:) = my_hist(allRho(step,:),rho);
            %     end

                for n = 1:squares
                    hist(n,:)...
                     = my_hist(allRho(startFrom:numberOfSteps,n)',rho);
                end

                rhoNorm = rho/(N/L^2);
                if save2data
                    obj.data.rhoNorm = rhoNorm;
                end

                if calcPL
                    %PL = mean(hist)/sq^2;
                    PL = mean(hist/(numberOfSteps-startFrom+1));
                    if save2data
                        obj.data.PL = PL;
                    end
                else
                    PL = [];
                    if save2data
                        obj.data.PL = [];
                    end
                end

       end

       
       function showStep(obj,step)
           
          if isnumeric(step)
                  ind = find(obj.data.stepInd == step);
          
                  if isempty(ind) %find the closest index saved
                      [~,ind] = min(abs(obj.data.stepInd - step));
                      step = obj.data.stepInd(1,ind);
                  end
          else
              
              switch step
                  case 'first'
                      ind = 1;
                      step = 0;
                  case 'mid'
                     
                      ind = ceil(obj.indIndata/2);
                      step = obj.data.stepInd(1,ind);
                  case 'last'
                      ind  = obj.indIndata;
                      step = obj.data.stepInd(1,ind);
              end
              
          end 
          
          plotParticles(obj.data.allCoords(:,:,ind),obj.simulationParam.L...
              ,obj.simulationParam.r,obj.currentAngs);
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
           addOptional(p, 'keepFigOpen', false);
           parse(p, varargin{:});
           Results = p.Results;
           plotVarVsStep = Results.plotVarVsStep;
           plotVarVarVsStep = Results.plotVarVarVsStep;
           saveFig = Results.saveFig;
           keepFigOpen = Results.keepFigOpen;
           startFromStep = Results.startFromStep;
           
           if startFromStep == 0
               minInd = 1;
           else
               % find the relevent initial index
               [~, minInd] = min(abs(startFromStep - obj.data.stepInd));
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
           
           steps = obj.data.stepInd(1,minInd:obj.indIndata);
           
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
               saveas(gcf,[name '.fig']);
               saveas(gcf,[name '.jpg']);
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
               saveas(gcf,[name '.fig']);
               saveas(gcf,[name '.jpg']);
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
           indIndata = obj.indIndata;
           
           meanU = obj.data.meanUlrcEq;
           U = obj.data.allUlrc(1,firstSteps2ignore:indIndata);
           sigU = mean((U - meanU).^2);

           obj.data.sigU = sigU;
           obj.data.cv = N + (sigU/T^2);
       end
        
    end    

end


        function [dist,particlesPosition] = ...
            createInitialConfig(L,N,r,initialConfig)

            possibleInitialConfigs = {'random','hex','auto'};
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
                if sepDist < 2*r
                    % Display an error message
                    disp('density is too high');
                    particlesPosition = [];
                    return
                end

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
                    % make d a row vector
                    d = reshape(allDists,1,[]);
                    % ignore zeros
                    d = nonzeros(d); 
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
        
