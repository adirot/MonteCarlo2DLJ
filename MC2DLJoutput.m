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
%       simulationParam.rCutoff - the cutoff distance for the energy
%       simulationParam.rl - rl for verelet algorithm. if empty - verelet
%                               algorithm will not be used.
%       simulationParam.pressure - true or false, calculate the pressure or
%                                   not.

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

% data inclides:
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
%       histo - histogram of the radial distribution function (for every
%                   step saved)
%       bins - bins of the radial distribution function


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
          currentVir,Ulrc,Plrc,currentStep;
      fileName,data,indIndata;
      RDFhisto, RDFbins;
      
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
                    if ~isempty(obj.data.allV)
                        obj.currentVir = obj.data.allV(1,obj.indIndata);
                        obj.currentPressure = obj.data.allP(1,obj.indIndata);
                    end
                    obj.currentStep = obj.data.stepInd(1,obj.indIndata);

                end
                    
                if nargin >= 7
                    
                    
                    p = inputParser();
                    addOptional(p, 'verelet', []);
                    addOptional(p, 'pressure', false);
                    addOptional(p, 'fileNameInit', '');
                    parse(p, varargin{8:end});
                    Results = p.Results;
                    rl = Results.verelet;
                    pressure = Results.pressure;
                    fileNameInit = Rsults.fileNameInit;
                    
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
                    
                    obj.currentmaxdr = obj.simulationParam.initialmaxdr;
                    obj.moveCount = 0;
                    obj.indIndata = 0;
                    
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
                    
                    obj.fileName = [fileNameInit 'N' num2str(N)...
                        'T' my_num2str(T)...
                        'rho' my_num2str(rho) 'initialmaxdr'...
                        my_num2str(initialmaxdr) 'initialconfig_'...
                        initialConfig 'rCutoff'...
                         my_num2str(rCutoff) vereletstr '_' pressurestr...
                         'date'...
                        nowdatetimestr()];
                        
                    % create initial configuration
                    allCoords = zeros(2,N,2);
                    allDists = zeros(N,N,2);
                    [allDists(:,:,1),allCoords(:,:,1)] = ...
                            createInitialConfig(obj.simulationParam.L,N,r...
                            ,initialConfig);
                        
                    % calculate initial energy
                    d = reshapeDist(allDists);
                    allU = pairU(d,rCutoff);
                    allUlrc = allU - pi*rho*N/rCutoff^4;
                    
                    % calculate initial virial and pressere
                    if pressure
                        allV = ...
                            calcVirial(d,rho,12,6,N,rCutoff);
                        allP = T*rho + allV;
                        allPlrc = allP - 3*pi*rho^2/rCutoff^4;
                    else
                        allV = [];
                        allP = [];
                        allPlrc = [];
                    end
                    clear d;
                    
                    % save to data file
                    stepInd = 0;
                    moveCount = 0;
                    indIndata = 1;
                    currentmaxdr = initialmaxdr;
                    simulationParam = obj.simulationParam;
                    save(obj.fileName, 'allDists','allCoords',...
                            'allU','allUlrc','allV','allP','allPlrc','stepInd',...
                            'moveCount','indIndata'...
                            ,'currentmaxdr','simulationParam','-v7.3');
                    obj.data = matfile(obj.fileName);
                    obj.data = matfile(obj.fileName,'Writable',true);
                        
                    obj.currentCoords = allCoords;
                    obj.currentDists = allDists;
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
            N = obj.simulationParam.N; rho = obj.simulationParam.rho;
                rCutoff = obj.simulationParam.rCutoff;        
                    
            stepCount = 0;
            while(stepCount < Nsteps)
                [finalU,finalV,finalPressure,...
                    finalConfiguration,finalDistances,...
                    currentmoveCount] = ...
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
                    'virial',obj.currentVir);
                
                stepCount = stepCount + obj.simulationParam.N*saveEvery;
                obj.currentStep = obj.currentStep...
                    + obj.simulationParam.N*saveEvery;
                obj.moveCount = obj.moveCount + currentmoveCount;
                obj.currentU = finalU;
                obj.currentVir = finalV;
                obj.currentPressure = finalPressure;
                obj.currentCoords = finalConfiguration;
                obj.currentDists = finalDistances;
                
            
                % long range correction
                obj.Ulrc = obj.currentU - pi*N*rho/rCutoff^4;
                if obj.simulationParam.pressure
                    obj.Plrc = obj.currentPressure - 3*pi*rho^2/rCutoff^4;
                end
            
                obj = obj.addStep2data(obj.currentStep,finalConfiguration,...
                    finalDistances,finalU,finalV,finalPressure,...
                    obj.moveCount,obj.currentmaxdr,obj.Ulrc,obj.Plrc);
                
                clear finalU finalV finalConfiguration finalDistances...
                    currentmoveCount
            end
            
            
            
            
        end
        
       function obj = addStep2data(obj,newInd,newCoords,newDists,newU,newV...
                ,newP,moveCount,currentmaxdr,newUlrc,newpressurelrc)
            
            obj.data.stepInd(1,obj.indIndata+1) = newInd;
            s = zeros(2,obj.simulationParam.N,2);
            s(:,:,1) = newCoords(:,:,1);
            obj.data.allCoords(:,:,obj.indIndata+1) = s(:,:,1);
            s = zeros(obj.simulationParam.N,obj.simulationParam.N,2);
            s(:,:,1) = newDists(:,:,1);
            obj.data.allDists(:,:,obj.indIndata+1) = s(:,:,1);
            clear s;
            obj.data.allU(1,obj.indIndata+1) = newU;
            obj.data.allUlrc(1,obj.indIndata+1) = newUlrc;
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
       
       function obj = calcRDF(obj,maxDist,numOfBins)
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
            
              N = obj.simulationParam.N;
              rho = obj.simulationParam.rho;
              
              
              obj.data.RDFhisto = zeros(obj.indIndata,numOfBins);      
              bins = linspace(0,maxDist,numOfBins); 
              obj.data.RDFbins = bins;
              obj.RDFhisto = zeros(1,numOfBins);
              
              for step = 1:obj.indIndata

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
                        obj.data.RDFhisto(step,:) = histo;
                        obj.RDFhisto = obj.RDFhisto + histo;
              end
            obj.RDFbins = bins;
            obj.RDFhisto = obj.RDFhisto/obj.indIndata;
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
              ,obj.simulationParam.r);
          title(['snapshot of step: ' num2str(step)]);
       end
       
       function obj = meanProp(obj)
            % makes mean vectors for U,P.
            % meanProp(i) is the mean of all the values of property on
            % steps 1 to i.
            
            obj.data.meanU = my_mean(obj.data.allU);
            obj.data.meanUlrc = my_mean(obj.data.allUlrc);
            obj.data.meanP = my_mean(obj.data.allP);
            obj.data.meanPlrc = my_mean(obj.data.allPlrc);
            
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
        
    end
    
        
        

end

        function [dist,particlesPosition] = ...
            createInitialConfig(L,N,r,initialConfig)

        possibleInitialConfigs = {'random','hex'};
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
        
