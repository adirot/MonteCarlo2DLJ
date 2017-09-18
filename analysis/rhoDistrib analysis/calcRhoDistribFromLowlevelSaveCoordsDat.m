function [histxnumOfPartInSquare, histnumOfPartInSquare, meanhistnumOfPartInSquare] = ...
    calcRhoDistribFromLowlevelSaveCoordsDat(rhoDistribDatName,...
    CoordsDatName, N, rho, numOfSquares, varargin)
    
    p = inputParser();
    addOptional(p, 'talk', false);
    addOptional(p, 'plotHist', false); 
    addOptional(p, 'plotHist_times_partNum', false);
    parse(p, varargin{:});
    Results = p.Results;
    talk = Results.talk;
    plotHist = Results.plotHist;
    plotHist_times_partNum = Results.plotHist_times_partNum;
    
    CoordsFid = fopen(CoordsDatName,'r');
    rhoDistribFid = fopen(rhoDistribDatName,'w');
    
    L = sqrt(N/rho);
    histxnumOfPartInSquare = 0:N;
        
    stepCount = 0;
    coords = fread(CoordsFid,[2 N],'double');
    figure;
    hold on;
    while ~isempty(coords)
        stepCount = stepCount + 1;
        if talk
            disp(stepCount);
        end
        
        numOfPartInSquare = numOfCellsDistribution(coords, L, numOfSquares);
        if stepCount == 1
            histnumOfPartInSquare = histcounts(numOfPartInSquare,0:N+1);
            meanhistnumOfPartInSquare = histnumOfPartInSquare;
        else
            histnumOfPartInSquare = histnumOfPartInSquare +...
            histcounts(numOfPartInSquare,0:N+1);
        end
        s = sum(histxnumOfPartInSquare.*histnumOfPartInSquare);
        fwrite(rhoDistribFid,histnumOfPartInSquare/s,'double');
        if stepCount > 1
            meanhistnumOfPartInSquare = meanhistnumOfPartInSquare + histnumOfPartInSquare/s;
        end
        coords = fread(CoordsFid,[2 N],'double');
        
        if mod(stepCount,10000) == 0
            plot(histxnumOfPartInSquare, meanhistnumOfPartInSquare/stepCount);
        end
    end
    
    meanhistnumOfPartInSquare = meanhistnumOfPartInSquare/stepCount;
    
    fclose(CoordsFid);
    fclose(rhoDistribFid);
    
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
        m = sqrt(numOfSquares);
        if isinteger(m)
        error('The square root of numOfSquares must be an intiger');
        end


        % Get the density in each square  
        squareArea = L^2/numOfSquares;
        squareSide = sqrt(squareArea);
        Ni = zeros(m,m);
        indi = 0;

        for ii = (-L/2):squareSide:(L/2 - squareSide)
        indi = indi + 1;
        indj = 0;
        for jj = (-L/2):squareSide:(L/2 - squareSide)
            indj = indj + 1;
            Ni(indi,indj) =...
                countParticlesInSquare(ii,jj,coords(1,:),coords(2,:),squareSide);
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

    function coords = freadCoords(fid,firstStep,lastStep,stepPointers,N)
        % reads coordinates from 'firstStep' to 'lastStep' from a dat file in 'fid'.
        % 'fid' is a fid with a dat file. The fid sould be openned with an
        % 'r' flag. The output fid will have a pointer with the location after the
        % data read. 'coords' is a 2 by N*number of steps read.
        fseek(fid,stepPointers(firstStep),'bof');
        numOfSteps2read = lastStep - firstStep;
        coords = fread(fid,[2 N*numOfSteps2read],'double');
    end
end