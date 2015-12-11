function [finalU,finalConfiguration] = MonteCarlo2DLJ(N,T,rho,Nsteps,maxdr,initialConfig,rCutoff)


%% Monte-Carlo in NVT ensemble for Lennard-Jonse potantioal in 2D %%

L = sqrt(N/rho);
r = 2^(1/6)/2; % particle radius in reduced units 

%% create initial configuration, find all pair distances %%
[dist,particlesPosition] = createInitialConfig(L,N,r,initialConfig);

plotParticles(particlesPosition,L,2^(1/6)/2);
title('start');


%% monte carlo %%

% calculate initial energy 
d = reshape(dist,1,[]); 
d = nonzeros(d); % make d a row vector
initialU = pairU(d,rCutoff);

% perform Nsteps steps of monte carlo
[finalU,finalConfiguration,finalDistances,moveCount] = ...
    MonteCarlo2DLJHeart(N,T,rho,Nsteps,maxdr,particlesPosition...
    ,dist,initialU,rCutoff);

%% functions used in the main code %%

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
                yPos = (sqrt(3)/2)*xPos;

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
                        distPBC(x,y,particlesPosition(:,(par+1):N));
                end
                
        end
end