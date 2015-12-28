function [finalU,finalConfiguration,finalDistances,moveCount] = ...
    MonteCarlo2DLJHeart(N,T,rho,Nsteps,maxdr,initialConfig,rCutoff...
    ,initialDistances,initialU)

%% Monte-Carlo in NVT ensemble for Lennard-Jonse potantioal in 2D %%

% calculate the energy after 'Nstep' Monte-Carlo runs (Metropolis algorithm)
% for a system of 'N' particles in 2D interacting via Lennard-Jonse pair
% potantial. periodic boundary conditions (PBC) are apllied.

% inputs:
% N - number of particles
% T - reduced temperature
% Nsteps - number of steps
% maxdr - maximum particle displacement
% initialConfig - initial configuration of particles (2 by N matrix)
% initialU - initial energy of the configuration
% rCutoff - the cutoff distance for the energy

% outputs:
% finalU - the energy in each step (1 by Nstep matrix)
% finalConfigurations - the coordinates of all particles in each step (2 by N
%                     by Nstep matrix)
% finalDistances - the pair distances of all particles in each step (N by N
%                by Nstep matrix)
% moveCount - counts accepted moves

% the potantial of Lennard-Jonse in reduced units:
%   U = 4*[(1/r)^12 - (1/r)^6]

% the potantial of Lennard-Jonse in non-reduced units:
%   U = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]

% reduced units:
% T(reduced) = kT/epsilon | r(reduced) = r/sigma | U(reduced) = U/epsilon


dist = initialDistances;
particlesPosition = initialConfig;
U = initialU;
L = sqrt(N/rho); % board length in reduced units
moveCount = 0;
movedParticle = 0;

for step = 1:Nsteps
    
        % choose particle to move 
        movedParticle = movedParticle + 1;
        if movedParticle == N + 1
            movedParticle = 1;
        end
                            
        % choose displacement:
        displacex = maxdr*rand - (maxdr/2);
        displacey = maxdr*rand - (maxdr/2);

        % move particle
        newParticlesPosition = movePBC(particlesPosition,movedParticle,...
            displacex,displacey,L);

        % calculate new distances
        newDist = reCalcDist(dist,movedParticle,...
            particlesPosition,newParticlesPosition,N,L);
        
        % calculate the change in energy
        dU = Uchange(movedParticle,dist,newDist,N,rCutoff);

        % if dU < 0 eccept move
        %(if (1/T)*dU > 75, we are sure the move
        % will not be excepted, so we don't calculate exp(1/T)*dU to 
        % save calculation time)
            
        if (1/T)*dU < 75
            if dU < 0  
                U = U + dU;
                dist = newDist;
                particlesPosition = newParticlesPosition;
                moveCount = moveCount + 1; 
            else
                %% otherwise,
                % keep the new state with a probability corresponding to the
                % Boltzmann factor. if the new state is rejected, recount the
                % old configuration. 

                if rand < exp(-(1/T)*dU)
                    U = U + dU;
                    dist = newDist;
                    particlesPosition = newParticlesPosition;
                    moveCount = moveCount + 1;
                end
            end
        end
end

finalU = U;
finalConfiguration = particlesPosition;
finalDistances = dist;

%% functions used in main code %%

        function newparticlesPosition = movePBC(particlesPosition,...
                movedParticle,displacex,displacey,L)
            
                % if the particle gets out of the board, apply PBC
                newparticlesPosition = particlesPosition;
                x = particlesPosition(1,movedParticle)+ displacex;
                y = particlesPosition(2,movedParticle)+ displacey;
                
                if x > L/2
                    x = x - L;
                end
                
                if x < (-L/2) 
                    x =  x + L;
                end
                
                if y > L/2 
                    y = y - L;
                end
                
                if y < (-L/2) 
                    y = y + L;
                end
                
                newparticlesPosition(1,movedParticle) = x;
                newparticlesPosition(2,movedParticle) = y;
        end
    
        function newdist = reCalcDist(dist,movedParticle,...
                particlesPosition,newParticlesPosition,N,L)
            
                % recalculates pair distances after moving a particle
            
                xi = newParticlesPosition(1,movedParticle);
                yi = newParticlesPosition(2,movedParticle);
                newdist = dist;
                
                % recalculate the relevent row elements in dist matrix
                
                if movedParticle > 1
                    newdist(movedParticle,1:(movedParticle-1)) =...
                        distPBC(xi,yi,...
                            particlesPosition(:,1:(movedParticle-1)),L);
                end
                
                % recalculate the relevent column elements in dist matrix
                
                if movedParticle < N
                    newdist((movedParticle + 1):N,movedParticle) =...
                        distPBC(xi,yi,...
                            particlesPosition(:,(movedParticle+1):N),L);
                end
                
        end
    
        function dU = Uchange(movedParticle,dist,newDist,N,rCutoff)
        % calculates the change in energy after a particle has moved
                
                % calculate the old energy for the relevant particle pairs
                
                if movedParticle > 1
                    oldUrow = ...
                        pairU(dist(movedParticle,1:(movedParticle - 1)),rCutoff);
                else 
                    oldUrow = 0;
                end
                
                if movedParticle < N
                    oldUcol = ...
                        pairU(dist((movedParticle + 1):N,movedParticle),rCutoff);
                else 
                    oldUcol = 0;
                end

                oldU = oldUrow + oldUcol;
                
                % calculate the new energy for the relevant particle pairs
                
                if movedParticle > 1
                    newUrow = pairU(newDist...
                            (movedParticle,1:(movedParticle - 1)),rCutoff);
                else 
                    newUrow = 0;
                end
                
                if movedParticle < N
                    newUcol = pairU(newDist...
                        ((movedParticle + 1):N,movedParticle),rCutoff);
                else 
                    newUcol = 0;
                end
                
                newU = newUrow + newUcol;
                
                % clculate the change in energy
                
                dU = newU - oldU;
        end
    
end
