function [finalU,finalVirial,finalPressure,finalConfiguration,...
    finalDistances,moveCount,finalAngs,finalAlphas,finalThetas] = ...
    MonteCarlo2DLJHeart(N,T,rho,Nsteps,maxdr,initialConfig,rCutoff...
    ,initialDistances,initialU,varargin)

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
% optional: 'verelet', rl : use verelet neighbor algorithm with rCutoff
%           (rc) inner radius and rl outer radius. see section 5.3.1 in the
%           book computer simulation of liquids for more information.
% optional: the 'm' in the intetraction (U = 4((1/r)^12 - (1/r)^-m)),
%           default is 6.


% outputs:
% finalU - the energy in each step (1 by Nstep matrix)
% finalConfigurations - the coordinates of all particles in each step (2 by N
%                     by Nstep matrix)
% finalDistances - the pair distances of all particles in each step (N by N
%                by Nstep matrix)
% moveCount - counts accepted moves

% the potantial of Lennard-Jonse in reduced units:
%   U = 4*[(1/r)^12 - (1/r)^6] 
%   (you can change the 6 with the optional input 'm')

% the potantial of Lennard-Jonse in non-reduced units:
%   U = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]

% reduced units:
% T(reduced) = kT/epsilon | r(reduced) = r/sigma | U(reduced) = U/epsilon

% usage examples: 

% with verelet algorithm
% [finalU,finalConfiguration,finalDistances,moveCount] = ...
%     MonteCarlo2DLJHeart(N,T,rho,Nsteps,maxdr,initialConfig,rCutoff...
%     ,initialDistances,initialU,'verelet',2.7)

% without verelet algorithm
% [finalU,finalConfiguration,finalDistances,moveCount] = ...
%     MonteCarlo2DLJHeart(N,T,rho,Nsteps,maxdr,initialConfig,rCutoff...
%     ,initialDistances,initialU)


% parse input parameters
p = inputParser();
addOptional(p, 'verelet', []); 
addOptional(p, 'virial', []);
addOptional(p, 'm', 6);
addOptional(p, 'angleDependent',false);
addOptional(p, 'angleDependence',false);
addOptional(p, 'initialAng', []);
addOptional(p, 'initialAlphas', []);
addOptional(p, 'initialThetas', []);
addOptional(p, 'maxdAng', []);
addOptional(p, 'ufunc', []);
parse(p, varargin{:});
Results = p.Results;
rl = Results.verelet;
virial = Results.virial;
m = Results.m;
angleDependent = Results.angleDependent;
angleDependence = Results.angleDependence;
initialAngs = Results.initialAng;
initialAlphas = Results.initialAlphas;
initialThetas = Results.initialThetas;
maxdAng = Results.maxdAng;
ufunc = Results.ufunc;

if isempty(ufunc)
    ufunc = @(r) 4*(((1./r).^12)-((1./r).^m));
end

% initiate virables
dist = initialDistances;
particlesPosition = initialConfig;
if angleDependent
    particlesAngs = initialAngs;
    particlesAlphas = initialAlphas;
    particlesThetas = initialThetas;
end
U = initialU;
if ~isempty(virial)
    V = virial;
end
L = sqrt(N/rho); % board length in reduced units
moveCount = 0;
movedParticle = 0;

% calculate nieghbour list
if ~isempty(rl)
    
    % construct neighbors list object
    nlist = verelet(dist,rl,N);
else 
    nlist = [];

end

for step = 1:Nsteps
    
        % choose particle to move 
        movedParticle = movedParticle + 1;
        if movedParticle == N + 1
            movedParticle = 1;
        end
        
        % update verelet neighbors list
        if (~isempty(nlist))&&(nlist.summaxdisplace > (rl-rCutoff))
              
            % we need to update the distances
            movedParticles = find(nlist.dispacements > 0);
            dist = reCalcDist(dist,movedParticles,particlesPosition,N,L,[]);
            
            % and the verlet list
            nlist = verelet(dist,rl,N);
        end
       
                            
        % choose displacement:
        displacex = maxdr*rand - (maxdr/2);
        displacey = maxdr*rand - (maxdr/2);
        displace = sqrt(displacex^2 + displacey^2);
        
        % choose rotation:
        if angleDependent
            dAng = maxdAng*rand - (maxdAng/2);
        end
        
        % move particle
        newParticlesPosition = movePBC(particlesPosition,movedParticle,...
            displacex,displacey,L);
        
        % rotate particle
        if angleDependent
            newParticlesAngs = particlesAngs;
            newParticlesAngs(movedParticle) =...
                particlesAngs(movedParticle) + dAng; 
        end

        % calculate new distances
        newDist = reCalcDist(dist,movedParticle,...
            newParticlesPosition,N,L,nlist);
        
        % calculate new relative angles
        if angleDependent
            newAlphas = reCalcAlphas(particlesAlphas,movedParticle,N,dAng);
            newThetas = reCalcThetas(movedParticle,...
                newParticlesPosition,N,particlesThetas,newAlphas);
        else
            newAlphas = [];
            newThetas = [];
        end
        
        % calculate the change in energy
        dU = Uchange(movedParticle,dist,newDist,N,rCutoff,m,...
            particlesAlphas,newAlphas,particlesThetas,newThetas,...
            angleDependence,ufunc);
        
        % calculate the change in the virial 
        if ~isempty(virial)
            dV = Vchange(movedParticle,dist,newDist,N,rCutoff,rho,m);
        end
        
        % if dU < 0 eccept move
        %(if (1/T)*dU > 75, we are sure the move
        % will not be excepted, so we don't calculate exp(1/T)*dU to 
        % save calculation time)
            
        if (1/T)*dU < 75
            if dU < 0  
                U = U + dU;
                if ~isempty(virial)
                    V = V + dV;
                end
                dist = newDist;
                particlesPosition = newParticlesPosition;
                moveCount = moveCount + 1; 
                
                if angleDependent
                    particlesAngs = newParticlesAngs;
                    particlesAlphas = newAlphas;
                    particlesThetas = newThetas;
                end
                    
                % update verelet displacement count
                if ~isempty(nlist)
                    nlist = nlist.updateDisplace(movedParticle,displace);
                end
                
            else
                %% otherwise,
                % keep the new state with a probability corresponding to the
                % Boltzmann factor. if the new state is rejected, recount the
                % old configuration. 

                if rand < exp(-(1/T)*dU)
                    U = U + dU;
                    if ~isempty(virial)
                        V = V + dV;
                    end
                    dist = newDist;
                    particlesPosition = newParticlesPosition;
                    moveCount = moveCount + 1;
                    
                    if angleDependent
                        particlesAngs = newParticlesAngs;
                        particlesAlphas = newAlphas;
                        particlesThetas = newThetas;
                    end
                    
                    % update verelet displacement count
                    if ~isempty(nlist)
                        nlist = nlist.updateDisplace(movedParticle,displace);
                    end
                    
                end
            end
        end
end

finalU = U;
if ~isempty(virial)
    finalVirial = V;
    finalPressure = T*rho + V;
else
    finalVirial = [];
    finalPressure = [];
end
finalConfiguration = particlesPosition;
finalDistances = dist;

if angleDependent
    finalAngs = particlesAngs;
    finalAlphas = particlesAlphas;
    finalThetas = particlesThetas;
else 
    finalAngs = [];
    finalAlphas = [];
    finalThetas = [];
end

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
    
        function newdist = reCalcDist(dist,movedParticles,...
                newParticlesPosition,N,L,nlist)
                
                for i = 1:length(movedParticles)
                    movedP = movedParticles(i);
                    % recalculates pair distances after moving a particle    

                    xi = newParticlesPosition(1,movedP);
                    yi = newParticlesPosition(2,movedP);
                    newdist = dist;

                    % recalculate the relevent row elements in dist matrix

                    if movedP > 1
                        if ~isempty(nlist)
                            neiInd =...
                                nlist.neighborsindy(nlist.neighborsindx == movedP);
                            newdist(movedP,neiInd) =...
                            distPBC(xi,yi,...
                            newParticlesPosition(:,neiInd),...
                            L);

                        else

                            newdist(movedP,1:(movedP-1)) =...
                                distPBC(xi,yi,...
                                newParticlesPosition(:,1:(movedP-1)),L);
                        end
                    end

                    % recalculate the relevent column elements in dist matrix

                    if movedP < N
                        if ~isempty(nlist)
                            neiInd =...
                                nlist.neighborsindx(nlist.neighborsindy == movedP);
                            newdist(neiInd,movedP) =...
                            distPBC(xi,yi,...
                                newParticlesPosition(:,neiInd),L);

                        else
                            newdist((movedP + 1):N,movedP) =...
                            distPBC(xi,yi,...
                                newParticlesPosition(:,(movedP+1):N),L);
                        end
                    end
                end   
        end
        
        function newAlphas = reCalcAlphas(alphas,movedParticles,N,dAng)
                
                for i = 1:length(movedParticles)
                    movedP = movedParticles(i);
                    % recalculates alphas after moving a particle    

                    newAlphas = alphas;

                    % recalculate the relevent row elements in alpha matrix

                    if movedP > 1                        
                        newAlphas(movedP,1:(movedP-1)) =...
                            alphas(movedP,1:(movedP-1)) + dAng;
                    end

                    % recalculate the relevent column elements in alpha matrix

                    if movedP < N
                        newAlphas((movedP + 1):N,movedP) =...
                            alphas((movedP + 1):N,movedP) + dAng;
                    end
                end   
        end
    
        function newThetas = reCalcThetas(movedParticles,...
                newParticlesPosition,N,thetas,newAlphas)
                
                for i = 1:length(movedParticles)
                    movedP = movedParticles(i);
                    % recalculates thetas after moving a particle    

                    xi = newParticlesPosition(1,movedP);
                    yi = newParticlesPosition(2,movedP);

                    newThetas = thetas;
                    
                    % recalculate the relevent row elements in theta matrix

                    if movedP > 1
                        newThetas(movedP,1:(movedP-1)) =...
                            newAlphas(movedP,1:(movedP-1)) - atan(yi/xi);
                    end

                    % recalculate the relevent column elements in theta matrix

                    if movedP < N
                        newThetas((movedP + 1):N,movedP) =...
                            newAlphas((movedP + 1):N,movedP) - atan(yi/xi);
                    end
                end   
        end
        
        function dU = Uchange(movedParticle,dist,newDist,N,rCutoff,m,...
                alphas,newAlphas,thetas,newThetas,angleDependence,ufunc)
        % calculates the change in energy after a particle has moved
        
                % calculate the old energy for the relevant particle pairs
                
                if movedParticle > 1
                    
                    if ~isempty(angleDependence)
                        relAng = [];
                        relAng(:,1) = alphas(movedParticle,1:(movedParticle - 1));
                        relAng(:,2) = thetas(movedParticle,1:(movedParticle - 1));
                    else
                        relAng = [];
                    end
                    distrow = dist(movedParticle,1:(movedParticle - 1))';
                    oldUrow = ...
                        pairU(distrow,rCutoff,m,...
                        'angleDependence',angleDependence,...
                        'relativeCellAngles',relAng,'ufunc',ufunc);
                else 
                    oldUrow = 0;
                end
                
                if movedParticle < N
                    if ~isempty(angleDependence)
                        relAng = [];
                        relAng(:,1) =...
                            alphas((movedParticle + 1):N,movedParticle);
                        relAng(:,2) =...
                            thetas((movedParticle + 1):N,movedParticle);
                    else
                        relAng = [];
                    end
                    distcol = dist((movedParticle + 1):N,movedParticle);
                    oldUcol = ...
                        pairU(distcol,rCutoff,m,...
                        'angleDependence',angleDependence,...
                        'relativeCellAngles',relAng,'ufunc',ufunc);
                else 
                    oldUcol = 0;
                end

                oldU = oldUrow + oldUcol;
                
                % calculate the new energy for the relevant particle pairs
                
                if movedParticle > 1
                    if ~isempty(angleDependence)
                        relAng = [];
                        relAng(:,1) =...
                            newAlphas(movedParticle,1:(movedParticle - 1));
                        relAng(:,2) =...
                            newThetas(movedParticle,1:(movedParticle - 1));
                    else
                        relAng = [];
                    end
                    distrow = newDist(movedParticle,1:(movedParticle - 1))';
                    newUrow = pairU(distrow,rCutoff,m,...
                        'angleDependence',angleDependence,...
                        'relativeCellAngles',relAng,'ufunc',ufunc);
                else 
                    newUrow = 0;
                end
                
                if movedParticle < N
                    if ~isempty(angleDependence)
                        relAng = [];
                        relAng(:,1) =...
                            newAlphas((movedParticle + 1):N,movedParticle);
                        relAng(:,2) =...
                            newThetas((movedParticle + 1):N,movedParticle);
                    else
                        relAng = [];
                    end
                    distcol = newDist((movedParticle + 1):N,movedParticle);
                    newUcol = pairU(distcol,rCutoff,m,...
                        'angleDependence',angleDependence,...
                        'relativeCellAngles',relAng,'ufunc',ufunc);
                else 
                    newUcol = 0;
                end
                
                newU = newUrow + newUcol;
                
                % clculate the change in energy
                
                dU = newU - oldU;
        end
    
        function dV = Vchange(movedParticle,dist,newDist,N,rCutoff,rho,m)
        % calculates the change in the virial after a particle has moved
                
                % calculate the old virial for the relevant particle pairs
                
                if movedParticle > 1
                    oldVrow = ...
                        calcVirial(dist(movedParticle,1:(movedParticle - 1))...
                        ,rho,12,m,N,rCutoff);
                else 
                    oldVrow = 0;
                end
                
                if movedParticle < N
                    oldVcol = ...
                        calcVirial(dist((movedParticle + 1):N,movedParticle)...
                        ,rho,12,m,N,rCutoff);
                else 
                    oldVcol = 0;
                end

                oldV = oldVrow + oldVcol;
                
                % calculate the new virial for the relevant particle pairs
                
                if movedParticle > 1
                    newVrow = calcVirial(newDist...
                            (movedParticle,1:(movedParticle - 1))...
                            ,rho,12,m,N,rCutoff);
                else 
                    newVrow = 0;
                end
                
                if movedParticle < N
                    newVcol = calcVirial(newDist...
                        ((movedParticle + 1):N,movedParticle)...
                        ,rho,12,m,N,rCutoff);
                else 
                    newVcol = 0;
                end
                
                newV = newVrow + newVcol;
                
                % clculate the change in the virial
                
                dV = newV - oldV;
        end
        
    
            
end