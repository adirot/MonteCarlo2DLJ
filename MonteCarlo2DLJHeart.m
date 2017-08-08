function [finalU,finalVirial,finalPressure,finalConfiguration,...
    finalDistances,moveCount,finalAngs,finalBettas] = ...
    MonteCarlo2DLJHeart(N,T,rho,Nsteps,maxdr,initialConfig,rc...
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
% rc - the cutoff distance for the energy
% optional: 'vereletRadius', rl : use verelet neighbor algorithm with rc
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
%   U = (1/r)^12 - (1/r)^6 
%   (you can change the 6 with the optional input 'm')

% the potantial of Lennard-Jonse in non-reduced units:
%   U = epsilon*[(sigma/r)^12 - (sigma/r)^6]

% reduced units:
% T(reduced) = kT/epsilon | r(reduced) = r/sigma | U(reduced) = U/epsilon

% usage example: 

% [finalU,finalConfiguration,finalDistances,moveCount] = ...
%     MonteCarlo2DLJHeart(N,T,rho,Nsteps,maxdr,initialConfig,rc...
%     ,initialDistances,initialU)


% parse input parameters
p = inputParser();
addOptional(p, 'vereletRadius', []);
addOptional(p, 'energyCutoffRadius', []); 
addOptional(p, 'virial', []);
addOptional(p, 'm', 6);
addOptional(p, 'angleDependent',false);
addOptional(p, 'angleDependence',[]);
addOptional(p, 'initialAng', []);
%addOptional(p, 'initialAlphas', []);
%addOptional(p, 'initialThetas', []);
addOptional(p, 'initialBettas', []);
addOptional(p, 'maxdAng', []);
addOptional(p, 'ufunc', []);
addOptional(p, 'hardCoreRepRad', 0);
addOptional(p, 'TalkEvery', []);
addOptional(p, 'dipoleStrength', []);
parse(p, varargin{:});
Results = p.Results;
rl = Results.vereletRadius;
rc = Results.energyCutoffRadius;
virial = Results.virial;
m = Results.m;
angleDependent = Results.angleDependent;
angleDependence = Results.angleDependence;
initialAngs = Results.initialAng;
% initialAlphas = Results.initialAlphas;
% initialThetas = Results.initialThetas;
initialBettas = Results.initialBettas;
maxdAng = Results.maxdAng;
ufunc = Results.ufunc;
hardCoreRepRad = Results.hardCoreRepRad;
TalkEvery = Results.TalkEvery;
dipoleStrength = Results.dipoleStrength;

if isempty(ufunc)
    ufunc = @(r,m) (((1./r).^12)-((1./r).^m));
end

% initiate virables
dist = initialDistances;
particlesPosition = initialConfig;
if angleDependent
    particlesAngs = initialAngs;
%    particlesAlphas = initialAlphas;
%    particlesThetas = initialThetas;
    particlesBettas = initialBettas;
else 
%    particlesAlphas = [];
%    particlesThetas = [];
    particlesBettas = [];
    particlesAngs = [];
end
U = initialU;
if ~isempty(virial)
    V = virial;
end
L = sqrt(N/rho); % board length in reduced units
moveCount = 0;
movedParticle = 0;

if isempty(dipoleStrength)
    dipoleStrength = ones(1,N);
end

% calculate nieghbour list
if ~isempty(rl)
    if isempty(rc)
        error('Cutoff radius input is needed for verelet algoritm');
    else
        % construct neighbors list object (see verelet.m)
        nlist = verelet(dist,rl,rc,N);
    end
else 
    nlist = [];
end

for step = 1:Nsteps
        
        % talk
        if ~isempty(TalkEvery)
            if mod(step,TalkEvery) == 0
                disp([num2str(step) 'steps of MC done']);
            end
        end
    
        % choose particle to move 
        movedParticle = movedParticle + 1;
        if movedParticle == N + 1
            movedParticle = 1;
        end
        
        % update verelet neighbors list
        if ~isempty(nlist)
            if nlist.needs2beInitialized
              
                % Update the neighbors list
                
                    movedParticles = find(nlist.dispacements > 0);
                    dist = reCalcDist(dist,movedParticles,particlesPosition,N,L,[]);
            
                % and the verlet list
                nlist = verelet(dist,rl,N);
            end
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
            if newParticlesAngs(movedParticle) > 2*pi
                newParticlesAngs(movedParticle) = ...
                    newParticlesAngs(movedParticle) - 2*pi;
            end
        end

        % calculate new distances
        newDist = reCalcDist(dist,movedParticle,...
            newParticlesPosition,N,L,nlist);
        
        % plot particles (for debug)
        % plotParticles(newParticlesPosition,L,hardCoreRepRad,'angles',newParticlesAngs)
          
        if ~check_if_cores_touch(newDist,hardCoreRepRad)% no particles are inside each other 
        
            % calculate new relative angles
            if angleDependent
                %newAlphas = reCalcAlphas(particlesAlphas,movedParticle,N,dAng);
                %newThetas = reCalcThetas(movedParticle,...
                %    newParticlesPosition,N,particlesThetas,newAlphas);
                newBettas =...
                    reCalcBetta(movedParticle,particlesBettas,newParticlesPosition,N);
            else
%                 newAlphas = [];
%                 newThetas = [];
                newBettas = [];
                newParticlesAngs = [];
            end

            % calculate the change in energy
            dU = Uchange(movedParticle,dist,newDist,N,rc,m,...
                particlesBettas,newBettas,particlesAngs,newParticlesAngs,...
                angleDependence,ufunc,dipoleStrength);

            % calculate the change in the virial 
            if ~isempty(virial)
                dV = Vchange(movedParticle,dist,newDist,N,rc,rho,m);
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
%                         particlesAlphas = newAlphas;
%                         particlesThetas = newThetas;
                        particlesBettas = newBettas;
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
%                         particlesAlphas = newAlphas;
%                         particlesThetas = newThetas;
                            particlesBettas = newBettas;
                        end

                        % update verelet displacement count
                        if ~isempty(nlist)
                            nlist = nlist.updateDisplace(movedParticle,displace);
                        end

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
%    finalAlphas = particlesAlphas;
%    finalThetas = particlesThetas;
    finalBettas = particlesBettas;
else 
    finalAngs = [];
%    finalAlphas = [];
%    finalThetas = [];
    finalBettas = [];
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
                            largerThan2pi =...
                                newAlphas(movedP,1:(movedP-1)) > 2*pi;
                            newAlphas(movedP,1:(movedP-1)) = ...
                                newAlphas(movedP,1:(movedP-1)) -...
                                2*pi*largerThan2pi;
                        end
                    

                    % recalculate the relevent column elements in alpha matrix

                    if movedP < N
                        newAlphas((movedP + 1):N,movedP) =...
                            alphas((movedP + 1):N,movedP) + dAng;
                        largerThan2pi =...
                                newAlphas((movedP + 1):N,movedP) > 2*pi;
                            newAlphas((movedP + 1):N,movedP) = ...
                                newAlphas((movedP + 1):N,movedP) -...
                                2*pi*largerThan2pi;
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
                        largerThan2pi =...
                                newThetas(movedP,1:(movedP-1)) > 2*pi;
                        newThetas(movedP,1:(movedP-1)) = ...
                                newThetas(movedP,1:(movedP-1)) -...
                                2*pi*largerThan2pi;
                    end

                    % recalculate the relevent column elements in theta matrix

                    if movedP < N
                        newThetas((movedP + 1):N,movedP) =...
                            newAlphas((movedP + 1):N,movedP) - atan(yi/xi);
                        largerThan2pi =...
                                newThetas((movedP + 1):N,movedP) > 2*pi;
                        newThetas((movedP + 1):N,movedP) = ...
                                newThetas((movedP + 1):N,movedP) -...
                                2*pi*largerThan2pi;
                    end
                end   
        end
    
        function newBettas =...
                reCalcBetta(movedParticles,bettas,newParticlesPosition,N)
                
                for i = 1:length(movedParticles)
                    movedP = movedParticles(i);

                    xi = newParticlesPosition(1,movedP);
                    yi = newParticlesPosition(2,movedP);
                    
                    newBettas = bettas;
                    
                    % recalculate the relevent row elements in dist matrix

                    if movedP > 1
                        newBettas(movedP,1:(movedP-1)) =...
                            atan((yi-newParticlesPosition(2,1:(movedP-1)))...
                            /(xi-newParticlesPosition(1,1:(movedP-1))));
                    end

                    % recalculate the relevent column elements in dist matrix

                    if movedP < N
                        newBettas((movedP + 1):N,movedP) =...
                            atan((yi-newParticlesPosition(2,(movedP + 1):N))...
                            /(xi-newParticlesPosition(1,(movedP + 1):N)));
                       
                    end
                end   
        end
        
        function dU = Uchange(movedParticle,dist,newDist,N,rc,m,...
                bettas,newBettas,angs,newAngs,angleDependence,ufunc,...
                dipoleStrength)
        % calculates the change in energy after a particle has moved
        
                % calculate the old energy for the relevant particle pairs
                
                % the relevent row:
                if movedParticle > 1
                    
                    if ~isempty(angleDependence)
%                         relAng = [];
%                         relAng(:,1) = alphas(movedParticle,1:(movedParticle - 1));
%                         relAng(:,2) = thetas(movedParticle,1:(movedParticle - 1));
                          relAng{1,1} =...
                              bettas(movedParticle,1:(movedParticle - 1));
                          relAng{1,2} = angs(1,1:(movedParticle - 1));
                          relAng{1,3} = ...
                              angs(1,movedParticle)*ones(1,movedParticle - 1);
                    else
                        relAng = [];
                    end
                    distrow = dist(movedParticle,1:(movedParticle - 1));
%                     oldUrow = ...
%                         pairU(distrow,rc,m,...
%                         'angleDependence',angleDependence,...
%                         'relativeCellAngles',relAng,'ufunc',ufunc);

                    dipolePairStrength =...
                        dipoleStrength(movedParticle)*dipoleStrength(1:(movedParticle-1));
                    oldUrow = ...
                        pairU(distrow,rc,m,...
                        'angleDependence',angleDependence,...
                        'relativeCellAngles',relAng,'ufunc',ufunc,...
                        'numOfrelativeCellAngles',3,...
                        'dipolePairStrength',dipolePairStrength);
                                         
                    
                else 
                    oldUrow = 0;
                end
                
                % relevent column:
                if movedParticle < N
                    if ~isempty(angleDependence)
                        relAng = [];
%                         relAng(:,1) =...
%                             alphas((movedParticle + 1):N,movedParticle);
%                         relAng(:,2) =...
%                             thetas((movedParticle + 1):N,movedParticle);
                        relAng{1,1} =...
                            bettas((movedParticle + 1):N,movedParticle)';
                        relAng{1,2} = angs(1,(movedParticle + 1):N);
                        relAng{1,3} = ...
                              angs(1,movedParticle)*ones(1,N - movedParticle);
                    else
                        relAng = [];
                    end
                    distcol = dist((movedParticle + 1):N,movedParticle)';
%                     oldUcol = ...
%                         pairU(distcol,rc,m,...
%                         'angleDependence',angleDependence,...
%                         'relativeCellAngles',relAng,'ufunc',ufunc);

                    dipolePairStrength =...
                        dipoleStrength(movedParticle)*dipoleStrength((movedParticle+1):N);
                    

                    oldUcol = ...
                        pairU(distcol,rc,m,...
                        'angleDependence',angleDependence,...
                        'relativeCellAngles',relAng,'ufunc',ufunc,...
                        'numOfrelativeCellAngles',3,...
                        'dipolePairStrength',dipolePairStrength);
                else 
                    oldUcol = 0;
                end

                oldU = oldUrow + oldUcol;
                
                % calculate the new energy for the relevant particle pairs
                
                if movedParticle > 1
                    if ~isempty(angleDependence)
                        relAng = [];
%                         relAng(:,1) =...
%                             newAlphas(movedParticle,1:(movedParticle - 1));
%                         relAng(:,2) =...
%                             newThetas(movedParticle,1:(movedParticle - 1));
                        relAng{1,1} =...
                            newBettas(movedParticle,1:(movedParticle - 1));
                        relAng{1,2} = newAngs(1,1:(movedParticle - 1));
                        relAng{1,3} = ...
                              newAngs(1,movedParticle)*ones(1,movedParticle - 1);
                    else
                        relAng = [];
                    end
                    distrow = newDist(movedParticle,1:(movedParticle - 1));
%                     newUrow = pairU(distrow,rc,m,...
%                         'angleDependence',angleDependence,...
%                         'relativeCellAngles',relAng,'ufunc',ufunc);

                    dipolePairStrength = ...
                        dipoleStrength(movedParticle)*dipoleStrength(1:(movedParticle-1));
                    newUrow = ...
                        pairU(distrow,rc,m,...
                        'angleDependence',angleDependence,...
                        'relativeCellAngles',relAng,'ufunc',ufunc,...
                        'numOfrelativeCellAngles',3,...
                        'dipolePairStrength',dipolePairStrength);
                else 
                    newUrow = 0;
                end
                
                if movedParticle < N
                    if ~isempty(angleDependence)
                        relAng = [];
%                         relAng(:,1) =...
%                             newAlphas((movedParticle + 1):N,movedParticle);
%                         relAng(:,2) =...
%                             newThetas((movedParticle + 1):N,movedParticle);
                        relAng{1,1} =...
                            newBettas((movedParticle + 1):N,movedParticle)';
                        relAng{1,2} = newAngs(1,(movedParticle + 1):N);
                        relAng{1,3} = ...
                              newAngs(1,movedParticle)*ones(1,N - movedParticle);
                    else
                        relAng = [];
                    end
                    distcol = newDist((movedParticle + 1):N,movedParticle)';
%                     newUcol = pairU(distcol,rc,m,...
%                         'angleDependence',angleDependence,...
%                         'relativeCellAngles',relAng,'ufunc',ufunc);

                    dipolePairStrength = ...
                        dipoleStrength(movedParticle)*dipoleStrength((movedParticle+1):N);
                    newUcol = ...
                        pairU(distcol,rc,m,...
                        'angleDependence',angleDependence,...
                        'relativeCellAngles',relAng,'ufunc',ufunc,...
                        'numOfrelativeCellAngles',3,...
                        'dipolePairStrength',dipolePairStrength);
                else 
                    newUcol = 0;
                end
                
                newU = newUrow + newUcol;
                
                % clculate the change in energy
                
                dU = newU - oldU;
        end
    
        function dV = Vchange(movedParticle,dist,newDist,N,rc,rho,m)
        % calculates the change in the virial after a particle has moved
                
                % calculate the old virial for the relevant particle pairs
                
                if movedParticle > 1
                    oldVrow = ...
                        calcVirial(dist(movedParticle,1:(movedParticle - 1))...
                        ,rho,12,m,N,rc);
                else 
                    oldVrow = 0;
                end
                
                if movedParticle < N
                    oldVcol = ...
                        calcVirial(dist((movedParticle + 1):N,movedParticle)...
                        ,rho,12,m,N,rc);
                else 
                    oldVcol = 0;
                end

                oldV = oldVrow + oldVcol;
                
                % calculate the new virial for the relevant particle pairs
                
                if movedParticle > 1
                    newVrow = calcVirial(newDist...
                            (movedParticle,1:(movedParticle - 1))...
                            ,rho,12,m,N,rc);
                else 
                    newVrow = 0;
                end
                
                if movedParticle < N
                    newVcol = calcVirial(newDist...
                        ((movedParticle + 1):N,movedParticle)...
                        ,rho,12,m,N,rc);
                else 
                    newVcol = 0;
                end
                
                newV = newVrow + newVcol;
                
                % clculate the change in the virial
                
                dV = newV - oldV;
        end
        
    function hardCoretouchs = check_if_cores_touch(dists,r)
        hardCoretouchs = ~isempty(find(tril(dists < 2*r,-1),1));
    end
            
end
