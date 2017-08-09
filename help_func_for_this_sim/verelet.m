classdef verelet
% verelet object: saves a nieghbor list for each particle, and monitors
% the displacement of each particle.

% properties:
% N - number of particles.
% rc - cutoff distance for energy calculation of Monte Carlo.
% rl - cutoff distance for the neighbors list.
% neighbors - N by N logical: if i and j are neighbors (thier distance is
%             less than rl) that the i,j place in the matrix is 1.
% dispacements - counts the total displacement of each particle during the
%                MC run.
% sumMaxDisplacement - the sum of the two maximal displacements.

% functions:

%   constructor:
%       Creates a verelet object. Input: distances matrix of all particle
%       cuples, rl, rc, N. Finds the neibors list for each particle, and
%       creates the 'neighbors' matrix. Initializes all displacements to
%       zero.

%   updateDisplacement: 
%        Adds a displacement 'displacement' to a particle
%        'movedParticle'.

%   clearDisplacements: 
%       Initializes all displacements to zero.

%   needs2beInitialized:
%       True if sumMaxDisplacement is larger than rl-rc (the neighbors list
%       may need to be initialized if two particles have gone more than
%       rl-rc from each other).


   properties
      neighbors, neighborsindx, neighborsindy, ...
          dispacements , sumMaxDisplacement, rl, rc, N;
   end
   methods
       
       % constructor
      function obj = verelet(dists,rl,rc,N)
            obj.neighbors = logical((dists <= rl).*(dists > 0));
            [obj.neighborsindx,obj.neighborsindy] = find(obj.neighbors);
            obj.dispacements = zeros(1,N);
            obj.sumMaxDisplacement = 0;
            obj.rl = rl;
            obj.rc = rc;
            obj.N = N;
      end
      
      function obj = updateDisplacement(obj, movedParticleInd, diplacements)
            obj.dispacements(1, movedParticleInd) = ...
                obj.dispacements(1, movedParticleInd) + diplacements;
            [max1, Ind] = max(obj.dispacements);
            
            % get the second maximal value
            nomax1 = [obj.dispacements(1:(Ind-1))...
                obj.dispacements((Ind+1):end)];
            max2 = max(nomax1);
            
            obj.sumMaxDisplacement = max1 + max2;
      end
      
      function obj = clearDisplacements(obj)
          obj.dispacements = zeros(1,obj.N);
      end
      
      function answer = needs2beInitialized(obj)
          answer = obj.sumMaxDisplacement > (obj.rl - obj.rc);
      end
   end
end