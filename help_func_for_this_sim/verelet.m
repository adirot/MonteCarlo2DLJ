classdef verelet
   properties
      neighbors, neighborsindx, neighborsindy, ...
          dispacements , summaxdisplace;
   end
   methods
       
       % constructor
      function obj = verelet(dists,rl,N)
            obj.neighbors = logical((dists <= rl).*(dists > 0));
            [obj.neighborsindx,obj.neighborsindy] = find(obj.neighbors);
            obj.dispacements = zeros(1,N);
            obj.summaxdisplace = 0;
      end
      
      function obj = updateDisplace(obj,movedParticle,diplacement)
            obj.dispacements(1,movedParticle) = ...
                obj.dispacements(1,movedParticle) + diplacement;
            [max1,Ind] = max(obj.dispacements);
            
            % get the second maximal value
            nomax1 = [obj.dispacements(1:(Ind-1))...
                obj.dispacements((Ind+1):end)];
            max2 = max(nomax1);
            
            obj.summaxdisplace = max1 + max2;
      end
      
   end
end