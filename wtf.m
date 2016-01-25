classdef wtf
    
    properties
        wtff,derrr;
    end
    
    methods
        % constructor
        function obj = wtf()
            obj.wtff = 0;
        end
        
        % add a function
        function obj = add(obj,a)
            obj.wtff = obj.wtff+a; 
        end
        
        function obj = addit(obj,a)
            obj = obj.add(a); % this changes obj
            obj.add(a); % this doesn't change obj, but why??
        end
    end
end