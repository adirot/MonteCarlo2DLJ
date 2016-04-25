classdef wtf
    
    properties
        wtff;
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
        
    end
end
