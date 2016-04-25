classdef myclass
        properties
            prop;
        end
        methods
            % constructor
            function obj = myclass()
                 obj.prop = 0;
            end
            % add function 
            function obj = add(obj,a)
                obj.prop = obj.prop + a;
            end
        end
    end
