classdef chin
    %CHIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m
        n
    end
    
    methods
        function obj = chin(m, n)
            obj.m = m;
            obj.n = n;
        end
        
        function F = Qed(obj, u)
            F = u/(obj.m + obj.n*u);
        end
        
        function w = cedimento(obj, F)
            w = obj.m*F/(1-obj.n*F);
        end
    end
    
end

