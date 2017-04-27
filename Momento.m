classdef Momento < handle
    %MOMENTO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sup
        inf
        mx
        dx = .001;  % risoluzione di mx
        L = 1;  % [m] lunghezza della trave
    end
    
    methods
        
        function obj = Momento(sup, inf, L)
            obj.sup = sup;
            obj.inf = inf;
            obj.L = L;
            obj = obj.setMx;
        end
        
        function obj = setMx(obj)
            [a, f] = parabola3P([0, obj.sup], [obj.L/2, obj.inf], [obj.L, obj.sup]);
            obj.mx = f(0:obj.dx:obj.L);
        end
        
        function obj = setVal(obj, val)
            if val > 0
                obj.inf = val;
            else
                obj.sup = val;
            end
            obj = obj.setMx;
        end
    end
    
end

