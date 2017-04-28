classdef sezione
    %SEZIONE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        b
        h
        x0
        y0
        section
        % valori di default
        numRetX = 1
        numRetY = 100
        area_sezione
        peso_sezione
        Jx
        Jy
    end
    
    properties (SetAccess = protected)
        peso_cls = 25;
    end
    
    methods
        
        function obj = sezione(b, h, x0, y0, numRetX, numRetY)
            if ~(length(b) == length(h)) || nargin > 2 && ~(length(b) == length(x0) && length(h) == length(y0))
                error('b, h, x0, y0 must be of the same length')
            end
            obj.b = b;
            obj.h = h;
            
            % se il numero di argomenti == 2 allora posiziona i rettangoli
            % uno sopra l'altro, lungo l'asse y
            if nargin == 2


                y = 0;
                for i = 1:length(h)
                    y = y + h(i)/2;
                    obj.x0(i) = 0;
                    obj.y0(i) = y;
                    y = h(i);
                end
                
            elseif nargin == 4
                    obj.x0 = x0;
                    obj.y0 = y0;
                
            elseif nargin == 6
                obj.x0 = x0;
                obj.y0 = y0;
                obj.numRetX = numRetX;
                obj.numRetY = numRetY;
            else
                error('invalid number of arguments: %d \n(2, 4 or 6 arguments are allowed)', nargin)
            end
            obj = obj.discretizza;
            obj = obj.area;
            obj = obj.peso;
            obj = obj.inerzia;
        end
               
        function obj = discretizza(obj)
            %ritorna la matrice [xm, ym, dx, dy]
            obj.section = rettangolo(obj.b, obj.h, obj.x0, obj.y0, obj.numRetX, obj.numRetY);
        end
        
        function obj = area(obj)
            obj.area_sezione = sum(obj.b .* obj.h);
        end
        
        function obj = peso(obj)
            obj.peso_sezione =  obj.area.area_sezione * 1e-6 * obj.peso_cls;
        end
        
        function obj = inerzia(obj)
            obj.Jx = obj.b^3 * obj.h / 12;  % direzione x (lungo b)
            obj.Jy = obj.b * obj.h^3 / 12;  % direzione y (lungo y)
        end
    end
    
end

