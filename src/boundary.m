classdef boundary
    % This class contain the information of the measurement boundary
    
    properties (SetAccess = public)
        % center of domain
        xb
        yb
        % radius of domain
        rb
    end
     
    methods
        % constructor
        function BD = boundary(xb,yb,rb)
            BD.xb = xb;
            BD.yb = yb;
            BD.rb = rb;
        end
    end
end