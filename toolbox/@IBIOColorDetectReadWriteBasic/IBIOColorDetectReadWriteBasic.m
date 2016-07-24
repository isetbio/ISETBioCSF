classdef IBIOColorDetectReadWriteBasic < IBIOColorDetectReadWrite
% A basic class for reading and writing IBIOColorDetect
% results.
% 
% 6/26/16  dhb  Started in on this

    % Public read/write properties
    properties
    end
    
    % Dependent properties, computed from other parameters
    properties
    end
        
    % Public, read-only properties.  These can be set by methods of the
    % parent class (that is, this class) but not by methods of subclasses.
    properties (SetAccess = private, GetAccess = public)  
    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected, GetAccess = public)
    end
    
    % Private properties. Only methods of the parent class can set or read these
    properties (Access = private)      
    end
    
    % Public methods
    %
    % Methods defined in separate files are public by default, so we don't
    % explicitly decare them.  But, if you wanted to write the whole body
    % of some short public method here, you could do so.
    methods (Access=public)
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    methods (Abstract, Access=public)

    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
        dirname = paramsToResponseGenerationDirName(obj,rparams);        
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)

    end
    
end
