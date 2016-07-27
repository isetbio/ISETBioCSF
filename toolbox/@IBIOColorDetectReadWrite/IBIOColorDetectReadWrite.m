classdef IBIOColorDetectReadWrite < handle
% The abstract parent class for reading and writing IBIOColorDetect
% results.
% 
% 6/26/16  dhb  Started in on this

    % Public read/write properties
    properties
%         % Name of programs that is are charge of writing/reading data
%         writeProgram = '';
%         readProgram = '';
%         
%         % Cell array list of parameters of parent output tree.  The current data should
%         % be associated with this parent tree.  Can be empty, in which
%         % case the current data is taken to be at the top level.
%         parentParamsList = {};
%         
%         % Cell array list of parameters of the current input or output.
%         currentParamsList = {};
%         
%         % Extra information that we want to store, expressed as a struct.
%         extraParams = [];
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
        % Parameters describing the system and code that led to the data,
        % obtained via a method of the class.
        systemParams = [];
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
        % Write data
        write(obj,name,data,varargin);
        
        % Get unique identifier for data
        fileid = getid(obj,name,varagin);
        
        % Get unique identifier for data
        data = read(obj,fileid,varargin);
        
        % Get me a scratch directory
        tempDir = tempdir(obj,varargin);
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end
