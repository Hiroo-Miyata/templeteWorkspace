% Class that contains the properties of a reach target (though you can also
% use this for a center target, just ignore reward probs). Classes are
% pretty much used just like structures in terms of calling their
% properties and stuff.
%
% If certain properties don't work or appear missing, ensure that this file
% is the one being used for the class instantiation, and not another older
% copy elsewhere.
%
% Adam Smoulder, 3/1/21 (edit 10/21/21 to add CatchTarget)

classdef Target
    properties
        Name
        Radius
        Location
        Color
        Reward
        InscribedImage
        CatchTarget
    end
    methods
        % Constructor
        function obj = Target(varargin)
            % Defaults
            obj.Name = 'Target';
            obj.Radius = 10;                % in mm
            obj.Location = [0 0];           % in mm (w.r.t. center)
            obj.Color = [0.5 0.5 0.5];      % RGB on [0 1]
            obj.Reward = 100;               % in ms
            obj.InscribedImage = [];        % string (file name, no .png needed)
            obj.CatchTarget = 0;
            
            % Arguments parsed in order
            if nargin > 0
                obj.Name = varargin{1};
            end
            if nargin > 1
                obj.Radius = varargin{2};
            end
            if nargin > 2
                obj.Location = varargin{3};
            end
            if nargin > 3
                obj.Color = varargin{4};
            end
            if nargin > 4
                obj.Reward = varargin{5};
            end
            if nargin > 5
                obj.InscribedImage = varargin{6};
            end
            if nargin > 6
                obj.CatchTarget = varargin{7};
            end
        end
    end
end