% Class that contains the properties of a "change target", specifically 
% made for the discrimination go cue center out task. The idea is that the
% target will change it's color/image from some trigger, so this object 
% adds "FinalColor" "FinalImage" separate fields
%
% If certain properties don't work or appear missing, ensure that this file
% is the one being used for the class instantiation, and not another older
% copy elsewhere.
%
% Adam Smoulder, 5/11/22 

classdef ChangeTarget < Target
    properties
        FinalColor
        FinalImage
    end
    methods
        % Constructor
        function obj = ChangeTarget(varargin)
            % Defaults
            obj.FinalColor = [0 0 0];          % RGB on [0 1]
            obj.FinalImage = [];               % string (file name, no .png needed)
            
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
            if nargin > 7
                obj.FinalColor = varargin{8};
            end
            if nargin > 8
                obj.FinalImage = varargin{9};
            end
        end
    end
end