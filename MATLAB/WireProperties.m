%------ WireProperties -----%
%                           %
%     Cameron Robertson     %
%       April 23, 2012      %
%                           %
%---------------------------%
%
% This code outputs the material properties of wire options used in the HPH
% Project. 
% 
% Input Variables:
%   	type_flag - flag indicating wire material to be used
%           1 - Pianowire (High-Tensile Steel)
%           2 - Vectran
%
% Output Variables:
%
%       The properties contained in this function are:
%       RHO - Material density, calculated for woven material if necessary [Km/m^3]
%       E - Elastic modulus [Pa]
%       ULTIMATE - Failure strength [Pa]
%  
% Also included is D, an array of the available diameters of each material
% [m]

function [RHO, E, ULTIMATE] = WireProperties(type_flag)


if type_flag == 1
    % 1 - Pianowire (High-Tensile Steel)
    RHO = 7.85E+3; % From ASTM228 Standard, Accessed Online at MatWeb
    E = 2.10E+11; % From ASTM228 Standard, Accessed Online at MatWeb
    ULTIMATE = 2.62E+9; % FortePiano.com
    D = [.35 .40 .45 .50 .55 .60 .65 .72 ]; % FortePiano.com

elseif type_flag == 2
    % 2 - Vectran
    RHO = 1.1065E+3;
    E = 3.921E+10;  % Effective E of braided Vectran, averaged from data supplied by Cortland Cable
    ULTIMATE = 9.828E8; % Based on minimum tensile strength supplied by Cortland Cable
    D = [1.5 2 ];

    
    
end
