%------ SparProperties -----%
%                           %
%     Cameron Robertson     %
%        Feb 28, 2011       %
%                           %
%---------------------------%

% Computes the structural properties of a CFRP spar given the diameter, d,
% wrap angle, theta, number of tube layers, nTube, and number of cap
% strips, nCap. Properties are computed for each element along the span,
% where the node locations are given by yN.

% To check (notes from Todd)
% - So we use MTM28-M46J for everything now right? (high mod, high strength stuff)
% - To make nCap a continuous variable, I linearly interpolated
% - mSpar(s) is the mass of the element (linear density)*dy(s)
% - check that these values agree with measured masses in HPO spreadsheet

function [EIx, EIz, EA, GJ, mSpar] = SparProperties(yN, d, theta, nTube, nCap, lBiscuit, flags)

% MATERIAL PROPERTIES
[RHO_TUBE, T_PLY_TUBE, E_11_TUBE, E_22_TUBE, G_12_TUBE, V_12_TUBE,...
 ULTIMATE_11_TENS, ULTIMATE_11_COMP,ULTIMATE_22_TENS, ULTIMATE_22_COMP, ULTIMATE_12]...
 = PrepregProperties(flags.CFRPType);

[RHO_CAP, T_PLY_CAP, E_11_CAP, E_22_CAP, G_12_CAP, V_12_CAP,...
 ULTIMATE_11_TENS, ULTIMATE_11_COMP,ULTIMATE_22_TENS, ULTIMATE_22_COMP, ULTIMATE_12]...
 = PrepregProperties(flags.CFRPType);

% Other Material Properties
RHO_BALSA = 1.60E2; % Materials Database, Balsa
RHO_STRUCTURAL_FOAM = 3.1E1;   % Materials Database, Rohacell 31 IG-F


% Spar Geometry Constants
ALPHA_BASE = 45*(pi/180); % Half-angle of base cap ply width, based on 90-degree rule
PLY_TAPER_RATIO = 0.05; % Based on HPO wing, 0.05 for 3" spar, ~0.01 for smaller diameters

% Sizing Factors
AF_biscuit = 1.133;
thickness_biscuit_core = (1/4)*(0.0254); % 1/4"
thickness_biscuit_plate = 2*2*(3/32)*(0.0254); % 2 3/32" plies per face
biscuit_face_fraction = 0.70; % Percent of biscuit surface area occupied by plate

% Determine the length of each spar element
Ns = length(yN)-1; %number of elements
dy = zeros(Ns,1);
for s=1:Ns
    dy(s) = yN(s+1) - yN(s); %length of each element
end

% Pre-allocate EIx, EIz, EA, GJ, mSpar
EIx = zeros(Ns,1);
EIz = zeros(Ns,1);
EA = zeros(Ns,1);
GJ = zeros(Ns,1);
mSpar = zeros(Ns,1);

% Pre-populate Q matrix for composite angle transformation
Q = zeros(3);   % Preallocate Q-matrix
    
% Matrix of elastic constants
Q(1,1) = E_11_TUBE;
Q(2,2) = E_22_TUBE;
Q(1,2) = E_22_TUBE*V_12_TUBE;
Q(2,1) = Q(1,2);
Q(3,3) = G_12_TUBE;

% Compute EIx, EIz, EA, GJ, mSpar for each element 1:Ns
for s = 1:Ns
    % Determine E_xx_tube and G_xy_tube using Q_Transform Code
    x = theta(s); % Composite angle (in radians)

    % Transformation matrix
    T = [cos(x)^2 sin(x)^2  2*sin(x)*cos(x);
        sin(x)^2  cos(x)^2  -2*sin(x)*cos(x);
        -sin(x)*cos(x)  sin(x)*cos(x)   (cos(x)^2)-(sin(x)^2)];

    % Transform the elastic constants using the transformation matrix to obtain the
    % elastic constants at the composite angle.
    Qbar = (T\Q)/T';

    % Breakout tube elastic constants at the transformed angle
    E_xx_tube = Qbar(1,1);
    G_xy_tube = Qbar(3,3);

    % I Tube
    I_tube = pi*(((d(s)/2) + (1/2)*nTube(s)*T_PLY_TUBE)^3)*(nTube(s)*T_PLY_TUBE); % HPO2:194, I_tube = pi*r_avg^3*t
    % A Tube
    A_tube = pi*((((d(s)/2) + nTube(s)*T_PLY_TUBE)^2)-((d(s)/2)^2)); % HPO2:139

    % Initialize Ix_cap(s) (In-plane), Iz_cap(s) (Out-of-plane), A_cap
    Ix_cap_0 = 0;
    Iz_cap_0 = 0;
    A_cap_0 = 0;
    Ix_cap_1 = 0;
    Iz_cap_1 = 0;
    A_cap_1 = 0;

    % Linearly interpolate between discrete values of nCap
    if nCap(s) < 0
        nCap(s) = 0;
    end
    nCap_0 = floor(nCap(s));
    nCap_1 = nCap_0 + 1;
    
    if nCap_0 ~=0
        width_ply_0 = zeros(1,nCap_0);
    else
        width_ply_0 = zeros(1,1);
    end
    width_ply_1 = zeros(1,nCap_1);
    
    for i = 1:nCap_0
        % width_ply_0 = [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032]; % Spar 2009-1 cap widths
        % width_ply_0 = [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032 0.028 0.024]; % Spar 2009-2 & 2009-3 cap widths
        width_ply_0(i) = d(s)*ALPHA_BASE - (i-1)*(PLY_TAPER_RATIO*d(s)); % Sets base cap width from 90 degree rule, tapers each subsequent layer based on tube diameter
        alpha_ply = width_ply_0(i)/(2*((d(s)/2) + nTube(s)*T_PLY_TUBE + (i-(1/2))*T_PLY_CAP)); % HPO2:195, alpha_ply = (width_ply)/(2*r_ply)
        Ix_cap_0 = Ix_cap_0 + 2*(alpha_ply*(((d(s)/2) + (i - (1/2))*(T_PLY_CAP))^3)*T_PLY_CAP); % HPO2:194
        Iz_cap_0 = Iz_cap_0 + 2*((alpha_ply + (1/2)*sin(2*alpha_ply))*(((d(s)/2) + (i - (1/2))*(T_PLY_CAP))^3)*T_PLY_CAP); % HPO2:194
        A_cap_0 = A_cap_0 + 2*(2*alpha_ply*((d(s)/2) + (i - (1/2))*(T_PLY_CAP))*T_PLY_CAP); % HPO2:194
    end
    for i = 1:nCap_1
        % width_ply_1 = [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032]; % Spar 2009-1 cap widths
        % width_ply_1 = [0.060 0.056 0.052 0.048 0.044 0.040 0.036 0.032 0.028 0.024]; % Spar 2009-2 & 2009-3 cap widths    
        width_ply_1(i) = d(s)*ALPHA_BASE - (i-1)*(PLY_TAPER_RATIO*d(s)); % Sets base cap width from 90 degree rule, tapers each subsequent layer based on tube diameter
        alpha_ply = width_ply_1(i)/(2*((d(s)/2) + nTube(s)*T_PLY_TUBE + (i-(1/2))*T_PLY_CAP)); % HPO2:195, alpha_ply = (width_ply)/(2*r_ply)
        Ix_cap_1 = Ix_cap_1 + 2*(alpha_ply*(((d(s)/2) + (i - (1/2))*(T_PLY_CAP))^3)*T_PLY_CAP); % HPO2:194
        Iz_cap_1 = Iz_cap_1 + 2*((alpha_ply + (1/2)*sin(2*alpha_ply))*(((d(s)/2) + (i - (1/2))*(T_PLY_CAP))^3)*T_PLY_CAP); % HPO2:194
        A_cap_1 = A_cap_1 + 2*(2*alpha_ply*((d(s)/2) + (i - (1/2))*(T_PLY_CAP))*T_PLY_CAP); % HPO2:194
    end
    Ix_cap = Ix_cap_0 + (Ix_cap_1-Ix_cap_0)*(nCap(s)-nCap_0);
    Iz_cap = Iz_cap_0 + (Iz_cap_1-Iz_cap_0)*(nCap(s)-nCap_0);
    A_cap = A_cap_0 + (A_cap_1-A_cap_0)*(nCap(s)-nCap_0);
    
    % GJ_spar
    r_tube_avg = (d(s)/2) + (1/2)*nTube(s)*T_PLY_TUBE;
    r_cap_avg = ((d(s)/2) + nTube(s)*T_PLY_TUBE) + (1/2)*nCap(s)*T_PLY_CAP;
    width_cap_avg = (1/2)*(mean(width_ply_0) + mean(width_ply_1));
    t_cap_avg = nCap(s)*T_PLY_CAP; % Don't know if the average is actually necessary here, and this value is close regardless
    r_spar_avg = ((pi*r_tube_avg - width_cap_avg)/(pi*r_tube_avg))*(r_tube_avg) + ((width_cap_avg)/(pi*r_tube_avg))*(r_cap_avg); %HPV1:71
    GJ(s) = (4*(pi^2)*(r_spar_avg^4))/(((2*pi*r_tube_avg - 2*width_cap_avg)/(G_xy_tube*nTube(s)*T_PLY_TUBE)) + ((2*width_cap_avg)/(G_xy_tube*nTube(s)*T_PLY_TUBE + G_12_CAP*t_cap_avg))); %HPV1:71

    % Biscuit Mass
    mass_biscuit =  (AF_biscuit)*(dy(s)/lBiscuit(s))*((pi*(d(s)/2)^2)*(RHO_BALSA*thickness_biscuit_plate*biscuit_face_fraction + RHO_STRUCTURAL_FOAM*thickness_biscuit_core));
    
    % Determine Total Spar Properties
    EIx(s) = E_xx_tube*I_tube + E_11_CAP*Ix_cap;
    EIz(s) = E_xx_tube*I_tube + E_11_CAP*Iz_cap;
    EA(s) = E_xx_tube*A_tube + E_11_CAP*A_cap;
    mSpar(s) = (A_tube*RHO_TUBE + A_cap*RHO_CAP)*dy(s) + mass_biscuit;

end