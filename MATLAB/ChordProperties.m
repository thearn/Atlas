%----- ChordProperties -----%
%                           %
%     Cameron Robertson     %
%        Feb 28, 2011       %
%                           %
%---------------------------%

% Computes the mass of the ribs, trailing edge, LE-sheeting and covering
% given the chord. The mass is computed for each spanwise element. If
% flagTESpar is 1 then the mass of a trailing edge spar and in-plane truss
% is added as well.


function [mChord, xCGChord] = ChordProperties(yN, c, d, flagGWing, xtU)

% Determine the length of each spar element
Ns = length(yN)-1; %number of elements
dy = zeros(Ns,1);
for s=1:Ns
    dy(s) = yN(s+1) - yN(s); %length of each element
end

% Pre-allocate mChord
mChord = zeros(Ns,1);
xCGChord = zeros(Ns,1);

% Material Properties
RHO_BALSA = 1.60E2; % Materials Database, Balsa
RHO_BASSWOOD = 3.87E2;   % Materials Database, Basswood
RHO_CARBON = 1.6106E3; % Materials Database, MTM28-M46J 140 37 %RW 12"
T_PLY_CARBON = 1.47E-4; % Materials Database, MTM28-M46J 140 37 %RW 12"
RHO_EPS = 1.6E1; % Materials Database, CURRENTLY TYPE I EPS FROM PICOTTE
RHO_KEVLAR = 1.3E3;  % Materials Database, 1.8oz Kevlar Weave
T_PLY_KEVLAR = 1.27E-4; % Materials Database, 1.8oz Kevlar Weave
RHO_MYLAR = 1.4E3;  % Materials Database, Melinex S Type, 48 gauge
T_MYLAR = 1.2E-5;   % Materials Database, Melinex S Type, 48 gauge 
RHO_STRUCTURAL_FOAM = 3.1E1;   % Materials Database, Rohacell 31 IG-F
RHO_XPS = 2.467E1; % OC Data, Celfort 200 4"
RHO_STEEL_WIRE = 7.85E3; % Density of Piano Wire (Spring Steel)

% Airfoil Geometry Properties (DT1068156, HPO root airfoil)
AREA_AIRFOIL = 0.08233; % Area of airfoil, in percentage chord^2, from HPO Wing Design Spreadsheet, cell BU262 (multiply by chord^2 for area in m^2)
PERIMETER_AIRFOIL = 2.06814; % Perimeter of airfoil, in percentage chord, from HPO Wing Design Spreadsheet, cell BU263 (method assumes top and bottom perimeter are same, not exactly correct)
XCG_AIRFOIL = 0.39003; % Xcg of airfoil

% Wing Design Variables
thickness_rib = 0.005; %(1/4)*(0.0254); % 1/4"
rib_spacing = 12*(0.0254); %12"
thickness_LE_sheeting = 0.003; % 3mm
thickness_rib_caps = (1/32)*(0.0254); % 1/32"
percent_rib_caps = 0.97;    % Percent chord of rib caps, assumed same as HPO, see HPVDT1:83
thickness_spar_plate = (1/32)*(0.0254); % 1/32"
percent_radius_spar_plate = 0.0632; % Rib plate radius as a percent of chord (assumes full depth of a 15% thick airfoil)
percent_LE_sheeting_top = 0.60; % Percent chord of LE sheeting on top
percent_LE_sheeting_top_GWing = xtU; % Percent chord of LE sheeting on top for Gossamer Wing
percent_LE_sheeting_bottom = 0.10; % Percent chord of LE sheeting on bottom
d_TE_spar = (3/4)*(0.0254); % 3/4"
nTube_TE_spar = 4;  % Number of layers in TE spar
d_comp_member = (3/4)*(0.0254); % 3/4"
nTube_comp_member = 4;  % Number of layers in compression members
spacing_comp_member = 3;    % Spacing of In-Plane Truss compression members, 3m
percent_length_comp_member = 0.782-0.25; % Length, in percent chord, of in-plane-truss compression member spacing, see HPVDT1:83
d_cross_bracing = 0.001; % Cross-bracing wire diameter, same as for HPO
height_TE_foam = 7.0E-3; % Height of triangular TE foam, assuming same as HPO
length_TE_foam = 5.05E-2; % Length of triangular TE foam, assuming same as HPO
percent_height_TE_plate = 0.05;    % Height of triangular TE balsa plate, percent chord, see HPVDT1:83
percent_length_TE_plate = 0.24;    % Length of triangular TE balsa plate, percent chord, see HPVDT1:83
thickness_TE_plate = (1/32)*(0.0254); % 1/32"
diameter_piano_wire = (0.022)*(0.0254); %0.022", same as Gossamer Albatross
spar_location = 0.25; % Location of main spar, in percent chord, assumed at 25%

% Adjustment Factors
AF_ribs = 1.569;
AF_TE = 1.602;
AF_leading_edge_sheeting = 1.302;
AF_leading_edge_sheeting_GWing = 1.210; % Determined from HPO, but without XPS riblet mass
AF_covering = 0.972;
AF_TE_spar = 0.529;
AF_comp_member = 1.000;
AF_cross_bracing = 4.603;

% Compute mChord for each element 1:Ns
for s = 1:Ns
    % Rib Mass & Xcg
    mass_rib_foam = RHO_EPS*(thickness_rib*((c(s)^2)*AREA_AIRFOIL));  
    mass_rib_caps = RHO_BASSWOOD*(thickness_rib*thickness_rib_caps*(percent_rib_caps*PERIMETER_AIRFOIL*c(s)));
    mass_rib_plate_spar = RHO_BALSA*(thickness_spar_plate*(pi*(((c(s)*percent_radius_spar_plate)^2)-((d(s)/2)^2)))); % Assumes circle with interior cutout for spar, with radius as percent of chord
    if flagGWing == 0
        mass_rib_plate_TE = RHO_BALSA*(thickness_TE_plate*((1/2)*(c(s)^2)*percent_height_TE_plate*percent_length_TE_plate));
        mass_rib = AF_ribs*((dy(s)/rib_spacing)*(mass_rib_foam + mass_rib_caps + 2*mass_rib_plate_spar + 2*mass_rib_plate_TE));
        Xcg_rib = XCG_AIRFOIL;
    else
        mass_rib = AF_ribs*(0.66)*((dy(s)/rib_spacing)*(mass_rib_foam + mass_rib_caps + 2*mass_rib_plate_spar)); % Assumes 1 1/3 Gossamer ribs for every 2 Daedalus ribs
        Xcg_rib = ((XCG_AIRFOIL + (2/3)*(spar_location))/2); % Assumes riblet is a triangle with spar at 1/4 chord
    end
    
    % Trailing Edge Mass & Xcg
    if flagGWing == 0
        mass_TE = AF_TE*(RHO_STRUCTURAL_FOAM*dy(s)*((1/2)*(length_TE_foam*height_TE_foam)) +...
            RHO_KEVLAR*dy(s)*T_PLY_KEVLAR*(length_TE_foam + height_TE_foam + sqrt(length_TE_foam^2 + height_TE_foam^2))); % (Constant profile, thickness, Rohacell 31 IG-F assumed)
        Xcg_TE = (c(s)-(2/3)*length_TE_foam)/c(s);
    else
        mass_TE = RHO_STEEL_WIRE*dy(s)*pi*((diameter_piano_wire/2)^2);
        Xcg_TE = 1;
    end
    
    % Leading Edge Sheeting Mass 
    if flagGWing == 0
        mass_LE_sheeting = AF_leading_edge_sheeting*(RHO_XPS*dy(s)*c(s)*(percent_LE_sheeting_top + percent_LE_sheeting_bottom)*PERIMETER_AIRFOIL*thickness_LE_sheeting);
        Xcg_LE_sheeting = (1/2)*((1/2)*(percent_LE_sheeting_top)+(1/2)*(percent_LE_sheeting_bottom));
    else
        mass_LE_sheeting = AF_leading_edge_sheeting_GWing*(RHO_EPS*dy(s)*c(s)*(percent_LE_sheeting_top_GWing(s) + percent_LE_sheeting_bottom)*PERIMETER_AIRFOIL*thickness_LE_sheeting) + ...
            RHO_STEEL_WIRE*dy(s)*pi*((diameter_piano_wire/2)^2); % Includes LE wire
        Xcg_LE_sheeting = (1/2)*((1/2)*(percent_LE_sheeting_top_GWing(s))+(1/2)*(percent_LE_sheeting_bottom));
    end

    % Covering Mass (Mylar, use perimeter estimates from HPO airfoils)
    mass_covering = AF_covering*(RHO_MYLAR*dy(s)*c(s)*PERIMETER_AIRFOIL*T_MYLAR);
    Xcg_covering = 0.5;

    % Trailing Edge Spar and In-Plane Truss Mass & Xcg (TE Spar, Kevlar cross-bracing, compression members)
    if flagGWing == 0
        mass_TE_spar = AF_TE_spar*(RHO_CARBON*dy(s)*(pi*(((d_TE_spar/2) + nTube_TE_spar*T_PLY_CARBON)^2 - (d_TE_spar/2)^2)));
        Xcg_TE_spar = 0.9;
        mass_comp_member = AF_comp_member*(dy(s)/spacing_comp_member)*RHO_CARBON*(c(s)*percent_length_comp_member)*(pi*(((d_comp_member/2) + nTube_comp_member*T_PLY_CARBON)^2 - (d_comp_member/2)^2));
        Xcg_comp_member = 0.25 + (1/2)*(percent_length_comp_member);
        mass_cross_bracing = AF_cross_bracing*RHO_KEVLAR*(dy(s)/spacing_comp_member)*2*(pi*(d_cross_bracing/2)^2)*((((c(s)*percent_length_comp_member)^2) + ((spacing_comp_member)^2))^(1/2));
        Xcg_cross_bracing = 0.25 + (1/2)*(percent_length_comp_member);
    else
        mass_TE_spar = 0;
        Xcg_TE_spar = 0;
        mass_comp_member = 0;
        Xcg_comp_member = 0;
        mass_cross_bracing = 0;
        Xcg_cross_bracing = 0;
    end
    
    % Total Chord Mass
    mChord(s) = mass_rib + mass_TE + (mass_LE_sheeting) + (mass_covering)...      % Rib mass + TE mass + LE sheeting mass + covering mass...
        + (mass_TE_spar + mass_comp_member + mass_cross_bracing);  % + TE spar mass
    
    xCGChord(s) = ((mass_rib*Xcg_rib) + (mass_TE*Xcg_TE) + (mass_LE_sheeting*Xcg_LE_sheeting) + ...
        (mass_covering*Xcg_covering) + (mass_TE_spar*Xcg_TE_spar) + (mass_comp_member*Xcg_comp_member) + ...
        (mass_cross_bracing*Xcg_cross_bracing)) / (mass_rib + mass_TE + mass_LE_sheeting + ...
        mass_covering + mass_TE_spar + mass_comp_member + mass_cross_bracing);
end