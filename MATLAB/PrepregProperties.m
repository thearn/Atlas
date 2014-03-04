%---- PrepregProperties ----%
%                           %
%     Cameron Robertson     %
%       April 21, 2012      %
%                           %
%---------------------------%
%
% This code outputs the material properties of the CFRP prepregs used in
% the HPH Project. 
% 
% Input Variables:
%   	type_flag - flag indicating prepreg to be used
%           1 - NCT301-1X HS40 G150 33 +/-2%RW
%           2 - HexPly 6376 HTS-12K 268gsm 35%RW
%           3 - MTM28-1/IM7-GP-145gsm 32 +/-2%RW
%           4 - HexPly 8552 IM7 160gsm 35%RW
%           5 - NCT304-1 HR40 G80 40 +/-2%RW
%           6 - MTM28-1B/M46J-140-37%RW
%
% Output Variables:
%
%       The properties contained in this function are:
%       RHO - Prepreg laminate density [Km/m^3]
%       T_PLY - Ply thickness [m]
%       E_11 - Elastic modulus in the material 11-axis [Pa]
%       E_22 - Elastic modulus in the material 22-axis [Pa]
%       G_12 - Shear modulus in the material 12-axis [Pa]
%       ULTIMATE_11_TENS - Failure strength in the material 11-axis, tension
%       ULTIMATE_11_COMP - Failure strength in the material 11-axis,
%           compression
%       ULTIMATE_22_TENS - Failure strength in the material 22-axis, tension
%       ULTIMATE_22_COMP - Failure strength in the material 22-axis,
%           compression
%       ULTIMATE_12 - Failure strength in shear in the material
%           12-axis
%       V_12 - Poisson's ratio of the material in the 12-axis [N/A]

function [RHO, T_PLY, E_11, E_22, G_12, V_12, ULTIMATE_11_TENS, ULTIMATE_11_COMP, ULTIMATE_22_TENS, ULTIMATE_22_COMP, ULTIMATE_12] = PrepregProperties(type_flag)


if type_flag == 1
    % 1 - NCT301-1X HS40 G150 33 +/-2%RW
    RHO = 1.5806E+03;
    T_PLY = 1.4173E-04;
    E_11 = 2.1066E+11;
    E_22 = 9.0101E+09;
    G_12 = 3.3426E+09;
    V_12 = 0.27;
    ULTIMATE_11_TENS = 2.0469E+09;
    ULTIMATE_11_COMP = 8.4593E+08;
    ULTIMATE_22_TENS = 3.4903E+07;
    ULTIMATE_22_COMP = 2.0806E+08;
    ULTIMATE_12 = 1.1974E+08;

elseif type_flag == 2
    % 2 - HexPly 6376 HTS-12K 268gsm 35%RW
    RHO = 1.5763E+03;
    T_PLY = 2.6157E-04;
    E_11 = 1.2329E+11;
    E_22 = 9.2393E+09;
    G_12 = 2.8958E+09;
    V_12 = 0.36;
    ULTIMATE_11_TENS = 2.0342E+09;
    ULTIMATE_11_COMP = 1.2689E+09;
    ULTIMATE_22_TENS = 5.7517E+07;
    ULTIMATE_22_COMP = 2.0884E+08;
    ULTIMATE_12 = 5.4882E+07;

elseif type_flag == 3
    % 3 - MTM28-1/IM7-GP-145gsm 32 +/-2%Rw
    RHO = 1530.90379;
    T_PLY = 1.3970E-04;
    E_11 = 1.3607E+11;
    E_22 = 7.4808E+09;
    G_12 = 2.8958E+09;
    V_12 = 0.36;
    ULTIMATE_11_TENS = 2.2732E+09;
    ULTIMATE_11_COMP = 1.0378E+09;
    ULTIMATE_22_TENS = 3.7301E+07;
    ULTIMATE_22_COMP = 1.8181E+08;
    ULTIMATE_12 = 5.4882E+07;
    
elseif type_flag == 4
    % 4 - HexPly 8552 IM7 160gsm 35%RW
    RHO = 1575.7808;
    T_PLY = 1.5669E-04;
    E_11 = 1.3410E+11;
    E_22 = 8.4806E+09;
    G_12 = 4.4816E+09;
    V_12 = 0.356;
    ULTIMATE_11_TENS = 1.8805E+09;
    ULTIMATE_11_COMP = 1.5406E+09;
    ULTIMATE_22_TENS = 5.1021E+07;
    ULTIMATE_22_COMP = 2.6745E+08;
    ULTIMATE_12 = 8.8598E+07;
    
elseif type_flag == 5
    % 5 - NCT304-1 HR40 G80 40 +/-2%RW
    RHO = 1508.2873;
    T_PLY = 9.5250E-05;
    E_11 = 1.9010E+11;
    E_22 = 9.0101E+09;
    G_12 = 3.3426E+09;
    V_12 = 0.27;
    ULTIMATE_11_TENS = 2.1614E+09;
    ULTIMATE_11_COMP = 1.0003E+09;
    ULTIMATE_22_TENS = 3.4903E+07;
    ULTIMATE_22_COMP = 2.0806E+08;
    ULTIMATE_12 = 1.1974E+08;
    
elseif type_flag == 6
    % 6 - MTM28-1B/M46J-140-37%RW
    RHO = 1548.7788;
    T_PLY = 1.4700E-04;
    E_11 = 2.1713E+11;
    E_22 = 6.7341E+09;
    G_12 = 3.3426E+09;
    V_12 = 0.27;
    ULTIMATE_11_TENS = 1.9436E+09;
    ULTIMATE_11_COMP = 6.0103E+08;
    ULTIMATE_22_TENS = 1.6470E+07;
    ULTIMATE_22_COMP = 1.4617E+08;
    ULTIMATE_12 = 1.0255E+08;
    
end