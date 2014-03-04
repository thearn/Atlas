% Validation run file

% Validation case
flagCase = 4;

% Flags
flagOpt = 0;       % 0 - single run, 1 - optimization
flagQuad = 1;      % 0 - prop drive, 1 - quad rotor
flagFreeWake = 0;  % 0 - momentum theory, 2 - free vortex ring wake
flagLoad = 2;      % 0 - normal run, 1 - gravity forces only, 2 - prescribed load from presLoad
multiCore = 1; %set number of cores


% General Propteries
Ns = 25;
rho = 1.18;
visc = 1.78*10^(-5);
vc = 0;

Cl_  = [1.2  1.2  1.2];
Cm_  = [-0.12 -0.12 -0.12];
t_   = [0.14 0.14 0.14];
xtU_ = [0.15  0.15  0.15];
xtL_ = [1    1    1  ];
xEA_ = [0.25 0.25 0.25]; % percent chord

flagTESpar = 0;
flagGWing = 0;
mPilot = 80; %kg

ycmax = 2.9;
R = 10;
b = 2;
h = 3.5;
Omega = 0.16*2*pi; %rad/s
c_   = [0 -1 1.45 0 0.2];

d_   =   [2.1    2.2    1   ]*2.54/100; %in -> m
theta_ = [30   30   30  ]*pi/180; %deg -> rad
nTube_ = [4    4    4   ];
nCap_ =  [0    0    0   ];
lBiscuit_ = [12 12  12  ]*2.54/100; %in -> m

mElse = 6.670;
dQuad = 3*2.54/100; %in -> m
thetaQuad = 20*pi/180;
nTubeQuad = 4;
hQuad = 3;
etaP = 0;
yWire = 4.1;
zWire = 0.5;
tWire = 0.002;
TWire = 0;

vrOpt = [];
maxDelta = 0;

switch flagCase
    
    % FEA: Bending, point load at tip
    case 1
        
        Ns = 25;
        presLoad.y = 9.9999; % Point load location
        presLoad.pointZ = 50; % N
        presLoad.pointM = 0; % Nm
        presLoad.distributedX = 0; %N/m
        presLoad.distributedZ = 0; %N/m
        presLoad.distributedM = 0; %Nm/m
        
        R = 10;
        d_   =   [3    3    3   ]*2.54/100; %in -> m
        nTube_ = [4    4    4   ];
        nCap_ =  [0    0    0   ];
        
    % FEA: Bending, distributed load
    case 2
        
        Ns = 25;
        presLoad.y = 9.9999; % Point load location
        presLoad.pointZ = 0; % N
        presLoad.pointM = 0; % Nm
        presLoad.distributedX = 0; %N/m
        presLoad.distributedZ = 25; %N/m
        presLoad.distributedM = 0; %Nm/m
        
        R = 10;
        d_   =   [3    3    3   ]*2.54/100; %in -> m
        nTube_ = [4    4    4   ];
        nCap_ =  [0    0    0   ];
        
    % FEA: Torsion, point torque at tip
    case 3
        
        Ns = 25;
        presLoad.y = 9.9999; % Point load location
        presLoad.pointZ = 0; % N
        presLoad.pointM = 25; % Nm
        presLoad.distributedX = 0; %N/m
        presLoad.distributedZ = 0; %N/m
        presLoad.distributedM = 0; %Nm/m
        
        R = 10;
        d_   =   [3    3    3   ]*2.54/100; %in -> m
        nTube_ = [4    4    4   ];
        nCap_ =  [0    0    0   ];
        
    % FEA: Torsion, distributed torque
    case 4
        
        Ns = 100;
        presLoad.y = 9.9999; % Point load location
        presLoad.pointZ = 0; % N
        presLoad.pointM = 0; % Nm
        presLoad.distributedX = 0; %N/m
        presLoad.distributedZ = 0; %N/m
        presLoad.distributedM = 10; %Nm/m
        
        R = 10;
        d_   =   [3    3    3   ]*2.54/100; %in -> m
        nTube_ = [4    4    4   ];
        nCap_ =  [0    0    0   ];
        
    % Properties: Bending & Torsion, Spar 2009-1
    case 5
        
        presLoad.y = 10; % Point load location
        presLoad.pointZ = 0; % N
        presLoad.pointM = 0; % Nm
        presLoad.distributedX = 0; %N/m
        presLoad.distributedZ = 0; %N/m
        presLoad.distributedM = 0; %Nm/m
        
        R = 10;
        d_   =   [3    3    3   ]*2.54/100; %in -> m
        theta_ = [30   30   30  ]*pi/180; %deg -> rad
        nTube_ = [4    4    4   ];
        nCap_ =  [8    8    8   ];
        
    % Properties: Bending & Torsion, Spar 2009-2
    case 6
        
        presLoad.y = 2.09; % Point load location
        presLoad.pointZ = 1078.09938; % N
        presLoad.pointM = 0; % Nm
        presLoad.distributedX = 0; %N/m
        presLoad.distributedZ = 0; %N/m
        presLoad.distributedM = 0; %Nm/m
        
        R = 10;
        d_   =   [3    3    3   ]*2.54/100; %in -> m
        theta_ = [30   30   30  ]*pi/180; %deg -> rad
        nTube_ = [4    4    4   ];
        nCap_ =  [10    10    10   ];
            
    % Properties: Bending & Torsion, Spar 2009-3
    case 7
        
        presLoad.y = 1.688; % Point load location
        presLoad.pointZ = 1193.877; % N
        presLoad.pointM = 0; % Nm
        presLoad.distributedX = 0; %N/m
        presLoad.distributedZ = 0; %N/m
        presLoad.distributedM = 0; %Nm/m
        
        R = 10;
        d_   =   [3    3    3   ]*2.54/100; %in -> m
        theta_ = [30   30   30  ]*pi/180; %deg -> rad
        nTube_ = [4    4    4   ];
        nCap_ =  [10    10    10   ];
        
end