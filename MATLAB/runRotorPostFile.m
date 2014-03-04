% HeliOpt run file

% Flags
flags.Opt = 0;       % 0 - single run, 1 - optimization
flags.ConFail = 0;
flags.ConWireCont = 0;
flags.ConJigCont = 0;
flags.ConDef = 0;
flags.MultiPoint = 0;% 0 - single point optimization, 1 - multipoint (h=0.5, h=3, wind case, gravity load)
flags.Quad = 1;      % 0 - prop drive, 1 - quad rotor
flags.FreeWake = 0;  % 0 - momentum theory, 1 - free vortex ring wake
flags.Cover = 0;     % 0 - no cover, 1 - cover
flags.Load = 2;      % 0 - normal run, 1 - gravity forces only, 2 - prescribed load from pLoad
flags.Cdfit = 1;     % 0 - analytic model, 1 - curve fit on BE airfoils
flags.GWing = 1;     % 0 - Daedalus style wing, 1 - Gossamer style wing
flags.AeroStr = 1;   % 0 - Assume flat wing, 1 - take deformation into account
flags.Movie = 0;     % 0 - don't save animation, 1 - save animation
flags.wingWarp = 0;  % 0 - no twist constraint, >0 - twist constraint at flags.wingWarp
flags.CFRPType = 1;  % 1 - high mod
flags.WireType = 1;  % 1 - Steel, 2 - Vectran

multiCore = 4; %set number of cores

% General Propteries
Ns = 10;
rho = 1.18;
visc = 1.78*10^(-5);
vw = 0;
vc = 0;

% Cl_  = [1.4 1.3 1.1]; % 0.5 m No canard lift distribution
% Omega = 0.18*2*pi; %rad/s

% Cl_  = [1.4 1.25 0.8]; % 0.5 m Canard at negative aoa (causes twist of 2.25deg)
% Omega = 0.195*2*pi; %rad/s

Cl_  = [1.4 1.35 1.4]; % 0.5 m Canard at positive aoa (causes twist of 2.25deg)
Omega = 0.165*2*pi; %rad/s

Cm_  = [-0.15 -0.12 -0.12 ];
t_   = [0.14 0.14 0.14 ];
xEA_ = [0.27 0.33 0.24]; % percent chord

% Cl_  = [1.5  1.3  1.3];
% Cm_  = [-0.15 -0.12 -0.12];
% t_   = [0.14 0.14 0.14];
% xEA_ = [0.25 0.25 0.25]; % percent chord

TEtension = 50; %N

mPilot = 72; %kg
presLoad = [];


% Quad-rotor properties
if flags.Quad
        
    
    R = 1.2;
    b = 2;
    H = 0.5; %Height of aircraft
    mElseRotor = 5.11;
    mElseCentre = 6.487+3;
    mElseR = 0.032;
    collective = 0*pi/180; %angle in rad
    
    
    
    % Optimization Constraints
    vrCon.MaxDelta = -0.1;
    vrCon.MinDelta = 0.1;
    vrCon.Wind = 0/3.6;
    vrCon.OmegaRatio = 2;
    vrCon.ClMax = [1.4 1.25 0.85]; % min control
    vrCon.ClMax = [1.4 1.35 1.55]; % max control
    %vrCon.ClMax = [1.4 1.3  1.2]; % max control
    vrCon.AltRatio = 35/60; %proportion of time near ground
    vrCon.Alt = [0.5 3.5];
    vrCon.FOSmat = 0.5; %1.3;
    vrCon.FOSbuck = 0.5; %1.3;
    vrCon.FOSquadbuck = 5;
    vrCon.FOStorbuck = 0.5; %1.5;
    vrCon.FOSwire = 0.5; %2;
    

        ycmax = [1.4656 3.2944];
        


        %c_ = [0 1 1.4 0 0.35]; %linear   [root 1/r start curve end]
        c_ = [0 0.8 1.4 0.4 0.36]; %ideal curve
        
        
        xtU_ = [0.15  7  0.15];
        xtL_ = [0.30  7  0.30];

        d_   =   [22.2    24.8   27.3  ]/1000; %mm -> m
        theta_ = [20   20   20  ]*pi/180; %deg -> rad
        nTube_ = [23    14    4   ];
        nCap_ =  [0    0    0   ];
        lBiscuit_ = [12 12 6  ]*2.54/100; %in -> m


        dQuad = 4*2.54/100; %in -> m
        thetaQuad = 35*pi/180;
        nTubeQuad = 4;
        lBiscuitQuad = 12*2.54/100;
        hQuad = 3;
        etaP = 0;
        yWire = 0.6; %actual spars
        zWire = 1;
        tWire = 0.0028; %Vectran
        tWire = 0.0016; %Steel
        TWire = 1100;
        TWireMulti = [900 2100 110];
        %anhedral = 3.3*pi/180; %deg -> rad
        anhedral = 0;

        h = H + zWire; % Height of rotor


        % Prescribed load
        presLoad.y = 0.6; % Point load location
        %presLoad.y = 7; % Point load location
        presLoad.pointZ = 335; % N
        %presLoad.pointZ = 380; % N
        presLoad.pointM = 0; % Nm
        presLoad.distributedX = 0; %N/m
        presLoad.distributedZ = 0; %N/m
        presLoad.distributedM = 0; %Nm/m
        
        %Spar test 1
%         vrCon.Wind = 10/3.6;
%         vrCon.OmegaRatio = 1;
%         Omega = 0.186*2*pi; %rad/s
%         vrCon.ClMax = [1.4 1.23 0.99];
%         TWireMulti = [1200 2100 200];
%         yWire = 5.7328;
%         anhedral = 0;
%         
        %Spar test 2
%         vrCon.Wind = 10/3.6;
%         vrCon.OmegaRatio = 0.45;
%         Omega = 0.18*2*pi; %rad/s
%         vrCon.ClMax = [1.4 1.3 1.1];
%         TWireMulti = [1200 980 200];
%         yWire = 5.8852; %actual spars
%         anhedral = 0;
        
        
        % Optimization ranges Quad
        if flags.Opt

            i = 1;

            vrOpt(i).name = 'Omega';
            vrOpt(i).u = 0.19*2*pi; %25
            vrOpt(i).l = 0.15*2*pi; %15
            i = i+1;
            
%             vrOpt(i).name = 'TWire';
%             vrOpt(i).u = 1300;
%             vrOpt(i).l = 800;
%             i = i+1;
            
%             vrOpt(i).name = 'TWireMulti(1:2)';
%             vrOpt(i).u = [1300 2500];
%             vrOpt(i).l = [800 1500];
%             i = i+1;
            
%             vrOpt(i).name = 'ycmax(2)';
%             vrOpt(i).u = 6;
%             vrOpt(i).l = 4;
%             i = i+1;
    

%             vrOpt(i).name = 'c_(2:3)';
%             vrOpt(i).u = [1 0.45];
%             vrOpt(i).l = [0.7 0.25];
%             i = i+1;

%             vrOpt(i).name = 'd_(2)';
%             vrOpt(i).u = [2.5]*2.54/100;
%             vrOpt(i).l = [2]*2.54/100;
%             i = i+1;


%             vrOpt(i).name = 'dQuad';
%             vrOpt(i).u = 4*2.54/100;
%             vrOpt(i).l = 3.5*2.54/100;
%             i = i+1;


%             vrOpt(i).name = 'yWire';
%             vrOpt(i).u = 7.5;
%             vrOpt(i).l = 5.5;
%             i = i+1;


            
%             vrOpt(i).name = 'xtU_(2)';
%             vrOpt(i).u = 8;
%             vrOpt(i).l = 5;
%             i = i+1;
            
           

%             vrOpt(i).name = 'anhedral';
%             vrOpt(i).u = 3*pi/180;
%             vrOpt(i).l = 0*pi/180;
%             i = i+1;
% 
%             vrOpt(i).name = 'tWire';
%             vrOpt(i).u = 0.003;
%             vrOpt(i).l = 0.001;
%             i = i+1;
            

%             vrOpt(i).name = 'xtU_(2)';
%             vrOpt(i).u = 8;
%             vrOpt(i).l = 3;
%             i = i+1;

    %         vrOpt(i).name = 'nCap_(2)';
    %         vrOpt(i).u = 5;
    %         vrOpt(i).l = 0;
    %         i = i+1;

    %         vrOpt(i).name = 'nTube_(2)';
    %         vrOpt(i).u = 7;
    %         vrOpt(i).l = 4;
    %         i = i+1;


            
%             vrOpt(i).name = 'collective';
%             vrOpt(i).u = [5]*pi/180;
%             vrOpt(i).l = [0];
%             i = i+1;

    %         vrOpt(i).name = 'zWire';
    %         vrOpt(i).u = 1.2;
    %         vrOpt(i).l = 0.8;
    %         i = i+1;


    %         vrOpt(i).name = 'R';
    %         vrOpt(i).u = 25;
    %         vrOpt(i).l = 14;
    %         i = i+1;

    %         vrOpt(i).name = 'Cl_(2:3)';
    %         vrOpt(i).u = [1.5 1.3];
    %         vrOpt(i).l = [1.2 1];
    %         i = i+1;


    
    
            % Linked parameters
            i = 1;
%             vrLink(i).equation = 'Cl_(3) = Cl_(2);';
%             i = (i+1);
               vrLink(i).equation = 'blah = 0;';
               i = (i+1);
                

            % Multipoint Optimization

            i = 1;
            vrOptH(i).name = 'Omega';
            vrOptH(i).u = 0.19*2*pi; %25
            vrOptH(i).l = 0.15*2*pi; %18
            i = i+1;

            vrOptH(i).name = 'Cl_([1 2])';
            vrOptH(i).u = [1.4 1.3];
            vrOptH(i).l = [0.8 0.8];
            i = i+1;
            
%             vrOptH(i).name = 'Cl_';
%             vrOptH(i).u = [1.5 1.3 1.3];
%             vrOptH(i).l = [0.8 0.8 0.8];
%             i = i+1;


        else
            vrOpt = [];
            vrOptH = [];
            vrLink = [];
        end
        
        
        
        
        
        
        
        if R == 10
    elseif R == 8
    
        ycmax = [0.1 5];
        Omega = 0.2*2*pi; %rad/s
        c_   = [1.4 0 0.3]; %c_2 0:1/r 1:curve.    c_4 curvature

        d_   =   [2.5    2    1.5   ]*2.54/100; %in -> m
        theta_ = [30   30   30  ]*pi/180; %deg -> rad
        nTube_ = [4    4    2   ];
        nCap_ =  [0    0    0   ];
        lBiscuit_ = [12 12  12  ]*2.54/100; %in -> m


        dQuad = 3.5*2.54/100; %in -> m
        thetaQuad = 20*pi/180;
        nTubeQuad = 4;
        lBiscuitQuad = 12*2.54/100;
        hQuad = 3.2;
        etaP = 0;
        yWire = 6.5;
        zWire = 1;
        tWire = 0.002;
        TWire = 1000; %2135;
        TWireMulti = [1000 2500 50];
        anhedral = 2*pi/180; %deg -> rad

        h = H + zWire; % Height of rotor

        

        % Optimization ranges Quad
        if flags.Opt

            i = 1;

    %         vrOpt(i).name = 'ycmax';
    %         vrOpt(i).u = 3;
    %         vrOpt(i).l = 1;
    %         i = i+1;
    % 
    %         vrOpt(i).name = 'c_(2:5)';
    %         vrOpt(i).u = [1 1.5 1 0.5];
    %         vrOpt(i).l = [0 1.3 -1 0.2];
    %         i = i+1;

            vrOpt(i).name = 'c_(5)';
            vrOpt(i).u = [0.5];
            vrOpt(i).l = [0.3];
            i = i+1;

            vrOpt(i).name = 'd_(2)';
            vrOpt(i).u = [3]*2.54/100;
            vrOpt(i).l = [1.5]*2.54/100;
            i = i+1;

            vrOpt(i).name = 'Omega';
            vrOpt(i).u = 0.25*2*pi;
            vrOpt(i).l = 0.17*2*pi;
            i = i+1;

            vrOpt(i).name = 'dQuad';
            vrOpt(i).u = 4*2.54/100;
            vrOpt(i).l = 2.5*2.54/100;
            i = i+1;


            vrOpt(i).name = 'yWire';
            vrOpt(i).u = 6.5;
            vrOpt(i).l = 4;
            i = i+1;

            vrOpt(i).name = 'TWire';
            vrOpt(i).u = 1200;
            vrOpt(i).l = 800;
            i = i+1;


    %         vrOpt(i).name = 'nCap_(2)';
    %         vrOpt(i).u = 5;
    %         vrOpt(i).l = 0;
    %         i = i+1;

    %         vrOpt(i).name = 'nTube_(2)';
    %         vrOpt(i).u = 7;
    %         vrOpt(i).l = 4;
    %         i = i+1;

            vrOpt(i).name = 'TWireMulti(1:2)';
            vrOpt(i).u = [1200 3500];
            vrOpt(i).l = [800 2000];
            i = i+1;

            vrOpt(i).name = 'anhedral';
            vrOpt(i).u = 4*pi/180;
            vrOpt(i).l = 0*pi/180;
            i = i+1;

            vrOpt(i).name = 'tWire';
            vrOpt(i).u = 0.003;
            vrOpt(i).l = 0.001;
            i = i+1;

    %         vrOpt(i).name = 'zWire';
    %         vrOpt(i).u = 1.2;
    %         vrOpt(i).l = 0.8;
    %         i = i+1;


    %         vrOpt(i).name = 'R';
    %         vrOpt(i).u = 25;
    %         vrOpt(i).l = 14;
    %         i = i+1;

    %         vrOpt(i).name = 'Cl_(2:3)';
    %         vrOpt(i).u = [1.5 1.3];
    %         vrOpt(i).l = [1.2 1];
    %         i = i+1;


            % Linked parameters
            i = 1;
            vrLink(i).equation = 'd_(1) = d_(2) + 0.5*2.54/100;';
            i = (i+1);

            % Multipoint Optimization

            i = 1;
            vrOptH(i).name = 'Omega';
            vrOptH(i).u = 0.26*2*pi;
            vrOptH(i).l = 0.18*2*pi;
            i = i+1;

            vrOptH(i).name = 'Cl_(2:3)';
            vrOptH(i).u = [1.5 1.3];
            vrOptH(i).l = [0.9 0.9];
            i = i+1;


        else
            vrOpt = [];
            vrOptH = [];
        end
        
    
    elseif R == 12
        ycmax = 3;
        Omega = 0.15*2*pi; %rad/s
        c_   = [0 0 1.4 0 0.3]; %c_2 0:1/r 1:curve.    c_4 curvature

        d_   =   [3.5    3    1.5   ]*2.54/100; %in -> m
        theta_ = [30   30   30  ]*pi/180; %deg -> rad
        nTube_ = [4    4    2   ];
        nCap_ =  [0    0    0   ];
        lBiscuit_ = [12 12  12  ]*2.54/100; %in -> m


        dQuad = 5*2.54/100; %in -> m
        thetaQuad = 20*pi/180;
        nTubeQuad = 4;
        lBiscuitQuad = 12*2.54/100;
        hQuad = 3.2;
        etaP = 0;
        yWire = 7;
        zWire = 1;
        tWire = 0.002;
        TWire = 1300; %2135;
        TWireMulti = [1300 3000 50];
        anhedral = 2*pi/180; %deg -> rad

        h = H + zWire; % Height of rotor



        % Optimization ranges Quad
        if flags.Opt

            i = 1;

    %         vrOpt(i).name = 'ycmax';
    %         vrOpt(i).u = 3;
    %         vrOpt(i).l = 1;
    %         i = i+1;
    % 
    %         vrOpt(i).name = 'c_(2:5)';
    %         vrOpt(i).u = [1 1.5 1 0.5];
    %         vrOpt(i).l = [0 1.3 -1 0.2];
    %         i = i+1;

            vrOpt(i).name = 'c_(5)';
            vrOpt(i).u = [0.5];
            vrOpt(i).l = [0.3];
            i = i+1;

            vrOpt(i).name = 'd_(2)';
            vrOpt(i).u = [4]*2.54/100;
            vrOpt(i).l = [3]*2.54/100;
            i = i+1;

            vrOpt(i).name = 'Omega';
            vrOpt(i).u = 0.15*2*pi;
            vrOpt(i).l = 0.08*2*pi;
            i = i+1;

            vrOpt(i).name = 'dQuad';
            vrOpt(i).u = 7*2.54/100;
            vrOpt(i).l = 4*2.54/100;
            i = i+1;


            vrOpt(i).name = 'yWire';
            vrOpt(i).u = 9;
            vrOpt(i).l = 6;
            i = i+1;

            vrOpt(i).name = 'TWire';
            vrOpt(i).u = 1700;
            vrOpt(i).l = 1000;
            i = i+1;


    %         vrOpt(i).name = 'nCap_(2)';
    %         vrOpt(i).u = 5;
    %         vrOpt(i).l = 0;
    %         i = i+1;

    %         vrOpt(i).name = 'nTube_(2)';
    %         vrOpt(i).u = 7;
    %         vrOpt(i).l = 4;
    %         i = i+1;

            vrOpt(i).name = 'TWireMulti(1:2)';
            vrOpt(i).u = [1700 5500];
            vrOpt(i).l = [1000 3500];
            i = i+1;

            vrOpt(i).name = 'anhedral';
            vrOpt(i).u = 5*pi/180;
            vrOpt(i).l = 0*pi/180;
            i = i+1;

            vrOpt(i).name = 'tWire';
            vrOpt(i).u = 0.004;
            vrOpt(i).l = 0.001;
            i = i+1;

    %         vrOpt(i).name = 'zWire';
    %         vrOpt(i).u = 1.2;
    %         vrOpt(i).l = 0.8;
    %         i = i+1;


    %         vrOpt(i).name = 'R';
    %         vrOpt(i).u = 25;
    %         vrOpt(i).l = 14;
    %         i = i+1;

    %         vrOpt(i).name = 'Cl_(2:3)';
    %         vrOpt(i).u = [1.5 1.3];
    %         vrOpt(i).l = [1.2 1];
    %         i = i+1;


            % Linked parameters
            i = 1;
            vrLink(i).equation = 'd_(1) = d_(2) + 0.5*2.54/100;';
            i = (i+1);

            % Multipoint Optimization

            i = 1;
            vrOptH(i).name = 'Omega';
            vrOptH(i).u = 0.18*2*pi;
            vrOptH(i).l = 0.1*2*pi;
            i = i+1;

            vrOptH(i).name = 'Cl_(2:3)';
            vrOptH(i).u = [1.5 1.3];
            vrOptH(i).l = [0.9 0.9];
            i = i+1;


        else
            vrOpt = [];
            vrOptH = [];
        end
    end

% Prop-driven properties
else
    
    ycmax = 5;
    R = 20;
    b = 3;
    h = 5;
    Omega = 0.1*2*pi; %rad/s
    c_   = [0 -1 1.4 0 0.3];
    
    d_   =   [4.1    4.5    1   ]*2.54/100; %in -> m
    theta_ = [30   30   30  ]*pi/180; %deg -> rad
    nTube_ = [4    4    4   ];
    nCap_ =  [0    0    0   ];
    %nCap_ =  [5    5    5   ];
    lBiscuit_ = [12 12  12  ]*2.54/100; %in -> m
    
    mElse = 6.370;
    dQuad = 0;
    thetaQuad = 0;
    nTubeQuad = 0;
    RQuad = 0;
    hQuad = 0;
    etaP = 0.9;
    yWire = 10.1;
    zWire = 2;
    tWire = 0.001;
    TWire = 2500;
    
    % Optimization ranges Prop-Driven
    if flags.Opt
        maxDelta = 2;
        
        i = 1;

%         vrOpt(i).name = 'ycmax';
%         vrOpt(i).u = 8;
%         vrOpt(i).l = 5;
%         i = i+1;
% 
        vrOpt(i).name = 'c_(3:5)';
        vrOpt(i).u = [1.5 0.5 0.3];
        vrOpt(i).l = [0.8 -0.5 0.1];
        i = i+1;
% 
% %         vrOpt(i).name = 'nCap_';
% %         vrOpt(i).u = [10 10 3];
% %         vrOpt(i).l = [0 0 0];
% %         i = i+1;
% 
        vrOpt(i).name = 'd_';
        vrOpt(i).u = [4.5 5 1.5]*2.54/100;
        vrOpt(i).l = [3.5 4 1]*2.54/100;
        i = i+1;

        vrOpt(i).name = 'Omega';
        vrOpt(i).u = 0.12*2*pi;
        vrOpt(i).l = 0.06*2*pi;
        i = i+1;

       
%         vrOpt(i).name = 'yWire';
%         vrOpt(i).u = 12;
%         vrOpt(i).l = 8;
%         i = i+1;

        vrOpt(i).name = 'TWire';
        vrOpt(i).u = 3500;
        vrOpt(i).l = 1500;
        i = i+1;
        
%         vrOpt(i).name = 'R';
%         vrOpt(i).u = 25;
%         vrOpt(i).l = 15;
%         i = i+1;

    else
        vrOpt = [];
    end
    
end