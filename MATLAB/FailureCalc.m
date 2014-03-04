%------- FailureCalc -------%
%                           %
%     Cameron Robertson     %
%       March 11, 2011      %
%                           %
%---------------------------%

% FailureCalc computes material failure, euler buckling failure and
% torsional buckling failure. All failure modes are computed at the nodes.
% For material failure it looks at a sample
% laminate on the top, bottom, back and front of the spar. For each side of
% the spar it computes the material failure for the cap, the plus angle
% plys and the minus angle plys. For each ply it computes the fracture of
% failure in the fibre direction, matrix direction and shear. Thus there
% are 4x3x3=36 (side x lamina x direction) possible failure modes in the
% tube.
% ex. fail.top.cap = 3x(Ns+1) vector
%     fail.top.plus = 3x(Ns+1) vector
%     fail.top.minus = 3x(Ns+1) vector
%
% Stresses and strain are given as a 3x(Ns+1) vector
% ex. sigma_11, sigma_22, sigma_12
% Positive sign denotes tensile stresses, negative denotes compressive

function fail = FailureCalc(yN, Finternal, strain, d, theta, nTube, nCap, yWire, zWire, EIxJ, EIzJ, lBiscuit,...
                            dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, RQuad, hQuad, TQuad, EIQuad, GJQuad, tWire, TWire, TEtension,flags)
                        
                        

Ns = length(yN)-1; %number of elements
dy = zeros(Ns,1);
for s=1:Ns
    dy(s) = yN(s+1) - yN(s); %length of each element
end

% Material failure 
%-----------------

fail.top = MaterialFailure(Ns,strain.top,theta,nCap,flags);
fail.bottom = MaterialFailure(Ns,strain.bottom,theta,nCap,flags);
fail.back = MaterialFailure(Ns,strain.back,theta,0,flags);
fail.front = MaterialFailure(Ns,strain.front,theta,0,flags);


% Euler Buckling failure in main spar from wire force
%----------------------------------------------------

k = 0.7; % pinned-pinned = 1, fixed-pinned = 0.7 with correction factor
%kk = 1.42; % correction factor
kk = 1; % correction factor was in error, never saw buckling failure
thetaWire = atan2(zWire,yWire);
L = yWire; %wire attachment provides pinned end
F = TWire*cos(thetaWire)+TEtension;
for s = 1:Ns
    if yN(s) <= yWire
        critical_load_x = pi^2*EIxJ/(k*L)^2;
        critical_load_z = pi^2*EIzJ/(k*L)^2;
        fail.buckling.x(s) = kk*F/critical_load_x;
        fail.buckling.z(s) = kk*F/critical_load_z;
    else
        fail.buckling.x(s) = 0;
        fail.buckling.z(s) = 0;
    end
end
fail.buckling.x(Ns+1) = 0; %no buckling at tip
fail.buckling.z(Ns+1) = 0; %no buckling at tip

% k = 0.7; % pinned-pinned = 1, fixed-pinned = 0.7
% thetaWire = atan2(zWire,yWire);
% for w = 1:length(yWire)
%     L = yWire(w); %wire attachment provides pinned end
%     F = sum(TWire(w:end).*cos(thetaWire(w:end)))+TEtension;
%     for s = 1:Ns
%         if yN(s) <= yWire(w)
%             critical_load_x = pi^2*EIx(s)/(k*L)^2;
%             critical_load_z = pi^2*EIz(s)/(k*L)^2;
%             fail.buckling.x(w,s) = F/critical_load_x;
%             fail.buckling.z(w,s) = F/critical_load_z;
%         else
%             fail.buckling.x(w,s) = 0;
%             fail.buckling.z(w,s) = 0;
%         end
%     end
%     fail.buckling.x(w,Ns+1) = 0; %no buckling at tip
%     fail.buckling.z(w,Ns+1) = 0; %no buckling at tip
% end



% k = 0.7; %one fixed end, one pinned end
% L = yWire; %wire attachment provides pinned end
% for s = 1:Ns
%     critical_load_x = pi^2*EIx(s)/(k*L^2);
%     critical_load_z = pi^2*EIz(s)/(k*L^2);
%     fail.buckling.x(s) = -Finternal(2,s)/critical_load_x;
%     fail.buckling.z(s) = -Finternal(2,s)/critical_load_z;
% end
% fail.buckling.x(Ns+1) = 0; %no buckling at tip
% fail.buckling.z(Ns+1) = 0; %no buckling at tip




% Torsional Buckling failure
%---------------------------
fail.buckling.torsion = torsionalBucklingFailure(Ns,Finternal,d,theta,nTube,nCap,lBiscuit,flags);


% Quad Buckling failure
%-----------------------
if EIQuad ~= 0
    k = 1; %pinned-pinned
    L = sqrt(RQuad^2 + hQuad^2);
    alpha = atan2(hQuad,RQuad);
    P = TQuad/sin(alpha);
    critical_load = pi^2*EIQuad/(k*L)^2;
    fail.quad.buckling = P/critical_load;
else
    fail.quad.buckling = 0;
end


% Quad bending moment failure (does not include torsion since they don't occur at the same time)
%------------------------------
RotorMoment = 1400;
if EIQuad ~= 0
    TbottomWire = TQuad/tan(alpha);
    BM = TbottomWire*zWire + RotorMoment;
    strainQuad = -[BM*(dQuad/2)/EIQuad 0 0]';  %strain on compression side
    fail.quad.bend = MaterialFailure(1,strainQuad,thetaQuad,0,flags);
    fail.quad.bend = abs(fail.quad.bend.plus(1,1)); %use only compressive failure in fibre direction
else
    fail.quad.bend = 0;
end


% Quad torsional material failure
%--------------------------------
if GJQuad ~= 0
    strainQuad = [0 0 dQuad/2*RotorMoment/GJQuad]';
    fail.quad.torsion = MaterialFailure(1,strainQuad,thetaQuad,0,flags);
    fail.quad.torsion = abs(fail.quad.torsion.plus(1,1));
else
    fail.quad.torsion = 0;
end

% Quad torsional buckling failure
%--------------------------------
FRotorMoment(5,1) = RotorMoment;
fail.quad.torbuck = torsionalBucklingFailure(1,FRotorMoment,dQuad,thetaQuad,nTubeQuad,0,lBiscuitQuad,flags);
fail.quad.torbuck = fail.quad.torbuck(1);

% Wire tensile failure
%----------------------
for i = 1:length(yWire)
    stress_wire = TWire(i)/(pi*(tWire(i)/2)^2);
    [RHOWire, EWire, ULTIMATEWire] = WireProperties(flags.WireType);
    fail.wire(i) = stress_wire/ULTIMATEWire;
end

end

function failure = MaterialFailure(Ns,strain,theta,nCap,flags)

failure.cap = zeros(3,Ns+1);
failure.plus = zeros(3,Ns+1);
failure.minus = zeros(3,Ns+1);

% MATERIAL PROPERTIES
[RHO_TUBE, T_PLY_TUBE, E_11_TUBE, E_22_TUBE, G_12_TUBE, V_12_TUBE,...
 ULTIMATE_11_TENS_TUBE, ULTIMATE_11_COMP_TUBE,ULTIMATE_22_TENS_TUBE, ULTIMATE_22_COMP_TUBE, ULTIMATE_12_TUBE]...
 = PrepregProperties(flags.CFRPType);

% Cap Prepreg Properties (MTM28-M46J 140 37 %RW 12")
[RHO_CAP, T_PLY_CAP, E_11_CAP, E_22_CAP, G_12_CAP, V_12_CAP,...
 ULTIMATE_11_TENS_CAP, ULTIMATE_11_COMP_CAP,ULTIMATE_22_TENS_CAP, ULTIMATE_22_COMP_CAP, ULTIMATE_12_CAP]...
 = PrepregProperties(flags.CFRPType);

% Populate Q matrix for tube
Q_TUBE = zeros(3);   % Preallocate Q-matrix
Q_TUBE(1,1) = E_11_TUBE;
Q_TUBE(2,2) = E_22_TUBE;
Q_TUBE(1,2) = E_22_TUBE*V_12_TUBE;
Q_TUBE(2,1) = Q_TUBE(1,2);
Q_TUBE(3,3) = G_12_TUBE;

% Populate Q matrix for caps
Q_CAP = zeros(3);   % Preallocate Q-matrix
Q_CAP(1,1) = E_11_CAP;
Q_CAP(2,2) = E_22_CAP;
Q_CAP(1,2) = E_22_CAP*V_12_CAP;
Q_CAP(2,1) = Q_CAP(1,2);
Q_CAP(3,3) = G_12_CAP;

for s = 1:Ns
    
    % Compute stresses in structural axes for each lamina angle
    % Q is the matrix of elastic constants in the material axis
    % Q_bar is the matrix of elastic constants in the structural axes for a
    % lamina at a ply angle theta.
    
    % Failure is computed at each node, but using the ply angle of the
    % element (this shouldn't cause a large discrepency). Failure at tip is
    % zero, since stresses/strains at tip are zero.

    % Transform the elastic constants for the plus angle ply
    x = theta(s); % Composite angle (in radians)
    T_PLUS = [cos(x)^2 sin(x)^2  2*sin(x)*cos(x);
        sin(x)^2  cos(x)^2  -2*sin(x)*cos(x);
        -sin(x)*cos(x)  sin(x)*cos(x)   (cos(x)^2)-(sin(x)^2)];
    Qbar_TUBE_PLUS = (T_PLUS\Q_TUBE)/T_PLUS'; % transform elastic constants
    
    % Transform the elastic constants for the minus angle ply
    x = -theta(s); % Composite angle (in radians)
    T_MINUS = [cos(x)^2 sin(x)^2  2*sin(x)*cos(x);
        sin(x)^2  cos(x)^2  -2*sin(x)*cos(x);
        -sin(x)*cos(x)  sin(x)*cos(x)   (cos(x)^2)-(sin(x)^2)];
    Qbar_TUBE_MINUS = (T_MINUS\Q_TUBE)/T_MINUS'; % transform elastic constants
    
    % Compute stresses in structural coordinates
    stress.cap(:,s) = Q_CAP*strain(:,s); %using Q for the cap
    stress.plus(:,s) = Qbar_TUBE_PLUS*strain(:,s); %using Q_bar for the + angle plys
    stress.minus(:,s) = Qbar_TUBE_MINUS*strain(:,s); %using Q_bar for the - angle plys

    % Transform stresses to material axes for each lamina angle
    stress.plus(:,s) = T_PLUS*stress.plus(:,s);
    stress.minus(:,s) = T_MINUS*stress.minus(:,s);


    %Determine fraction of failure for each lamina angle

    %ULTIMATE_11_TENS and ULTIMATE_!!_COMP are both positive values
    %indicating the maximum tensile and compressive stress before failure.
    %failure.cap(1,s) will be positive for tensile failures and negative for
    %compressive failures.

    % Cap failure
    if nCap == 0
        failure.cap(1,s) = 0;
        failure.cap(2,s) = 0;
        failure.cap(3,s) = 0;
    else
        if stress.cap(1,s) > 0 %tensile stress in fibre
            failure.cap(1,s) = stress.cap(1,s)/ULTIMATE_11_TENS_CAP;
        else
            failure.cap(1,s) = stress.cap(1,s)/ULTIMATE_11_COMP_CAP;
        end
        if stress.cap(2,s) > 0 %tensile stress in matrix
            failure.cap(2,s) = stress.cap(2,s)/ULTIMATE_22_TENS_CAP;
        else
            failure.cap(2,s) = stress.cap(2,s)/ULTIMATE_22_COMP_CAP;
        end
        failure.cap(3,s) = stress.cap(3,s)/ULTIMATE_12_CAP;
    end

    % Plus angle ply failure
    if stress.plus(1,s) > 0 %tensile stress in fibre
        failure.plus(1,s) = stress.plus(1,s)/ULTIMATE_11_TENS_TUBE;
    else
        failure.plus(1,s) = stress.plus(1,s)/ULTIMATE_11_COMP_TUBE;
    end
    if stress.plus(2,s) > 0 %tensile stress in matrix
        failure.plus(2,s) = stress.plus(2,s)/ULTIMATE_22_TENS_TUBE;
    else
        failure.plus(2,s) = stress.plus(2,s)/ULTIMATE_22_COMP_TUBE;
    end
    failure.plus(3,s) = stress.plus(3,s)/ULTIMATE_12_TUBE;

    % Minus angle ply failure
    if stress.minus(1,s) > 0 %tensile stress in fibre
        failure.minus(1,s) = stress.minus(1,s)/ULTIMATE_11_TENS_TUBE;
    else
        failure.minus(1,s) = stress.minus(1,s)/ULTIMATE_11_COMP_TUBE;
    end
    if stress.minus(2,s) > 0 %tensile stress in matrix
        failure.minus(2,s) = stress.minus(2,s)/ULTIMATE_22_TENS_TUBE;
    else
        failure.minus(2,s) = stress.minus(2,s)/ULTIMATE_22_COMP_TUBE;
    end
    failure.minus(3,s) = stress.minus(3,s)/ULTIMATE_12_TUBE;

end
end

function failure = torsionalBucklingFailure(Ns,Finternal,d,theta,nTube,nCap,lBiscuit,flags)

    % MATERIAL PROPERTIES
    [RHO_TUBE, T_PLY_TUBE, E_11_TUBE, E_22_TUBE, G_12_TUBE, V_12_TUBE,...
     ULTIMATE_11_TENS_TUBE, ULTIMATE_11_COMP_TUBE,ULTIMATE_22_TENS_TUBE, ULTIMATE_22_COMP_TUBE, ULTIMATE_12_TUBE]...
     = PrepregProperties(flags.CFRPType);

    V_21_TUBE = V_12_TUBE*(E_22_TUBE/E_11_TUBE);

    % Coordinate system: x is axial direction, theta is circumferential direction
    mu_prime_x = V_12_TUBE;
    mu_prime_theta = V_21_TUBE;

    Q = zeros(3);   % Preallocate Q-matrix

    % Matrix of elastic constants
    Q(1,1) = E_11_TUBE/(1 - V_12_TUBE*V_21_TUBE); % AER1401, Composite Lamina, slide 8
    Q(2,2) = E_22_TUBE/(1 - V_12_TUBE*V_21_TUBE); % AER1401, Composite Lamina, slide 8
    Q(1,2) = E_22_TUBE*V_12_TUBE/(1 - V_12_TUBE*V_21_TUBE); % AER1401, Composite Lamina, slide 8
    Q(2,1) = Q(1,2); % AER1401, Composite Lamina, slide 8
    Q(3,3) = G_12_TUBE; % AER1401, Composite Lamina, slide 8

    failure = zeros(1,Ns+1);
    for s = 1:Ns

        if nCap(s) ~= 0
            AF_torsional_buckling = 1.25; % See "Validation - Torsional Buckling.xlsx"
        else
            AF_torsional_buckling = 1;
        end   

        % Determine elastic properties of rotated tube laminate
        x = theta(s); % Composite angle (in radians)

        % Transformation matrix
        T = [cos(x)^2 sin(x)^2  2*sin(x)*cos(x);
            sin(x)^2  cos(x)^2  -2*sin(x)*cos(x);
            -sin(x)*cos(x)  sin(x)*cos(x)   (cos(x)^2)-(sin(x)^2)];

        % Transform the elastic constants using the transformation matrix to obtain the
        % elastic constants at the composite angle.
        Qbar = (T\Q)/T';

        % Breakout tube elastic constants at the transformed angle
        E_x = Qbar(1,1);
        E_theta = Qbar(2,2);

        % Calculate tube geometric properties
        t_tube = nTube(s)*T_PLY_TUBE;  % Shell thickness, CR-912 p.xxxvi
        R = (d(s) + t_tube)/2; % Radius from axis of rotation to centroidal surface of cylinder wall, CR-912 p.xxxiv
        L = lBiscuit(s); % Unsupported length of cylinder, CR-912 p.xxx

        % Calculate tube elastic properties

        D_x = E_x*(1/12)*(t_tube^3); % CR-912 p.576
        D_theta = E_theta*(1/12)*(t_tube^3); % CR-912 p.576
        B_x = E_x*t_tube; % CR-912 p.576
        B_theta = E_theta*t_tube;% CR-912 p.576

        % Calculate Gamma
        rho = ((D_x*D_theta)/(B_x*B_theta))^(1/4); % CR-912, p.602
        Gamma = (0.00000036125)*((R/(rho*1000))^6) + ...
            (-0.000019724)*((R/(rho*1000))^5) + ...
            (0.0004283)*((R/(rho*1000))^4) + ...
            (-0.0048315)*((R/(rho*1000))^3) + ...
            (0.031801)*((R/(rho*1000))^2) + ...
            (-0.12975)*(R/(rho*1000)) + ...
            0.88309; % CR-912 P.602 for plot, fit from Excel

        % Calculate factors required in critical torque equation
        Z = ((B_theta*(1-mu_prime_x*mu_prime_theta)*(L^4))/(12*D_x*(R^2)))^(1/2);   % CR-912 p. 600
        Z_s = ((D_theta/D_x)^(5/6))*((B_x/B_theta)^(1/2))*Z;    % CR-912 p.600
        K_s = 0.89*(Z_s^(3/4)); % CR-912 p.599
        N_x_theta = (Gamma*K_s*(pi^2)*D_x)/(L^2); % CR-912 p.599

        % Calculate critical torque
        critical_torque = AF_torsional_buckling*N_x_theta*2*pi*(R^2); % CR-912 p.599

        failure(s) = abs(Finternal(5,s)/critical_torque);

    end
    failure(Ns+1) = 0; %no torsion at tip    

end