%--------- StrCalc ---------%
%                           %
%       Todd Reichert       %
%        Feb 28, 2011       %
%                           %
%---------------------------%

% StrCalc performs three tasks:
% 1. Computes the structural properties of the spars (EI, GJ, m) based on
% the spar properties, and returns the total mass of the helicopter
% 2. Computes the deformation of the spar
% 3. Computes the factor of safety for each of the failure modes of the
% spar.

% Need to add
% - proper calculation of where xAC is
% - proper calculation of where xCG is

function [Mtot, mSpar, mChord, mQuad, mCover, mWire, EIx, EIz, EA, GJ, q, EIQuad, GJQuad, Finternal, strain, fail] ...
    = StrCalc(yN, R, b, cE, xEA, xtU, d, theta, nTube, nCap, Jprop, yWire, zWire, tWire, TWire, TEtension, ycmax, lBiscuit, ...
    dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, RQuad, hQuad, mElseRotor, mElseCentre, mElseR, mPilot, Fblade, presLoad, flags)

Ns = length(yN)-1; %number of elements
dy = zeros(Ns,1);
for s=1:Ns
    dy(s) = yN(s+1) - yN(s); %length of each element
end

% Compute structural properties and mass of each element
%--------------------------------------------------------
[EIx, EIz, EA, GJ, mSpar] = SparProperties(yN, d, theta, nTube, nCap, lBiscuit, flags);
[EIxJ, EIzJ, EAJ, GJJ, mSparJ] = SparProperties([0 1], Jprop.d, Jprop.theta, Jprop.nTube, Jprop.nCap, Jprop.lBiscuit, flags);
[mChord, xCGChord] = ChordProperties(yN, cE, d, flags.GWing,xtU);
xCG = (xCGChord.*mChord + xEA.*mSpar)./(mChord+mSpar);
if flags.Quad % Quad rotor configuration
    lQuad = sqrt(RQuad^2 + hQuad^2);
    [EIxQuad, EIzQuad, EAQuad, GJQuad, mQuad] = SparProperties([0 lQuad], dQuad, thetaQuad, nTubeQuad, [0 0], lBiscuitQuad, flags); 
    EIQuad = EIxQuad;
else % Prop driven configuration
    mQuad = 0;
    EIQuad = 0;
end

if flags.Cover
    mCover = (ycmax^2*0.0528 + ycmax*0.605/4)*1.15;
else
    mCover = 0;
end

[RHOWire, EWire, ULTIMATEWire] = WireProperties(flags.WireType);
LWire = sqrt(zWire^2 + yWire^2);
mWire = pi*(tWire/2)^2*RHOWire*LWire;

if flags.Quad
    Mtot = (sum(mSpar)*b + sum(mChord)*b + mWire*b + mQuad + mCover)*4 + mElseRotor + mElseCentre + mElseR*R + mPilot;
else
    Mtot = sum(mSpar)*b + sum(mChord)*b + mWire*b + mCover + mElseRotor + mElseCentre + mElseR*R + mPilot;
end


% FEM computation for structural deformations
%--------------------------------------------

% Initialize global stiffness matrix 
K = zeros((Ns+1)*6,(Ns+1)*6);  %global stiffness
F = zeros((Ns+1)*6,1); %global force vector

% Create global stiffness maxtrix and force vector
k = zeros(12,12,Ns);
for s = 1:Ns
    
    % Local elastic stiffness matrix 
    k(1,1,s) = 12*EIx(s)/(dy(s)*dy(s)*dy(s));
    k(1,6,s) = -6*EIx(s)/(dy(s)*dy(s));              k(6,1,s)=k(1,6,s);
    k(1,7,s) = -12*EIx(s)/(dy(s)*dy(s)*dy(s));        k(7,1,s)=k(1,7,s);
    k(1,12,s) = -6*EIx(s)/(dy(s)*dy(s));             k(12,1,s)=k(1,12,s); 
    k(2,2,s) = EA(s)/dy(s);
    k(2,8,s) = -EA(s)/dy(s);                      k(8,2,s)=k(2,8,s);
    k(3,3,s) = 12*EIz(s)/(dy(s)*dy(s)*dy(s));
    k(3,4,s) = 6*EIz(s)/(dy(s)*dy(s));             k(4,3,s)=k(3,4,s);
    k(3,9,s) = -12*EIz(s)/(dy(s)*dy(s)*dy(s));        k(9,3,s)=k(3,9,s);
    k(3,10,s) = 6*EIz(s)/(dy(s)*dy(s));            k(10,3,s)=k(3,10,s);
    k(4,4,s) = 4*EIz(s)/dy(s);
    k(4,9,s) = -6*EIz(s)/(dy(s)*dy(s));              k(9,4,s)=k(4,9,s);
    k(4,10,s) = 2*EIz(s)/dy(s);                   k(10,4,s)=k(4,10,s);
    k(5,5,s) = GJ(s)/dy(s);
    k(5,11,s) = -GJ(s)/dy(s);                     k(11,5,s)=k(5,11,s);
    k(6,6,s) = 4*EIx(s)/dy(s);
    k(6,7,s) = 6*EIx(s)/(dy(s)*dy(s));             k(7,6,s)=k(6,7,s);
    k(6,12,s) = 2*EIx(s)/dy(s);                   k(12,6,s)=k(6,12,s);
    k(7,7,s) = 12*EIx(s)/(dy(s)*dy(s)*dy(s));
    k(7,12,s) = 6*EIx(s)/(dy(s)*dy(s));            k(12,7,s)=k(7,12,s);
    k(8,8,s) = EA(s)/dy(s);
    k(9,9,s) = 12*EIz(s)/(dy(s)*dy(s)*dy(s));
    k(9,10,s) = -6*EIz(s)/(dy(s)*dy(s));             k(10,9,s)=k(9,10,s);
    k(10,10,s) = 4*EIz(s)/dy(s);
    k(11,11,s) = GJ(s)/dy(s);
    k(12,12,s) = 4*EIx(s)/dy(s);

    % Perform dihedral and sweep rotations here if needed
    
    % Assemble global stiffness matrix
    K(6*s-5:6*s+6,6*s-5:6*s+6) = K(6*s-5:6*s+6,6*s-5:6*s+6) + k(:,:,s);
       
    if flags.Load == 0 %Include aero forces

        % Aerodynamic forces
        xAC(s) = 0.25;
        Faero(1) = Fblade.Fx(s)/2;
        Faero(2) = 0;
        Faero(3) = Fblade.Fz(s)/2;
        Faero(4) = Fblade.Fz(s)*dy(s)/12;
        Faero(5) = Fblade.My(s)/2 + (xEA(s)-xAC(s))*cE(s)*Fblade.Fz(s)/2;
        Faero(6) = -Fblade.Fx(s)*dy(s)/12;
    
    else
        Faero = zeros(6,1);
    end
    
    if (flags.Load == 0) || (flags.Load == 1) %Include gravity and wire forces
        % Gravitational forces
        g = 9.81;
        Fg(1) = 0;
        Fg(2) = 0;
        Fg(3) = -(mSpar(s) + mChord(s))*g/2;
        Fg(4) = -(mSpar(s) + mChord(s))*g*dy(s)/12;
        Fg(5) = (xCG(s)-xEA(s))*cE(s)*(mSpar(s)+mChord(s))*g/2;
        Fg(6) = 0;
        
        % Wire forces (using consistent force vector)
        for w = 1:length(yWire)
            if (yWire(w) >= yN(s)) && (yWire(w) < yN(s+1))
                thetaWire = atan2(zWire,yWire(w));
                a = yWire(w)-yN(s);
                L = dy(s);
                FxWire = -cos(thetaWire)*TWire(w);
                FzWire = -sin(thetaWire)*TWire(w);
                Fwire(1) = 0;
                Fwire(2) = FxWire*(1-a/L);
                Fwire(3) = FzWire*(2*(a/L)^3 - 3*(a/L)^2 + 1);
                Fwire(4) = FzWire*a*((a/L)^2 - 2*(a/L) + 1);
                Fwire(5) = 0;
                Fwire(6) = 0;
                Fwire(7) = 0;
                Fwire(8) = FxWire*(a/L);
                Fwire(9) = FzWire*(-2*(a/L)^3 + 3*(a/L)^2);
                Fwire(10) = FzWire*a*((a/L)^2 - (a/L));
                Fwire(11) = 0;
                Fwire(12) = 0;
            else
                Fwire = zeros(12,1);
            end
        end
    else
        Fg = zeros(6,1);
        Fwire = zeros(12,1);
    end
    
    if flags.Load == 2
        % Prescribed point load (using consistent force vector)
        if (presLoad.y >= yN(s)) && (presLoad.y < yN(s+1))
            a = presLoad.y-yN(s);
            L = dy(s);
            Fpres(1) = 0;
            Fpres(2) = 0;
            Fpres(3) = presLoad.pointZ*(2*(a/L)^3 - 3*(a/L)^2 + 1);
            Fpres(4) = presLoad.pointZ*a*((a/L)^2 - 2*(a/L) + 1);
            Fpres(5) = presLoad.pointM*(1-a/L);
            Fpres(6) = 0;
            Fpres(7) = 0;
            Fpres(8) = 0;
            Fpres(9) = presLoad.pointZ*(-2*(a/L)^3 + 3*(a/L)^2);
            Fpres(10) = presLoad.pointZ*a*((a/L)^2 - (a/L));
            Fpres(11) = presLoad.pointM*(a/L);
            Fpres(12) = 0;
        else
            Fpres = zeros(12,1);
        end
    
        % Prescribed distributed load
        Fpres(1) = Fpres(1) + presLoad.distributedX*dy(s)/2;
        Fpres(2) = Fpres(2) + 0;
        Fpres(3) = Fpres(3) + presLoad.distributedZ*dy(s)/2;
        Fpres(4) = Fpres(4) + presLoad.distributedZ*dy(s)*dy(s)/12;
        Fpres(5) = Fpres(5) + presLoad.distributedM*dy(s)/2;
        Fpres(6) = Fpres(6) - presLoad.distributedX*dy(s)*dy(s)/12;
        Fpres(7) = Fpres(7) + presLoad.distributedX*dy(s)/2;
        Fpres(8) = Fpres(8) + 0;
        Fpres(9) = Fpres(9) + presLoad.distributedZ*dy(s)/2;
        Fpres(10) = Fpres(10) - presLoad.distributedZ*dy(s)*dy(s)/12;
        Fpres(11) = Fpres(11) + presLoad.distributedM*dy(s)/2;
        Fpres(12) = Fpres(12) + presLoad.distributedX*dy(s)*dy(s)/12;
    else
        Fpres = zeros(12,1);
    end
    
    
    % Assemble global force vector
    F((s-1)*6+1) = F((s-1)*6+1)  +  Fpres(1) + Fwire(1) + Fg(1) + Faero(1);  % x force     
    F((s-1)*6+2) = F((s-1)*6+2)  +  Fpres(2) + Fwire(2) + Fg(2) + Faero(2);              % y force
    F((s-1)*6+3) = F((s-1)*6+3)  +  Fpres(3) + Fwire(3) + Fg(3) + Faero(3);  % z force
    F((s-1)*6+4) = F((s-1)*6+4)  +  Fpres(4) + Fwire(4) + Fg(4) + Faero(4); % x moment
    F((s-1)*6+5) = F((s-1)*6+5)  +  Fpres(5) + Fwire(5) + Fg(5) + Faero(5); %y moment
    F((s-1)*6+6) = F((s-1)*6+6)  +  Fpres(6) + Fwire(6) + Fg(6)  + Faero(6);% z moment
    F((s-1)*6+7) = F((s-1)*6+7)  +  Fpres(7) + Fwire(7) +  Fg(1) + Faero(1);  % x force     
    F((s-1)*6+8) = F((s-1)*6+8)  +  Fpres(8) + Fwire(8) +  Fg(2) + Faero(2);              % y force
    F((s-1)*6+9) = F((s-1)*6+9)  +  Fpres(9) + Fwire(9) +  Fg(3) + Faero(3);  % z force
    F((s-1)*6+10) = F((s-1)*6+10)+  Fpres(10) + Fwire(10) - Fg(4) - Faero(4); % x moment
    F((s-1)*6+11) = F((s-1)*6+11)+  Fpres(11) + Fwire(11) + Fg(5) + Faero(5); %y moment
    F((s-1)*6+12) = F((s-1)*6+12)+  Fpres(12) + Fwire(12) - Fg(6) - Faero(6); % z moment
end

% Add constraints to all 6 dof at root

if flags.wingWarp > 0 % Also add wingWarping constraint
    ii = [];
    for ss = 1:(Ns+1)*6
        if (ss > 6) && (ss ~= flags.wingWarp*6+5)
            ii = [ii ss];
        end
    end
    %F(5*6+5) = sum(K(:,(flags.wingWarp)*6+5)*2*pi/180);
    Fc = F(ii);
    Kc = K(ii,ii);
    %Fc = Fc - Kc(:,(flags.wingWarp-1)*6+5)*2*pi/180; % use K not Kc (column was eliminated)
else
    Fc = F(7:end);
    Kc = K(7:end,7:end);
end
    

% Solve constrained system
qc = Kc\Fc;

%Add constrained displacements back in
if flags.wingWarp > 0
    q(ii,1) = qc;
    %q(flags.wingWarp*6+5,1) = 2*pi/180;
    
%     sss = 1;
%     for ss = 1:(Ns+1)*6
%         if (ss > 6) %&& (ss ~= flags.wingWarp*6+5)
%             q(ss,1) = qc(sss,1);
%             sss = sss+1;
%         else
%             q(ss,1) = 0;
%         end
%     end
else
    q = [0; 0; 0; 0; 0; 0; qc];
end
    

% Compute internal forces and strains
%---------------------------------------

Ftemp = zeros(12,Ns);
Finternal = zeros(6,Ns+1);

strain.top = zeros(3,Ns+1);
strain.bottom = zeros(3,Ns+1);
strain.back = zeros(3,Ns+1);
strain.front = zeros(3,Ns+1);

strain.bending_x = zeros(1,Ns+1);
strain.bending_z = zeros(1,Ns+1);
strain.axial_y   = zeros(1,Ns+1);
strain.torsion_y = zeros(1,Ns+1);

for s = 1:Ns
    
    % Determine internal forces acting at the nodes of each element
    Ftemp(:,s) = -(k(:,:,s)*q((s-1)*6+1:(s-1)*6+12) - F((s-1)*6+1:(s-1)*6+12));
    Finternal(:,s) = Ftemp(1:6,s);
    % Finternal(1,s) = x-shear load
    % Finternal(2,s) = y-axial load
    % Finternal(3,s) = z-shear load
    % Finternal(4,s) = x-bending moment
    % Finternal(5,s) = y-torsional load
    % Finternal(6,s) = z-bending moment
    
    % Determine strains at each node
    x_hat = d(s)/2; % distance of max strain from bending centroid
    z_hat = d(s)/2; % distance of max strain from bending centroid
    r_hat = d(s)/2; % distance of max strain from torsional centroid

    % Break out displacement vector for element
    qq = q((s-1)*6+1 : (s-1)*6+12);
    
    strain.bending_x(1,s) = -[(-(6*x_hat)/(dy(s)^2)) ((4*x_hat)/dy(s)) ((6*x_hat)/(dy(s)^2)) ((2*x_hat)/dy(s))]*[qq(1) qq(6) qq(7) qq(12)]'; 
    strain.bending_z(1,s) = -[(-(6*z_hat)/(dy(s)^2)) ((-4*z_hat)/dy(s)) ((6*z_hat)/(dy(s)^2)) ((-2*z_hat)/dy(s))]*[qq(3) qq(4) qq(9) qq(10)]'; 
    strain.axial_y(1,s) = [(-1/dy(s)) (1/dy(s))]*[qq(2) qq(8)]'; 
    strain.torsion_y(1,s) = r_hat*[(-1/dy(s)) (1/dy(s))]*[qq(5) qq(11)]';
    
    strain.top(:,s) = [strain.bending_z(1,s)+strain.axial_y(1,s) 0 strain.torsion_y(1,s)]';
    strain.bottom(:,s) = [-strain.bending_z(1,s)+strain.axial_y(1,s) 0 strain.torsion_y(1,s)]';
    strain.back(:,s) = [strain.bending_x(1,s)+strain.axial_y(1,s) 0 strain.torsion_y(1,s)]';
    strain.front(:,s) = [-strain.bending_x(1,s)+strain.axial_y(1,s) 0 strain.torsion_y(1,s)]';
    
end

% Loads at the tip are zero
Finternal(:,Ns+1) = [0 0 0 0 0 0]';

% Strains at tip are zero
strain.top(:,Ns+1) = [0 0 0]';
strain.bottom(:,Ns+1) = [0 0 0]'; 
strain.back(:,Ns+1) = [0 0 0]'; 
strain.front(:,Ns+1) = [0 0 0]'; 

% Strains at tip are zero
strain.bending_x(1,Ns+1) = 0;
strain.bending_z(1,Ns+1) = 0; 
strain.axial_y(1,Ns+1) = 0; 
strain.torsion_y(1,Ns+1) = 0; 


% Compute factor of safety for each failure mode
%-----------------------------------------------
TQuad = sum(Fblade.Fz)*b - (sum(mSpar + mChord)*b + mElseRotor/4)*9.81;

fail = FailureCalc(yN, Finternal, strain, d, theta, nTube, nCap, yWire, zWire, EIxJ, EIzJ, lBiscuit, ...
                   dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, RQuad, hQuad, TQuad, EIQuad, GJQuad, tWire, TWire, TEtension,flags);
