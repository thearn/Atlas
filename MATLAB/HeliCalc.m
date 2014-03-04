%--------- HeliCalc --------%
%                           %
%       Todd Reichert       %
%        Feb 28, 2011       %
%                           %
%---------------------------%

% HeliCalc performs an aerodynamic and structural computation on a single
% configuration given the full set of design parameters. The aerodynamic
% computation returns the thrust, torque, power, and induced velocity. The
% structural computation first computes the mass of the helicopter based on
% the structural description of the spars and chord lengths. It then
% computes the deformation of the spars, the strains, and the resulting
% factor of safety for each of the failure modes.

% HeliCalc Upgrades
% - Figure out why convergence is an issue when constraining lift
% coefficient consistency
% - Figure out dynamic climb model
% - Combine with structural optimization of FEM truss
% - Make universal 6-dof FEM code to use for all programs


function out = HeliCalc(Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
            yWire, zWire, tWire, TWire, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective, ...
            dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags)

        
        
% Discretize properties
%----------------------
% Takes parameterized properties, and returns property for each element
% c - chord
% cE - chord of each element
% cN - chord at each node
% Cl - lift coefficient
% t - airfoil thickness
% d - spar diameter
% theta - CFRP wrap angle
% nTube - number of tube layers
% nCap - number of cap strips
% yN - node locations
% Ns - number of elements
[cE, cN, c100, Cl, Cm, t, xtU, xtL, xEA, d, theta, nTube, nCap, lBiscuit, yN, yE] ...
    = DiscretizeProperties(Ns, ycmax, R, c_, Cl_, Cm_, t_, xtU_, xtL_, xEA_, yWire, d_, theta_, nTube_, nCap_, lBiscuit_);

% Quad-rotor configuration
RQuad = sqrt(2*R^2)+0.05;

% Properties at joint location for buckling analysis
Jprop.d = d_(2);
Jprop.theta = theta_(2);
Jprop.nTube = nTube_(2);
Jprop.nCap = nCap_(2);
Jprop.lBiscuit = lBiscuit_(2);

% Aerodynamic Calculation
%------------------------
% yN - node locations
% T - thrust distribution
% P - power distribution
% Q - torque distribution
% vi - induced downwash distribution
% rho - air density
% visc - air viscosity
% vw - wind velocity
% vc - vertical velocity
% R - rotor radius
% b - number of blades
% h - height of rotor
% Omega - rotor angular velocity
% c - chord distribution
% Cl - lift coefficient distribution
% yWire - location of wire attachment along span
% zWire - depth of wire attachement
% tWire - thickness of wire
% flags.Quad - 0 for prop-driven, 1 for quad-rotor configuration
if flags.Load == 0
    q = zeros(1,(Ns+1)*6);
    if flags.AeroStr && flags.FreeWake
        flagFreeWakeTemp = flags.FreeWake;
        flags.FreeWake = 0;
        %First run AeroCalc with simple blade-element model. Lift is
        %accurate with this simple model since Cl is pre-defined
        [Fblade, vi, phi, Re, Cd, ring] = AeroCalc(yN, rho, visc, vw, vc, b, h, Omega, cE, Cl, Cm, t, xtU, xtL, yWire, zWire, tWire, d, q, anhedral, ycmax(1), flags);

        %Then run StrCalc (simply to determine the spar deflection for accurate
        %ground effect computation)
        [Mtot, mSpar, mChord, mQuad, mCover, mWire, EIx, EIz, EA, GJ, q, EIQuad, GJQuad, Finternal, strain, fail] ...
        = StrCalc(yN, R, b, cE, xEA, xtU, d, theta, nTube, nCap, Jprop, yWire, zWire, tWire, TWire, TEtension, ycmax(1), lBiscuit, ...
        dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, RQuad, hQuad, mElseRotor, mElseCentre, mElseR, mPilot, Fblade, presLoad, flags);

        %Then run more accurate Vortex method (if the FreeWake flag is set)
        flags.FreeWake = flagFreeWakeTemp;
        [Fblade, vi, phi, Re, Cd, ring] = AeroCalc(yN, rho, visc, vw, vc, b, h, Omega, cE, Cl, Cm, t, xtU, xtL, yWire, zWire, tWire, d, q, anhedral, ycmax(1), flags);
    else
        [Fblade, vi, phi, Re, Cd, ring] = AeroCalc(yN, rho, visc, vw, vc, b, h, Omega, cE, Cl, Cm, t, xtU, xtL, yWire, zWire, tWire, d, q, anhedral, ycmax(1), flags);
    end
else
    Fblade.Fz = zeros(Ns,1);
    Fblade.Fx = zeros(Ns,1);
    Fblade.My = zeros(Ns,1);
    Fblade.Q = zeros(Ns,1);
    Fblade.P = zeros(Ns,1);
    Fblade.Pi = zeros(Ns,1);
    Fblade.Pp = zeros(Ns,1);
    vi = zeros(Ns,1);
    phi = zeros(Ns,1);
    Re = zeros(Ns,1);
    Cd = zeros(Ns,1);
    Cl = zeros(Ns,1);
    ring = [];
end

% Structural Calculation
%-----------------------
% yN - node locations
% Mtot - total helicopter mass
% mSpar - mass of each spar element
% mChord - mass ribs, skin, trailing edge, LES for each spanwise element
% mQuad - mass of quad rotor struts
% q - deformation
% fosOutPlane - fos on compressive fibre failure in the caps or tube
% fosInPlane - fos on....*don't know the structure yet
% fosShear - fos on fibre shear failure in tube from torsion and vertical loads
% fosBuckOutPlane - fos on Euler buckling in out-of-plane direction
% fosBuckInPlane - fos on Euler buckling in plane
% fosBuckShear - fos on torsional buckling
% fosBuckQuad - fos on Euler buckling of quad rotor struts
% cE - chord of each element
% d - spar diameter distribution
% theta - CFRP wrap angle distribution
% nTube - number of tube layers distribution
% nCap - number of cap strips distribution
% yWire - location of wire attachment along span
% zWire - depth of wire attachement
% tWire - thickness of wire
% lBiscuit - unsupported buiscuit length
% dQuad - diameter of quad rotor struts
% thetaQuad - wrap angle of quad rotor struts
% nTubeQuad - number of CFRP layers in quad rotor struts
% RQuad - distance from centre of helicopter to centre of quad rotors
% hQuad - height of quad-rotor truss

% Perform structural calculation once more with more accurate idea of drag
[Mtot, mSpar, mChord, mQuad, mCover, mWire, EIx, EIz, EA, GJ, q, EIQuad, GJQuad, Finternal, strain, fail] ...
    = StrCalc(yN, R, b, cE, xEA, xtU, d, theta, nTube, nCap, Jprop, yWire, zWire, tWire, TWire, TEtension, ycmax(1), lBiscuit, ...
    dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, RQuad, hQuad, mElseRotor, mElseCentre, mElseR, mPilot, Fblade, presLoad, flags);

% Compute aerodynamic jig angle
alphaJig = zeros(size(cE));
qq(:,1) = [0 0 0 0 0 0];
for s = 2:Ns+1
    qq(:,s) = q((s-1)*6+1:(s-1)*6+6);
end
Clalpha = 2*pi;
for s = 1:length(yN)-1
    alphaJig(s) = Cl(s)/Clalpha - (qq(5,s)+qq(5,s+1))/2 + phi(s) - collective;
end

% Compute dihedral angle
di = zeros(Ns,1);
for s = 1:Ns
    di(s) = atan2(qq(3,s+1)-qq(3,s), yN(s+1)-yN(s));
end
    
    
% Compute totals
if flags.Quad % Quad-rotor configuration
    Ttot = sum(Fblade.Fz.*cos(di))*b*4;
    Qtot = sum(Fblade.Q)*b*4;
    MomRot = sum(Fblade.Fz.*yE);
    Pitot = sum(Fblade.Pi)*b*4;
    Pptot = sum(Fblade.Pp)*b*4;
    Ptot = Pptot + Pitot; %Non-covered centre
    
else % Prop-driven configuration
    Ttot = sum(Fblade.Fz.*cos(di))*b;
    Qtot = sum(Fblade.Q)*b;
    MomRot = sum(Fblade.Fz.*yE);
    Ptot = sum(Fblade.P)*b*(1/etaP);
    Pitot = sum(Fblade.Pi)*b*(1/etaP);
    Pptot = sum(Fblade.Pp)*b*(1/etaP);
end

out.Ptot = Ptot;
out.Pitot = Pitot;
out.Pptot = Pptot;
out.Ttot = Ttot;
out.Mtot = Mtot;
out.MomRot = MomRot;
out.cE = cE;
out.cN = cN;
out.c100 = c100;
out.Cl = Cl;
out.Cm = Cm;
out.t = t;
out.xtU = xtU;
out.xtL = xtL;
out.alphaJig = alphaJig;
out.xEA = xEA;
out.d = d;
out.theta = theta;
out.nTube = nTube;
out.nCap = nCap;
out.lBiscuit = lBiscuit;
out.yN = yN;
out.Fblade = Fblade;
out.vi = vi;
out.phi = phi;
out.Re = Re;
out.Cd = Cd;
out.ring = ring;
out.mSpar = mSpar;
out.mChord = mChord;
out.mQuad = mQuad;
out.mCover = mCover;
out.mWire = mWire;
out.EIx = EIx;
out.EIz = EIz;
out.EA = EA;
out.GJ = GJ;
out.q = q;
out.EIQuad = EIQuad;
out.GJQuad = GJQuad;
out.Finternal = Finternal;
out.strain = strain;
out.fail = fail;

