%-------- AeroCalc ---------%
%                           %
%       Todd Reichert       %
%        March 5, 2011      %
%                           %
%---------------------------%

% AeroCalc computes the forces on a single rotor blade.

function [Fblade, vi, phi, Re, Cd, ring] = ...
    AeroCalc(yN, rho, visc, vw, vc, b, h, Omega, c, Cl, Cm, t, xtU, xtL, yWire, zWire, tWire, d, q, anhedral, ycmax, flags)

Ns = length(yN)-1; %number of elements
dr = zeros(Ns,1);
r = zeros(Ns,1);
for s=1:Ns
    dr(s) = yN(s+1) - yN(s); %length of each element
    r(s) = 0.5*(yN(s) + yN(s+1)); %radial location of each element
end
R = yN(Ns+1);

Fblade.Fx = zeros(Ns,1);
Fblade.Fz = zeros(Ns,1);
Fblade.My = zeros(Ns,1);
Fblade.Q = zeros(Ns,1);
Fblade.P = zeros(Ns,1);
Fblade.Pi = zeros(Ns,1);
Fblade.Pp = zeros(Ns,1);
phi = zeros(Ns,1);
vi = zeros(Ns,1);
dT = zeros(Ns,1);
Re = zeros(Ns,1);
Cd = zeros(Ns,1);
chordFrac = ones(Ns,1);

% Compute multiplyer for partial element
for s = 1:Ns+1
    if yN(s) < ycmax
        sTrans = s;   %determine transitional partial element
    end
end
chordFrac(sTrans) = (yN(sTrans+1)-ycmax) / (yN(sTrans+1)-yN(sTrans)); %Fraction of element that has a chord

% Compute thrust assuming small angles
for s = 1:Ns
    dT(s) = chordFrac(s)*0.5*rho*(Omega*r(s))^2*Cl(s)*c(s)*dr(s);
end

if flags.FreeWake
    
    % Set fidelity
    if Ns == 15
        Nw = 15;
        Ntt = 5;
        Ntheta = 40;
    elseif Ns == 10
        Nw = 8;
        Ntt = 3;
        Ntheta = 20;
    else
        Nw = 5;
        Ntt = 1;
        Ntheta = 15;
    end
    
    % Break out deformations
    qq(:,1) = [0 0 0 0 0 0];
    for s = 2:Ns+1
        qq(:,s) = q((s-1)*6+1:(s-1)*6+6);
    end
    qh = qq(3,:)-yN'*anhedral;
    
    [vi, ring] = VortexWakeCover(yN, rho, dT, vc, Omega, b, h, Nw, Ntt, Ntheta, qh, ycmax, flags.Cover, flags.PlotWake, flags.DynamicClimb);
    
    %vi = vi+0.3; %add vertical velocity
    
else
    % Compute induced velocity using annual-ring actuator disk theory
    for s = 1:Ns
        vi(s) = -0.5*vc + sqrt(0.25*vc^2 + 0.25*b*dT(s)/(pi*rho*r(s)*dr(s))); %where vc is vertical velocity
        %vi(s) = sqrt(0.25*b*dT(s)/(pi*rho*r(s)*dr(s)));
    end

    % Add ground effect (Bramwell)
%     A = pi*R^2;
%     solidity = b*mean(c)/(pi*R);
%     tc = sum(dT)/(rho*solidity*A*Omega^2*R^2);
%     if h >= 0
%         vi = vi*(1 - 0.3925/((h/R)^2+0.1714*(h/R)+0.6497) / (1-0.553) * (0.3362*(tc/solidity)+0.8338)/((tc/solidity)+1.492));
%     end

    % Add ground effect Cheesemen & Benett's
    vi = vi/(1+(R/h/4)^2);

    % Add ground effect (vi fit)
%     hR = h/R;
%     rR = r/R;
%     if hR < 0.05
%         hR = 0.05;
%         disp('h out of bounds for ground effect model')
%     end
%     if hR > 0.5
%         hR = 0.05;
%         disp('h out of bounds for ground effect model')
%     end
%     vieps = 0.004467*exp(1.109*rR)+2.4466E-8*exp(18.04*rR);
%     vislope = -148.8*rR.^6 + 385.7*rR.^5 - 380.5*rR.^4 + 178.5*rR.^3 - 38.92*rR.^2 + 3.704*rR + 0.6167;
%     vi = vi.*(vieps + vislope*(hR-0.05));
    
    ring = [];
end

% Compute lift and drag using full angles
for s = 1:Ns
    U = sqrt( (Omega*r(s)+vw)^2 + (vc+vi(s))^2 ); %where vw is wind, vc is vertical velocity
    
    % Wing section
    if c(s) > 0.001
        Re(s) = rho*U*c(s)/visc;
        Cd(s) = dragCoefficientFit(Re(s),t(s),xtU(s),xtL(s));
        dL = 0.5*rho*U^2*Cl(s)*c(s)*dr(s);
        dD = 0.5*rho*U^2*Cd(s)*c(s)*dr(s);
    % Root spar section
    else
        Re(s) = rho*U*d(s)/visc;
        if Re < 3500
            Cd(s) = -1E-10*Re(s)^3 + 7E-7*Re(s)^2 - 0.0013*Re(s) + 1.7397;
        else
            Cd(s) = 1;
        end
        dL = 0;
        dD = 0.5*rho*U^2*Cd(s)*d(s)*dr(s);
    end
    
    % Add wire drag
    for w = 1:length(yWire)
        if yN(s) < yWire(w)
            if yN(s+1) < yWire(w)
                L = dr(s)*sqrt(zWire^2+yWire(w)^2)/yWire(w);
            else
                L = (yWire(w) - yN(s))*sqrt(zWire^2+yWire(w)^2)/yWire(w);
            end
            ReWire = rho*U*tWire/visc;
            CdWire = -1E-10*ReWire^3 + 7E-7*ReWire^2 - 0.0013*ReWire + 1.7397;
            dD = dD + 0.5*rho*U^2*CdWire*tWire*L;
        end
    end
    
    phi(s) = atan2(vc+vi(s),vw+Omega*r(s)); %where vw is wind, vc is vertical velocity
    Fblade.Fz(s) = chordFrac(s)*(dL*cos(phi(s)) - dD*sin(phi(s)));
    Fblade.Fx(s) = chordFrac(s)*(dD*cos(phi(s)) + dL*sin(phi(s)));
    Fblade.My(s) = chordFrac(s)*(0.5*rho*U^2*Cm(s)*c(s)*c(s)*dr(s));
    Fblade.Q(s) = Fblade.Fx(s)*r(s);
    Fblade.P(s) = Fblade.Q(s)*Omega;
    Fblade.Pi(s) = chordFrac(s)*(dL*sin(phi(s))*r(s)*Omega);
    Fblade.Pp(s) = chordFrac(s)*(dD*cos(phi(s))*r(s)*Omega);
end

end
    
    
