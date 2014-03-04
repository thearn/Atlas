%--- DiscretizeProperties --%
%                           %
%       Todd Reichert       %
%        March 5, 2011      %
%                           %
%---------------------------%

% Discretize properties along rotor blade. Y defines the locations at which
% the properties are defined. Properties are linearly interpolated between
% Y locations.

function [cE, cN, c100, Cl, Cm, t, xtU, xtL, xEA, d, theta, nTube, nCap, lBiscuit, yN, yE] ...
    = DiscretizeProperties(Ns, ycmax, R, c_, Cl_, Cm_, t_, xtU_, xtL_, xEA_, yWire, d_, theta_, nTube_, nCap_, lBiscuit_)

yN = zeros(Ns+1,1);
cN = zeros(Ns+1,1);

yE = zeros(Ns,1);
cR = zeros(Ns,1);
cZ = zeros(Ns,1);
cE = zeros(Ns,1);
Cl = zeros(Ns,1);
Cm = zeros(Ns,1);
t = zeros(Ns,1);
xtU = zeros(Ns,1);
xtL = zeros(Ns,1);
xEA = zeros(Ns,1);
d = zeros(Ns,1);
theta = zeros(Ns,1);
nTube = zeros(Ns,1);
nCap = zeros(Ns,1);
lBiscuit = zeros(Ns,1);

Ns100 = 100;
yN100 = zeros(Ns100+1,1);
yE100 = zeros(Ns100,1);
cR100 = zeros(100,1);
cZ100 = zeros(100,1);
c100 = zeros(100,1);

% Compute node locations
for s = 1:Ns+1
    yN(s) = R/Ns*(s-1); %linear distribution of nodes
    if yN(s) < ycmax(1)
        sTrans(1) = s;
    end
    if yN(s) < ycmax(2)
        sTrans(2) = s;
    end
end

% Compute element locations
for s=1:Ns
    yE(s) = 0.5*(yN(s) + yN(s+1));
end

% Compute chord length at elements
for s = 1:Ns
    
    % If root section
    if s < sTrans(1)
        x = yE(s) / ycmax(1);
        cE(s) = c_(1) + x*(c_(2)-c_(1));
        
    % If chord section
    else
        % Compute curve component
        x = (yE(s)-ycmax(1)) / (R-ycmax(1));
        pStart = c_(3);
        pCurve = c_(4);
        pEnd = c_(5);
        xx = x*(1-pCurve) + sin(x*pi/2)*pCurve;   
        cZ(s)= pStart+(pEnd-pStart)*xx;
    
        % Compute 1/r component
        c3 = c_(3)/(c_(5)*R/ycmax(1));
        cR(s) = c_(5)*R/yE(s)*(c3+(1-c3)*x);
        
        % Average based on c_(2)
        cE(s) = cR(s) + (cZ(s)-cR(s))*c_(2);
    end
    
    if cE(s) == 0
        cE(s) = 0.001;
    end
end

% Compute chord for display purposes
for s = 1:Ns100+1
    yN100(s) = R/Ns100*(s-1); %linear distribution of nodes
    if yN100(s) < ycmax
        sTrans100 = s;
    end
end

for s=1:Ns100
    yE100(s) = 0.5*(yN100(s) + yN100(s+1));
end

for s = 1:Ns100
    
    % If root section
    if s < sTrans100
        x = yE100(s) / ycmax(1);
        c100(s) = c_(1) + x*(c_(2)-c_(1));
        
    % If chord section
    else
        % Compute curve component
        x = (yE100(s)-ycmax(1)) / (R-ycmax(1));
        pStart = c_(3);
        pCurve = c_(4);
        pEnd = c_(5);
        xx = x*(1-pCurve) + sin(x*pi/2)*pCurve;   
        cZ100(s)= pStart+(pEnd-pStart)*xx;
    
        % Compute 1/r component
        c3 = c_(3)/(c_(5)*R/ycmax(1));
        cR100(s) = c_(5)*R/yE100(s)*(c3+(1-c3)*x);
        
        % Average based on c_(2)
        c100(s) = cR100(s) + (cZ100(s)-cR100(s))*c_(2);
    end
    
    if c100(s) == 0
        c100(s) = 0.001;
    end
end


% Compute aero properties for each element
Y = [ycmax(1) ycmax(2) R];
for s = 1:Ns
    
    % If root section
    if s < sTrans(1)
        Cl(s) = Cl_(1);
        Cm(s) = Cm_(1);
        t(s) =  t_(1);
        xEA(s) = xEA_(1);
    
    else

        % Check which segment the element is on
        for j = 1:length(Y)
            if Y(j) > yE(s)
                break
            end
        end
        if s == sTrans(1)
            j = 2;
        end
        x = (yE(s)-Y(j-1)) / (Y(j)-Y(j-1)); 

        % Linearly interpolate between Y(j) and Y(j-1)
        Cl(s) = Cl_(j-1) + x*(Cl_(j) - Cl_(j-1));
        Cm(s) = Cm_(j-1) + x*(Cm_(j) - Cm_(j-1));
        t(s) = t_(j-1) + x*(t_(j) - t_(j-1));
        xEA(s) = xEA_(j-1) + x*(xEA_(j) - xEA_(j-1));
    end
end
Cl(Ns) = Cl(Ns)*2/3;
        
% Compute xtU and xtL for each element
% Changes instantly from xtU(1) to xtU(3) at point xtU(2)
for s = 1:Ns+1
    if yN(s) < ycmax(1)
        sTrans(1) = s;
    end
    if yN(s) < xtU_(2)
        sTrans(2) = s;
    end
end
for s = 1:Ns
    
    % If root section
    if s < sTrans(1)
        xtU(s) = 0.05;
        xtL(s) = 0.05;
    
    elseif s < sTrans(2)
        xtU(s) = xtU_(1);
        xtL(s) = xtL_(1);
        
    elseif s == sTrans(2);
        x = (xtU_(2)-yN(s))/(yN(s)-yN(s-1));
        xtU(s) = (1-x)*xtU_(1) + x*xtU_(3);
        xtL(s) = (1-x)*xtL_(1) + x*xtL_(3);
        
    elseif s > sTrans(2)
        xtU(s) = xtU_(3);
        xtL(s) = xtL_(3);
        
    end
end


% Compute str properties for each element
Y = [0 yWire(1) R];
for s=1:Ns
    
    % Check which segment the element is on
    for j = 1:length(Y)
        if Y(j) > yE(s)
            break
        end
    end
    x = (yE(s)-Y(j-1)) / (Y(j)-Y(j-1)); 
    
    % Linearly interpolate between Y(j) and Y(j-1)
    d(s) = d_(j-1) + x*(d_(j) - d_(j-1));
    theta(s) = theta_(j-1) + x*(theta_(j) - theta_(j-1));
    nTube(s) = nTube_(j-1) + x*(nTube_(j) - nTube_(j-1));
    nCap(s) = nCap_(j-1) + x*(nCap_(j) - nCap_(j-1));
    lBiscuit(s) = lBiscuit_(j-1) + x*(lBiscuit_(j) - lBiscuit_(j-1));
end


end
