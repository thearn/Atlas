%----- dragCoefficient -----%
%                           %
%       Todd Reichert       %
%        Feb 1, 2011        %
%                           %
%---------------------------%

% dragCoefficient.m returns the drag coefficient of an airfoil at Reynolds 
% number Re, with thickness to chord ratio tc, and with xtcU and xtcL
% fraction of laminar flow on the upper and lower surfaces respectively.


function Cd = dragCoefficient(Re,tc,xtcU,xtcL)
    CfU = frictionCoefficient(Re,xtcU); 
    CfL = frictionCoefficient(Re,xtcL);
    Cfflat = (CfU + CfL)/2;

    Cd = 2*Cfflat*(1 + 2*tc + 60*(tc)^4);

end

% Computes Cf of a flat plate at Re, with xtc fraction of laminar flow
function Cfflat = frictionCoefficient(Re,xtc)
    Re,xtc
    % Fully turbulent
    if xtc == 0 
        Cfflat = 0.072/(Re)^0.2;
    
    % Fully laminar
    elseif xtc == 1
        Cfflat = (1.328/sqrt(Re));
    
    % Partially laminar
    else
        Cflam = (1.328/sqrt(Re))*xtc^(-0.5); %Cf of laminar part
        deltalamc = (5/sqrt(Re))*sqrt(xtc);  %boundary layer thickness, delta/c, of laminar part
        deltaturbc = (0.13/0.097)*deltalamc; %boundary layer thickness, delta/c, of turbulent part
        x0c = xtc - (Re^0.2*deltaturbc/0.375)^(1/0.8); %imaginary start point of turbulent BL
        CfturbFull = 0.072/((1-x0c)*Re)^0.2; %Cf of flat plate of length c-x0
        CfturbStart = 0.072/((xtc-x0c)*Re)^0.2; %Cf of imaginary part of turb BL
        Cfturb = (CfturbFull*(1-x0c) - CfturbStart*(xtc-x0c))/(1-xtc); %Cf of turbulent part
        Cfflat = Cflam*xtc + Cfturb*(1-xtc);
    end
    Cfflat
end