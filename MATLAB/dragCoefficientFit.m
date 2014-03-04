%--- dragCoefficientFit ----%
%                           %
%       Todd Reichert       %
%        Feb 1, 2011        %
%                           %
%---------------------------%

% dragCoefficient.m returns the drag coefficient of an airfoil at Reynolds 
% number Re, with thickness to chord ratio tc, and with xtcU and xtcL
% fraction of laminar flow on the upper and lower surfaces respectively.
% The result if a fit on existing HPH airfoils.


function Cd = dragCoefficientFit(Re,tc,xtcU,xtcL)


    Cf15_15 = 0.6798*Re^(-0.283);
    Cf60_100 = 22.09*Re^(-0.604);
    
    xtc = xtcU+xtcL;
    Cd = Cf15_15 + (Cf60_100-Cf15_15)*(xtc-0.3)/(1.6-0.3);
    
end