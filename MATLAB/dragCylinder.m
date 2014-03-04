%------- dragCylinder ------%
%                           %
%       Todd Reichert       %
%       April 11, 2011      %
%                           %
%---------------------------%


function Cd = dragCylinder(Re,d)

    if Re < 3500
        Cd = -1E-10*Re^3 + 7E-7*Re^3 - 0.0013*Re + 1.7397;
    else
        Cd = 1;
    end

end