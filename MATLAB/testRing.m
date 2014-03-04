% testRing is used to check the discretized integral of the Biot-Savard law
% against the analytic integral given by Yoon and Heister.

% At about Ntheta = 20 discrete steps the discrete method converges on the
% analytic result.

function testRing

    close all
    clear all

    zi = 0;
    ri = 10;
    z = 0.2;
    r = 9.8;
    Ntheta = [5 10 15 20 25 30 35 40];
    Gamma = 1;
    
    for i = 1:length(Ntheta)
        [urD(i) uzD(i)] = discreteMethod(Gamma,zi,ri,z,r,Ntheta(i));
        [urA uzA] = analyticMethod(Gamma,zi,ri,z,r);
    end
    
    subplot(1,2,1)
    plot(Ntheta,urD,0,urA,'o');
    title('Ur');
    subplot(1,2,2)
    plot(Ntheta,uzD,0,uzA,'o');
    title('Uz');

end

function [ur uz] = discreteMethod(Gamma,zi,ri,z,r,Ntheta)

    cr = 0.5;

    dtheta = pi/Ntheta;
    thetaArray = linspace(0+dtheta/2,pi-dtheta/2,Ntheta);

    M = Gamma*ri*dtheta/(2*pi);
    X2 = (-ri*sin(thetaArray)).^2;
    Y2 = (r-ri*cos(thetaArray)).^2;
    Z2 = (z-zi).^2;
    Normal = sqrt(X2 + Y2 + Z2);
    for iii = 1:length(thetaArray)
        if Normal(iii) < cr
            Normal(iii) = cr;
        end
    end
    Norm3 = Normal.^3;

    ur =  sum(-cos(thetaArray)*(z-zi)./Norm3)*M;
    uz =  sum((cos(thetaArray)*r-ri)./Norm3)*M;

end

function [ur uz] = analyticMethod(Gamma,zi,ri,z,r)

    a = sqrt((r+ri)^2 + (z-zi)^2);
    m = 4*r*ri/a^2;

    a0 = 1.3862944;
    a1 = 0.1119723;
    a2 = 0.0725296;
    b0 = 0.5;
    b1 = 0.1213478;
    b2 = 0.0288729;
    K = (a0+a1*m+a2*m^2) + (b0+b1*m+b2*m^2)*log(1/m);

    a1 = 0.4630151;
    a2 = 0.1077812;
    b1 = 0.2452727;
    b2 = 0.0412496;
    E = (1+a1*m+a2*m^2) + (b1*m+b2*m^2)*log(1/m);
    
    [K,E] = ellipke(m);

    I1 = 4/a*K;
    I2 = 4/a^3*E/(1-m);

    A = (z-zi)^2 + r^2 + ri^2;
    B = -2*r*ri;

    uz = -Gamma/(4*pi)*ri*((ri+r*A/B)*I2 - r/B*I1);
    ur = -Gamma/(4*pi)*ri*(z-zi)/B*(I1-A*I2);
end