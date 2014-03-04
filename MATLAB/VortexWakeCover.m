%-------- VortexWake -------%
%                           %
%       Todd Reichert       %
%        March 27, 2011     %
%                           %
%---------------------------%

% VortexWake computes the induced velocity on the rotor blades given the
% thrust distribution on the rotor.

function [vi, ring] = VortexWakeCover(yN, rho, dT, vc, Omega, b, h, Nw, Ntt, Ntheta, qh, ycmax, flagCover, flagPlot, flagDynamicClimb)
%tic

flagHelix = 1;
flagMovie = 0;

Ns = length(yN)-1; %number of elements
dy = zeros(Ns,1);
yE = zeros(Ns,1);
for s=1:Ns
    dy(s) = yN(s+1) - yN(s); %length of each element
    yE(s) = 0.5*(yN(s) + yN(s+1)); %radial location of each element
end
R = yN(Ns+1);

cr = 0.5*mean(dy);
dtheta = pi/Ntheta;
thetaArray = linspace(0+dtheta/2,pi-dtheta/2,Ntheta);

if flagPlot
    scrsz = get(0,'ScreenSize');
    figsz = round(scrsz*0.8);
    if flagMovie
        figure('Position',[(scrsz(3)-figsz(3))/2 (scrsz(4)-figsz(4))/2 figsz(3) figsz(4)/2])
    else
        figure('Position',[scrsz(3)*0.1 scrsz(4)*0.25 scrsz(3)*0.8 scrsz(4)/2])
    end  
    set(gcf,'color',[1 1 1]);
end



% Pre-allocate memory
ring.Gamma = zeros(Nw+1,Ns+1);
ring.z = zeros(Nw+1,Ns+1);
ring.r = zeros(Nw+1,Ns+1);
ring.vz = zeros(Nw+1,Ns+1);
ring.vr = zeros(Nw+1,Ns+1);
vi = zeros(Ns,1);

if flagCover
    for s = 1:Ns+1
        if yN(s) < ycmax
            sCover = s;
        end
    end
    C = zeros(sCover,sCover);
    cover.Gamma = zeros(sCover,1);
    cover.z = qh(2:sCover+1);
    cover.r = yN(2:sCover+1);
    cover.vi = zeros(sCover,1);
end

% Create nacent vortex rings
GammaBound = dT./(rho*(Omega*yE).*dy);
ring.Gamma(1,1) = -GammaBound(1);
for s = 2:Ns
    ring.Gamma(1,s) = GammaBound(s-1) - GammaBound(s);
end
ring.Gamma(1,Ns+1) = GammaBound(Ns);
ring.r(1,:) = yN;
ring.z(1,:) = qh(:);



% Free-wake time stepping
for t = 1:Nw
    
    % Proceed with substeps
    for tt = 1:Ntt

        %% Compute cover circulation
        if flagCover
            
            % Compute the induced velocity on callocation point s from all iix(Ns+1) rings
            cover.vi = zeros(sCover,1);
            for s = 1:sCover
                for ii = 1:t %add the velocity induced from each disk
                    for ss = 2:Ns+1 %and each ring on each disk (inner ring cancels itself out)

                        zr = ring.z(ii,ss);
                        r  = ring.r(ii,ss);
                        zp = (qh(s)+qh(s+1))/2;
                        yp = yE(s);

                        M = ring.Gamma(ii,ss)*r*dtheta/(2*pi);
                        X2 = (-r*sin(thetaArray)).^2;
                        Y2 = (yp-r*cos(thetaArray)).^2;
                        Z2 = (zp-zr).^2;
                        Normal = sqrt(X2 + Y2 + Z2);
                        for iii = 1:length(thetaArray)
                            if Normal(iii) < cr
                                Normal(iii) = cr;
                            end
                        end
                        Norm3 = Normal.^3;

                        cover.vi(s) = cover.vi(s) + sum((cos(thetaArray)*yp-r)./Norm3)*M;

                        % Ground effect ring
                        zr = -2*h - zr;
                        Z2 = (zp-zr).^2;
                        Normal = sqrt(X2 + Y2 + Z2);
                        for iii = 1:length(thetaArray)
                            if Normal(iii) < cr
                                Normal(iii) = cr;
                            end
                        end
                        Norm3 = Normal.^3;

                        cover.vi(s) = cover.vi(s) - sum((cos(thetaArray)*yp-r)./Norm3)*M;

                    end
                end
            end
            
            % Compute the induced velocity coefficient on point s from all cover vortices
            for s = 1:sCover %infuence on s
                for ss = 1:sCover %from ss

                    zr = cover.z(ss);
                    r  = cover.r(ss);
                    zp = (qh(s)+qh(s+1))/2;
                    yp = yE(s);

                    M = 1*r*dtheta/(2*pi);
                    X2 = (-r*sin(thetaArray)).^2;
                    Y2 = (yp-r*cos(thetaArray)).^2;
                    Z2 = (zp-zr).^2;
                    Normal = sqrt(X2 + Y2 + Z2);
                    for iii = 1:length(thetaArray)
                        if Normal(iii) < cr
                            Normal(iii) = cr;
                        end
                    end
                    Norm3 = Normal.^3;

                    C(s,ss) = sum((cos(thetaArray)*yp-r)./Norm3)*M;

                    % Ground effect ring
                    zr = -2*h - zr;
                    Z2 = (zp-zr).^2;
                    Normal = sqrt(X2 + Y2 + Z2);
                    for iii = 1:length(thetaArray)
                        if Normal(iii) < cr
                            Normal(iii) = cr;
                        end
                    end
                    Norm3 = Normal.^3;

                    C(s,ss) = C(s,ss) - sum((cos(thetaArray)*yp-r)./Norm3)*M;

                end
            end
            
            % Solve for cover vortex strengths
            cover.Gamma = C\(-cover.vi);
        end
        

        %% Compute induced velocity on all ix(Ns+1) rings from all iix(Ns+1) rings
        ring.vz = zeros(Nw+1,Ns+1);
        ring.vr = zeros(Nw+1,Ns+1); 
        for i = 1:t %for each disk
            for s = 1:Ns+1 %and for each ring on each disk

                %v = [0 0 0];
                for ii = 1:t %add the velocity induced from each disk
                    for ss = 2:Ns+1 %and each ring on each disk (inner ring cancels itself out)

                        zr = ring.z(ii,ss);
                        r  = ring.r(ii,ss);
                        zp = ring.z(i,s);
                        yp = ring.r(i,s);

                        M = ring.Gamma(ii,ss)*r*dtheta/(2*pi);
                        X2 = (-r*sin(thetaArray)).^2;
                        Y2 = (yp-r*cos(thetaArray)).^2;
                        Z2 = (zp-zr).^2;
                        Normal = sqrt(X2 + Y2 + Z2);
                        for iii = 1:length(thetaArray)
                            if Normal(iii) < cr
                                Normal(iii) = cr;
                            end
                        end
                        Norm3 = Normal.^3;

                        ring.vr(i,s) = ring.vr(i,s) + sum(-cos(thetaArray)*(zp-zr)./Norm3)*M;
                        ring.vz(i,s) = ring.vz(i,s) + sum((cos(thetaArray)*yp-r)./Norm3)*M;

                        % Ground effect ring
                        zr = -2*h - ring.z(ii,ss);
                        Z2 = (zp-zr).^2;
                        Normal = sqrt(X2 + Y2 + Z2);
                        for iii = 1:length(thetaArray)
                            if Normal(iii) < cr
                                Normal(iii) = cr;
                            end
                        end
                        Norm3 = Normal.^3;

                        ring.vr(i,s) = ring.vr(i,s) - sum(-cos(thetaArray)*(zp-zr)./Norm3)*M;
                        ring.vz(i,s) = ring.vz(i,s) - sum((cos(thetaArray)*yp-r)./Norm3)*M;

                    end
                end
            end
        end
        
        %% Compute induced velocity on all ix(Ns+1) rings from all cover vortices
        if flagCover
            for i = 1:t %for each disk
                for s = 1:Ns+1 %and for each ring on each disk

                    for ss = 1:sCover %add the velocity from each cover vortice

                        zr = cover.z(ss);
                        r  = cover.r(ss);
                        zp = ring.z(i,s);
                        yp = ring.r(i,s);

                        M = cover.Gamma(ss)*r*dtheta/(2*pi);
                        X2 = (-r*sin(thetaArray)).^2;
                        Y2 = (yp-r*cos(thetaArray)).^2;
                        Z2 = (zp-zr).^2;
                        Normal = sqrt(X2 + Y2 + Z2);
                        for iii = 1:length(thetaArray)
                            if Normal(iii) < cr
                                Normal(iii) = cr;
                            end
                        end
                        Norm3 = Normal.^3;

                        ring.vr(i,s) = ring.vr(i,s) + sum(-cos(thetaArray)*(zp-zr)./Norm3)*M;
                        ring.vz(i,s) = ring.vz(i,s) + sum((cos(thetaArray)*yp-r)./Norm3)*M;

                        % Ground effect ring
                        zr = -2*h - zr;
                        Z2 = (zp-zr).^2;
                        Normal = sqrt(X2 + Y2 + Z2);
                        for iii = 1:length(thetaArray)
                            if Normal(iii) < cr
                                Normal(iii) = cr;
                            end
                        end
                        Norm3 = Normal.^3;

                        ring.vr(i,s) = ring.vr(i,s) - sum(-cos(thetaArray)*(zp-zr)./Norm3)*M;
                        ring.vz(i,s) = ring.vz(i,s) - sum((cos(thetaArray)*yp-r)./Norm3)*M;

                    end
                end
            end
        end
    
        % Compute altitude and time power approximation
        if tt == 1
            PiApprox = 8*sum(dT.*vi);
        end
        realtime = 2*pi/Omega/b*((t-1)*Ntt+tt)/Ntt;
        altitude = vc*realtime;
        
        
    
        %% Plot
        if flagPlot
            
            % Compute induced velocity on rotor from rings (rp = [0 r(s) 0])
            for s = 1:Ns %for each element
                vi(s) = 0;
                for ii = 1:t %add the velocity induced from each disk
                    for ss = 2:Ns+1 %and each ring on each disk (inner ring cancels itself out)
                        
                        if (ii == 1) && flagHelix
                            ringFrac = 0.675;
                        else
                            ringFrac = 1;
                        end
                        
                        zr = ring.z(ii,ss);
                        r  = ring.r(ii,ss);
                        zp = (qh(s)+qh(s+1))/2;
                        yp = yE(s);

                        M = ring.Gamma(ii,ss)*r*dtheta/(2*pi);
                        X2 = (-r*sin(thetaArray)).^2;
                        Y2 = (yp-r*cos(thetaArray)).^2;
                        Z2 = (zp-zr).^2;
                        Normal = sqrt(X2 + Y2 + Z2);
                        for iii = 1:length(thetaArray)
                            if Normal(iii) < cr
                                Normal(iii) = cr;
                            end
                        end
                        Norm3 = Normal.^3;

                        vi(s) = vi(s) + ringFrac*sum((cos(thetaArray)*yp-r)./Norm3)*M;

                        % Ground effect ring
                        zr = -2*h - ring.z(ii,ss);
                        Z2 = (zp-zr).^2;
                        Normal = sqrt(X2 + Y2 + Z2);
                        for iii = 1:length(thetaArray)
                            if Normal(iii) < cr
                                Normal(iii) = cr;
                            end
                        end
                        Norm3 = Normal.^3;

                        vi(s) = vi(s) - ringFrac*sum((cos(thetaArray)*yp-r)./Norm3)*M;

                    end
                end
            end
            
            % Compute induced velocity on from cover (rp = [0 r(s) 0])
            if flagCover
                for s = 1:Ns %for each element
                    for ss = 1:sCover %add the velocity induced from the cover

                        if flagHelix
                            ringFrac = 0.675;
                        else
                            ringFrac = 1;
                        end
                        
                        zr = cover.z(ss);
                        r  = cover.r(ss);
                        zp = (qh(s)+qh(s+1))/2;
                        yp = yE(s);

                        M = cover.Gamma(ss)*r*dtheta/(2*pi);
                        X2 = (-r*sin(thetaArray)).^2;
                        Y2 = (yp-r*cos(thetaArray)).^2;
                        Z2 = (zp-zr).^2;
                        Normal = sqrt(X2 + Y2 + Z2);
                        for iii = 1:length(thetaArray)
                            if Normal(iii) < cr
                                Normal(iii) = cr;
                            end
                        end
                        Norm3 = Normal.^3;

                        vi(s) = vi(s) + ringFrac*sum((cos(thetaArray)*yp-r)./Norm3)*M;

                        % Ground effect ring
                        zr = -2*h - zr;
                        Z2 = (zp-zr).^2;
                        Normal = sqrt(X2 + Y2 + Z2);
                        for iii = 1:length(thetaArray)
                            if Normal(iii) < cr
                                Normal(iii) = cr;
                            end
                        end
                        Norm3 = Normal.^3;

                        vi(s) = vi(s) - ringFrac*sum((cos(thetaArray)*yp-r)./Norm3)*M;

                    end
                end
            end
            
            vi = -vi; %vi is positive downwards
            
            if flagMovie
                rr = [yN'; ring.r(1:t,:)];
                zz = [qh(:); ring.z(1:t,:)];
                set(gca,'nextplot','replacechildren','visible','off');
                
%                 f = getframe;
%                 [im,map] = rgb2ind(f.cdata,256,'nodither');
%                 im(1,1,1,20) = 0;
                
                plot(rr,zz,'k-o',-rr,zz,'k-o', [-R R], [0 0], 'k', [-100 100], [-h -h],'b')
                axis([-30 30 -5 5])
                f = getframe(gca);
                if (t == 1) && (tt == 1)
                    [im(:,:,1,(t-1)*(Ntt) + tt),map] = rgb2ind(f.cdata,256,'nodither');
                else
                    im(:,:,1,(t-1)*(Ntt) + tt) = rgb2ind(f.cdata,map,'nodither');
                end
            else
                subplot(1,2,1);
                rr = [yN'; ring.r(1:t,:)];
                
                if flagDynamicClimb
                    zz = [qh+altitude; ring.z(1:t,:)];
                else
                    zz = [qh; ring.z(1:t,:)];
                end                
                zz = [qh+altitude; ring.z(1:t,:)];
                plot(rr,zz,'k-o',-rr,zz,'k-o', [-R R], [0 0], 'k', [-100 100], [-h -h],'b')
                axis([-30 30 -10 10])
                text(-25,-3,['Pi = ',num2str(PiApprox,'%0.1f'),' W']);
                text(-25,-5,['t = ',num2str(realtime,'%0.1f'),' s']);
                text(-25,-7,['H = ',num2str(altitude,'%0.1f'),' m']);
                subplot(1,2,2);
                plot(yE,vi,'k',[-100 100],[0 0],'b');
                axis([0 R -1 1])
                pause(2);
            end
            
        end
        
        %% Convect rings downstream
        dt = 2*pi/Omega/b/Ntt;
        for i = 1:t %for each disk
            for s = 1:Ns+1 %and for each ring on each disk
                ring.z(i,s) = ring.z(i,s) + ring.vz(i,s)*dt;
                ring.r(i,s) = ring.r(i,s) + ring.vr(i,s)*dt;
                ring.Gamma(i,s) = ring.Gamma(i,s);
                if (flagCover) && (tt == 1) && (i==1) && (s <= sCover+1) && (s>1)
                    ring.Gamma(i,s) = ring.Gamma(i,s) + cover.Gamma(s-1);
                end
            end
        end

    end
    
    
    %% Shift elements in ring array
    for i = t:-1:1
        ring.Gamma(i+1,:) = ring.Gamma(i,:);
        ring.r(i+1,:) = ring.r(i,:);
        ring.z(i+1,:) = ring.z(i,:);
    end
    
    %% Create nacent vortex rings
    GammaBound = dT./(rho*(Omega*yE).*dy);
    ring.Gamma(1,1) = -GammaBound(1);
    for s = 2:Ns
        ring.Gamma(1,s) = GammaBound(s-1) - GammaBound(s);
    end
    ring.Gamma(1,Ns+1) = GammaBound(Ns);
    ring.r(1,:) = yN;
    if flagDynamicClimb
        ring.z(1,:) = qh(:) + altitude; %Add altitude to displacement
    else
        ring.z(1,:) = qh(:);
    end
end


%% Compute cover circulation
if flagCover

    % Compute the induced velocity on callocation point s from all iix(Ns+1) rings
    cover.vi = zeros(sCover,1);
    for s = 1:sCover
        for ii = 1:t %add the velocity induced from each disk
            for ss = 2:Ns+1 %and each ring on each disk (inner ring cancels itself out)

                zr = ring.z(ii,ss);
                r  = ring.r(ii,ss);
                zp = (qh(s)+qh(s+1))/2;
                yp = yE(s);

                M = ring.Gamma(ii,ss)*r*dtheta/(2*pi);
                X2 = (-r*sin(thetaArray)).^2;
                Y2 = (yp-r*cos(thetaArray)).^2;
                Z2 = (zp-zr).^2;
                Normal = sqrt(X2 + Y2 + Z2);
                for iii = 1:length(thetaArray)
                    if Normal(iii) < cr
                        Normal(iii) = cr;
                    end
                end
                Norm3 = Normal.^3;

                cover.vi(s) = cover.vi(s) + sum((cos(thetaArray)*yp-r)./Norm3)*M;

                % Ground effect ring
                zr = -2*h - zr;
                Z2 = (zp-zr).^2;
                Normal = sqrt(X2 + Y2 + Z2);
                for iii = 1:length(thetaArray)
                    if Normal(iii) < cr
                        Normal(iii) = cr;
                    end
                end
                Norm3 = Normal.^3;

                cover.vi(s) = cover.vi(s) - sum((cos(thetaArray)*yp-r)./Norm3)*M;

            end
        end
    end

    % Compute the induced velocity coefficient on point s from all cover vortices
    for s = 1:sCover %infuence on s
        for ss = 1:sCover %from ss

            zr = cover.z(ss);
            r  = cover.r(ss);
            zp = (qh(s)+qh(s+1))/2;
            yp = yE(s);

            M = 1*r*dtheta/(2*pi);
            X2 = (-r*sin(thetaArray)).^2;
            Y2 = (yp-r*cos(thetaArray)).^2;
            Z2 = (zp-zr).^2;
            Normal = sqrt(X2 + Y2 + Z2);
            for iii = 1:length(thetaArray)
                if Normal(iii) < cr
                    Normal(iii) = cr;
                end
            end
            Norm3 = Normal.^3;

            C(s,ss) = sum((cos(thetaArray)*yp-r)./Norm3)*M;

            % Ground effect ring
            zr = -2*h - zr;
            Z2 = (zp-zr).^2;
            Normal = sqrt(X2 + Y2 + Z2);
            for iii = 1:length(thetaArray)
                if Normal(iii) < cr
                    Normal(iii) = cr;
                end
            end
            Norm3 = Normal.^3;

            C(s,ss) = C(s,ss) - sum((cos(thetaArray)*yp-r)./Norm3)*M;

        end
    end

    % Solve for cover vortex strengths
    cover.Gamma = C\(-cover.vi);
end

%% Compute induced velocity on rotor (rp = [0 r(s) 0])
for s = 1:Ns %for each element
    vi(s) = 0;
    for ii = 1:t %add the velocity induced from each disk
        for ss = 2:Ns+1 %and each ring on each disk (inner ring cancels itself out)

            if (ii == 1) && flagHelix
                ringFrac = 0.675;
            else
                ringFrac = 1;
            end
            
            zr = ring.z(ii,ss);
            r  = ring.r(ii,ss);
            zp = (qh(s)+qh(s+1))/2;
            yp = yE(s);

            M = ring.Gamma(ii,ss)*r*dtheta/(2*pi);
            X2 = (-r*sin(thetaArray)).^2;
            Y2 = (yp-r*cos(thetaArray)).^2;
            Z2 = (zp-zr).^2;
            Normal = sqrt(X2 + Y2 + Z2);
            for iii = 1:length(thetaArray)
                if Normal(iii) < cr
                    Normal(iii) = cr;
                end
            end
            Norm3 = Normal.^3;

            vi(s) = vi(s) + ringFrac*sum((cos(thetaArray)*yp-r)./Norm3)*M;

            % Ground effect ring
            zr = -2*h - ring.z(ii,ss);
            Z2 = (zp-zr).^2;
            Normal = sqrt(X2 + Y2 + Z2);
            for iii = 1:length(thetaArray)
                if Normal(iii) < cr
                    Normal(iii) = cr;
                end
            end
            Norm3 = Normal.^3;

            vi(s) = vi(s) - ringFrac*sum((cos(thetaArray)*yp-r)./Norm3)*M;

        end
    end
end

% Compute induced velocity from cover (rp = [0 r(s) 0])
if flagCover
    for s = 1:Ns %for each element
        for ss = 1:sCover %add the velocity induced from the cover

            if flagHelix
                ringFrac = 0.675;
            else
                ringFrac = 1;
            end
            
            zr = cover.z(ss);
            r  = cover.r(ss);
            zp = (qh(s)+qh(s+1))/2;
            yp = yE(s);

            M = cover.Gamma(ss)*r*dtheta/(2*pi);
            X2 = (-r*sin(thetaArray)).^2;
            Y2 = (yp-r*cos(thetaArray)).^2;
            Z2 = (zp-zr).^2;
            Normal = sqrt(X2 + Y2 + Z2);
            for iii = 1:length(thetaArray)
                if Normal(iii) < cr
                    Normal(iii) = cr;
                end
            end
            Norm3 = Normal.^3;

            vi(s) = vi(s) + ringFrac*sum((cos(thetaArray)*yp-r)./Norm3)*M;

            % Ground effect ring
            zr = -2*h - zr;
            Z2 = (zp-zr).^2;
            Normal = sqrt(X2 + Y2 + Z2);
            for iii = 1:length(thetaArray)
                if Normal(iii) < cr
                    Normal(iii) = cr;
                end
            end
            Norm3 = Normal.^3;

            vi(s) = vi(s) - ringFrac*sum((cos(thetaArray)*yp-r)./Norm3)*M;

        end
    end
end

vi = -vi; %vi is positive downwards

%toc
if flagPlot
    if flagMovie
        plot(ring.r(1:t,:), ring.z(1:t,:), 'k-o', -ring.r(1:t,:), ring.z(1:t,:), 'k-o', [-R R], [0 0], 'k', [-100 100], [-h -h],'b');
        axis([-30 30 -3 15])
        imwrite(im,map,'DancingPeaks.gif','DelayTime',0.1,'LoopCount',inf)
    else
        subplot(1,2,1);
        plot(ring.r(1:t,:), ring.z(1:t,:), 'k-o', -ring.r(1:t,:), ring.z(1:t,:), 'k-o', [-R R], [0 0], 'k', [-100 100], [-h -h],'b');
        axis([-30 30 -10 10])
        text(-25,-3,['Pi = ',num2str(PiApprox,'%0.1f'),' W']);
        text(-25,-5,['t = ',num2str(realtime,'%0.1f'),' s']);
        text(-25,-7,['H = ',num2str(altitude,'%0.1f'),' m']);
        subplot(1,2,2);
        plot(yE,vi,'k',[-100 100],[0 0],'b');
        axis([0 R -1 1])
        pause(0.5);
    end
end
                    

end
    
    
    

