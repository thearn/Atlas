function testWake()

close all
clear all

test = 4;

b = 2;
R = 10;
Omega = 0.16*2*pi; %rad/s
rho = 1.18;
Cl = 1.2;
ycmax = 2.5;
c_ = [0 1.3 0.2];



switch test
    case 1
        
        yN = [0 2 4 6 8 10]';
        %dT = [0.5 1.5 2.5 3.5 4.5]';
        %dT = [0.05 0.15 0.25 0.35 0.45]';
        dT = [0.1 1 2.2 3.5 4.2]';
        h = 1.5;
        Ns = length(yN)-1;
        q = zeros(1,(Ns+1)*6);

        VortexWake(yN, rho, dT, Omega, b, h, q, 1)
    
    % Validation Study    
    case 2
        
        h = [0.01 0.05 0.1 0.2 0.3 0.4 0.6 0.8 1 1.5 2 10000]*R;
        %h = [0.1 1 100];
        
        % Set fidelity
        Ns = 15;
        Nw = 15;
        Ntt = 5;
        Ntheta = 40;
        [Pi PiBram PiCheese PiHayden] = singlePoint(Ns, Nw, Ntt, Ntheta, h, 0);
        plot(h/R,Pi/Pi(end))
        hold on
        
        Ns = 10;
        Nw = 8;
        Ntt = 3;
        Ntheta = 20;
        [Pi PiBram PiCheese PiHayden] = singlePoint(Ns, Nw, Ntt, Ntheta, h, 0);
        plot(h/R,Pi/Pi(end),'--')

        Ns = 8;
        Nw = 5;
        Ntt = 1;
        Ntheta = 15;
        [Pi PiBram PiCheese PiHayden] = singlePoint(Ns, Nw, Ntt, Ntheta, h, 0);
        plot(h/R,Pi/Pi(end),':',h/R,PiBram, h/R,PiCheese, h/R,PiHayden)

        legend('Vortex HF','Vortex MF','Vortex LF','Bramwell','Cheeseman','Hayden');
        axis([0 2 0 1])
        
        
    % Convergence Study
    case 3
        h = [1.5 3.5 1000000];
        Ns = 20;
        Nw = 20;
        Ntt = 5;
        Ntheta = 15;
        
        NsArr = [5 10 15 20 25 30];
        for ii = 1:length(NsArr)
            [Pi(ii,:) PiBram(ii,:)] = singlePoint(NsArr(ii), Nw, Ntt, Ntheta, h);
        end
        for ii = 1:length(h)
            PinormNs(:,ii) = Pi(:,ii)/Pi(end,ii);
        end
        Pi
        clear Pi
        
        NwArr = [1 5 10 15 20 25 30];
        for ii = 1:length(NwArr)
            [Pi(ii,:) PiBram(ii,:)] = singlePoint(Ns, NwArr(ii), Ntt, Ntheta, h);
        end
        for ii = 1:length(h)
            PinormNw(:,ii) = Pi(:,ii)/Pi(end,ii);
        end
        Pi
        clear Pi
        
        NttArr = [1 2 3 5 10];
        for ii = 1:length(NttArr)
            [Pi(ii,:) PiBram(ii,:)] = singlePoint(Ns, Nw, NttArr(ii), Ntheta, h);
        end
        for ii = 1:length(h)
            PinormNtt(:,ii) = Pi(:,ii)/Pi(end,ii);
        end
        Pi
        clear Pi
        
        NthetaArr = [5 10 15 20 30 40 50];
        for ii = 1:length(NthetaArr)
            [Pi(ii,:) PiBram(ii,:)] = singlePoint(Ns, Nw, Ntt, NthetaArr(ii), h);
        end
        for ii = 1:length(h)
            PinormNtheta(:,ii) = Pi(:,ii)/Pi(end,ii);
        end
        Pi
        clear Pi
        
        subplot(2,2,1)
        plot(NsArr,PinormNs);
        title('Ns')
        legend('1 m','5 m','inf');
        subplot(2,2,2)
        plot(NwArr,PinormNw);
        title('Nw')
        legend('1 m','5 m','inf');
        subplot(2,2,3)
        plot(NttArr,PinormNtt);
        title('Ntt')
        legend('1 m','5 m','inf');
        subplot(2,2,4)
        plot(NthetaArr,PinormNtheta);
        title('Ntheta')
        legend('1 m','5 m','inf');
    
    %vi vs. r study
    case 4
        h = [0.05 0.1 0.25 0.5 1 1000000]*R;
        h=3;
        
        % Set High fidelity
        Ns = 15;
        Nw = 15;
        Ntt = 5;
        Ntheta = 40;
        
        % Set Med fidelity
        Ns = 10;
        Nw = 20;
        Ntt = 10;
        Ntheta = 15;
        
        [Pi PiBram PiCheese PiHayden r vi] = singlePoint(Ns, Nw, Ntt, Ntheta, h, 1);
        for ii = 1:length(h)
            vivi(ii,:) = vi(ii,:)./vi(end,:);
        end
        plot(r,vi)
        xlabel('r/R')
        ylabel('v_i / v_{i \infty}')
        title('v_i / v_{i \infty} in Ground Effect for Actual Chord Distribution')
        legend('h/R = 0.05','h/R = 0.1','h/R = 0.25','h/R = 0.5','h/R = 1')
        
%         Ns = 10;
%         Nw = 8;
%         Ntt = 3;
%         Ntheta = 20;
%         [Pi PiBram PiCheese PiHayden r vi] = singlePoint(Ns, Nw, Ntt, Ntheta, h, 1);
%         for ii = 1:length(h)
%             vivi(ii,:) = vi(ii,:)./vi(end,:);
%         end
%         plot(r,vivi)

%         Ns = 8;
%         Nw = 5;
%         Ntt = 1;
%         Ntheta = 15;
%         [Pi PiBram PiCheese PiHayden r vi] = singlePoint(Ns, Nw, Ntt, Ntheta, h, 1);
%         for ii = 1:length(h)
%             vivi(ii,:) = vi(ii,:)./vi(end,:);
%         end
%         plot(r,vivi)

        %legend('Vortex HF','Vortex MF','Vortex LF','Bramwell','Cheeseman','Hayden');
        %axis([0 2 0 1])
        
        
end

%Subfunction
function [Pi,PiBram,PiCheese,PiHayden, r, vi] = singlePoint(Ns, Nw, Ntt, Ntheta, h, cIdeal)

        yN = zeros(Ns+1,1);
        yE = zeros(Ns,1);
        dy = zeros(Ns,1);
        c = zeros(Ns,1);
        dT = zeros(Ns,1);
        q = zeros(1,(Ns+1)*6);
        
        % Discretize space
        yN = linspace(0,R,Ns+1);
        for s = 1:Ns
            r(s) = (yN(s)+yN(s+1))/2;
            dr(s) = yN(s+1)-yN(s);
        end
        
        % Compute chord length at nodes
        if cIdeal
            for s = 1:Ns
                c(s) = 1/r(s)*c_(3)*R;
            end
            %c(1) = 0.5;
        else
        
            Y = [0 ycmax R];
            for s = 1:Ns+1

                % Check which segment the node is on
                for j = 1:length(Y)
                    if Y(j) > yN(s)
                        break
                    end
                end
                x = (yN(s)-Y(j-1)) / (Y(j)-Y(j-1)); 

                if j == 2
                    pStart = c_(1);
                    pCurve = 0;
                    pEnd = c_(2);
                else
                    pStart = c_(2);
                    pCurve = 0;
                    pEnd = c_(3);
                end
                xx = x*(1-pCurve) + sin(x*pi/2)*pCurve;   
                cN(s)= pStart+(pEnd-pStart)*xx;
                
                if j == 2
                    cN(s) = 0.001;
                end
            end
            
            for s = 1:Ns
                c(s) = cN(s) + cN(s+1);
            end
        end
       
        
        % Compute thrust assuming small angles
        for s = 1:Ns
            dT(s) = 0.5*rho*Omega^2*r(s)^2*Cl*c(s)*dr(s);
        end
        
        for i = 1:length(h)
            [vi(i,:), ring] = VortexWakeCover(yN', rho, dT, Omega, b, h(i), Nw, Ntt, Ntheta, q, ycmax, 0, 1);
            % Compute lift and drag using full angles
            for s = 1:Ns
                U = sqrt( (Omega*r(s))^2 + (vi(i,s))^2 );
                dL = 0.5*rho*U^2*Cl*c(s)*dr(s);
                phi(s) = atan2(vi(i,s),Omega*r(s));
                Fblade.Pi(s) = dL*sin(phi(s))*r(s)*Omega;
            end
            Pi(i) = sum(Fblade.Pi);
        end
        
        % Branwell's ground effect model
        A = pi*R^2;
        solidity = b*mean(c)/(pi*R);
        tc = sum(dT)/(rho*solidity*A*Omega^2*R^2);
        for i = 1:length(h)
            PiBram(i) = (1 - 0.3925/((h(i)/R)^2+0.1714*(h(i)/R)+0.6497) / (1-0.553) * (0.3362*(tc/solidity)+0.8338)/((tc/solidity)+1.492));
        end
        
        % Cheesemen & Benett's ground effect model
        for i = 1:length(h)
            PiCheese(i) = 1/(1+(R/h(i)/4)^2);
        end
        
        % Hayden's ground effect model
        for i = 1:length(h)
            PiHayden(i) = 1/(0.9926 + 0.03794*(2*R/h(i))^2);
        end
        

end

end

