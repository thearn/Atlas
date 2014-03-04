%---------- HeliOpt --------%
%                           %
%       Todd Reichert       %
%        March 5, 2011      %
%                           %
%---------------------------%

% HeliOpt performs an aero-structural optimization of a helicopter
% configuration.

function HeliOpt

close all
clear all

%runFile
runFileValidation

global Ptot Pitot Pptot Ttot Mtot cE cN c100 Cl Cm t xtU xtL xEA d theta nTube nCap lBiscuit yN
global Fblade vi phi Re Cd ring mSpar mChord mQuad EIx EIz EA GJ q EIQuad Finternal strain fail figsz

scrsz = get(0,'ScreenSize');
figsz = round(scrsz*0.9);
figure('Position',[(scrsz(3)-figsz(3))/2 (scrsz(4)-figsz(4))/2 figsz(3) figsz(4)])
%set(gcf,'Position',[0 0 100 100])

if isempty(vrOpt)
    X = [];
    optimValues = [];
    state = 'iter';
    maxDelta = 0;
    
    objectiveFcn(X, Ns, ycmax, rho, visc, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, flagFreeWake, ...
    yWire, zWire, tWire, TWire, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, flagGWing,...
    flagQuad, dQuad, thetaQuad, nTubeQuad, hQuad, etaP, mElse, mPilot, presLoad, flagLoad, vrOpt, maxDelta);
    
    outputFcn(X,optimValues,state,Ns, ycmax, rho, visc, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, flagFreeWake,...
    yWire, zWire, tWire, TWire, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, flagGWing,...
    flagQuad, dQuad, thetaQuad, nTubeQuad, hQuad, etaP, mElse, mPilot, presLoad, flagLoad, vrOpt, maxDelta)

    %keyboard

else
    
    % Set number of cores
    %if multiCore > 1
    %    matlabpool('open',multiCore);
    %end

    % Set linear constraints
    AA = [];
    bb = [];
    Aeq = [];
    beq = [];

    % Set bounds
    ii = 1;
    for i = 1:length(vrOpt)
        for j = 1:length(vrOpt(i).u)
            ub(ii) = vrOpt(i).u(j);
            lb(ii) = vrOpt(i).l(j);
            ii = ii+1;
        end
    end

    % Set initial condition (random)
    for ii = 1:length(ub)
        X0(ii) = rand()*(ub(ii) - lb(ii)) + lb(ii);
    end

    %options = optimset('LargeScale','off','display','off','diagnostics','off','TolFun',1e-6,'TolCon',1e-6);
    if multiCore > 1
        options = optimset('OutputFcn',@outputFcn,'LargeScale','off','display','iter','diagnostics','off','TolFun',1e-4,'TolCon',1e-4,'MaxIter',300,'MaxFunEvals',300*(length(ub)+1),'UseParallel','always'); 
    else
        options = optimset('OutputFcn',@outputFcn,'LargeScale','off','display','iter','diagnostics','off','TolFun',1e-4,'TolCon',1e-4,'MaxIter',300,'MaxFunEvals',300*(length(ub)+1)); 
    end

    [X,P,exitflag] = ...
        fmincon(@objectiveFcn,X0,AA,bb,Aeq,beq,lb,ub,@constraintFcn,options, ...
                Ns, ycmax, rho, visc, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, flagFreeWake, ...
                yWire, zWire, tWire, TWire, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, flagGWing,...
                flagQuad, dQuad, thetaQuad, nTubeQuad, hQuad, etaP, mElse, mPilot, presLoad, flagLoad, vrOpt, maxDelta);

 
    % Set design variables
    ii = 1;
    for i = 1:length(vrOpt)
        eval([vrOpt(i).name,' = X(ii:ii+length(vrOpt(i).u)-1)']);
        ii = ii + length(vrOpt(i).u);
    end
    
    %clean up parallel jobs
    if multiCore > 1
        matlabpool close
    end

end
end


%---------- Objective Function ----------%
function P = ...
    objectiveFcn(X, Ns, ycmax, rho, visc, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, flagFreeWake, ...
    yWire, zWire, tWire, TWire, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, flagGWing,...
    flagQuad, dQuad, thetaQuad, nTubeQuad, hQuad, etaP, mElse, mPilot, presLoad, flagLoad, vrOpt, maxDelta)
    
    global Ptot Pitot Pptot Ttot Mtot cE cN c100 Cl Cm t xtU xtL xEA d theta nTube nCap lBiscuit yN
    global Fblade vi phi Re Cd ring mSpar mChord mQuad EIx EIz EA GJ q EIQuad Finternal strain fail

 
    % Set design variables
    ii = 1;
    for i = 1:length(vrOpt)
        eval([vrOpt(i).name,' = X(ii:ii+length(vrOpt(i).u)-1);']);
        ii = ii + length(vrOpt(i).u);
    end
    
    %         for j = 1:length(vrOpt(i).u)
%             eval([vrOpt(i).name,'(',num2str(j),')',' = X(ii);']);
%             ii = ii+1;
%         end

    % Call HeliCalc
    [Ptot, Pitot, Pptot, Ttot, Mtot, cE, cN, c100, Cl, Cm, t, xtU, xtL, xEA, d, theta, nTube, nCap, lBiscuit, yN...
        Fblade, vi, phi, Re, Cd, ring, mSpar, mChord, mQuad, EIx, EIz, EA, GJ, q, EIQuad, Finternal, strain, fail]...
        = HeliCalc(Ns, ycmax, rho, visc, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, flagFreeWake, ...
            yWire, zWire, tWire, TWire, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, flagGWing,...
            flagQuad, dQuad, thetaQuad, nTubeQuad, hQuad, etaP, mElse, mPilot, presLoad, flagLoad);
       
    P = Ptot;
end

%---------- constraintFcn -----------%
function [Con, ConEq] = constraintFcn(X, Ns, ycmax, rho, visc, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, flagFreeWake, ...
    yWire, zWire, tWire, TWire, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, flagGWing,...
    flagQuad, dQuad, thetaQuad, nTubeQuad, hQuad, etaP, mElse, mPilot, presLoad, flagLoad, vrOpt, maxDelta)
    
    global Ttot Mtot fail q
    
    ConLW = Mtot*9.8 - Ttot;
    
    FOS = 2;
    for s = 1:Ns+1
        i = (s-1)*(36+1+length(yWire)) + 1;
        ConFail(i:i+2) = -1+FOS*abs(fail.top.cap(:,s));
        ConFail(i+3:i+5) = -1+FOS*abs(fail.top.plus(:,s));
        ConFail(i+6:i+8) = -1+FOS*abs(fail.top.minus(:,s));
        ConFail(i+9:i+11) = -1+FOS*abs(fail.bottom.cap(:,s));
        ConFail(i+12:i+14) = -1+FOS*abs(fail.bottom.plus(:,s));
        ConFail(i+15:i+17) = -1+FOS*abs(fail.bottom.minus(:,s));
        ConFail(i+18:i+20) = -1+FOS*abs(fail.front.cap(:,s));
        ConFail(i+21:i+23) = -1+FOS*abs(fail.front.plus(:,s));
        ConFail(i+24:i+26) = -1+FOS*abs(fail.front.minus(:,s));
        ConFail(i+27:i+29) = -1+FOS*abs(fail.back.cap(:,s));
        ConFail(i+30:i+32) = -1+FOS*abs(fail.back.plus(:,s));
        ConFail(i+33:i+35) = -1+FOS*abs(fail.back.minus(:,s));
        ConFail(i+36) = -1+FOS*abs(fail.buckling.torsion(s));
        ConFail(i+37:i+36+length(yWire)) = -1+FOS*abs(fail.buckling.x(:,s));
        ConFail(i+37+length(yWire):i+36+length(yWire)*2) = -1+FOS*abs(fail.buckling.z(:,s));
    end
    i = (s)*(36+1+length(yWire)) + 1;
    ConFail(i) = -1+FOS*abs(fail.buckling.quad);
    
    % Break out deformations
    qq(:,1) = [0 0 0 0 0 0];
    for s = 2:Ns+1
        qq(:,s) = q((s-1)*6+1:(s-1)*6+6);
    end
    ConDelta = qq(3,end) - maxDelta;
    
    Con = [ConLW ConFail ConDelta];
    ConEq = [];
end

%------------ outputFcn --------------%
function stop = ...
    outputFcn(X,optimValues,state,Ns, ycmax, rho, visc, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, flagFreeWake, ...
    yWire, zWire, tWire, TWire, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, flagGWing,...
    flagQuad, dQuad, thetaQuad, nTubeQuad, hQuad, etaP, mElse, mPilot, presLoad, flagLoad, vrOpt, maxDelta)
    
    global Ptot Pitot Pptot Ttot Mtot cE cN c100 Cl Cm t xtU xtL xEA d theta nTube nCap lBiscuit yN
    global Fblade vi phi Re Cd ring mSpar mChord mQuad EIx EIz EA GJ q EIQuad Finternal strain fail figsz


    stop = false; 

    switch state
       case 'init'

       case 'iter'
           % Concatenate current point and objective function
           % value with history. x must be a row vector.
           %history.fval = [history.fval; optimValues.fval];
           %plot(x(1),x(2),'o');
           % Label points with iteration number.
           % Add .15 to x(1) to separate label from plotted 'o'
           %text(x(1)+.15,x(2),num2str(optimValues.iteration));
           
            % Set design variables
            ii = 1;
            for i = 1:length(vrOpt)
                eval([vrOpt(i).name,' = X(ii:ii+length(vrOpt(i).u)-1);']);
                ii = ii + length(vrOpt(i).u);
            end
            
            % Determine the length of each spar element
            Ns = length(yN)-1; %number of elements
            dy = zeros(Ns,1);
            for s=1:Ns
                dy(s) = yN(s+1) - yN(s); %length of each element
            end
            
            yE = zeros(Ns,1);
            for s=1:Ns
                yE(s) = 0.5*(yN(s) + yN(s+1)); %radial location of each element
            end

            % Break out deformations
            qq(:,1) = [0 0 0 0 0 0];
            for s = 2:Ns+1
                qq(:,s) = q((s-1)*6+1:(s-1)*6+6);
            end
            
            % Rotor mass
            if flagQuad
               MRotor = sum(mSpar+mChord)*b*4;
            else
               MRotor = sum(mSpar+mChord)*b;
            end
            MQuad = mQuad*4;
            
            % Plot colors
            k = [0 0 0];
            b = [0 0 1];
            r = [1 0 0];
            g = [0 0.5 0];
            c = [0 0.75 0.75];
            m = [0.75 0 0.75];
            y = [0.75 0.75 0];
            
           hh = subplot(3,4,1);
           set(hh,'position',[0.03 0.71 0.2 0.25]);
           plot([0 80],[0 80],'.');
           title('General Definition');
           text(10,75,['Power = ',num2str(Ptot,'%0.1f'),' W']);
           text(10,70,['Pi = ',num2str(Pitot,'%0.1f'),' W']);
           text(10,65,['Pp = ',num2str(Pptot,'%0.1f'),' W']);
           text(10,55,['Weight = ',num2str(Mtot*9.8,'%0.1f'),' N']);
           text(10,50,['Thrust = ',num2str(Ttot,'%0.1f'),' N']);
           text(10,45,['Radius = ',num2str(R,'%0.2f'),' m']);
           text(10,40,['Omega = ',num2str(Omega/2/pi,'%0.3f'),' Hz']);
           text(10,30,['T wire = ',num2str(TWire,'%0.1f'),' N']);
           
           text(50,75,['M total = ',num2str(Mtot,'%0.1f'),' kg']);
           text(50,70,['M pilot = ',num2str(mPilot,'%0.1f'),' kg']);
           text(50,65,['M HPH = ',num2str(Mtot-mPilot,'%0.1f'),' kg']);
           text(50,60,['M rotors = ',num2str(MRotor,'%0.2f'),' kg']);
           text(50,55,['M quad = ',num2str(MQuad,'%0.2f'),' kg']);
           text(50,50,['M other = ',num2str(mElse,'%0.2f'),' kg']);
           
         
           y100 = linspace(0,R,100);
           Y = [0 ycmax R];
           hh = subplot(3,4,2);
           set(hh,'position',[0.28 0.71 0.2 0.25]);
           set(gca,'colororder',[k; b; r; c],'NextPlot','replacechildren');
           plot(y100,c100,yE,Cl,yE,Cd*100,yE,Re/1000000); %add aerodynamic twist
           title('Aerodynamic Definition');
           xlabel('r (m)')
           legend('Chord (m)','C_l','C_d (10^{-2})','Re (10^6)')
           v = axis;
           axis([0 v(2) 0 3]);
           
           hh = subplot(3,4,3);
           set(hh,'position',[0.53 0.71 0.2 0.25]);
           set(gca,'colororder',[k; g; b; c; m; y; k; g; b; y],'NextPlot','replacechildren');
           plot(yE,d*100/2.54,yE,theta*180/pi/10,yE,nTube,yE,nCap,yE,lBiscuit*100/2.54/12,yWire,zWire,'o',...
                0,dQuad*100/2.54,'s',0,thetaQuad*180/pi/10,'s',0,nTubeQuad,'s',0,hQuad,'s');
           title('Structural Definition');
           xlabel('r (m)')
           legend('Spar diameter (in)','Wrap angle (10^1 deg)','Tube layers','Cap layers','Biscuit spacing (ft)','Wire location (m)','Quad definition');
           v = axis;
           axis([0 v(2) 0 v(4)]);
           
           hh = subplot(3,4,4);
           set(hh,'position',[0.78 0.71 0.2 0.25]);
           set(gca,'colororder',[k; b; r; g; m; m; b],'NextPlot','replacechildren');
           plot(yE,EA/(10^7),yE,EIz/(10^4),yE,EIx/(10^4),yE,GJ/(10^4),yE,mSpar./dy,yE,mChord./dy,'--',...
                0,EIQuad/(10^4),'s');
           title('Structural Properties');
           xlabel('r (m)')
           legend('EA (10^7)','EIz (10^4)','EIx (10^4)','GJ (10^4)','mSpar (kg/m)','mChord (kg/m)', 'EIquad (10^4)');
           
           hh = subplot(3,4,5);
           set(hh,'position',[0.03 0.38 0.2 0.25]);
           set(gca,'colororder',[k; b; r; c; m],'NextPlot','replacechildren');
           plot(yE,Fblade.P./dy,yE,Fblade.Pi./dy,yE,Fblade.Pp./dy,yE,vi*10^2,yE,phi*10*180/pi);
           title('Aerodynamic Power');
           xlabel('r (m)')
           legend('Power (W/m)','P_i (W/m)','P_p (W/m)','v_i (10^{-2} m/s)','\phi  (10^{-1} deg)', 'location','NorthWest')
           v = axis;
           axis([v(1) v(2) 0 v(4)]);
           
           hh = subplot(3,4,6);
           set(hh,'position',[0.28 0.38 0.2 0.25]);
           set(gca,'colororder',[b; r; c; g],'NextPlot','replacechildren');
           plot(yE,Fblade.Fz./dy,yE,Fblade.Fx*10./dy,yE,Fblade.Q./dy,yE,Fblade.My./dy);
           title('Aerodynamic Forces');
           xlabel('r (m)')
           legend('Thrust (N/m)','Drag (10^{-1} N/m)','Torque (Nm/m)','Moment (Nm/m)', 'location','NorthWest')
           
           hh = subplot(3,4,7);
           set(hh,'position',[0.53 0.38 0.2 0.25]);
           set(gca,'colororder',[b; b; r; r; k; g],'NextPlot','replacechildren');
           plot(yN,Finternal(3,:),'--',yN,Finternal(4,:),yN,Finternal(1,:)*10,'--',yN,-Finternal(6,:)*10,yN,Finternal(2,:),yN,Finternal(5,:));
           title('Structural Loads');
           xlabel('r (m)')
           legend('Z (N)','Z'' (Nm)','X (10^{-1} N)','X'' (10^{-1} Nm)','Y (N)','Torque (Nm)','location','southeast');
           
           hh = subplot(3,4,8);
           set(hh,'position',[0.78 0.38 0.2 0.25]);
           set(gca,'colororder',[b; r; g],'NextPlot','replacechildren');
           plot(yN,qq(3,:),yN,qq(1,:),yN,qq(5,:)*180/pi);
           title('Structural Deformation');
           xlabel('r (m)')
           legend('Z Deflection','X Deflection','Twist');
           
           hh = subplot(3,4,9);
           set(hh,'position',[0.03 0.05 0.2 0.25]);
           if flagFreeWake
            plot(ring.r, ring.z, 'k-o', -ring.r, ring.z, 'k-o',[-R R],[-h -h],'k');
           end
           title('Free Wake');
           xlabel('r (m)')
%            set(gca,'colororder',[b; r; k; g],'NextPlot','replacechildren');
%            plot(yN,strain.bending_z,yN,strain.bending_x,yN,strain.axial_y,yN,strain.torsion_y);
%            title('Structural Strain');
%            xlabel('r (m)')
%            legend('Z Bending','X Bending','Axial','Torsion');
           
           hh = subplot(3,4,10);
           set(hh,'position',[0.28 0.05 0.2 0.25]);
           set(gca,'colororder',[k; k; k; b; b; b; r; r; r],'NextPlot','replacechildren');
           plot(yN,fail.top.cap(1,:),'-',yN,fail.top.cap(2,:),'--',yN,fail.top.cap(3,:),'-.', ...
           	    yN,fail.top.plus(1,:),'-',yN,fail.top.plus(2,:),'--',yN,fail.top.plus(3,:),'-.', ...
                yN,fail.top.minus(1,:),'-',yN,fail.top.minus(2,:),'--',yN,fail.top.minus(3,:),'-.');
           title('Out-of-Plane Failure');
           xlabel('r (m)')
           legend('Top_c 11','Top_c 22','Top_c 12',...
                  'Top_+ 11','Top_+ 22','Top_+ 12',...
                  'Top_- 11','Top_- 22','Top_- 12');
              
           hh = subplot(3,4,11);
           set(hh,'position',[0.53 0.05 0.2 0.25]);
           set(gca,'colororder',[b; b; b; r; r; r],'NextPlot','replacechildren');
           plot(yN,fail.back.plus(1,:),'-',yN,fail.back.plus(2,:),'--',yN,fail.back.plus(3,:),'-.', ...
                yN,fail.back.minus(1,:),'-',yN,fail.back.minus(2,:),'--',yN,fail.back.minus(3,:),'-.');
           title('In-Plane Failure');
           xlabel('r (m)')
           legend('Back_+ 11','Back_+ 22','Back_+ 12',...
                  'Back_- 11','Back_- 22','Back_- 12');
              
           hh = subplot(3,4,12);
           set(hh,'position',[0.78 0.05 0.2 0.25]);
           if length(yWire) == 1
                set(gca,'colororder',[b; r; g; k; y],'NextPlot','replacechildren');
                plot(yN,fail.buckling.z,yN,fail.buckling.x,yN,fail.buckling.torsion,...
                    0,fail.buckling.quad,'s',yWire,fail.wire,'o')
                legend('Buckle_Z','Buckle_X','Buckle_{Tor}','Buckle_{Quad}','Wire');
           else
               set(gca,'colororder',[b; b; r; r; g; k; y],'NextPlot','replacechildren');
                plot(yN,fail.buckling.z(1,:),yN,fail.buckling.z(2,:),'--',...
                    yN,fail.buckling.x(1,:),yN,fail.buckling.x(2,:),'--',yN,fail.buckling.torsion,...
                    0,fail.buckling.quad,'s',yWire,fail.wire,'o')
                legend('Buckle_Z 1','Buckle_Z 2','Buckle_X 1','Buckle_X 2','Buckle_{Tor}','Buckle_{Quad}','Wire');
           end
               
           title('Buckling Failure');
           xlabel('r (m)')
           
 
           pause(0.01);
           
       case 'done'
           
       otherwise
   end
end

