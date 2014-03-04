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

runFile
%runFileValidation
%runFileDriveTube
%runRotorPostFile

global out

flagMovie = flags.Movie;
%set(gcf,'Position',[0 0 100 100])

if isempty(vrOpt)
    X = [];
    optimValues = [];
    maxDelta = 0;
    
    objectiveFcn(X, Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
    yWire, zWire, tWire, TWire, TWireMulti, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective, ...
    dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags,...
    vrOpt, vrOptH, vrCon, vrLink);

    state = 'init';
    outputFcn(X,optimValues,state,Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_,...
    yWire, zWire, tWire, TWire, TWireMulti, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective,...
    dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags,...
    vrOpt, vrOptH, vrCon, vrLink)

    state = 'iter';
    outputFcn(X,optimValues,state,Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_,...
    yWire, zWire, tWire, TWire, TWireMulti, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective,...
    dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags,...
    vrOpt, vrOptH, vrCon, vrLink)

    %keyboard

else
    
    % Set number of cores
    if multiCore > 1
        matlabpool('open',multiCore);
    end

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
    
    if flags.MultiPoint
        for i = 1:length(vrOptH)
            for j = 1:length(vrOptH(i).u)
                ub(ii) = vrOptH(i).u(j);
                lb(ii) = vrOptH(i).l(j);
                ii = ii+1;
            end
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
                Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
                yWire, zWire, tWire, TWire, TWireMulti, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective,...
                dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags,...
                vrOpt, vrOptH, vrCon, vrLink);

 
    
            
    % Output results
    fprintf(1,'Variable          lb       ub       x \n');
    fprintf(1,'------------------------------------- \n');
    ii = 1;
    for i = 1:length(vrOpt)
        for j = 1:length(vrOpt(i).u)
            if (abs(X(ii)-vrOpt(i).u(j)) < options.TolFun) || (abs(X(ii)-vrOpt(i).l(j)) < options.TolFun)
                xx = 'x';
            else
                xx = ' ';
            end
            fprintf(1,'%10s: %s %8.2f %8.2f %8.2f \n',vrOpt(i).name,xx,vrOpt(i).l(j),vrOpt(i).u(j),X(ii));
            ii = ii + 1;
        end
    end
    
    if flags.MultiPoint
        fprintf(1,'------------------------------------- \n');
        for i = 1:length(vrOptH)
            for j = 1:length(vrOptH(i).u)
                if (abs(X(ii)-vrOptH(i).u(j)) < options.TolFun) || (abs(X(ii)-vrOptH(i).l(j)) < options.TolFun)
                    xx = 'x';
                else
                    xx = ' ';
                end
                fprintf(1,'%10s: %s %8.2f %8.2f %8.2f \n',vrOptH(i).name,xx,vrOptH(i).l(j),vrOptH(i).u(j),X(ii));
                ii = ii + 1;
            end
        end
    end
    
    % Set design variables
    ii = 1;
    for i = 1:length(vrOpt)
        eval([vrOpt(i).name,' = X(ii:ii+length(vrOpt(i).u)-1)'])
        ii = ii + length(vrOpt(i).u);
    end
    
    % Save movie
    if flagMovie > 0
        imwrite(im,map,'DancingPeaks.gif','DelayTime',0.1,'LoopCount',inf)
    end
    
    % Clean up parallel jobs
    if multiCore > 1
        matlabpool close
    end
    

end
end


%---------- Objective Function ----------%
function P = ...
    objectiveFcn(X, Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
    yWire, zWire, tWire, TWire, TWireMulti, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective,...
    dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags,...
    vrOpt, vrOptH, vrCon, vrLink)
    
    global out

 
    % Set design variables
    ii = 1;
    for i = 1:length(vrOpt)
        eval([vrOpt(i).name,' = X(ii:ii+length(vrOpt(i).u)-1);']);
        ii = ii + length(vrOpt(i).u);
    end
    
    % Set linked variables
    for i = 1:length(vrLink)
        eval(vrLink(i).equation);
    end

    % Call HeliCalc
    if ~flags.MultiPoint
        out = HeliCalc(Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
            yWire, zWire, tWire, TWire, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective,...
            dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags);
        
        P = out(1).Ptot;
    
    else
        % Set design variables for low altitude
        h = vrCon.Alt(1) + zWire;
        collectiveTemp = 0;
        out = HeliCalc(Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
            yWire, zWire, tWire, TWire, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collectiveTemp,...
            dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags);
        Omega1 = Omega;
        
        % Set design variables and linked variables for high altitude
        for i = 1:length(vrOptH)
            eval([vrOptH(i).name,' = X(ii:ii+length(vrOptH(i).u)-1);']);
            ii = ii + length(vrOptH(i).u);
        end
        for i = 1:length(vrLink)
            eval(vrLink(i).equation);
        end
        h = vrCon.Alt(2) + zWire;
        TWire = TWireMulti(1);
        out(2) = HeliCalc(Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
            yWire, zWire, tWire, TWire, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective,...
            dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags);
        Omega2 = Omega;
        
        % perform wind case
        TWire = TWireMulti(2);
        Cl_ = vrCon.ClMax;
        Omega = nthroot(Omega2^3*vrCon.OmegaRatio,3);
        vw = vrCon.Wind;
        collectiveTemp = 0;
        flags.FreeWake = 0;
        flags.AeroStr = 0;
%         out(3) = HeliCalc(Ns, ycmax, rho, visc, -vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
%            yWire, zWire, tWire, TWire, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collectiveTemp,...
%            dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags);
%         MomRotneg = out(3).MomRot;
        out(3) = HeliCalc(Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
           yWire, zWire, tWire, TWire, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collectiveTemp,...
           dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags);
        %out(3).MomRot = out(3).MomRot - MomRotneg;
        
        
        % perform gravity case
        TWire = TWireMulti(3);
        flags.Load = 1; %gravity and wire forces only
        flags.FreeWake = 0;
        flags.AeroStr = 0;
        out(4) = HeliCalc(Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
           yWire, zWire, tWire, TWire, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective,...
           dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags);

        P = (vrCon.AltRatio*out(1).Ptot + (1-vrCon.AltRatio)*out(2).Ptot);
    end
        
end

%---------- constraintFcn -----------%
function [Con, ConEq] = constraintFcn(X, Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
    yWire, zWire, tWire, TWire, TWireMulti, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective,...
    dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags,...
    vrOpt, vrOptH, vrCon, vrLink)
    
    global out
    
    Con = [];
    ConEq = [];
    
    % Set design variables
    ii = 1;
    for i = 1:length(vrOpt)
        eval([vrOpt(i).name,' = X(ii:ii+length(vrOpt(i).u)-1);']);
        ii = ii + length(vrOpt(i).u);
    end
    
    % Set linked variables
    for i = 1:length(vrLink)
        eval(vrLink(i).equation);
    end
    
    % Lift equals Weight
    ConLW = out(1).Mtot*9.8 - out(1).Ttot;
    Con = [Con ConLW];
    if flags.MultiPoint
        ConLW = out(2).Mtot*9.8 - out(2).Ttot;
        Con = [Con ConLW];
    end
        
    % Structural Failure in Rotor Spar
    if flags.ConFail
        if flags.MultiPoint
            fail = out(3).fail;
        else
            fail = out(1).fail;
        end    
        FOS = vrCon.FOSmat;
        FOSbuck = vrCon.FOSbuck;
        FOStorbuck = vrCon.FOStorbuck;
        FOSquadbuck = vrCon.FOSquadbuck;
        FOSwire = vrCon.FOSwire;
        for s = 1:Ns+1
            i = (s-1)*37 + 1;
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
            ConFail(i+36) = -1+FOStorbuck*abs(fail.buckling.torsion(s));
        end
        Con = [Con ConFail];

        % Buckling failure of spar
        ConFailBuck(1) = -1+FOSbuck*abs(fail.buckling.x(1));
        ConFailBuck(2) = -1+FOSbuck*abs(fail.buckling.z(1));
        Con = [Con ConFailBuck];
        
        % Tensile Failure in Wire
        if flags.MultiPoint
            fail = out(3).fail;
        else
            fail = out(1).fail;
        end   
        ConFailWire = -1+FOSwire*abs(fail.wire);
        Con = [Con ConFailWire];
        
        % Structural Failure in Quadboom (NOT USED)
        fail = out(1).fail;
        ConFailQuad = -1+FOSquadbuck*abs(fail.quad.buckling);
        %Con = [Con ConFailQuad]; 
    end
     
    
     % Rotor wind torque (NOT USED)
     if flags.MultiPoint
        ConTorque = out(3).MomRot - 1500;
        %Con = [Con ConTorque];
     end
     

    % Break out deformation
    yN = linspace(0,R,Ns+1);
    qq1(:,1) = [0 0 0 0 0 0];
    for s = 2:Ns+1
        qq1(:,s) = out(1).q((s-1)*6+1:(s-1)*6+6);
    end
    qh1 = qq1(3,:)-yN*anhedral;
    if flags.MultiPoint
        qq2(:,1) = [0 0 0 0 0 0];
        qq3(:,1) = [0 0 0 0 0 0];
        qq4(:,1) = [0 0 0 0 0 0];
        for s = 2:Ns+1
            qq2(:,s) = out(2).q((s-1)*6+1:(s-1)*6+6);
            qq3(:,s) = out(3).q((s-1)*6+1:(s-1)*6+6);
            qq4(:,s) = out(4).q((s-1)*6+1:(s-1)*6+6);
        end
        qh2 = qq2(3,:)-yN*anhedral;
        qh3 = qq3(3,:)-yN*anhedral;
        qh4 = qq4(3,:)-yN*anhedral;
    end
    
    % Constraints on Maximum Deformation
    if flags.ConDef
        ConDelta = qh1(end) - vrCon.MaxDelta;
        Con = [Con ConDelta];
    end

    % Consitent jig twist
    if flags.MultiPoint && flags.ConJigCont
        for s = 1:Ns+1
            if yN(s) < ycmax(1)
                sTrans(1) = s;
            end
            if yN(s) < ycmax(2)
                sTrans(2) = s;
            end
%             if yN(s) < ycmax(3)
%                 sTrans(3) = s;
%             end
        end
        ConAlphaJig(1) = out(2).alphaJig(sTrans(1)+1) - out(1).alphaJig(sTrans(1)+1);
        
        %ConAlphaJig(2) = out(2).alphaJig(sTrans(2)) - out(1).alphaJig(sTrans(2));
%         ConAlphaJig(3) = out(1).alphaJig(sTrans(3)) - out(3).alphaJig(sTrans(3));
        
%         ConAlphaJig(1) = out(1).alphaJig(sTrans(1)+1) - out(2).alphaJig(sTrans(1)+1);
        ConAlphaJig(2) = out(2).alphaJig(Ns-1) - out(1).alphaJig(Ns-1);
        
        
        ConEq = [ConEq ConAlphaJig];
    end
    
    %Ground clearance on gravity load
    if flags.MultiPoint
        minDelta = -zWire + vrCon.MinDelta;
        ConClearance = minDelta - qh4(end);
        %Con = [Con ConClearance];
    end

    % Wire stretch consistency
    if flags.MultiPoint && flags.ConWireCont
        for s = 1:Ns
            if out(1).yN(s) < yWire(1)
                y = s;
            end
        end
        x = (yWire(1)-out(1).yN(y))/(out(1).yN(y+1)-out(1).yN(y));
        hh1 = qh1(y)*(1-x) + qh1(y+1)*x;
        hh2 = qh2(y)*(1-x) + qh2(y+1)*x;
        hh3 = qh3(y)*(1-x) + qh3(y+1)*x;
        hh4 = qh4(y)*(1-x) + qh4(y+1)*x;
        [RHOWire, EWire, ULTIMATEWire] = WireProperties(flags.WireType);
        EAWire = pi*(tWire/2)^2*EWire;
        L0 = sqrt(yWire(1)^2 + zWire^2);
        L1 = L0 + TWire*L0/EAWire;
        L2 = L0 + TWireMulti(1)*L0/EAWire;
        L3 = L0 + TWireMulti(2)*L0/EAWire;
        L4 = L0 + TWireMulti(3)*L0/EAWire;
        delh1 = sqrt(L1^2-yWire(1)^2) - zWire;
        delh2 = sqrt(L2^2-yWire(1)^2) - zWire;
        delh3 = sqrt(L3^2-yWire(1)^2) - zWire;
        delh4 = sqrt(L4^2-yWire(1)^2) - zWire;
        conWire(1) = (hh1-delh1) - (hh2-delh2);
        conWire(2) = (hh1-delh1) - (hh3-delh3);
        conWire(3) = (hh1-delh1) - (hh4-delh4);
        ConEq = [ConEq conWire];
    end

end

%------------ outputFcn --------------%
function stop = ...
    outputFcn(X,optimValues,state,Ns, ycmax, rho, visc, vw, vc, R, b, h, Omega, c_, Cl_, Cm_, t_, xtU_, xtL_, ...
    yWire, zWire, tWire, TWire, TWireMulti, TEtension, xEA_, d_, theta_, nTube_, nCap_, lBiscuit_, anhedral, collective,...
    dQuad, thetaQuad, nTubeQuad, lBiscuitQuad, hQuad, etaP, mElseRotor, mElseCentre, mElseR, mPilot, presLoad, flags,...
    vrOpt, vrOptH, vrCon, vrLink)
    
    global out flagMovie map im hfig


    stop = false; 

    switch state
       case 'init'
           scrsz = get(0,'ScreenSize');
           figsz = round(scrsz*0.9);
           
           if ~flags.MultiPoint
               hfig = figure('Position',[(scrsz(3)-figsz(3))/2 (scrsz(4)-figsz(4))/2 figsz(3) figsz(4)]);
               set(gcf,'color',[1 1 1]);
           else
               for j = 1:4
                    hfig(j) = figure('Position',[(scrsz(3)-figsz(3))/2 (scrsz(4)-figsz(4))/2 figsz(3) figsz(4)]);
                    set(gcf,'color',[1 1 1]);
               end
           end

       case 'iter'
           % Concatenate current point and objective function
           % value with history. x must be a row vector.
           %history.fval = [history.fval; optimValues.fval];
           %plot(x(1),x(2),'o');
           % Label points with iteration number.
           % Add .15 to x(1) to separate label from plotted 'o'
           %text(x(1)+.15,x(2),num2str(optimValues.iteration));
           for j = 1:length(out)
           
               Ptot = out(j).Ptot;
               Pitot = out(j).Pitot;
               Pptot = out(j).Pptot;
               Ttot = out(j).Ttot;
               Mtot = out(j).Mtot;
               MomRot = out(j).MomRot;
               cE = out(j).cE;
               cN = out(j).cN;
               c100 = out(j).c100;
               Cl = out(j).Cl;
               Cm = out(j).Cm;
               t = out(j).t;
               xtU = out(j).xtU;
               xtL = out(j).xtL;
               alphaJig = out(j).alphaJig;
               xEA = out(j).xEA;
               d = out(j).d;
               theta = out(j).theta;
               nTube = out(j).nTube;
               nCap = out(j).nCap;
               lBiscuit = out(j).lBiscuit;
               yN = out(j).yN;
               Fblade = out(j).Fblade;
               vi = out(j).vi;
               phi = out(j).phi;
               Re = out(j).Re;
               Cd = out(j).Cd;
               ring = out(j).ring;
               mSpar = out(j).mSpar;
               mChord = out(j).mChord;
               mQuad = out(j).mQuad;
               mCover = out(j).mCover;
               mWire = out(j).mWire;
               EIx = out(j).EIx;
               EIz = out(j).EIz;
               EA = out(j).EA;
               GJ = out(j).GJ;
               q = out(j).q;
               EIQuad = out(j).EIQuad;
               GJQuad = out(j).GJQuad;
               Finternal = out(j).Finternal;
               strain = out(j).strain;
               fail = out(j).fail;

                % Set design variables
                ii = 1;
                for i = 1:length(vrOpt)
                    eval([vrOpt(i).name,' = X(ii:ii+length(vrOpt(i).u)-1);']);
                    ii = ii + length(vrOpt(i).u);
                end
                
                % Set linked variables
                for i = 1:length(vrLink)
                    eval(vrLink(i).equation);
                end
                
                % Set multipoint variables
                if j > 1
                    for i = 1:length(vrOptH)
                        eval([vrOptH(i).name,' = X(ii:ii+length(vrOptH(i).u)-1);']);
                        ii = ii + length(vrOptH(i).u);
                    end
                end
                
                % Set worst case load variables
                if j == 3
                    Omega = nthroot(Omega^3*vrCon.OmegaRatio,3);
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
                qh = qq(3,:)-yN'*anhedral;

                % Rotor mass
                if flags.Quad
                   MRotor = sum(mSpar+mChord)*b*4;
                   MSpar = sum(mSpar)*b*4;
                   MChord = sum(mChord)*b*4;
                   MCover = mCover*4;
                   MWire = mWire*b*4;
                else
                   MRotor = sum(mSpar+mChord)*b;
                   MSpar = sum(mSpar)*b;
                   MChord = sum(mChord)*b;
                   MCover = mCover;
                   MWire = mWire*b;
                end
                MQuad = mQuad*4;

                % Plot colors
                k = [0 0 0];
                bb = [0 0 1];
                r = [1 0 0];
                g = [0 0.5 0];
                c = [0 0.75 0.75];
                m = [0.75 0 0.75];
                y = [0.75 0.75 0];
                gr = [0.5 0.5 0.5];

               set(0,'CurrentFigure',hfig(j))
                
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
               text(10,35,['Tip Vel = ',num2str(Omega*R,'%0.3f'),' m/s']);

               text(10,25,['T wire_1 = ',num2str(TWire,'%0.1f'),' N']);
               text(10,20,['T wire_2 = ',num2str(TWireMulti(1),'%0.1f'),' N']);
               text(10,15,['T wire_{wind} = ',num2str(TWireMulti(2),'%0.1f'),' N']);
               text(10,10,['T wire_{grav} = ',num2str(TWireMulti(3),'%0.1f'),' N']);
               text(10,5,['anhedral = ',num2str(anhedral*180/pi,'%0.1f'),' deg']);
               
               text(50,75,['M total = ',num2str(Mtot,'%0.1f'),' kg']);
               text(50,70,['M pilot = ',num2str(mPilot,'%0.1f'),' kg']);
               text(50,65,['M HPH = ',num2str(Mtot-mPilot,'%0.1f'),' kg']);
               text(50,60,['M rotors = ',num2str(MRotor,'%0.2f'),' kg']);
               text(50,50,['M chord = ',num2str(MChord,'%0.2f'),' kg']);
               text(50,45,['M spars = ',num2str(MSpar,'%0.2f'),' kg']);
               text(50,40,['M quad = ',num2str(MQuad,'%0.2f'),' kg']);
               text(50,35,['M cover = ',num2str(MCover,'%0.2f'),' kg']);
               text(50,30,['M wire = ',num2str(MWire,'%0.2f'),' kg']);
               
               text(50,25,['M else rotor = ',num2str(mElseRotor,'%0.2f'),' kg']);
               text(50,20,['M else centre = ',num2str(mElseCentre,'%0.2f'),' kg']);
               text(50,15,['M else R = ',num2str(mElseR*R,'%0.2f'),' kg']);
               
               text(50,5,['Rotor Moment = ',num2str(MomRot,'%0.2f'),' Nm']);


               y100 = linspace(0,R,100);
               hh = subplot(3,4,2);
               set(hh,'position',[0.28 0.71 0.2 0.25]);
               set(gca,'colororder',[k; bb; r; g; c; gr],'NextPlot','replacechildren');
               plot(yE,cE,yE,Cl,yE,Cd*100,yE,Cl./Cd/100,yE,Re/1000000,y100,c100); %add aerodynamic twist
               title('Aerodynamic Definition');
               xlabel('r (m)')
               legend('Chord (m)','C_l','C_d (10^{-2})','L/D (10^{2})','Re (10^6)')
               v = axis;
               axis([0 v(2) 0 3]);

               hh = subplot(3,4,3);
               set(hh,'position',[0.53 0.71 0.2 0.25]);
               set(gca,'colororder',[k; g; bb; c; m; y; k; g; bb; m; y],'NextPlot','replacechildren');
               plot(yE,d*100/2.54,yE,theta*180/pi/10,yE,nTube,yE,nCap,yE,lBiscuit*100/2.54/12,yWire,zWire,'o',...
                    0,dQuad*100/2.54,'s',0,thetaQuad*180/pi/10,'s',0,nTubeQuad,'s',0,lBiscuitQuad,'s',0,hQuad,'s');
               title('Structural Definition');
               xlabel('r (m)')
               legend('Spar diameter (in)','Wrap angle (10^1 deg)','Tube layers','Cap layers','Biscuit spacing (ft)','Wire location (m)','Quad definition');
               v = axis;
               axis([0 v(2) 0 v(4)]);

               hh = subplot(3,4,4);
               set(hh,'position',[0.78 0.71 0.2 0.25]);
               set(gca,'colororder',[k; bb; r; g; m; m; bb; g],'NextPlot','replacechildren');
               plot(yE,EA/(10^7),yE,EIz/(10^4),yE,EIx/(10^4),yE,GJ/(10^4),yE,mSpar./dy,yE,mChord./dy,'--',...
                    0,EIQuad/(10^4),'s',0,GJQuad/(10^4),'s');
               title('Structural Properties');
               xlabel('r (m)')
               legend('EA (10^7)','EIz (10^4)','EIx (10^4)','GJ (10^4)','mSpar (kg/m)','mChord (kg/m)', 'EIquad (10^4)', 'GJquad (10^4)');

               hh = subplot(3,4,5);
               set(hh,'position',[0.03 0.38 0.2 0.25]);
               set(gca,'colororder',[c; m; r],'NextPlot','replacechildren');
               plot(yE,vi*10^1,yE,phi*180/pi,yE,alphaJig*180/pi,[0 R],[0 0],'k');
               title('Aerodynamic Velocities and Angles');
               xlabel('r (m)')
               legend('v_i (10^{-1} m/s)','\phi  (deg)', '\alpha_{jig} (deg)', 'location','SouthEast')
               v = axis;
               axis([v(1) v(2) -20 v(4)]);

               hh = subplot(3,4,6);
               set(hh,'position',[0.28 0.38 0.2 0.25]);
               set(gca,'colororder',[bb; r; m; g; k; k; k],'NextPlot','replacechildren');
               plot(yE,Fblade.Fz./dy,yE,Fblade.Fx*10./dy,yE,Fblade.Q./dy,yE,Fblade.My./dy,yE,Fblade.P./dy,yE,Fblade.Pi./dy,'-.',yE,Fblade.Pp./dy,'--');
               title('Aerodynamic Forces and Power');
               xlabel('r (m)')
               legend('Thrust (N/m)','Drag (10^{-1} N/m)','Torque (Nm/m)','Moment (Nm/m)','Power (W/m)','P_i (W/m)','P_p (W/m)', 'location','NorthWest')

               hh = subplot(3,4,7);
               set(hh,'position',[0.53 0.38 0.2 0.25]);
               set(gca,'colororder',[bb; bb; r; r; k; g],'NextPlot','replacechildren');
               plot(yN,Finternal(3,:),'--',yN,Finternal(4,:),yN,Finternal(1,:)*10,'--',yN,-Finternal(6,:)*10,yN,Finternal(2,:),yN,Finternal(5,:));
               title('Structural Loads');
               xlabel('r (m)')
               legend('Z (N)','Z'' (Nm)','X (10^{-1} N)','X'' (10^{-1} Nm)','Y (N)','Torque (Nm)','location','southeast');

               hh = subplot(3,4,8);
               set(hh,'position',[0.78 0.38 0.2 0.25]);
               set(gca,'colororder',[bb; r; g],'NextPlot','replacechildren');
               plot(yN,qh,yN,qq(1,:),yN,qq(5,:)*180/pi);
               title('Structural Deformation');
               xlabel('r (m)')
               legend('Z Deflection','X Deflection','Twist');

               hh = subplot(3,4,9);
               set(hh,'position',[0.03 0.05 0.2 0.25]);
               if (flags.FreeWake) && ((j==1) || (j==2))
                %plot(ring.r, ring.z, 'k-o', -ring.r, ring.z, 'k-o',[-R R],[-h -h],'k');
                plot(ring.r, ring.z, 'k-o', -ring.r, ring.z, 'k-o');
               end
               title('Free Wake');
               xlabel('r (m)')
    %            set(gca,'colororder',[bb; r; k; g],'NextPlot','replacechildren');
    %            plot(yN,strain.bending_z,yN,strain.bending_x,yN,strain.axial_y,yN,strain.torsion_y);
    %            title('Structural Strain');
    %            xlabel('r (m)')
    %            legend('Z Bending','X Bending','Axial','Torsion');

               hh = subplot(3,4,10);
               set(hh,'position',[0.28 0.05 0.2 0.25]);
               set(gca,'colororder',[k; k; k; bb; bb; bb; r; r; r],'NextPlot','replacechildren');
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
               set(gca,'colororder',[bb; bb; bb; r; r; r],'NextPlot','replacechildren');
               plot(yN,fail.back.plus(1,:),'-',yN,fail.back.plus(2,:),'--',yN,fail.back.plus(3,:),'-.', ...
                    yN,fail.back.minus(1,:),'-',yN,fail.back.minus(2,:),'--',yN,fail.back.minus(3,:),'-.');
               title('In-Plane Failure');
               xlabel('r (m)')
               legend('Back_+ 11','Back_+ 22','Back_+ 12',...
                      'Back_- 11','Back_- 22','Back_- 12');

               hh = subplot(3,4,12);
               set(hh,'position',[0.78 0.05 0.2 0.25]);
               if length(yWire) == 1
%                     set(gca,'colororder',[bb; r; g; k; bb; g; g; y],'NextPlot','replacechildren');
%                     plot(yN,fail.buckling.z,yN,fail.buckling.x,yN,fail.buckling.torsion,...
%                         0,fail.quad.buckling,'s',0,fail.quad.bend,'s',0,fail.quad.torsion,'s',0,fail.quad.torbuck,'d',yWire,fail.wire,'o')
%                     legend('Buckle_Z','Buckle_X','Buckle_{Tor}','Quad_{Buckle}','Quad_{Bend}','Quad_{Tor}','Quad_{TorBuck}','Wire');
                    set(gca,'colororder',[bb; r; g; k; bb; g; g; y],'NextPlot','replacechildren');
                    plot(yN,fail.buckling.z,yN,fail.buckling.x,yN,fail.buckling.torsion,...
                        0,fail.quad.buckling,'s',yWire,fail.wire,'o')
                    legend('Buckle_Z','Buckle_X','Buckle_{Tor}','Quad_{Buckle}','Wire');
               else
                   set(gca,'colororder',[bb; bb; r; r; g; k; y],'NextPlot','replacechildren');
                    plot(yN,fail.buckling.z(1,:),yN,fail.buckling.z(2,:),'--',...
                        yN,fail.buckling.x(1,:),yN,fail.buckling.x(2,:),'--',yN,fail.buckling.torsion,...
                        0,fail.buckling.quad,'s',yWire,fail.wire,'o')
                    legend('Buckle_Z 1','Buckle_Z 2','Buckle_X 1','Buckle_X 2','Buckle_{Tor}','Buckle_{Quad}','Wire');
               end

               title('Buckling Failure');
               xlabel('r (m)')


               pause(0.01);
               
           end
           
           if flagMovie == 1
                f = getframe(gcf);
                [im(:,:,1,flagMovie),map] = rgb2ind(f.cdata,1000,'nodither');
                flagMovie = flagMovie + 1;
           elseif flagMovie > 1
                f = getframe(gcf);
                im(:,:,1,flagMovie) = rgb2ind(f.cdata,map,'nodither');
                flagMovie = flagMovie + 1;
           end
           
       case 'done'
           
       otherwise
   end
end

