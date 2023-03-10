%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Cello Bowed Stiff String Coupled with distributed bridge
%           Modal system - Iterative (Newton Raphson)
%               Riccardo Russo
%             University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear 
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
realTimeDraw = false;

play = true;

%if set true the system is solved with Sherman Morrison formula, otherwise with backslash
smSolver = true;

%sets if to use the improved friction model from desvages
desvagesFriction = false;

stringToPlay = 2;   %0=A3, 1=D3, 2=G2, 3=C2

%sets if to let the string free to vibrate at the end or to stop it
freeVib = true;

saveAudio = false;
if saveAudio
    cond = '_Stop';
    if freeVib
        cond = '_Free';
    end
    switch stringToPlay
    case 0
        string = 'A3';
    case 1
        string = 'D3';
    case 2
        string = 'G2';
    case 3
        string = 'C2';
end
    fileName = strcat('../Sounds/Cello/Notes/',string,cond,'.wav');
end

osFac = 1;          %Oversampling factor
SR = 44100*osFac;   %Sample rate [Hz]
T = 5;              %Time of Simulation [s]

h = 0.5e-2;

omegaMax = 14e3*2*pi;

%%%%% Time Parameters
k = 1/SR;
timeSamples = floor(T*SR);                       %Length of sim [samples]
timeVec = (1:floor(timeSamples))*k;              %Time vector [sec]
timeVecUS = (1:floor(timeSamples/osFac))*k*osFac;%Time vector after undersampling

%%%%% String Length & Fretting Position
baseLength = 0.69;  %String base length [m]
frettingPos = 1;
                %27/32; %3rdm
                %64/81; 3rdM
                %8/9;   %1 semitone
                %4/5;   %4th
                %2/3;   %5th
                %1/2^(nSemitone/12) %semitones

L = baseLength*frettingPos;           

Fb = zeros(1,timeSamples);

%%%%% Cello Strings
switch stringToPlay
    case 0
        % A3           
        radius = 3.75e-04;                  % string radius [m]
        rho = 3.7575e3;                     % string Density [kg/m] 
        T0 = 153;                           % tension [N] 
        E = 25e9;                          % young modulus [Pa]

        excitPos = 0.833*L;
        
        if freeVib
            startFb = 0.01; maxFb = 0.05; endFb = 0;
        else 
            maxFb = 0.05;
        end
    case 1
        % D3           
        radius = 4.4e-04;
        rho = 4.1104e3;                     
        T0 = 102.6;                         
        E = 25e9;                           

        excitPos = 0.833*L;

        if freeVib
            startFb = 0.01; maxFb = 0.05; endFb = 0;
        else 
            maxFb = 0.05;
        end
    case 2
        % G2           
        radius = 6.05e-04;
        rho = 5.3570e3;                          
        T0 = 112.67;                           
        E = 8.6e9;                          

        excitPos = 0.733*L; %G2 C2
        
        if freeVib
            startFb = 0.01; maxFb = 0.05; endFb = 0;
        else 
            maxFb = 0.05;
        end
    case 3
        % C2 
        radius = 7.2e-04;
        rho = 1.3017e4;       
        T0 = 172.74;
        E = 22.4e9;                          

        excitPos = 0.733*L; %G2 C2

        if freeVib
            startFb = 0.01; maxFb = 0.05; endFb = 0;
        else 
            maxFb = 0.05;
        end
end
A = pi*radius^2;
rA = rho*A;
Inertia = (pi*radius^4)/ 4;   
EI = E*Inertia;
K = sqrt(E*Inertia/(rA*L^4));
c = sqrt(T0/rA);

a = 100;            %Bow free parameter
muD = 0.3;          %Desvages friction parameter


%%%%% Bow Speed & Pressure
bowVel = zeros(1,timeSamples);

% %Linear ramp
frac = 3; %number of sections in which the whole bow speed vector will be divided
bowRampLength = floor(2*timeSamples/frac); %Time dedicated to speed variation. After this time the bow velocity will be zero
maxVb = 0.2; startVb = 0.1; endVb = 0.0;
timeFrac = 3; %number of sections in which the velocity ramp will be divided
bowVel(1:floor(1*bowRampLength/timeFrac)) = linspace(startVb,maxVb,floor(1*bowRampLength/timeFrac));
bowVel(floor(1*bowRampLength/timeFrac)+1:floor(2*bowRampLength/timeFrac)) = maxVb;
bowVel(floor(2*bowRampLength/timeFrac)+1:floor(3*bowRampLength/timeFrac))= linspace(maxVb,endVb,bowRampLength/timeFrac);

if freeVib
    Fb(1:floor(1*bowRampLength/5)) = linspace(startFb,maxFb,floor(1*bowRampLength/5));
    Fb(floor(1*bowRampLength/5)+1:floor(2*bowRampLength/5)) = maxFb;
    Fb(floor(2*bowRampLength/5)+1:floor(3*bowRampLength/5))= linspace(maxFb,endFb,bowRampLength/5);
else
    Fb(:) = maxFb;
end

%%%%% Bridge Parameters
LB = 7e-2;                          % string length [m]
radiusB = 5e-4;
rhoB =8e3;                          % string Density [kg/m^3] 
AB = pi*radiusB^2;
rAB = rhoB*AB;
EB = 2.5e15;                          % young modulus [Pa]
InertiaB = (pi*radiusB^4)/ 4;         % moment of inertia
EIB    = EB*InertiaB;

zContact = 0.43*LB;

outPos1 = floor(0.334*LB/h);
outPos2 = outPos1;% floor(0.773*LB/h);


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%% Initializing Discrete Operators
Ms = floor(L/h);
Mb = floor(LB/h);

Jvec = sparse(Mb - 1, 1); 
alphaInterp = zContact/ h - floor(zContact/h);
Jvec(floor(zContact/h)) = (1-alphaInterp)/h;
Jvec(floor(zContact/h) + 1) = alphaInterp/h;
JJT = Jvec*Jvec.';

vs = ones(Ms-1,1);
DxxS = spdiags([vs/h^2 -2*vs/h^2 vs/h^2], -1:1, Ms-1, Ms-1);
DxxxxS = DxxS*DxxS;
vb = ones(Mb-1,1);
DxxB = spdiags([vb/h^2 -2*vb/h^2 vb/h^2], -1:1, Mb-1, Mb-1);
DxxxxB = DxxB * DxxB;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%% Computing eigenfrequencies and modes
zeroVec = sparse(1,Mb - 1);
vec1 = sparse(-Jvec*EI/h^3);
vec2 = sparse(Jvec*(T0/h + 2*EI/h^3));
zeroMat = sparse(Ms-4, Mb-1);
angleMat = [zeroMat; zeroVec; vec1.'; vec2.'];

K = -[(T0*DxxS - EI*DxxxxS), angleMat;
     angleMat.', (-EIB*DxxxxB + (- T0 - EI/h^2)*JJT)];
     
fullK = full(K);
M1 = diag([rA * ones(Ms-1, 1); rAB * ones(Mb - 1, 1)]);
M2 = [zeros(Ms-1), zeros(Ms-1,Mb-1);
     zeros(Ms-1,Mb-1).', rA*h^2*JJT];
M = sparse(M1+M2);
R = full(M\K);
[V, D] = eig(R);

[Dsort, Dind] = sort(diag(D));
Vsort = V(:,Dind);
VsortM = Vsort^-1;

spaceSamplesTot = (Ms-1) + (Mb-1);

eigenFreqs = sqrt(Dsort);
eigenFreqs = eigenFreqs(eigenFreqs<=omegaMax & eigenFreqs>=20*2*pi);
modesNumber = length(eigenFreqs);
modesMat = Vsort(:,1:modesNumber);

inProd = VsortM*Vsort;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing States & Matrices

%Dirac Delta approx. Delta is computed onto the string length. Then delta
%is concatenated with an array of zeros representing the bridge, like what
%was done in case of the initial conditions
excitSample = excitPos/h;
excitPosFloor = floor(excitSample);
betaInterp = excitSample - excitPosFloor;
delta = zeros(Ms-1,1);
delta(excitPosFloor) = (1-betaInterp)/h;
delta(excitPosFloor+1) = betaInterp/h;

deltaTot = sparse([delta; zeros(Mb-1,1)]);

I = speye(modesNumber);

alphaVec = VsortM*M^-1*deltaTot;
betaVecTR = h*deltaTot.'*Vsort;

alphaVec = alphaVec(1:modesNumber);
betaVecTR = betaVecTR(1:modesNumber);

alphaBetaTR = alphaVec*betaVecTR;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing States
xNext = zeros(modesNumber,1); 
xPrev = zeros(modesNumber,1);

%The correct initialization  considers uPrev as u^0, u as u^1, therefore 
%u needs to be initialized taking into account the bow velocity.
x = -0.5*k^2*alphaVec*Fb(1)*sqrt(2*a)*(-bowVel(1))*exp(-a*(-bowVel(1))^2+ 0.5); 

tol = 0; tolThresh = 1e-9;

%%%%% Initializing outputs
Out = zeros(timeSamples,2);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
 fontSize = 14;
 lineWidth = 1.5;
 set(gca,'FontSize',fontSize);

itersAverage = 0;
tic
for i = 1:timeSamples
    coeff1 = Fb(i)*sqrt(2*a);
    %Newton-Raphson solver
    b = (-2 + k^2*eigenFreqs.^2).*x + xPrev*2;
    tol = 1; iters = 1;
    g = x - xPrev; %this because in this point of the code x = xNext
    while tol>tolThresh && iters < 1000
        iters   = iters + 1;
        eta = (0.5/k)*betaVecTR*g - bowVel(i);
        expCoeff = exp(-a*eta^2 + 0.5);
        f = g + b + alphaVec*k^2*coeff1*eta*expCoeff;
        fp = I + 0.5*k*coeff1*(1 - 2*a*eta^2)*expCoeff*alphaBetaTR;
        gNext = g - fp \ f;
        tol = max(abs(g - gNext));
        g = gNext;
    end

    itersAverage = itersAverage + iters;

    xNext = g + xPrev;
    
    xPrev = x;
    x = xNext;

    vNext = modesMat*xNext;
    
    Out(i,:)= [vNext(Ms-1 + outPos1), vNext(Ms-1 + outPos1)];

    if realTimeDraw
        plot(vNext,'LineWidth',lineWidth,'LineWidth',lineWidth,'LineWidth',lineWidth);
        ylabel("u(t,x_o) (m)");
        ylim([-2e-4,2e-4]);
        drawnow
    %     pause(0.05);
        hold off
    end
end
realTimeFrac = toc/T

itersAverage = itersAverage/timeSamples

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Retrieving undersampled output

%Calculating derivative for better sound output
OutDiff1 = diff(Out(:,1));
OutDiff2 = diff(Out(:,2));
OutDiff = [OutDiff1,OutDiff2];

finalOSFac = 1;
if finalOSFac>osFac disp("Undersampling Error."); return; end

OutPlay = zeros(floor(timeSamples/(osFac/finalOSFac)),2);

%lowpassing before downsampling for recording
lowpass(Out(:,1),20000,SR);
lowpass(Out(:,2),20000,SR);
for i=1:size(Out,1)
    if ~mod(i,osFac) || mod(i,osFac) == osFac/finalOSFac
        index = i/(osFac/finalOSFac);
        OutPlay(index,:) = Out(i,:);
    end
end

OutLP = OutPlay;

OutPlay1 = diff(OutPlay(:,1));
OutPlay2 = diff(OutPlay(:,2));
OutPlay = [OutPlay1/max(abs(OutPlay1)),OutPlay2/max(abs(OutPlay2))];

%if play soundsc(OutPlay,SR/osFac*finalOSFac); end
if play soundsc(OutDiff,SR); end
if saveAudio
    audiowrite(fileName,OutPlay,SR/osFac*finalOSFac);
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Plots
fontSize = 18;
lineWidth = 1.5;

figure
plot(timeVec*1000,Out(:,1))
hold on
plot(timeVec*1000,Out(:,2))
ylim([min(min(Out)),max(max(Out))]);
xlabel('Time [s]');
ylabel("u(t,x_o) [m]");

hold on
plot(timeVec*1000,bowVel*1e-4);
hold on
plot(timeVec*1000,Fb*1e-6);
