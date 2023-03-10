%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%     Cello Bowed Stiff String 
%    Modal system - Non Iterative
%         Riccardo Russo
%       University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear 
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
play = true;

%if set true the system is solved with Sherman Morrison formula, otherwise with backslash
smSolver = true;

%sets if to use the improved friction model from desvages
desvagesFriction = false;

stringToPlay = 0;   %0=A3, 1=D3, 2=G2, 3=C2

%sets if to let the string free to vibrate at the end or to stop it
freeVib = false;

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

osFac = 2;          %Oversampling factor
SR = 44100*osFac;   %Sample rate [Hz]
T = 5;              %Time of Simulation [s]

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
        A = pi*radius^2;
        rA = rho*A;
        E = 25e9;                          % young modulus [Pa]
        Inertia = (pi*radius^4)/ 4;        % moment of inertia
        EI = E*Inertia;
        K = sqrt(E*Inertia/(rA*L^4));
        c = sqrt(T0/rA);

        excitPos = 0.533;
        %double output for stereo sound
        outPos1 = 0.33*L;
        outPos2 = 0.77*L;
        
        if freeVib
            startFb = 20; maxFb = 20; endFb = 0;
        else 
            maxFb = 1;
        end
    case 1
        % D3           
        radius = 4.4e-04;
        rho = 4.1104e3;                     
        T0 = 102.6;                         
        A = pi*radius^2;
        rA = rho*A;
        E = 25e9;                           
        Inertia = (pi*radius^4)/ 4;    
        EI = E*Inertia;
        K = sqrt(E*Inertia/(rA*L^4));
        c = sqrt(T0/rA);

        excitPos = 0.533;
        %double output for stereo sound
        outPos1 = 0.33*L;
        outPos2 = 0.77*L;

        if freeVib
            startFb = 30; maxFb = 30; endFb = 0;
        else 
            maxFb = 20;
        end
    case 2
        % G2           
        radius = 6.05e-04;
        rho = 5.3570e3;                          
        T0 = 112.67;                           
        A = pi*radius^2;
        rA = rho*A;
        E = 8.6e9;                          
        Inertia = (pi*radius^4)/ 4;   
        EI = E*Inertia;
        K = sqrt(E*Inertia/(rA*L^4));
        c = sqrt(T0/rA);

        excitPos = 0.533*L; %G2 C2
        %double output for stereo sound
        outPos1 = 0.53*L;
        outPos2 = 0.77*L;
        
        if freeVib
            startFb = 10; maxFb = 10; endFb = 0;
        else 
            maxFb = 10;
        end
    case 3
        % C2 
        radius = 7.2e-04;
        rho = 1.3017e4;                          
        T0 = 172.74;                           
        A = pi*radius^2;
        rA = rho*A;
        E = 22.4e9;                          
        Inertia = (pi*radius^4)/ 4;   
        EI = E*Inertia;
        K = sqrt(E*Inertia/(rA*L^4));
        c = sqrt(T0/rA);

        excitPos = 0.533*L; %G2 C2
        %double output for stereo sound
        outPos1 = 0.33*L;
        outPos2 = 0.57*L;

        if freeVib
            startFb = 10; maxFb = 10; endFb = 0;
        else 
            maxFb = 10;
        end
end

a = 100;            %Bow free parameter
muD = 0.3;          %Desvages friction parameter

%%%%% Bridge Parameters
LB = 7e-2;                          % string length [m]
radiusB = 1e-3;
rhoB = 8e3;                          % string Density [kg/m^3] 
AB = pi*radiusB^2;
rAB = rhoB*AB;
EB = 3e10;                          % young modulus [Pa]
InertiaB = (pi*radiusB^4)/ 4;         % moment of inertia
EIB    = EB*InertiaB;

zContact = 0.13*LB;

%%%%% Bow Speed & Pressure
bowVel = zeros(1,timeSamples);

% %Linear ramp
bowRampLength = 2*timeSamples/T;
maxVb = 0.2; startVb = 0.1; endVb = 0.0;
timeFrac = 2;
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

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Computing eigenfrequencies and modes

%first guess omega vector
omegaMax = 14e3 * 2 * pi;
omegaFirst = (20 * 2 * pi : 1 : omegaMax)';

%calculating number of modes in isolation (can't be much much greater in the
%coupled case)
nModesFirst = 0;
eigenFreq = 0;
while eigenFreq < omegaMax
    nModesFirst = nModesFirst + 1;
    eigenFreq = ComputeEigenFreqIsolBridge(EB,InertiaB,LB,rAB,nModesFirst);
end

contactModesBridge = ComputeBridgeMode(LB, zContact, (1:nModesFirst));
stiffMatDiag = (-pi*(1:nModesFirst)/LB).^4';
stiffMat = diag(stiffMatDiag);
idVec = ones(nModesFirst,1);
id = eye(nModesFirst);
ce1 = zeros(size(omegaFirst,1), 1);

rP    = (T0 + sqrt(T0^2 + 4*EI*rA*omegaFirst.^2))/2/EI;
rM    = (sqrt(T0^2 + 4*EI*rA*omegaFirst.^2) - T0)/2/EI;

lambdaM   = sqrt(rM);
lambdaP   = sqrt(rP);

th    = tanh(lambdaP*L);
cs    = cos(lambdaM*L);
sn    = sin(lambdaM*L);

%characteristicEq
for i = 1:size(omegaFirst)
    omegaTerm = (-rAB*omegaFirst(i)^2*idVec + EIB*stiffMatDiag).^(-1);
    ce1(i) = (contactModesBridge.*omegaTerm)'*contactModesBridge;
end
ce0 = - sn.*(1 + (lambdaM.^2./lambdaP.^2));
ce2 = -T0*(lambdaM.*cs + (lambdaM.^2./lambdaP).*(sn./th));
ce3 = + EI*(-lambdaM.^3.*cs + (lambdaP.*lambdaM.^2).*(sn./th));
%charactEq = (1./ce1).*ce0 + (ce2 + ce3);
charactEq = ce0 + ce1.*(ce2 + ce3);

figure(1)
plot(omegaFirst,charactEq)
ylim([-10,10])
hold on
grid on;

%max of characteristic equation (for plotting purposes)
mm = max(abs(charactEq)) ;

%Search for zeros (first approximate search for NR 1st guesses)
zom = [] ;
for i = 1 : length(omegaFirst) - 1
    %It checks when there is a change in sign
    tt = charactEq(i) * charactEq(i+1) ;
    if tt <= 0
        if abs(charactEq(i))<1
            %check for singularities
            zom = [zom;omegaFirst(i)] ; 
        end
    end  
end

%-- initialise zero vector
eigenFreqs = zeros(length(zom),1) ;

%-- newton Raphson
for i = 1 : length(zom) 
    %first guess
    omega = zom(i) ;
    tol = 1;
    iter = 0;
    
    while tol > 1e-13 && iter < 20
        iter = iter + 1 ;
        rP = (T0 + sqrt(T0^2 + 4*EI*rA*omega.^2))/2/EI ;
        rM = (sqrt(T0^2 + 4*EI*rA*omega.^2) - T0)/2/EI ;
        lambdaM = sqrt(rM) ;
        lambdaP = sqrt(rP) ;
        th = tanh(lambdaP*L) ;
        cs = cos(lambdaM*L) ;
        sn = sin(lambdaM*L) ;
        
        ce0 = - sn.*(1 + (lambdaM.^2./lambdaP.^2));
        massStiffMatBridge = -omega^2*rAB*idVec + EIB*stiffMatDiag;
        ce1 = (contactModesBridge.*(massStiffMatBridge.^(-1)))'*contactModesBridge;
        ce2 = -T0*(lambdaM.*cs + (lambdaM.^2./lambdaP).*(sn./th));
        ce3 = + EI*(-lambdaM.^3.*cs + (lambdaP.*lambdaM.^2).*(sn./th));
        charactEq = ce0 + ce1 * (ce2 + ce3);
        
        dce0 = (L*omega*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/(T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2) - 1))/((-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) - sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((2*EI*omega*rA)/((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) + (2*EI*omega*rA*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)^2*(T0^2 + 4*EI*rA*omega^2)^(1/2)));
        dMassStiffMatBridge = -2*omega*rAB*idVec;
        dce1 = (-contactModesBridge.*(massStiffMatBridge.^(-1).*dMassStiffMatBridge.*massStiffMatBridge.^(-1)))'*contactModesBridge;
        dce2 = -T0*((omega*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)))/((-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) - (L*omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)))/(T0^2 + 4*EI*rA*omega^2)^(1/2) + (2*omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)))/(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) + (omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(3/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) - (L*omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))^2 - 1)*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))^2*(T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) - (L*omega*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)));
        dce3 = -EI*((3*omega*rA*cos((-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))/(T0^2 + 4*EI*rA*omega^2)^(1/2) + (omega*rA*sin((-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(EI*(T0^2 + 4*EI*rA*omega^2)^(1/2)) - (2*omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))/(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0^2 + 4*EI*rA*omega^2)^(1/2)) + (omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) + (L*omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))^2 - 1)*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))^2*(T0^2 + 4*EI*rA*omega^2)^(1/2)) + (L*omega*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))/(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)));
 
        dCharactEq = dce0 + dce1*(ce2 + ce3) + ce1*(dce2 + dce3);

        omegaNew = omega - charactEq/dCharactEq;
        
        tol = abs(omegaNew-omega);
        omega = omegaNew;  
    end
    
    if isnan(omega)
        disp("Newton Raphson Failed")
        eigenFreqs(i) = zom(i);
    else
        eigenFreqs(i) = omega;
    end
    
    line([omega,omega],[-mm,mm],'linestyle','--','color','k') ;
    
end
modesNumber = length(eigenFreqs);

rPlus    = (T0 + sqrt(T0^2 + 4*EI*rA*eigenFreqs.^2))/2/EI;
rMinus    = (sqrt(T0^2 + 4*EI*rA*eigenFreqs.^2) - T0)/2/EI;

lambdaM   = sqrt(rMinus);
lambdaP   = sqrt(rPlus);

sinhL    = sinh(lambdaP*L);
sinL    = sin(lambdaM*L);

const = (lambdaM.^2./lambdaP.^2).*(sinL./sinhL);

modesOut1 = (sin(lambdaM*outPos1) + const.*sinh(lambdaP*outPos1));
modesOut2 = (sin(lambdaM*outPos2) + const.*sinh(lambdaP*outPos2));
modesIn = (sin(lambdaM*excitPos) + const.*sinh(lambdaP*excitPos));

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%Calculating integrals for modal expansion (trapezoidal rule)
h = 1e-5 * L;
x = (0 : h : L).';

inProdInt = zeros(modesNumber);
initCondInt = zeros(modesNumber,1);

for m = 1 : modesNumber
    omegam = eigenFreqs(m);

    rPm   = (T0 + sqrt(T0^2 + 4*EI*rA*omegam^2))/2/EI;
    rMm   = (sqrt(T0^2 + 4*EI*rA*omegam^2) - T0)/2/EI;
    lambdaMm   = sqrt(rMm);
    lambdaPm   = sqrt(rPm);

    sinhLm  = sinh(lambdaPm*L);
    sinLm  = sin(lambdaMm*L);
    constm    = (lambdaMm^2/lambdaPm^2)*(sinLm/sinhLm);

    psim    = sin(lambdaMm*x) + constm*sinh(lambdaPm*x) ;
    psim(1) = 0.5*psim(1) ; psim(end) = 0.5*psim(end) ;

    for n = m : modesNumber
        omegan = eigenFreqs(n) ;
        
        rPn   = (T0 + sqrt(T0^2 + 4*EI*rA*omegan^2))/2/EI;
        rMn   = (sqrt(T0^2 + 4*EI*rA*omegan^2) - T0)/2/EI;
        lambdaMn   = sqrt(rMn);
        lambdaPn   = sqrt(rPn);
        
        sinhLn  = sinh(lambdaPn*L);
        sinLn  = sin(lambdaMn*L);
        constn    = (lambdaMn^2/lambdaPn^2)*(sinLn/sinhLn);
        
        psin    = sin(lambdaMn*x) + constn*sinh(lambdaPn*x) ;

        %-- trapezoid rule        
        %I have one "m" mode on each column
        inProdInt(n,m) = rA * psim.'*psin*h;
        if m ~= n
            %order of eigenfunction doesnt matter therefore the matrix is
            %symmetric
             inProdInt(m,n) = inProdInt(n,m);
        end
    end
end

invInProdInt = inProdInt^-1;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Computing Damping Parameters
sigma0 = zeros(modesNumber,1);
for i = 1:modesNumber
    sigma0(i) = - ComputeDamp(rho,radius,E,T0,eigenFreqs(i));
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing States & Matrices
I = speye(2*modesNumber);

x = [zeros(modesNumber,1) ; zeros(modesNumber,1)];

zeta = [zeros(modesNumber,1) ; inProdInt\modesIn];
zetaTR = zeta.';
zetaOutProd = zeta*zetaTR;

JOmega = [zeros(modesNumber),eye(modesNumber);-diag(eigenFreqs.^2),zeros(modesNumber)];

C = [zeros(modesNumber),zeros(modesNumber);zeros(modesNumber),diag(sigma0)];

%%%%% Initializing outputs
Out = zeros(timeSamples,2);

%%%%% Initializing Sherman Morrison Algorithm
%Offline computation of part of matrix A and extracting diagonal components
%into vectors for speeding up computation
A = full(I - 0.5*k*JOmega + 0.5*k*C);
A11 = diag(A(1:modesNumber,1:modesNumber));
A12 = diag(A(1:modesNumber,modesNumber + 1:end));
A21 = diag(A(modesNumber + 1:end,1:modesNumber));
A22 = diag(A(modesNumber + 1:end,modesNumber + 1:end));

%Computation of Shur Complement of A
invA11 = (1./A11); %Notice that A11=eye so invA11=A11
shurComp = A22 - A21.*(invA11.*A12);
invShurComp = 1./shurComp;

B1 = sparse(I + 0.5*k*JOmega - 0.5*k*C);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
tic
for i = 1:timeSamples
    
    %calculating bow input
    zeta1 = zetaTR*x;
    eta = zeta1 - bowVel(i);
    if desvagesFriction
        %Desvages friction
        d = sqrt(2*a)*exp(-a*eta^2 + 0.5) + 2*muD*atan(eta/0.02)/pi/eta;
        lambda = sqrt(2*a)*exp(-a*eta^2 + 0.5)*(1 - 2*a*eta^2) + 2*muD*50/pi/(2500*eta^2 + 1);
    else
        %Bilbao friction
        d = sqrt(2*a)*exp(-a*eta^2 + 0.5);
        lambda = sqrt(2*a)*exp(-a*eta^2 + 0.5)*(1 - 2*a*eta^2);
    end
    
    if smSolver
%     %Calculating known terms
        zeta2 = zeta*zeta1;
     
        B2 = B1*x;
        B = B2 + zeta2*0.5*k*Fb(i)*(lambda - 2*d) + k*Fb(i)*d*zeta*bowVel(i);
        
        b1 = B(1:modesNumber); b2 = B(modesNumber + 1:end);
        
        %Sherman morrison solver
        v = (0.5*k*Fb(i)*lambda)*zeta;
        v1 = v(1:modesNumber); v2 = v(modesNumber + 1:end);
        
        z1 = v2;
        invAv2 = invShurComp.*z1; invAv1 = - invA11.*A12.*invAv2;
        invAv = [invAv1;invAv2];
        
        y2 = A11.*b1; z2 = b2 - A21.*y2;
        invAb2 = invShurComp.*z2; invAb1 = y2 - invA11.*A12.*invAb2;
        invAb = [invAb1;invAb2];
        
        vt1 = zetaTR*invAv; vt2 = zetaTR*invAb;
        
        xNext = invAb - (1/(1+vt1))*invAv*vt2;
    else
        A = I + 0.5*k*Fb(i)*lambda*zetaOutProd - 0.5*k*JOmega;
        B = (I + 0.5*k*Fb(i)*(lambda - 2*d)*zetaOutProd + 0.5*k*JOmega - k*C)*x + k*Fb(i)*d*zeta*bowVel(i);

        xNext = A\B;
    end
    
    uNext1 = (modesOut1.')*xNext(1:modesNumber);
    uNext2 = (modesOut2.')*xNext(1:modesNumber);
    Out(i,:)= [uNext1, uNext2];


    x = xNext;
end
realTimeFrac = toc/T

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

OutPlay1 = diff(OutPlay(:,1));
OutPlay2 = diff(OutPlay(:,2));
OutPlay = [OutPlay1/max(abs(OutPlay1)),OutPlay2/max(abs(OutPlay2))];

%if play soundsc(OutPlay,SR/osFac*finalOSFac); end
% if play soundsc(OutDiff,SR); end
if play soundsc(OutDiff1,SR); end
if saveAudio
    audiowrite(fileName,OutPlay,SR/osFac*finalOSFac);
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Plots
fontSize = 18;
lineWidth = 1.5;

figure(2)
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

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Functions
function mode = ComputeBridgeMode(L, x, index)
    mode = (sqrt(2/L)*sin(index*pi*x/L))';
end

function omega = ComputeEigenFreqIsol(T0, E,I,L,rA,index)
    n = index*pi/L;
    omega = sqrt( ((T0/rA)*n.^2 + (E*I/rA)*n.^4) );
end

function omega = ComputeEigenFreqIsolBridge(E,I,L,rA,index)
    n = index*pi/L;
    omega = sqrt(((E*I/rA)*n.^4));
end


function sigma = ComputeDamp(rho,r,E,T0,omega)
    %Desvages
    rhoAir = 1.225;
    muAir = 1.619e-5;
    d0 = -2*rhoAir*muAir/(rho*r^2);
    d1 = -2*rhoAir*sqrt(2*muAir)/(rho*r);
    d2 = -1/18000;
    d3 = -0.003*E*rho*pi^2*r^6/(4*T0^2); 

    sigma = d0 + d1*sqrt(omega) + d2*omega + d3*omega^3;
    
    %Issanchou
%     rA = rho*pi*r^2;
%     nu = omega/2/pi;
%     R = 2*pi*rhoAir + 2*pi*2*r*sqrt(pi*rhoAir*muAir*nu);
%     Qa = R/(2*pi*rA*nu);
%     I = (pi*r)/ 4;
%     Qve = (4*pi^2*E*I*nu^2*4.5e-3)/T0^2;
%     Qdisl = 2.02e-4;
% 
%     sigma = pi*nu*(Qa + Qdisl + Qve);
end