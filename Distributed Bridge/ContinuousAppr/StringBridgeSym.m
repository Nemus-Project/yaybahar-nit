%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Eigenfrequencies of Modal Stiff String with distributed bridge
%                    Riccardo Russo
%                 University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear 
close all

%% Computing Eigenfrequencies
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
SR = 44100;                         % sample rate
k = 1 / SR;                         % time step
durSec = 5;                         % simulation length [s]

timeSamples = floor(durSec*SR);     % simulation length [samples]
timeVec = (1:timeSamples)*k;

%%%%% String Parameters
L = 0.69;                          % string length [m]
% radius = 4.55e-04;
% rho = 5535;                          % string Density [kg/m^3] 
% T0 = 147.7;                           % tension [N] 
% A = pi*radius^2;
% rA = rho*A;
% E = 2.5e10;                          % young modulus [Pa]
% Inertia = (pi*radius^4)/ 4;         % moment of inertia
% EI    = E*Inertia ;
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

outPoint = 0.5233*L;
inPoint = 0.35432*L;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% Initial Conditions: Raised Cosine
%{
    Position of excitation (at maximum) in percentage, 
    starting from left.
    0< excitPosRatio <100
%}
h = 1e-5 * L;
x = (0 : h : L).';
N = size(x,1);
excitPosRatio = 80;
excitPos = floor(N*excitPosRatio/100);
widthMeter = 0.1;
width = floor(widthMeter / h);
amplitude = 0.1;
initDistr = zeros(N,1);
for i = 1:N
    if abs(i-excitPos)<=width
        initDistr(i) = amplitude*0.5*(1+cos(pi*(i-excitPos)/width));
    else
        initDistr(i) = 0;
    end
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
stiffMatBridgeDiag = (-pi*(1:nModesFirst)/LB).^4';
% stiffMatBridge = diag(stiffMatBridgeDiag);

lambdaP = ComputeLambdaP(T0, EI, rA, omegaFirst);
lambdaM = ComputeLambdaM(T0, EI, rA, omegaFirst);

%characteristicEq
ce0 = - ComputeModeAtL(lambdaM, lambdaP, L);
ce1 = ComputeBridgeCoeff(rAB, EIB, stiffMatBridgeDiag, contactModesBridge, omegaFirst, nModesFirst);
ce2 = -T0*ComputeModeDiffAtL(lambdaM, lambdaP, L);
ce3 = + EI*ComputeModeDiff3AtL(lambdaM, lambdaP, L);
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

        lambdaP = ComputeLambdaP(T0, EI, rA, omega);
        lambdaM = ComputeLambdaM(T0, EI, rA, omega);
        
        ce0 = - ComputeModeAtL(lambdaM, lambdaP, L);
        ce1 = ComputeBridgeCoeff(rAB, EIB, stiffMatBridgeDiag, contactModesBridge, omega, nModesFirst);
        ce2 = -T0*ComputeModeDiffAtL(lambdaM, lambdaP, L);
        ce3 = + EI*ComputeModeDiff3AtL(lambdaM, lambdaP, L);
        charactEq = ce0 + ce1 * (ce2 + ce3);
        
        dce0 = ComputeModeAtLPrime(T0, EI, rA, omega, L);
        dce1 = ComputeBridgeCoeffPrime(rAB, EIB, stiffMatBridgeDiag, contactModesBridge, omega, nModesFirst);
        dce2 = -T0*ComputeModeDiffAtLPrime(T0, EI, rA, omega, L);
        dce3 = EI*ComputeModeDiff3AtLPrime(T0, EI, rA, omega, L);

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

lambdaP = ComputeLambdaP(T0, EI, rA, eigenFreqs);
lambdaM = ComputeLambdaM(T0, EI, rA, eigenFreqs);

modesOut = ComputeModeAtx(lambdaM, lambdaP, L, outPoint);
modesIn = ComputeModeAtx(lambdaM, lambdaP, L, inPoint);


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%Calculating integrals for modal expansion (trapezoidal rule)
MassMat = zeros(modesNumber);
initCondInt = zeros(modesNumber,1);

for m = 1 : modesNumber
   omegam = eigenFreqs(m);

    lambdaMm = ComputeLambdaM(T0, EI, rA, omegam);
    lambdaPm = ComputeLambdaP(T0, EI, rA, omegam);

    psim = ComputeModeAtx(lambdaMm, lambdaPm, L, x);
    psim(1) = 0.5*psim(1) ; psim(end) = 0.5*psim(end) ;

    for n = m : modesNumber
        omegan = eigenFreqs(n) ;

        lambdaMn   = ComputeLambdaM(T0, EI, rA, omegan);
        lambdaPn   = ComputeLambdaP(T0, EI, rA, omegan);
        
        psin = ComputeModeAtx(lambdaMn, lambdaPn, L, x);
        %-- trapezoid rule        
        %I have one "m" mode on each column
        MassMat(n,m) = rA * psim.'*psin*h;
        if m ~= n
            %order of eigenfunction doesnt matter therefore the matrix is
            %symmetric
             MassMat(m,n) = MassMat(n,m);
        end
    end
    initCondInt(m) = psim.'*initDistr*h;
end 
invMassMat = MassMat^-1;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Computing Damping Parameters
sigma0 = zeros(modesNumber,1);
for i = 1:modesNumber
    sigma0(i) = - ComputeDamp(rho,radius,E,T0,eigenFreqs(i));
end

%% Simulation
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing oscillators and coefficients

for i = 1:modesNumber
    if k>=(2/eigenFreqs(i))
        disp("stability condition violated");
        return;
    end
end

% q = zeros(modesNumber,1);
q = invMassMat*initCondInt;
qPrev = q;
qNext = zeros(modesNumber,1);

if isnan(q)
    disp("Error: Non Invertible Matrix");
    return
end

A = sigma0*k + 1;
B = 2 - eigenFreqs.^2*k^2;
C = sigma0*k - 1;
D = 2*k^2/L/rA;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
output = zeros(1,timeSamples);

tic
for n = 1:timeSamples
    qNext = q.*B./A + qPrev.*C./A;

    qPrev = q;
    q = qNext;

    for i = 1:modesNumber
        output(n) = output(n) + qNext(i)*modesOut(i);
    end
end
realTimeFrac = toc/durSec

figure(2)
plot(output);
soundsc(output,SR);
%audiowrite("StringBridge.wav",output,SR);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Functions
function mode = ComputeBridgeMode(L, x, index)
    mode = (sqrt(2/L)*sin(index*pi*x/L))';
end

function omega = ComputeEigenFreqIsolBridge(E,I,L,rA,index)
    n = index*pi/L;
    omega = sqrt(((E*I/rA)*n.^4));
end

function omega = ComputeEigenFreqIsol(T0, E,I,L,rA,index)
    n = index*pi/L;
    omega = sqrt( ((T0/rA)*n.^2 + (E*I/rA)*n.^4) );
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

function f = ComputeLambdaM(T0, EI, rA, omega)
    f = sqrt((sqrt(T0^2 + 4*EI*rA*omega.^2) - T0)/2/EI);
end
function f = ComputeLambdaP(T0, EI, rA, omega)
    f = sqrt((T0 + sqrt(T0^2 + 4*EI*rA*omega.^2))/2/EI);
end

function f = ComputeModeAtx(lambdaM, lambdaP, L, x)
    sinhL = sinh(lambdaP*L);
    sinL = sin(lambdaM*L);

    const = (lambdaM.^2./lambdaP.^2).*(sinL./sinhL);
    f = (sin(lambdaM*x) + const.*sinh(lambdaP*x));
end

function f = ComputeModeAtL(lambdaM, lambdaP, L)
    sn = sin(lambdaM*L);
    f = sn.*(1 + (lambdaM.^2./lambdaP.^2));
end
function f = ComputeModeDiffAtL(lambdaM, lambdaP, L)
    th = tanh(lambdaP*L);
    cs = cos(lambdaM*L);
    sn = sin(lambdaM*L);
    f = (lambdaM.*cs + (lambdaM.^2./lambdaP).*(sn./th));
end
function f = ComputeModeDiff3AtL(lambdaM, lambdaP, L)
    th = tanh(lambdaP*L);
    cs = cos(lambdaM*L);
    sn = sin(lambdaM*L);
    f = (-lambdaM.^3.*cs + (lambdaP.*lambdaM.^2).*(sn./th));
end
function f = ComputeBridgeCoeff(rAB, EIB, stiffMatDiag, contactModesBridge, omega, nModes)
    idVec = ones(nModes,1);
    f = zeros(size(omega,1), 1);
    for i = 1:size(omega)
        omegaTerm = (-rAB*omega(i)^2*idVec + EIB*stiffMatDiag).^(-1);
        f(i) = (contactModesBridge.*omegaTerm)'*contactModesBridge;
    end
end
function f = ComputeModeAtLPrime(T0, EI, rA, omega, L)
    f = (L*omega*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/(T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2) - 1))/((-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) - sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((2*EI*omega*rA)/((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) + (2*EI*omega*rA*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)^2*(T0^2 + 4*EI*rA*omega^2)^(1/2)));
end
function f = ComputeModeDiffAtLPrime(T0, EI, rA, omega, L)
    f = ((omega*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)))/((-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) - (L*omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)))/(T0^2 + 4*EI*rA*omega^2)^(1/2) + (2*omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)))/(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) + (omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(3/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) - (L*omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))^2 - 1)*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))^2*(T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) - (L*omega*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)));    
end
function f = ComputeModeDiff3AtLPrime(T0, EI, rA, omega, L)
    f = - ((3*omega*rA*cos((-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))/(T0^2 + 4*EI*rA*omega^2)^(1/2) + (omega*rA*sin((-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(EI*(T0^2 + 4*EI*rA*omega^2)^(1/2)) - (2*omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))/(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0^2 + 4*EI*rA*omega^2)^(1/2)) + (omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)) + (L*omega*rA*sin(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))^2 - 1)*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2))/(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))^2*(T0^2 + 4*EI*rA*omega^2)^(1/2)) + (L*omega*rA*cos(L*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))/(EI*tanh(L*((T0/2 + (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2))*(-(T0/2 - (T0^2 + 4*EI*rA*omega^2)^(1/2)/2)/EI)^(1/2)*(T0^2 + 4*EI*rA*omega^2)^(1/2)));
end
function f = ComputeBridgeCoeffPrime(rAB, EIB, stiffMatDiag, contactModesBridge, omega, nModes)
    idVec = ones(nModes,1);
    f = zeros(size(omega,1), 1);
        for i = 1:size(omega)
        omegaTerm = (-rAB*omega(i)^2*idVec + EIB*stiffMatDiag).^(-1);
        dOmegaTerm = -2*omega*rAB*idVec;
        f(i) = (-contactModesBridge.*(omegaTerm.*dOmegaTerm.*omegaTerm))'*contactModesBridge;
        end  
end