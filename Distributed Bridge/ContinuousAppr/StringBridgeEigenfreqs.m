%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Eigenfrequencies of Modal Stiff String with distributed bridge
%                    Riccardo Russo
%                 University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear 
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
SR = 44100;                         % sample rate
k = 1 / SR;                         % time step
durSec = 5;                         % simulation length [s]

timeSamples = floor(durSec*SR);     % simulation length [samples]
timeVec = (1:timeSamples)*k;

%%%%% String Parameters
T0    = 147.7 ;            %-- tension [N]
radius     = 4.55e-4 ;          %-- radius [m]
E     = 2.5e10 ;           %-- Young's mod [Pa]
rho   = 5535 ;             %-- density [kg/m^3]
L     = 0.69 ;             %-- length [m]
% L = 0.69;
% radius = 3.75e-04;                  % string radius [m]
%         rho = 3.7575e3;                     % string Density [kg/m] 
%         T0 = 153;                           % tension [N] 
%         E = 25e9;                          % young modulus [Pa]
        A = pi*radius^2;
        rA = rho*A;
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
%%%%% Computing eigenfrequencies and modes

%first guess omega vector
omegaMax = 14e3 * 2 * pi;
omegaFirst = (20 * 2 * pi : 1 : omegaMax)';

%calculating number of modes in isolation (can't be much much greater in the
%coupled case)
nModesFirst = 0;
eigenFreq = 0;
eigenFreqsBridge = 0;
while eigenFreq < omegaMax
    nModesFirst = nModesFirst + 1;
    eigenFreq = ComputeEigenFreqIsolBridge(EB,InertiaB,LB,rAB,nModesFirst);
    if eigenFreq < omegaMax
        eigenFreqsBridge(nModesFirst) = eigenFreq;
    end
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
charactEqFull = ce0 + ce1.*(ce2 + ce3);

figure(1)
plot(omegaFirst,charactEqFull)
ylim([-10,10])
hold on
grid on;

%max of characteristic equation (for plotting purposes)
mm = max(abs(charactEqFull)) ;

%Search for zeros (first approximate search for NR 1st guesses)
zom = [] ;
for i = 1 : length(omegaFirst) - 1
    %It checks when there is a change in sign
    tt = charactEqFull(i) * charactEqFull(i+1) ;
    if tt <= 0
        if abs(charactEqFull(i))<1
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

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%Calculating eigenfrequencies of stiff string in isolation for comparison
modesNumberIsol = 1;
while true
    if ComputeEigenFreqIsol(T0, E,Inertia,L,rA,modesNumberIsol)>omegaMax
        modesNumberIsol = modesNumberIsol-1;
        break; 
    end
    modesNumberIsol = modesNumberIsol+1;
end

eigenFreqsIsol = ComputeEigenFreqIsol(T0, E,Inertia,L,rA,(1:modesNumberIsol).');

modesNumberIsol = length(eigenFreqsIsol) ;

figure(2)
line([eigenFreqs,eigenFreqs],[-mm,mm],'linestyle','--','color','red','DisplayName','Eigenfreqs of Coupled System');
hold on
grid on
line([eigenFreqsIsol,eigenFreqsIsol],[-mm,mm],'linestyle','--','color','k','DisplayName','Eigenfreqs of String in Isolation');
title('Eigenfrequencies of Coupled System (red) and in Isolation (black)')

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