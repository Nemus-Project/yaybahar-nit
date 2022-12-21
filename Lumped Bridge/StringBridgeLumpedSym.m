%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Modal Stiff String with connected mass-spring
%         Riccardo Russo
%       University of Bologna
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

L = 0.69;                          % string length [m]
% radius = 4.55e-04;
% rho = 5535;                          % string Density [kg/m] 
% T0 = 147.7;                           % tension [N] 
% A = pi*radius^2;
% rA = rho*A;
% E = 2.5e10;                          % young modulus [Pa]
% Inertia = (pi*radius^4)/ 4;         % moment of inertia
% EI    = E*Inertia;
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

%-- bridge parameters
M     = 1e-2;             %-- mass [Kg]
K     = 1e10;             %-- stiffness [N/m]

outPoint = 0.5233*L;
inPoint = 0.35432*L;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% Initial Conditions: Raised Cosine
%{
    Position of excitation (at maximum) in percentage, 
    starting from left.
    0< excitPosRatio <100
%}
h = 1e-6 * L;
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

%-- newton-raphson routine for eigenvalues (eigenfrequencies)
omegaFirst    = 20 * 2 * pi : 0.1 : 10000 * 2 * pi;

rP    = (T0 + sqrt(T0^2 + 4*EI*rA*omegaFirst.^2))/2/EI;
rM    = (sqrt(T0^2 + 4*EI*rA*omegaFirst.^2) - T0)/2/EI;

lambdaM   = sqrt(rM);
lambdaP   = sqrt(rP);

th    = tanh(lambdaP*L);
cs    = cos(lambdaM*L);
sn    = sin(lambdaM*L);
g     = M*(-omegaFirst.^2) + K;

%-- eigenvale equation
charactEq   = -T0*(lambdaM.*th.*cs + lambdaM.^2./lambdaP.*sn) + EI*(-lambdaM.^3.*th.*cs + lambdaM.^2.*lambdaP.*sn) - g.*th.*sn.*(1+lambdaM.^2./lambdaP.^2);
%Qua plotta charactEq, cioè la equazione caratteristica. Non c''è un vero motivo
%in realtà ma è giusto per far vedere dove solo gli zeri
figure(1)
plot(omegaFirst,charactEq); 
hold on; 
grid on;

%max of characteristic equation (for plotting purposes)
mm = max(abs(charactEq));

%Search for zeros (first approximate search for NR 1st guesses)
zom = [];
for i = 1 : length(omegaFirst) - 1
    %It checks when there is a change in sign
    tt = charactEq(i) * charactEq(i+1);
    if tt <= 0
        zom = [zom;omegaFirst(i)];  
    end
    
end

%-- initialise zero vector
eigenFreqs = zeros(length(zom),1);

%-- newton Raphson
for n = 1 : length(zom)
    
    omega = zom(n);
    tol = 1;
    iter = 0;
    
    while tol > 1e-13 && iter < 20
        
        iter  = iter + 1;
        rP    = (T0 + sqrt(T0^2 + 4*EI*rA*omega.^2))/2/EI;
        rM    = (sqrt(T0^2 + 4*EI*rA*omega.^2) - T0)/2/EI;
        lambdaM   = sqrt(rM);
        lambdaP   = sqrt(rP);
        th    = tanh(lambdaP*L);
        cs    = cos(lambdaM*L);
        sn    = sin(lambdaM*L);
        g     = M*(-omega.^2) + K;

        charactEq     = -T0*(lambdaM.*th.*cs + lambdaM.^2./lambdaP.*sn) + EI*(-lambdaM.^3.*th.*cs + lambdaM.^2.*lambdaP.*sn) - g.*th.*sn.*(1+lambdaM.^2./lambdaP.^2);
        dCharactEq    = (L.*omega.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(- M.*omega.^2 + K).*((T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./(T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2) - 1))./((-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)) - T0.*((2.*omega.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)))./(((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)) + (omega.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)))./((-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)) - (L.*omega.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)))./(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2) + (omega.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2))./(EI.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(3./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)) - (L.*omega.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).^2 - 1).*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2))./(((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)) - (L.*omega.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2))./(EI.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2).*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2))) - sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(- M.*omega.^2 + K).*((2.*EI.*omega.*rA)./((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)) + (2.*EI.*omega.*rA.*(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2))./((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2).^2.*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2))) - 2.*M.*omega.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*((T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./(T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2) - 1) - (L.*omega.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).^2 - 1).*(- M.*omega.^2 + K).*((T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./(T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2) - 1))./(((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)) - EI.*((3.*omega.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2))./(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2) - (2.*omega.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2))./(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2) + (omega.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2))./(EI.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)) - (L.*omega.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).^2 - 1).*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(3./2))./(((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)) + (L.*omega.*rA.*sin(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*tanh(L.*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2))./(EI.*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)) + (L.*omega.*rA.*cos(L.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2)).*(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2).*((T0./2 + (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2))./(EI.*(-(T0./2 - (T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)./2)./EI).^(1./2).*(T0.^2 + 4.*EI.*rA.*omega.^2).^(1./2)));
        
        omegaNew     = omega - charactEq/dCharactEq;
        
        tol     = abs(omegaNew-omega);
        omega      = omegaNew;
        
    end
    
    eigenFreqs(n) = omega;
    
    line([omega,omega],[-mm,mm],'linestyle','--','color','k');
    
end
modesNumber = length(eigenFreqs);

rPlus    = (T0 + sqrt(T0^2 + 4*EI*rA*eigenFreqs.^2))/2/EI;
rMinus    = (sqrt(T0^2 + 4*EI*rA*eigenFreqs.^2) - T0)/2/EI;

lambdaM   = sqrt(rMinus);
lambdaP   = sqrt(rPlus);

sinhL    = sinh(lambdaP*L);
sinL    = sin(lambdaM*L);

const = (lambdaM.^2./lambdaP.^2).*(sinL./sinhL);

modesOut = (sin(lambdaM*outPoint) + const.*sinh(lambdaP*outPoint));
modesIn = (sin(lambdaM*inPoint) + const.*sinh(lambdaP*inPoint));

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%Calculating integrals for modal expansion (trapezoidal rule)
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

    for n = 1 : modesNumber
        
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
        inProdInt(n,m) = psim.'*psin*h;
    end
    initCondInt(m) = psim.'*initDistr*h;
end

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
q = inProdInt^(-1)*initCondInt;
qPrev = q;
qNext = zeros(modesNumber,1);

if isnan(q)
    disp("Non invertible Matrix");
    return;
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
audiowrite("StringBridgeLumped.wav",output,SR);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
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