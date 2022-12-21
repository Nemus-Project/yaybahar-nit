%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Modal Stiff String with connected mass-spring
%         Riccardo Russo
%       University of Bologna
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

L = 0.69;                          % string length [m]
radius = 4.55e-04;
rho = 5535;                          % string Density [kg/m] 
T0 = 147.7;                           % tension [N] 
A = pi*radius^2;
rA = rho*A;
E = 2.5e10;                          % young modulus [Pa]
Inertia = (pi*radius^4)/ 4;         % moment of inertia
EI    = E*Inertia;
sigma0 = 2;                       % damping parameter

%-- bridge parameters
M     = 1e-2;             %-- mass [Kg]
K     = 1e6;             %-- stiffness [N/m]

maxFreq = 10000*2*pi;

%where to calculate modes
modesPoints = 0:0.1:L;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Computing eigenfrequencies and modes
%--------------------------------------------------------

%-- newton-raphson routine for eigenvalues (eigenfrequencies)

omegaFirst    = 20 * 2 * pi : 0.1 : maxFreq;

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

modes = zeros(modesNumber,length(modesPoints));
for i=1:length(modesPoints)
    modes(:,i) = (sin(lambdaM*modesPoints(i)) + const.*sinh(lambdaP*modesPoints(i))) + i*2;
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%Calculating eigenfrequencies of stiff string in isolation for comparison
modesNumberIsol = 1;
while true
    if ComputeEigenFreqIsol(T0, E,Inertia,L,rA,modesNumberIsol)>(maxFreq)
        modesNumberIsol = modesNumberIsol-1;
        break; 
    end
    modesNumberIsol = modesNumberIsol+1;
end

eigenFreqsIsol = ComputeEigenFreqIsol(T0, E,Inertia,L,rA,(1:modesNumberIsol).');

modesNumberIsol = length(eigenFreqsIsol) ;

modesIsol = zeros(modesNumberIsol,length(modesPoints));
for i=1:length(modesPoints)
    modesIsol(:,i) = sin((1:modesNumberIsol)*pi*modesPoints(i)/L) + i*2;
end

figure(2)
line([eigenFreqs,eigenFreqs],[-mm,mm],'linestyle','--','color','red','DisplayName','Eigenfreqs of Coupled System');
hold on
grid on
line([eigenFreqsIsol,eigenFreqsIsol],[-mm,mm],'linestyle','--','color','k','DisplayName','Eigenfreqs of String in Isolation');
title('Eigenfrequencies of Coupled System (red) and in Isolation (black)')

figure(3)
plot(modes, 'color','k')
hold on
plot(modesIsol, 'color','red');
title('Modes of Coupled System (red) and in Isolation (black)')

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
function omega = ComputeEigenFreqIsol(T0, E,I,L,rA,index)
    n = index*pi/L;
    omega = sqrt( ((T0/rA)*n.^2 + (E*I/rA)*n.^4) );
end