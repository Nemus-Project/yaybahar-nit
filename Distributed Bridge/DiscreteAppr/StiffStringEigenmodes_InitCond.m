%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Stiff String with distributed bridge Initial Conditions
%                    Riccardo Russo
%                 University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear 
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
realTimeDraw = false;

SR = 44100;                         % sample rate
k = 1 / SR;                         % time step
durSec = 5;                         % simulation length [s]

h = 0.5e-2;

omegaMax = 14e3*2*pi;

timeSamples = floor(durSec*SR);     % simulation length [samples]
timeVec = (1:timeSamples)*k;

%%%%% String Parameters
T0    = 147.7 ;            %-- tension [N]
radius     = 4.55e-4 ;          %-- radius [m]
E  = 2.5e10 ;           %-- Young's mod [Pa]
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
rhoB =8e3;                          % string Density [kg/m^3] 
AB = pi*radiusB^2;
rAB = rhoB*AB;
EB = 3e15;                          % young modulus [Pa]
InertiaB = (pi*radiusB^4)/ 4;         % moment of inertia
EIB    = EB*InertiaB;

zContact = 0.43*LB;

outPoint = 0.5233*L;
inPoint = 0.35432*L;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%% Initializing Discrete Operators
Ms = floor(L/h);
Mb = floor(LB/h);

Jvec = sparse(Mb - 1, 1); 
a = zContact/ h - floor(zContact/h);
Jvec(floor(zContact/h)) = (1-a)/h;
Jvec(floor(zContact/h) + 1) = a/h;
JJT = Jvec*Jvec.';

vs = ones(Ms-1,1);
DxxS = spdiags([vs/h^2 -2*vs/h^2 vs/h^2], -1:1, Ms-1, Ms-1);
DxxxxS = DxxS*DxxS;
vb = ones(Mb-1,1);
DxxB = spdiags([vb/h^2 -2*vb/h^2 vb/h^2], -1:1, Mb-1, Mb-1);
DxxxxB = DxxB * DxxB;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%% Building Mass-Stiff Matrices
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

eigenFreqs = sqrt(Dsort);
eigenFreqs = eigenFreqs(eigenFreqs<=omegaMax & eigenFreqs>=20*2*pi);
modesNumber = length(eigenFreqs);
modesMat = Vsort(:,1:modesNumber);

inProd = Vsort^-1*Vsort;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% Initial Conditions: Raised Cosine
%{
    Position of excitation (at maximum) in percentage, 
    starting from left.
    0< excitPosRatio <100
%}
x = (0:Ms-2).';
excitPosRatio = 80;
excitPos = floor((Ms-1)*excitPosRatio/100);
widthMeter = 0.1;
width = floor(widthMeter / h);
amplitude = 0.1;
initDistr = zeros((Ms-1),1);
for i = 1:(Ms-1)
    if abs(i-excitPos)<=width
        initDistr(i) = amplitude*0.5*(1+cos(pi*(i-excitPos)/width));
    else
        initDistr(i) = 0;
    end
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

u = initDistr;
w = zeros(Mb-1,1);
v = [u;w];

q = Vsort^-1*v;
%updating only oscillators in audio band
q = q(1:modesNumber);
qPrev = q;
qNext = zeros(modesNumber,1);

A = sigma0*k + 1;
B = 2 - eigenFreqs.^2*k^2;
C = sigma0*k - 1;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
output = zeros(1,timeSamples);

tic
for n = 1:timeSamples
    qNext = q.*B./A + qPrev.*C./A;

    qPrev = q;
    q = qNext;

    v = modesMat*qNext;

    %%plot
    if realTimeDraw
        plot(v)
        ylim([-0.2,0.2])
        drawnow;
    %     pause(0.1)
    end

    output(n) = v(148);
end
realTimeFrac = toc/durSec

figure(2)
plot(output);
soundsc(output,SR);
%audiowrite("StringBridge.wav",output,SR);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Functions
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


