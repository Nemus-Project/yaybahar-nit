%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Eigenfrequencies and eigenmodes of Stiff String with distributed bridge
%                   FWD conditions
%          Contact on u_Ms: PASSIVE ENERGY BALANCE
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

zContact = 0.13*LB;

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

%Nota che bisogna usare V^-1 e non V.' perchè M^-1*K non è simmetrica!
%Infatti il teorema spettrale vale per matrici simmetriche. Il fatto è che
%una matrice è diagonalizzabile se esiste una base ortogonale di
%autovettori, e in questo caso c'è! Perchè il prodotto interno qua è
%sostanzialmente l'identità. Non capisco però come sia possibile che escano
%gli autovettori ortonormali, perchè la matrice non è simmetrica. Forse il
%motivo è esattamente questa "quasi" ortonormalità. Ma allora cosa succede
%per paramtri diversi di corda/ponte?
inProd = V^-1*V;

figure(1)
imagesc(inProd);

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
line([eigenFreqs/2/pi,eigenFreqs/2/pi],[-1,1],'linestyle','--','color','red','DisplayName','Eigenfreqs of Coupled System');
hold on
grid on
line([eigenFreqsIsol/2/pi,eigenFreqsIsol/2/pi],[-1,1],'linestyle','--','color','k','DisplayName','Eigenfreqs of String in Isolation');
title('Eigenfrequencies of Coupled System (red) and in Isolation (black)')

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Functions
function omega = ComputeEigenFreqIsol(T0, E,I,L,rA,index)
    n = index*pi/L;
    omega = sqrt( ((T0/rA)*n.^2 + (E*I/rA)*n.^4) );
end

