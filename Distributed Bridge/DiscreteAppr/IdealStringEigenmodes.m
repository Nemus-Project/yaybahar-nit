%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Eigenfrequencies and eigenmodes of Ideal String with distributed bridge
%                   FWD conditions
%               PASSIVE ENERGY BALANCE
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
rhoB =8e3;                          % string Density [kg/m^3] 
AB = pi*radiusB^2;
rAB = rhoB*AB;
EB = 3e10;                          % young modulus [Pa]
InertiaB = (pi*radiusB^4)/ 4;         % moment of inertia
EIB    = EB*InertiaB;

zContact = 0.43*LB;

outPoint = 0.5233*L;
inPoint = 0.35432*L;

Ms = floor(L/h);
Mb = floor(LB/h);

Jvec = sparse(Mb - 1, 1); Jvec(floor(zContact/h)) = 1/h;

vs = ones(Ms-1,1);
DxxS = T0*spdiags([vs/h^2 -2*vs/h^2 vs/h^2], -1:1, Ms-1, Ms-1);
vb = ones(Mb-1,1);
DxxB = spdiags([vb/h^2 -2*vb/h^2 vb/h^2], -1:1, Mb-1, Mb-1);
DxxxxB = EIB * DxxB * DxxB;

zero1 = sparse(Ms-1, 1); zero1(end) = T0/h^2;
zero2 = sparse(Ms-1, Mb -1);
vec1 = T0*Jvec.'/h;
vec2 = -DxxxxB - T0*(Jvec*Jvec.');
coeff =  -2*T0/h^2;

%nota che vec1, vec2 e coeff sono semplicemente delle continuazioni di
%DxxS, potevo semplicemente definire DxxS su Ms x Ms invece che su (Ms -
%1) x (Ms - 1), solo che poi era un casino concatenare tutti i pezzi

K = [-DxxS, -zero1, -zero2;
     -zero1.',-coeff, -vec1;
     -zero2.', -vec1.', -vec2];
M = [rA * ones(Ms, 1); rAB * ones(Mb - 1, 1)];
M = sparse(diag(M));
R = full(M\K);
[V, D] = eig(R);

[Dsort, Dind] = sort(diag(D));
Vsort = V(:,Dind);

%Nota che bisogna usare V^-1 e non V.' perchè M^-1*K non è simmetrica!
%Lo sarebbe se rAB fosse uguale a rA, a quel punto la matrice di massa
%sarebbe diagonale con la diagonale costante!
inProd = V^-1*V;

figure
imagesc(inProd);

