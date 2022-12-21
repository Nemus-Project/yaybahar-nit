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
StiffMatBridgeTot = -(EIB/rAB)*diag(stiffMatBridgeDiag) + contactModesBridge*(-T0*ComputeBridgeModeDiff(LB, zContact, (1:nModesFirst))' + EI*ComputeBridgeModeDiff3(LB, zContact, (1:nModesFirst))');

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Functions
function mode = ComputeBridgeMode(L, x, index)
    mode = (sqrt(2/L)*sin(index*pi*x/L))';
end
function mode = ComputeBridgeModeDiff(L, x, index)
    mode = (sqrt(2/L)*(index*pi/L).*cos(index*pi*x/L))';
end
function mode = ComputeBridgeModeDiff3(L, x, index)
    mode = -(sqrt(2/L)*(index*pi/L).^3.*cos(index*pi*x/L))';
end

function omega = ComputeEigenFreqIsolBridge(E,I,L,rA,index)
    n = index*pi/L;
    omega = sqrt(((E*I/rA)*n.^4));
end

function omega = ComputeEigenFreqIsol(T0, E,I,L,rA,index)
    n = index*pi/L;
    omega = sqrt( ((T0/rA)*n.^2 + (E*I/rA)*n.^4) );
end