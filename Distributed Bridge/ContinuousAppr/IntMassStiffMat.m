%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Calculating Integrals for Mass and Stiff matrices
%             Riccardo Russo
%          University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear 
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Computing eigenfrequencies and modes

syms x lambdaPm lambdaMm lambdaPn lambdaMn constn constm L

psim = sin(lambdaMm*x) + constm*sinh(lambdaPm*x);
psipm = diff(psim, x);
psipppm = diff(diff(psipm, x));

psin = sin(lambdaMn*x) + constn*sinh(lambdaPn*x);
psipn = diff(psin, x);
psipppn = diff(diff(psipn, x));

MassMatmn = int(psin*psim,x);
MassMatmm = int(psim*psim,x);

StiffMat1mn = int(psipn*psipm,x);
f1 = lambdaMm^2*cos(lambdaMm*x)^2 + 2*lambdaMm*cos(lambdaMm)*constm*lambdaPm*cosh(lambdaPm*x) + constm^2*lambdaPm^2*cosh(lambdaPm*x)*cosh(lambdaPm*x);
StiffMat1mm = int(f1,x);

StiffMat2mn = int(psipppn*psipppm,x);
f2 = lambdaMm^6*cos(lambdaMm*x)^2 - 2*lambdaMm^3*cos(lambdaMm*x)*constm*lambdaPm^3*cosh(lambdaPm*x) + constm^2*lambdaPm^6*cosh(lambdaPm*x)^2;
StiffMat2mm1 = int(lambdaMm^6*cos(lambdaMm*x)^2,x);
StiffMat2mm2 = int(- 2*lambdaMm^3*cos(lambdaMm*x)*constm*lambdaPm^3*cosh(lambdaPm*x),x);
StiffMat2mm3 = int(constm^2*lambdaPm^6*cosh(lambdaPm*x)^2,x);