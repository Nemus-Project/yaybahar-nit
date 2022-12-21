%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Calculating Derivative for Characteristic Equation
%             Riccardo Russo
%          University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear 
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Computing eigenfrequencies and modes

syms omega L EI T0 rA LB rAB EIB

rP    = (T0 + sqrt(T0^2 + 4*EI*rA*omega^2))/2/EI;
rM    = (sqrt(T0^2 + 4*EI*rA*omega^2) - T0)/2/EI;

lambdaM   = sqrt(rM);
lambdaP   = sqrt(rP);

f0 = - sin(lambdaM*L)*(1 + lambdaM^2/lambdaP^2);
f2 = - T0*(lambdaM*cos(lambdaM*L) + (lambdaM^2/lambdaP)*(sin(lambdaM*L)/tanh(lambdaP*L)));
f3 = EI*(-lambdaM^3*cos(lambdaM) + (sin(lambdaM*L)/tanh(lambdaP*L))*lambdaP*lambdaM^2);
df0 = diff(f0, omega);
df2 = diff(f2, omega);
df3 = diff(f3, omega);