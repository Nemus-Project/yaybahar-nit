%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%     Modal spring dispersion relation
%              Riccardo Russo
%           University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters

omegaMax = 20000*2*pi;

beta = (1:0.1:400).';

R = 5e-3;
r = 1e-3;
alpha = 8;

E = 2e11;
ni = 0.3;   
rho = 7.872e3;
L = 0.3;

kappaSq = (r/2)^2;

mu = tand(alpha);
l = R/(cosd(alpha)^2);
I = pi*r^4/4;
Iphi = 2*I;
G = E/(2*(1 + ni));

A = pi*r^2;

digits(5)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%% Dispersion Relation

[m1, m2, m3, m4] = ComputeDispRelMatElems(rho, A, l, mu, kappaSq, E, ni, beta);

delta = (m1+m4).^2 - 4*(m1.*m4 - m2.*m3);
checkDelta = (delta < 0);
if any(checkDelta) 
    disp("delta < 0");
    delta(checkDelta) = [];
    beta(checkDelta) = [];
    m1(checkDelta) = [];
    m4(checkDelta) = [];
end

omegaSqP = ((m1+m4) + sqrt(delta))/2;
omegaSqM = ((m1+m4) - sqrt(delta))/2;

dispRelP = sqrt(omegaSqP);
dispRelM = sqrt(omegaSqM);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%% Eigenfreqs & Eigenvecs

%%%%% Computing max modes number
m = 1;
S = pi*m/L;

omegaM = 0;
while true
    
    [m1P, m2P, m3P, m4P] = ComputeDispRelMatElems(rho, A, l, mu, kappaSq, E, ni, S);
    
    deltaP = (m1P+m4P)^2 - 4*(m1P*m4P - m2P*m3P);

    if delta >=0
        omegaP = sqrt(((m1P+m4P) + sqrt(deltaP))/2);
        omegaM = sqrt(((m1P+m4P) - sqrt(deltaP))/2);

        if omegaM >= omegaMax break; end

        m = m+1;
        S = m*pi/L;
    else
        disp("delta < 0\n");
        m;
    end
end
modesNumber = m - 1;

%%%%% Computing eigenfrequencies - METHOD 1
eigenFreqs = zeros(modesNumber,3);
eigenVecs = zeros(2*modesNumber,2);
for m = 1:modesNumber
    S = m*pi/L;

    [m1P, m2P, m3P, m4P] = ComputeDispRelMatElems(rho, A, l, mu, kappaSq, E, ni, S);
    
%     %Solving characteristic polynomial for each m
%     deltaP = (m1P+m4P)^2 - 4*(m1P*m4P - m2P*m3P);
%     omegaP = sqrt(((m1P+m4P) + sqrt(deltaP))/2);
%     omegaM = sqrt(((m1P+m4P) - sqrt(deltaP))/2);
    %Non ho capito come calcolare gli autovettori quindi uso semplicemente
    %eig

    [evcs,evls] = eig([m1P,m2P;m3P,m4P]);
    [evls, ind] = sort(diag(evls));
    evcs = evcs(:,ind);

    omegaM = sqrt(evls(1));
    omegaP = sqrt(evls(2));

    eigenFreqs(m,1) = omegaP;
    eigenFreqs(m,2) = omegaM;
    eigenFreqs(m,3) = S;
    eigenVecs(2*m - 1:2*m, :) = evcs;
end


%%%%% Computing eigenfrequencies and eigenvectors - METHOD 2
S = (1:modesNumber).'*pi/L;
[m1P, m2P, m3P, m4P] = ComputeDispRelMatElems(rho, A, l, mu, kappaSq, E, ni, S);

massStiffMat = [diag(m1P), diag(m2P); diag(m3P), diag(m4P)];
[V,D] = eig(massStiffMat);

[Dsort,Dind] = sort(diag(D));
Vsort = V(:,Dind);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%% Plotting
figure(1)
plot(beta,dispRelM/2/pi);
hold on
plot(beta,dispRelP/2/pi);
hold on
yline(20000);
hold on
plot(eigenFreqs(:,3), eigenFreqs(:,1)/2/pi, '*');
hold on
plot(eigenFreqs(:,3), eigenFreqs(:,2)/2/pi, '*');

eigenFreqsTot = [eigenFreqs(:,1); eigenFreqs(:,2)];
eigenFreqsTot = sort(eigenFreqsTot);

figure(2)
plot(eigenFreqsTot, '*');
hold on
plot(sort(sqrt(diag(D))), '+');


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%Functions

function [m1, m2, m3, m4] = ComputeDispRelMatElems(rho, A, l, mu, kappaSq, E, ni, beta)
    a4 = 1./(rho*A*(1+l^2*beta.^2));
    r2 = (1 - mu^2)/l - beta.^2*l;
    r4 = 2*mu*(1/l - beta.^2*l);
    d4 = E*kappaSq*A./(1+ni+l^2*beta.^2);
    
    m1 = beta.^2.*((2*mu/l)^2*E*kappaSq*A + r2.^2.*d4)/(rho*A);
    m2 = beta.^2.*(r2.*r4.*d4 - 2*mu*E*kappaSq*A*r2/l)/(rho*A);
    m3 = beta.^2.*(a4.*r4.*d4.*r2 - 2*E*kappaSq*A*a4.*r2*mu/l);
    m4 = beta.^2.*(E*kappaSq*A*r2.^2.*a4 + a4.*d4.*r4.^2);
end