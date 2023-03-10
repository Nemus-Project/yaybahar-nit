%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%           Modal spring simulation
%              Riccardo Russo
%           University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Excitation Signal
%'impulse', 'sin', 'square', 'sweep', 'audiofile'
inputType = 'square'; 
%in case inputType is 'audiofile', specify file name and path
audiofileName = 'dry_samples/DrumReference.wav';
%amplification factor
amp = 1; 
osFac = 1;
%if inputDur is set to 0 when audiofileName is 'audiofile', the entire file
%is played, otherwise only a portion of it
inputDur = 2;
%tail indicates how long sound should be played after the input signal is
%stopped. The value returned with durSec is inputDur + tail
tail = 1;
% [excit,SR,k,timeSamples,timeVec,durSec] = ExcitSignal(amp,osFac,durSec,tail,inputType,audiofileName);
[excit,SR,k,timeSamples,timeVec,durSec] = ExcitSignal(amp,osFac,inputDur,tail,inputType);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
play = true;
playDry = false;
saveAudio = false;

%Normalize output (before visualization)
normalizeOut = true;

% if 0 reads displacement, otherwise force at output
readForce = 0;

omegaLim = 2/k;

R = 1e-2;
r = 0.5e-3;
alpha = 2;

E = 2e11;
ni = 0.3;   
rho = 7.872e3;
L = 10;

kappaSq = (r/2)^2;

mu = tand(alpha);
l = R/(cosd(alpha)^2);
I = pi*r^4/4;
Iphi = 2*I;
G = E/(2*(1 + ni));

A = pi*r^2;

outPoint = 0.95*L;
inPoint = 0;%0.35432*L;

%sets how much excitation goes to one or the other polarization
alphaParam = 0.5;
normalAlphaParam = sqrt(alphaParam^2 - 2*alphaParam + 1);

sigma0 = 2;
sigma1 = 1e-6;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%% Computing eigenfrequencies and modes
%%%%% Computing max modes number
m = 1;
modeCounter = 1;
S = pi*m/L;

omegaM = 0;
while true
    
    [m1P, m2P, m3P, m4P] = ComputeDispRelMatElems(rho, A, l, mu, kappaSq, E, ni, S);
    
    deltaP = (m1P+m4P)^2 - 4*(m1P*m4P - m2P*m3P);

    if deltaP >=0
        omegaP = sqrt(((m1P+m4P) + sqrt(deltaP))/2);
        omegaM = sqrt(((m1P+m4P) - sqrt(deltaP))/2);

        if omegaP >= omegaLim break; end
        if omegaM >= 20*2*pi
            modesIndices(modeCounter) = m;
            modeCounter = modeCounter + 1;
        end
        m = m+1;
        S = m*pi/L;
    else
        disp("delta < 0\n");
        m;
    end
end
modesNumber = modeCounter - 1;
modesIndices = modesIndices.';

modesIn = sqrt(2/L)*cos(inPoint*pi*modesIndices/L);
modesOut = sqrt(2/L)*cos(outPoint*pi*modesIndices/L);
modesOutm = sqrt(2/L)*sin(outPoint*pi*modesIndices/L);

%%%%% Computing eigenfrequencies, eigenvectors and other vectors
eigenFreqs = zeros(2*modesNumber,1);
eigenVecs = zeros(2*modesNumber);
eigenVecsInv = zeros(2*modesNumber);
AmatsInv = zeros(2*modesNumber);
DinvRmats = zeros(2*modesNumber);
alphaVec = zeros(2*modesNumber,1);
betaVec = zeros(2*modesNumber,1);

for i = 1:modesNumber
    m = modesIndices(i);
    S = m*pi/L;

    [m1P, m2P, m3P, m4P] = ComputeDispRelMatElems(rho, A, l, mu, kappaSq, E, ni, S);

    [evcs,evls] = eig([m1P,m2P;m3P,m4P]);
    [evls, ind] = sort(diag(evls));
    evcs = evcs(:,ind);

    omegaM = sqrt(evls(1));
    omegaP = sqrt(evls(2));

    eigenFreqs(2*i - 1) = omegaM;
    eigenFreqs(2*i) = omegaP;
    eigenVecs(2*i - 1:2*i, 2*i - 1:2*i) = evcs;

    eigenVecsInv(2*i - 1:2*i, 2*i - 1:2*i) = inv(evcs);

    AmatsInv(2*i - 1,2*i - 1) = 1/(rho*A);
    AmatsInv(2*i,2*i) = 1/(rho*A*(1-l^2*S^2));

    alpha1 = alphaParam/normalAlphaParam;
    alpha2 = (1-alphaParam)/normalAlphaParam;

    betaVec(2*i-1:2*i) = modesIn(i)*eigenVecsInv(2*i-1:2*i,2*i - 1:2*i)*AmatsInv(2*i-1:2*i,2*i - 1:2*i)*[alpha1;alpha2];

    DinvRmats(2*i-1, 2*i - 1) = -2*mu*E*kappaSq*A/l;
    DinvRmats(2*i-1, 2*i) = E*kappaSq*A*((1-mu^2)/l - l*S^2);
    DinvRmats(2*i, 2*i - 1) = ((1-mu^2)/l - l*S^2)*(E*kappaSq*A)/(1+ni+l^2*S^2);
    DinvRmats(2*i-1, 2*i) = 2*mu*(1/l - l*S^2)*(E*kappaSq*A)/(1+ni+l^2*S^2);
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Damping Coefficients
%divided by 2 because here 2sigma0 is considered (Bilbao convention)
sigmaCoeffs = sigma0 + sigma1*((eigenFreqs*pi/L).^2); 

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%% Simulation

%%%% Initializing vectors
etau = zeros(modesNumber,1);
etaw = zeros(modesNumber,1);
etamu = zeros(modesNumber,1);
etamw = zeros(modesNumber,1);

zeta = zeros(2*modesNumber,1);
zetaNext = zeros(2*modesNumber,1);
zetaPrev = zeros(2*modesNumber,1);

output = zeros(1,timeSamples);

A = (2 - eigenFreqs(:,1).^2*k^2)./(1 + 0.5*sigmaCoeffs*k);
B = (0.5*sigmaCoeffs*k - 1)./(1 + 0.5*sigmaCoeffs*k);
C = betaVec./(1+0.5*sigmaCoeffs*k);

S = modesIndices*pi/L;

tic
for n = 1:timeSamples
    exc = excit(n);
    zetaNext = zeta.*A + zetaPrev.*B + C*exc;

    zetaPrev = zeta;
    zeta = zetaNext;

    eta = eigenVecs*zetaNext;

    if ~readForce
        u = modesOut.'*eta(1:2:end);
        w = modesOut.'*eta(2:2:end);
    
        output(n) = (u+w)/2;
    else
        etam = DinvRmats*eta;
        etamu = etam(1:2:end);
        etamw = etam(2:2:end);
        fout = (-2*mu/l)*modesOutm.'*etamu + ((1-mu^2)/l)*modesOutm.'*etamw - l*((pi*modesIndices/L).^2.*modesOutm).'*etamw;
        output(n) = fout;
    end

end
realTimeFrac = toc/durSec

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%% Plotting & Sound
if normalizeOut
    output = output/max(abs(output));
end

figure(1)
plot(output);

if ~strcmp(inputType,'audiofile')
    maxFreq = omegaLim/2/pi; %Hz
    fontSize = 15;
    
    figure(2)
    windowLength = 1024;
    [s,freqs,t,p] = spectrogram(output,blackmanharris(windowLength),floor(windowLength/2),windowLength,SR);
    colormap hot
    mesh([0,t(1:end-1)],freqs,20*log10(p));
    view(2)
    ylim([20,20000]);
    xlim([0,t(end)]);
    yticks([20 10000 20000])
    xticks([0 1 t(end)])
    set(gca,'FontSize',fontSize);
    ylabel("Frequency (Hz)");
    xlabel("Time (s)");
%     title(nLins(nLinType+1));
end

% figure(3)
% outPlot = output(floor(9*timeSamples/10):timeSamples);
% plot((1:length(outPlot))*SR/length(outPlot),abs(fft(outPlot)));
% xlim([2000,3000])

if play
    if playDry
        soundsc(excit,SR);
        pause(durSec)
    end
    diffOut = diff(output);
    soundsc(diffOut,SR);
end

if saveAudio
    outPlay = zeros(floor(timeSamples/(osFac/1)),2);
    lowpass(output,20000,SR);
    for i=1:length(output)
        if ~mod(i,osFac) || mod(i,osFac) == osFac
            index = i/(osFac);
            outPlay(index) = output(i);
        end
    end
    if physicalDamp cCoeff = 'Phys'; end    
    if strcmp(inputType, 'audiofile') 
        split1 = strsplit(audiofileName,'/');
        split2 = strsplit(string(split1(2)),'.');
        fileName = strcat('wet_samples/physicalDamp/','amp',string(amp),'_a',string(aCoeff),'_c',string(cCoeff),'_',string(split2(1)),'_',nLins(nLinType+1),'.wav');
    else
        fileName = strcat('wet_samples/Test/','amp',string(amp),'_a',string(aCoeff),'_c',string(cCoeff),'_',inputType,'_',nLins(nLinType+1),'.wav');
    end
    diffOutPlay = diff(outPlay);
    audiowrite(fileName,diffOutPlay/max(abs(diffOutPlay)),SR/osFac);
end

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

function [excit,SR,k,timeSamples,timeVec,durSec]=ExcitSignal(amp,OSFac,inputDur,tail,excitType,filename)
    switch nargin
        case 5
            if ~inputDur
                disp('Zero input duration');
                return
            end
            SR = OSFac*44100;
            durSec = inputDur + tail;
            timeSamples = durSec*SR;
            k = 1/SR;
            timeVec = (1:timeSamples)*k;
            if excitType == "impulse"
                excit = zeros(timeSamples,1);
                excit(5)=amp*1;
            elseif excitType == "sin"
                excit = zeros(timeSamples,1);
                excit(1:inputDur*SR) = amp*sin(200*2*pi*timeVec(1:SR*inputDur));
            elseif excitType == "square"
                excit = zeros(timeSamples,1);
                excit(1:inputDur*SR) = amp*square(200*2*pi*timeVec(1:SR*inputDur));
            elseif excitType == "sweep"
                startFreq = 200;
                endFreq = 5000;
                chirpLength = timeVec(1:SR*inputDur);
                
                envelope = [linspace(0,amp,SR*inputDur/10), amp*ones(1,SR*8*inputDur/10), linspace(amp,0,SR*inputDur/10)];
                
                forcing = envelope.*chirp(chirpLength,startFreq,inputDur,endFreq,'quadratic',90);
                excit = zeros(1,timeSamples);
                excit(1:length(forcing)) = forcing;
            else 
                disp('Wrong input type');
                return;
            end
        case 6
            if excitType == "audiofile"
                [excit,SR] = audioread(filename);
                if OSFac > 1
                    excit = resample(excit,OSFac,1);
                    SR = SR*OSFac;
                end
                k = 1/SR;
                if inputDur
                    excit = [amp*excit(1:floor(SR*inputDur),1),zeros(1,SR*tail)];
                else
                    excit = [amp*excit,zeros(1,SR*tail)];
                end
                timeSamples = length(excit);
                timeVec = (1:timeSamples)*k;
                durSec = timeSamples/SR;
            end
        otherwise
            disp("wrong input type");
            return
    end
end