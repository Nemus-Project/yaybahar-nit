%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%           Modal membrane simulation
%              Riccardo Russo
%           University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear
close all

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Excitation Signal
%'impulse', 'sin', 'square', 'sweep', 'audiofile'
inputType = 'audiofile'; 
%in case inputType is 'audiofile', specify file name and path
audiofileName = '../Spring/Sounds/Cello/C2_stop.wav';
%amplification factor
amp = 1; 
osFac = 1;
%if inputDur is set to 0 when audiofileName is 'audiofile', the entire file
%is played, otherwise only a portion of it
inputDur = 0;
%tail indicates how long sound should be played after the input signal is
%stopped. The value returned with durSec is inputDur + tail
tail = 4;
[excit,SR,k,timeSamples,timeVec,durSec] = ExcitSignal(amp,osFac,inputDur,tail,inputType,audiofileName);
% [excit,SR,k,timeSamples,timeVec,durSec] = ExcitSignal(amp,osFac,inputDur,tail,inputType);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Custom Parameters
play = true;
playDry = false;
saveAudio = true;

stereo = true;

%Normalize output (before visualization)
normalizeOut = true;

omegaLim = 2/k;

Lx = 0.5;                       %[m] Hor length
Ly = 0.5;                       %[m] Ver lentgh
Lz = 9e-4;%5e-5                    %[m] Thickness

rho = 1400; %Mylar density in kg/m^3

rhoH = rho*Lz;

T = 3000; %1000;                      %[N] Tension

sigma0 = 10; sigma1 = 5e-5;
% sigma0 = 2; sigma1 = 1e-5;

inPoint = [0.52*Lx,0.53*Ly];
outPoint1 = [0.47*Lx,0.62*Ly];
if stereo outPoint2 = [0.73*Lx,0.53*Ly]; end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Computing eigenfrequencies and eigenvectors
Nx = 0; Ny = 0; omegaCurr = 0;
while omegaCurr<omegaLim
    Nx = Nx + 1; 
    omegaCurr = ComputeOmega(Nx,0,T,rho,Lx,Ly,Lz);
end
omegaCurr = 0;
while omegaCurr<omegaLim
    Ny = Ny + 1;
    omegaCurr = ComputeOmega(0,Ny,T,rho,Lx,Ly,Lz);
end

eigenFreqs = zeros(Nx*Ny,3);

for i=1:Nx
    for j=1:Ny
        omegaCurr = ComputeOmega(i,j,T,rho,Lx,Ly,Lz);
        if omegaCurr<omegaLim && omegaCurr>(20*2*pi)
            
            eigenFreqs(j+(i-1)*Ny,1)=omegaCurr;
            eigenFreqs(j+(i-1)*Ny,2) = i;
            eigenFreqs(j+(i-1)*Ny,3) = j;
        end
    end
end
%removing invalid frequencies
eigenFreqs(~any(eigenFreqs,2),:) = [];
eigenFreqs = sortrows(eigenFreqs,1);
modesNumber = size(eigenFreqs,1);

modesIn = ComputeMode(inPoint(1),inPoint(2),eigenFreqs(:,2),...
    eigenFreqs(:,3),Lx,Ly);
modesOut = ComputeMode(outPoint1(1),outPoint1(2),eigenFreqs(:,2),...
    eigenFreqs(:,3),Lx,Ly);
if stereo modesOut2 = ComputeMode(outPoint2(1),outPoint2(2),...
            eigenFreqs(:,2),eigenFreqs(:,3),Lx,Ly); end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Damping Coefficients
%divided by 2 because here 2sigma0 is considered (Bilbao convention)
sigmaCoeffs = sigma0 + sigma1*((eigenFreqs(:,2)*pi/Lx).^2 + (eigenFreqs(:,3)*pi/Ly).^2); 

%% Simulation
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Initializing vectors
q = zeros(modesNumber,1);
qNext = zeros(modesNumber,1);
qPrev = zeros(modesNumber,1);

output = zeros(durSec,1);
outputRight = zeros(durSec,1);

A = (2 - eigenFreqs(:,1).^2*k^2)./(1 + 0.5*sigmaCoeffs*k);
B = (0.5*sigmaCoeffs*k - 1)./(1 + 0.5*sigmaCoeffs*k);
C = k^2*modesIn/(rho*Lz);

tic
for n = 1:timeSamples
    exc = excit(n);

    qNext = q.*A + qPrev.*B + C*exc;

    qPrev = q;
    q = qNext;

    output(n) = modesOut.'*qNext;
    if stereo outputRight(n) = modesOut2.'*qNext; end
end
realTimeFrac = toc/durSec

%% Plot & Sound
figure(1)
plot(output);
if stereo 
    hold on
    plot(outputRight);
end

if normalizeOut
    output = output/max(abs(output));
    if stereo
        outputRight = outputRight/max(abs(outputRight)); 
    end
end

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
end

if stereo 
    output = [output, outputRight];
end

if play
    if playDry
        soundsc(excit,SR);
        pause(durSec)
    end
    if stereo
        OutPlay1 = downsample(output(:,1),osFac);
        OutPlay2 = downsample(output(:,2),osFac);
        diffOutPlay1 = diff(OutPlay1);
        diffOutPlay2 = diff(OutPlay2);

        soundsc([diffOutPlay1,diffOutPlay2],SR);
    else
        OutPlay = downsample(output,osFac);
        diffOutPlay = diff(OutPlay);
        soundsc(diffOutPlay,SR);
    end
end

if saveAudio
%     outPlay = zeros(floor(timeSamples/(osFac/1)),2);
%     lowpass(output,20000,SR);
%     for i=1:length(output)
%         if ~mod(i,osFac) || mod(i,osFac) == osFac
%             index = i/(osFac);
%             outPlay(index) = output(i);
%         end
%     end    
    
    if strcmp(inputType, 'audiofile') 
        split1 = strsplit(audiofileName,'/');
        split2 = strsplit(string(split1(5)),'.');
        fileName = strcat('Sounds/',string(split2(1)),'.wav');
    else
        fileName = strcat('Sounds/Test/',inputType,'.wav');
    end
    if stereo
    %     audiowrite(fileName,diffOutPlay/max(abs(diffOutPlay)),SR/osFac);
       diffOutPlay1 = diffOutPlay1/max(abs(diffOutPlay1)+1e-4);
        diffOutPlay2 = diffOutPlay2/max(abs(diffOutPlay2)+1e-4);
        audiowrite(fileName,[diffOutPlay1,diffOutPlay2],SR/osFac);%,'BitsPerSample',64);
    else
        diffOutPlay = diffOutPlay/max(abs(diffOutPlay)+1e-4);
        audiowrite(fileName,diffOutPlay,SR/osFac);%,'BitsPerSample',64);
    end
end

%% Functions
function o = ComputeOmega(m1,m2,T,rho,Lx,Ly,Lz)
    coeff = (m1.^2*pi^2)/Lx^2 + (m2.^2*pi^2)/Ly^2;
    o = sqrt(T*coeff/(rho*Lz));
end
function m = ComputeMode(xp,yp,m1,m2,Lx,Ly)
    m = sqrt(4/Lx/Ly)*sin(m1*pi*xp/Lx).*sin(m2*pi*yp/Ly); 
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
                    excit = [amp*excit(1:floor(SR*inputDur),1);zeros(SR*tail,1)];
                else
                    excit = [amp*excit;zeros(SR*tail,1)];
                end
                timeSamples = length(excit);
                timeVec = (1:timeSamples)*k;
                durSec = floor(timeSamples/SR);
            end
        otherwise
            disp("wrong input type");
            return
    end
    if tail
        %smoothing out sound tail to avoid impulses
        fadeEnvLength = floor(SR/5);
        fadeEnv = linspace(1,0,fadeEnvLength).';
        fade = fadeEnv.*excit(end - (tail*SR) - fadeEnvLength + 1: end - (tail*SR));
        excit(end-(tail*SR)-fadeEnvLength+1: end - (tail*SR)) = fade;    
    end
end