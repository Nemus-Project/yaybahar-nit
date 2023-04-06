clear
close all

[data_bow,SR] = audioread('BowedString\Sounds\C2_Stop.wav');
data_spring = audioread('Spring\Sounds\Cello\C2_Stop.wav');
data_membr = audioread('Membrane\Sounds\C2_Stop_LD.wav');

timeVecBow = (1:length(data_bow))/SR;
timeVecSpring = (1:SR*8)/SR;
timeVecMembr = (1:SR*8)/SR;

data_spring = data_spring(1:SR*8);
data_membr = data_membr(1:SR*8,1);

fontSize = 11;
lineWidth = 2.5;

maxFreqHz = 10000;
maxFreq = maxFreqHz/2/pi; %Hz
windowLength = 2048;

figure(1)

% subplot(2,1,1)
% %subplot(2,3,1)
% plot(timeVecBow,data_bow); 
% set(gca,'FontSize',fontSize);
% xlabel('t(s)');
% 
% 
% subplot(2,1,2)
    subplot(3,1,1)
[s,freqs,t,p] = spectrogram(data_bow(42650:111850),parzenwin(windowLength),floor(windowLength*0.95),windowLength,SR);
colormap hot
mesh([0,t(1:end-1)],freqs,20*log10(p));
view(2)
ylim([20,maxFreqHz]);
xlim([0,t(end-1)]);
%     yticks([20 10000 20000])
%     xticks([0 1 2])
set(gca,'FontSize',fontSize);
xlabel('t (s)');
ylabel('f (Hz)');
title('(a) Bowed String');

% figure(2)
% subplot(2,1,1)
% %subplot(2,3,2)
% plot(timeVecSpring,data_spring); 
% set(gca,'FontSize',fontSize);
% xlabel('t(s)');
% title('Spring');
% 
% subplot(2,1,2)
subplot(3,1,2)
 [s,freqs,t,p] = spectrogram(data_spring(42650:111850),parzenwin(windowLength),floor(windowLength*0.95),windowLength,SR);
colormap hot
mesh([0,t(1:end-1)],freqs,20*log10(p));
view(2)
ylim([20,maxFreqHz]);
xlim([0,t(end-1)]);
%     yticks([20 10000 20000])
%     xticks([0 1 2])
set(gca,'FontSize',fontSize);
xlabel('t (s)');
ylabel('f (Hz)');
title('(b) Spring');

% figure(3)
% subplot(2,1,1)
% %subplot(2,3,3)
% plot(timeVecMembr,data_membr); 
% set(gca,'FontSize',fontSize);
% xlabel('t(s)');
% title('Membrane')
% 
% subplot(2,1,2)
subplot(3,1,3)
[s,freqs,t,p] = spectrogram(data_membr(42650:111850),parzenwin(windowLength),floor(windowLength*0.95),windowLength,SR);
colormap hot
mesh([0,t(1:end-1)],freqs,20*log10(p));
view(2)
ylim([20,maxFreqHz]);
xlim([0,t(end-1)]);
%     yticks([20 10000 20000])
%     xticks([0 1 2])
set(gca,'FontSize',fontSize);
xlabel('t (s)');
ylabel('f (Hz)');
title('(c) Membrane')





