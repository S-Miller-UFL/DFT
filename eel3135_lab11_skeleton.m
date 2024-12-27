%% QUESTION 2 
% DO NOT REMOVE THE LINE BELOW
% MAKE SURE 'eel3135_lab11_comment.m' IS IN THE SAME DIRECTORY AS THIS FILE
clear; close all; clc;
type('eel3135_lab11_comment.m')

%% Question 2
% This question will cover the creation of their own DFT function by
% modifying their DTFT function
% Keep in mind that the DFT can be thought of as a sampled version of the
% DTFT

n=0:59;
x = 0.75 + cos(pi*n/20) + cos(pi*n/15) + cos(pi*n + 2*pi/3);
w_DTFT = linspace(0, 2*pi-pi/5000, 10000);

% NOTE: USE THE FOLLOWING COMMENTED LINE FOR PLOTTING THE DFT ATOP THE DTFT
% (YOU NEED TO DEFINE w_DFT), THIS WILL MAKE THE PLOTS EASIER TO INTERPRET
%plot(w_DTFT,abs(X_DTFT)); 
%hold on; 
%plot(w_DFT,abs(X_DFT),'.', 'markersize', 10); 
%hold off;

%% Question 2a
% Take your DTFT function from previous labs, and modify it to output the
% DFT given an input X
%% Question 2b
% Given x and w_DTFT above, plot the magnitude of the DTFT. Be sure to 
% label axes with units and title the plot.
X_DTFT = DTFT(x,w_DTFT);
figure("Name","2b");
plot(w_DTFT,abs(X_DTFT));
title('Magnitude Response of x')
ylabel('Magnitude [rad]')
xlabel('Normalized Angular Frequency [rad/s]')
%% Question 2c
% Now use stem to plot the magnitude of the 60-length DFT of x. Use hold on
% and hold off to plot the DFT on top of the DTFT plot
npnt = 60;
w_DFT = (0:npnt-1)*(2*pi)/npnt;
X_DFT = DFT(x);
X_DTFT = DTFT(x,w_DTFT);
figure("Name","2c")
plot(w_DTFT,abs(X_DTFT)); 
hold on; 
plot(w_DFT,abs(X_DFT),'.', 'markersize', 10); 
hold off;
title('60 point DFT of x')
ylabel('Magnitude [rad]')
xlabel('Normalized Angular Frequency [rad/s]')
%% Question 2d
% Now use stem to plot the magnitude of the 55-length DFT (i.e. remove the
% last 5 values) of x. Use hold on and hold off to plot the DFT on top of
% the DTFT plot
npnt = 55;
w_DFT = (0:npnt-1)*(2*pi)/npnt;
X_DFT = DFT(x(1:55));
X_DTFT = DTFT(x(1:55),w_DTFT);
figure("Name","2d")
plot(w_DTFT,abs(X_DTFT)); 
hold on; 
plot(w_DFT,abs(X_DFT),'.', 'markersize', 10); 
hold off;
title('55 point DFT of x')
ylabel('Magnitude [rad]')
xlabel('Normalized Angular Frequency [rad/s]')
%% Question 2e
% Now plot the magnitude of the 65-length DFT (i.e. include 5 
% zeros values) of x. Use hold on and hold off to plot the DFT on top of
% the DTFT plot
npnt = 65;
x65 = zeros(1,65);
x65(1:length(x)) = x;
x = x65;
w_DFT = (0:npnt-1)*(2*pi)/npnt;
X_DFT = DFT(x);
X_DTFT = DTFT(x,w_DTFT);
figure("Name","2e")
plot(w_DTFT,abs(X_DTFT)); 
hold on; 
plot(w_DFT,abs(X_DFT),'.', 'markersize', 10); 
hold off;
title('65 point DFT of x')
ylabel('Magnitude [rad]')
xlabel('Normalized Angular Frequency [rad/s]')

%% Question 2f
% Now use stem to plot the magnitude of the 200-length DFT of x. Use hold 
% on and hold off to plot the DFT on top of the DTFT plot
npnt = 200;
x200 = zeros(1,200);
x200(1:length(x)) = x;
x = x200;
w_DFT = (0:npnt-1)*(2*pi)/npnt;
X_DFT = DFT(x);
X_DTFT = DTFT(x,w_DTFT);
figure("Name","2f")
plot(w_DTFT,abs(X_DTFT)); 
hold on; 
plot(w_DFT,abs(X_DFT),'.', 'markersize', 10); 
hold off;
title('200 point DFT of x')
ylabel('Magnitude [rad]')
xlabel('Normalized Angular Frequency [rad/s]')
%% Question 2g
% Answer in your comments:  Based on the last several questions, what is
% the relationship between the DTFT and the DFT? Under what conditions will
% the theoretical DTFT (i.e. when w_hat is continuous) and DFT have the
% same result?

% the DFT is a sampling of the DTFT
% the DFT and the DTFT will yield the same result when the DFT's sample
% number is high enough.

%% Question 3
% This question will cover the creation of your own DFT function by
% modifying your DTFT function. Keep in mind that the DFT can be thought of
% as a sampled version of the DTFT. In this problem. consider the following
% discrete-time signal: x[n] = u[n-10] - u[n-45]
clear all; clc;
x = [zeros(1,10) ones(1,35)];

%% Question 3a
% Modify the DFT question from Question 2 to create an inverse DFT function
% x = IDFT(x), which inputs a DFT-transformed signal X and outputs the
% time-domain signal x. Hint: the inverse DFT is defined in the lab
% document.

%% Question 3b
% Use stem to plot x[n] for n from 0 to 99, inclusive
n = 0:99;
xn = zeros(1,100);
xn(1:length(x)) = x;
figure("Name","3b")
stem(n,xn); 
title('x[n]')
ylabel('Magnitude')
xlabel('n')
%% Question 3c
% Use conv to compute y[n] = (x[n] convolved with x[n]) and then use stem
% to plot the result for n from 0 to 100, inclusive
n = 0:100;
xn = zeros(1,101);
xn(1:length(x)) = x;
y = conv(xn,xn);
y = y(1:length(xn));
figure("Name","3c")
stem(n,y); 
title('x[n] convoluted with x[n]')
ylabel('Magnitude')
xlabel('n')
%% Question 3d
% Use N = 100 length DFT function to compute y[n] = (x[n] 
% convolved with x[n]) in the frequency domain, then convert back to the
% time domain. Use stem to plot the result.
n = 0:99;
xn = zeros(1,100);
xn(1:length(x)) = x;

X_DFT = DFT(xn);
Y_DFT = X_DFT .* X_DFT;
Y_DFT = Y_DFT(1:length(n));
y = IDFT(Y_DFT);
figure("Name","3d")
stem(n,y); 
title('Impulse Response of conv(x[n],x[n]) through 100 point IDFT')
ylabel('Magnitude')
xlabel('n')
%% Question 3e
% Use N = 60, and repeat 3d
n = 0:59;
xn = zeros(1,60);
xn(1:length(x)) = x;

X_DFT = DFT(xn);
Y_DFT = X_DFT .* X_DFT;
Y_DFT = Y_DFT(1:length(n));
y = IDFT(Y_DFT);
figure("Name","3e")
stem(n,y); 
title('Impulse Response of conv(x[n],x[n]) through 60 point IDFT')
ylabel('Magnitude')
xlabel('n')
%% Question 3f
% Answer in your comments: What is the difference in the last two
% solutions? Why does this difference exist?

% the difference is that their shapes are completely different from
% eachother
% the reason for this difference is because of the number of samples taken.
% the more samples taken, the more accurate the DFT

%% Question 4
% This question will focus on computation time differences between the DFT
% and the FFT. Choose a song at least 3 minutes long to use in this
% problem, and include it in your submission. Load it into MATLAB using
% audioread. Note that mose audio files will be stereo, so you need to make
% Sure that you only use one column of audio data for this part of the lab
% This site has a large archive of free music that you can choose from:
% https://freemusicarchive.org/static
clear all; clc;
%song name is "feint operation" from virtua cop 3 by SEGA corporation
[x,fs] = audioread('vc3.wav');
x = x(:,1);
%soundsc(x,fs);
%% Question 4a
% Use DFT to plot the magnitude of the DFT of only the first 10000 samples
% of the audio. Use tic and toc to measure the length of time it takes to 
% compute the DFT. Display the result with the disp function.
npnt = 10000;
w_DFT = (0:npnt-1)*(2*pi)/npnt;
x10k = x(1:10000);
tic;
X_DFT = DFT(x10k);
toc;
figure("Name","4a")
plot(w_DFT,abs(X_DFT)); 
disp(toc-tic);
title('Magnitude response of audio using DFT')
ylabel('Magnitude')
xlabel('Normalized Angular Frequency [rad/s]')

%% Question 4b
% Use fft (read help fft) to plot the magnitude of the DFT of only the
% first 10000 samles of the audio. Use tic and toc to measure the length of
% time it takes to compute the DFT. Display the result with the disp 
% function.
npnt = 10000;
w_DFT = (0:npnt-1)*(2*pi)/npnt;
x10k = x(1:10000);
tic;
X_FFT = fft(x10k,10000);
toc;
figure("Name","4b")
plot(w_DFT,abs(X_FFT)); 
disp(toc-tic);
title('Magnitude response of audio using FFT')
ylabel('Magnitude')
xlabel('Normalized Angular Frequency [rad/s]')
%% Question 4c
% Answer in your comments: Are there any differences in the results? If so,
% why?

%there are no differences

%% Question 4d
% Answer in your comments: How much faster is the FFT algorithm compared
% with the DFT in this scenario?

%DFT took 1.344768 seconds to compute
%FFT took .004799 seconds to compute
%FFT is 28,021% faster than DFT

%% Question 4e
% Now use fft to plot the magnitude of the DFT of the ENTIRE audio signal.
npnt = length(x);
w_DFT = (0:npnt-1)*(2*pi)/npnt;
tic;
X_FFT = fft(x,length(x));
toc;
figure("Name","4e")
plot(w_DFT,abs(X_FFT)); 
title('Magnitude response of audio')
ylabel('Magnitude')
xlabel('Normalized Angular Frequency [rad/s]')
disp(toc-tic);
%% Functions provided for the lab
function H = DTFT(x,w)
% DTFT(X,W)  compute the Discrete-time Fourier Transform of signal X
% acroess frequencies defined by W. 

    H = zeros(1, length(w));
    for nn = 1:length(x)
        H = H + x(nn).*exp(-1j*w.*(nn-1));
    end
    
end

function X = DFT(x)
% DFT(x)  compute the N-point Discrete Fourier Transform of signal x  
% Where N is the length of signal x
    N = length(x);
    w = ((2*pi)/N).*(0:N-1); %w = (2pi/N) * k where k goes from 0 to length(x)
    X = zeros(1, length(w)); 
    for nn = 1:length(x)    
        X = X + x(nn).*exp(-j*w.*(nn-1));    
    end
end

function x = IDFT(X)
% IDFT(x)  compute the N-point Inverse Discrete Fourier Transform of signal
% X where N is the length of signal X
    N = length(X);
    w = ((2*pi)/N).*(0:N-1); %w = (2pi/N) * k where k goes from 0 to length(x)
    x = zeros(1, length(w));    
    for nn = 1:length(x)    
        x = x + X(nn).*exp(j*w.*(nn-1));
    end
    x=x/N;
end