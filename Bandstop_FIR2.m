f_samp = 260e3;

%Band Edge speifications
fs1 = 32.9e3;
fp1 = 28.9e3;
fp2 = 56.9e3;
fs2 = 52.9e3;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-8) / (2.285*0.03077*pi));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min +16;
%Ideal bandstop impulse response of length "n"
disp(n);
bs_ideal =  ideal_lp(pi,n) -ideal_lp(0.4223*pi,n) + ideal_lp(0.23765*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response
disp(bs_ideal);
%magnitude response
[H,f] = freqz(FIR_BandStop,1,2048, f_samp);
plot(f,abs(H))
line([0;14e4],[1.15;1.15],'Color', 'black');
line([28.9e3;28.9e3],[0;1.2],'Color', 'magenta');
line([32.9e3;32.9e3],[0;1.2],'Color', 'magenta');
line([56.9e3;56.9e3],[0;1.2],'Color', 'magenta');
line([52.9e3;52.9e3],[0;1.2],'Color', 'magenta');
line([0;14e4],[0.85;0.85],'Color', 'black');
line([0;14e4],[0.15;0.15],'Color', 'black');
xlabel('Frequency(Hz)');
ylabel('|H(f)');
grid
