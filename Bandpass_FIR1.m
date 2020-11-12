f_samp = 330e3;

%Band Edge speifications
fs1 = 30.5e3;
fp1 = 34.5e3;
fp2 = 54.5e3;
fs2 = 58.5e3;

Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end
disp(beta);
N_min = ceil((A-8) / (2.285*0.02424*pi));           %empirical formula for N_min

%Window length for Kaiser Window
n=N_min+18 ;
disp(n);
%Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(0.3424*pi,n) - ideal_lp(0.1969*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass);         %frequency response
disp(bp_ideal);
%magnitude response
[H,f] = freqz(FIR_BandPass,1,2048, f_samp);
plot(f,abs(H))
line([0;18e4],[1.15;1.15],'Color', 'black');
line([30.5e3;30.5e3],[0;1.2],'Color', 'magenta');
line([34.5e3;34.5e3],[0;1.2],'Color', 'magenta');
line([54.5e3;54.5e3],[0;1.2],'Color', 'magenta');
line([58.5e3;58.5e3],[0;1.2],'Color', 'magenta');
line([0;18e4],[0.85;0.85],'Color', 'black');
line([0;18e4],[0.15;0.15],'Color', 'black');
xlabel('Frequency(Hz)');
ylabel('|H(f)|');
grid