%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 4;

% Open CLHP Poles of the Chebyshev Polynomial of order 4
p1 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p2 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p3 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);
p4 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);        
disp(p1);


%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
den = conv(n1,n2);          %multiply n1 and n2, which are the two quadratic factors in the denominator
num = [den(5)*sqrt(1/(1+epsilon*epsilon))];        % even order, DC Gain set as 1/(1+ epsilon^2)^0.5
%disp(num);
%disp(den);
%Band Edge speifications
fp1 = 28.9;
fs1 = 32.9;
fs2 = 52.9;
fp2 = 56.9;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 260;
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);
%disp(ws1);
%disp(ws2);
%disp(wp1);
%disp(wp2);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;
%disp(W0);
%disp(B);
%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bsf(s) = analog_lpf((B*s)/(s*s +W0*W0));     %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BSF
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;
%disp(ds);
%disp(ns);
%coeffs of discrete BPF
[nz, dz] = numden(discrete_bsf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                       %frequency response in dB

%disp(nz);
%disp(dz);

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 260e3);
plot(f,abs(H))
line([0;14e4],[0.85;0.85],'Color', 'green','linestyle','--');
line([28.9e3;28.9e3],[0;1.2],'Color', 'black');
line([32.9e3;32.9e3],[0;1.2],'Color', 'black');
line([52.9e3;52.9e3],[0;1.2],'Color', 'black');
line([56.9e3;56.9e3],[0;1.2],'Color', 'black');
xlabel('Frequency(Hz)');
ylabel('|H(f)|');
grid