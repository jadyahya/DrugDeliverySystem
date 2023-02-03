%% SECTION 2
%% EXPERIMENTAL BODE
%Graphical Values:
%These were determined by hand
format shortE
f = 10^-6*[1, 10, 30, 50, 80, 100, 300, 350, 400, 500, 800, 1000, 2000, 5000, 8000, 10000, 50000, 100000, 500000, 700000, 800000, 900000, 1000000];
M = [3.469, 3.677, 4.881, 6.154, 7.336, 7.723, 6.469, 5.836, 5.216, 4.061, 1.416, 0.214, -2.704, -5.352, -7.373, -8.707, -26.384, -37.464, -64.86, -70.845, -73.191, -75.269, -77.07];
phi = [0, 5.869, 11.344, 11.84, 7.261, 0, -29.840, -33.449 -36.165, -40.594, -45.063, -45.837, -45, -57.037, -72.896, -87.186, -168.470, -209.531, -419.072, -522.554, -566.172, -619.997, -667.449];

%Construction of bode plot:
figure (1)
hold on 
subplot(2,1,1)
semilogx(f,M)
hold on
title('Bode Plot Experimental')
xlabel('frequency (rad/s)')
ylabel('Magnitude (dB)')
grid on
subplot(2,1,2)
semilogx(f,phi)
hold on
xlabel('frequency (rad/s)')
ylabel('Phase (degrees)')
grid on

%Construction of Assymptotes to bode plot:
figure(2)
semilogx(f,M)
hold on
title('Bode Plot Experimental: Magnitude')
xlabel('frequency (rad/s)')
ylabel('Magnitude (dB)')
semilogx([1*10^-6, 3.6*10^-5],[3.469, 3.469],'k-*')
semilogx([3.6*10^-5, 8*10^-5],[3.469,10.136],'k-*')
semilogx([8*10^-5, 2.85*10^-4],[10.136, 10.136],'k-*')
semilogx([2.85*10^-4, 1.5*10^-3],[10.136, -4.56],'k-*')
semilogx([1.5*10^-3, 7.75*10^-3],[-4.56, -4.56],'k-*')
semilogx([7.75*10^-3, 0.0287],[-4.56, -15.56],'k-*')
semilogx([0.0287,0.574],[-15.56, -67],'k-*')
grid on
hold off

%Transport Lag:
s = tf('s');
num2 = 1.4*10^-4*(s+1.5*10^-3)*(s+3.6*10^-5);
den2 = (s+7.75*10^-3)*(s+8*10^-5)*(s+2.87*10^-2)*(s+2.85*10^-4);
hs2 = num2/den2;
[mag_hs2,phase_hs2,freq_hs2] = bode(hs2);
mag_hs2 = 20*log10(mag_hs2);
figure(3)
semilogx(freq_hs2(:), mag_hs2(:))
hold on 
semilogx(f,M)
title('Bode Plot: Magnitude')
xlabel('Frequency (rad/s)')
ylabel('Magnitude (db)')
grid on
legend('Bode Experimental','Bode determined')
hold off

%% Pade Function

%Data from Frequency Response Estimation
magp = [ 3.4770e+00   3.4770e+00   3.4782e+00   3.4840e+00   3.4861e+00   3.5120e+00   3.5206e+00   3.6106e+00   3.6382e+00  3.8888e+00   3.9594e+00   4.5176e+00   4.6603e+00   5.6064e+00   5.8145e+00   6.9078e+00   7.0992e+00   7.8173e+00 7.8901e+00   7.8337e+00   7.7313e+00   6.7484e+00   6.4576e+00   4.6432e+00   4.2111e+00   1.9824e+00   1.5178e+00 -5.5873e-01  -9.4594e-01  -2.5017e+00  -2.7755e+00  -3.8972e+00  -4.1221e+00  -5.2643e+00  -5.5576e+00  -7.2573e+00 -7.7272e+00  -1.0344e+01  -1.1045e+01  -1.4693e+01  -1.5641e+01  -2.0323e+01  -2.1512e+01  -2.7089e+01  -2.8468e+01 -3.4631e+01  -3.6127e+01  -4.2583e+01  -4.4143e+01  -5.0719e+01  -5.2316e+01  -5.8931e+01  -6.0551e+01  -6.7173e+01 -6.8811e+01  -7.5426e+01  -7.7081e+01  -8.3684e+01  -9.1943e+01];
phasep = [-4.5000e+00  -4.4949e+00  -4.4003e+00  -3.9842e+00  -3.8482e+00  -2.5798e+00  -2.2466e+00   2.1218e-01   7.6794e-01  4.2115e+00   4.8958e+00   8.3844e+00   8.9302e+00   1.0420e+01   1.0332e+01   7.3298e+00   6.2346e+00  -1.7356e+00 -3.6776e+00  -1.4689e+01  -1.7004e+01  -2.8376e+01  -3.0504e+01  -3.9349e+01  -4.0721e+01  -4.5037e+01  -4.5444e+01 -4.5712e+01  -4.5538e+01  -4.4863e+01  -4.4931e+01  -4.7278e+01  -4.8343e+01  -5.6181e+01  -5.8548e+01  -7.2009e+01 -7.5506e+01  -9.3012e+01  -9.7249e+01  -1.1713e+02  -1.2182e+02  -1.4291e+02  -1.4781e+02  -1.6920e+02  -1.7420e+02  -1.9615e+02  -2.0148e+02  -2.2520e+02  -2.3111e+02  -2.5648e+02  -2.6259e+02  -2.8660e+02  -2.9198e+02  -3.1120e+02  -3.1523e+02  -3.2875e+02  -3.3149e+02  -3.4033e+02  -3.4772e+02];
freqp = [1.0000e-06   1.6103e-06   2.5929e-06   3.8857e-06   4.1753e-06   6.2513e-06   6.7234e-06   1.0057e-05   1.0826e-05 1.6180e-05   1.7433e-05   2.6031e-05   2.8072e-05   4.1879e-05   4.5204e-05   6.7375e-05   7.2790e-05   1.0839e-04 1.1721e-04   1.7439e-04   1.8874e-04   2.8055e-04   3.0392e-04   4.5136e-04   4.8939e-04   7.2615e-04   7.8805e-04 1.1682e-03   1.2690e-03   1.8795e-03   2.0434e-03   3.0237e-03   3.2903e-03   4.8646e-03   5.2983e-03   7.8263e-03 8.5317e-03   1.2591e-02   1.3738e-02   2.0257e-02   2.2122e-02   3.2589e-02   3.5622e-02   5.2430e-02   5.7362e-02 8.4350e-02   9.2367e-02   1.3570e-01   1.4874e-01   2.1832e-01   2.3950e-01   3.5124e-01   3.8566e-01   5.6507e-01 6.2102e-01   9.0910e-01   1.0000e+00   1.4626e+00   2.3530e+00];

%Construction of bode plot:
figure (4)
hold on 
subplot(2,1,1)
semilogx(freqp,magp)
hold on
title('Bode Plot Experimental')
xlabel('frequency (rad/s)')
ylabel('Magnitude (dB)')
subplot(2,1,2)
semilogx(freqp,phasep)
hold on
xlabel('frequency (rad/s)')
ylabel('Phase (degrees)')

%Assympote Construction:
figure(5)
semilogx(freqp,magp)
hold on
title('Bode Plot Experimental (Pade): Magnitude')
xlabel('frequency (rad/s)')
ylabel('Magnitude (dB)')
semilogx([1*10^-6, 3.6*10^-5],[3.469, 3.469],'k-*')
semilogx([3.6*10^-5, 8*10^-5],[3.469,10.136],'k-*')
semilogx([8*10^-5, 2.85*10^-4],[10.136, 10.136],'k-*')
semilogx([2.85*10^-4, 1.5*10^-3],[10.136, -4.56],'k-*')
semilogx([1.5*10^-3, 7.75*10^-3],[-4.56, -4.56],'k-*')
semilogx([7.75*10^-3, 0.0287],[-4.56, -15.56],'k-*')
semilogx([0.0287,  0.574],[-15.56, -67],'k-*')
grid on
hold off

%Non Minimum Term:
figure (6)
semilogx(freqp,phasep)
hold on
xlabel('frequency (rad/s)')
ylabel('Phase (degrees)')
title('Bode for Pade Approximation')
semilogx(freq_hs2(:),phase_hs2(:))
grid on
legend('Bode Experimental','Bode determined')
hold off

%Determining wnm:
figure(7)
semilogx(freqp,phasep)
hold on
xlabel('frequency (rad/s)')
ylabel('Phase (degrees)')
title('Bode for Pade Approximation')
semilogx([10^-3,10],[-180,-180],'r--')
semilogx([10^-3,10],[-360,-360],'r--')
semilogx([10^-3,10],[-270,-270],'r--')
grid on
hold off

%Checking results
gpade2_exp = num2*(-s+0.27)/den2/(s+0.27);
[m,p,f] = bode(gpade2_exp);
figure(8)
semilogx(f(:),p(:)-360)
hold on
semilogx(freqp,phasep)
xlabel('frequency (rad/s)')
ylabel('Phase (degrees)')
title('Bode for Pade Approximation')
legend('Bode Experimental','Bode determined')
grid on

%% Pade Approximation vs Exact
sys = exp(-8.5*s);
sysx = pade(sys,1);
gexact = num2*sys/den2;
gpade2 = num2*sysx/den2;
figure(9)
bode(gpade2)
hold on
bode(gexact)
bode(hs2)
legend('Bode Pade','Bode Exact','Bode Without Transport Lag')
grid on
hold off

%Divergence Frequency
w = 0.01:0.0001:0.1;
g = 0.95-2*atan(w./0.2353)./(8.5*w);
figure (10)
plot(w,g)
hold on 
xlabel('Frequency (rad/sec)')
ylabel('g(w)')
title('g vs Frequency')
grid on

%Higher orders
%Bode
figure(11)
bode(gexact)
hold on
for l =1:2:9
    sysx = pade(sys,l);
    gpadeHO = num2*sysx/den2;
    bode(gpadeHO)
end
legend('Bode Exact','Bode Pade Order 1','Bode Pade Order 3','Bode Pade Order 5','Bode Pade Order 7', 'Bode Pade Order 9')
grid on
hold off

%Root Locus
j =[12,13,14]
k=1;
for l =[1,4,9]
    sysx = pade(sys,l);
    gpadeHO = sysx*num2/den2;
    figure(j(k))
    rlocus(gpadeHO)
    k=k+1;
end

%Linearize
[A,B,C,D]=linmod('SimulinkModelLM');

%% Section III
%Function Gpade3
s = tf('s');
num3 = exp(-8.5*s)*0.17*(s+0.05)*(s+0.08);
den3 = (s+0.1)*(s+0.03)*(s+0.004)*(s+5*10^-4);
g3 = num3/den3;

%Simulating the system to a step response
g3_feedback = feedback(g3,1);
figure(15)
step(g3_feedback)

% LAG LEAD NETWORK
% _D is used to indicate values at the desired locations
% _B is used for Bisection method

%Finding Angle of Deficiency
s3 = -3.33*10^-3+2.9864*10^-3*i;
num3_D = (-s3+0.2353)*0.17*(s3+0.05)*(s3+0.08);
den3_D = (s3+0.2353)*(s3+0.1)*(s3+0.03)*(s3+0.004)*(s3+5*10^-4);
gpade3_D = num3_D/den3_D;
phase_gpade3_D=angle(gpade3_D);
deficiency = (pi-phase_gpade3_D)*180/pi;

% METHOD 1: Bisection

%Lead Network Design
glead3 = (s+3.52*10^(-3))/(s+5.68*10^(-3)); %without compensator gain
glead3_D = evalfr(glead3,s3);
phase_glead3_B_D = angle(glead3_D)*180/pi;
Kc = real(-1/(glead3_D*gpade3_D));
glead3 = Kc*glead3;

%Lag Network Design
T2 = 1:1:11000;
beta_B = 10.55;
mag_lag_B =[];
phase_lag_B = [];
for l = 1:length(T2)
    glag_D = (s3+1./T2(l))./(s3+1./(beta_B*T2(l)));
    mag_lag_B = [mag_lag_B abs(glag_D)];
    phase_lag_B = [phase_lag_B angle(glag_D)];
end

%We look over the matricies mag_lag_B and phase_lag_B to determine suitable
%values of T2 to satisfy the angle and magnitude conditions

T2_B1=2000;

%Test: Bisection
glag3 =(s+1/T2_B1)/(s+1/(beta_B*T2_B1));
gcompensated_3 = glag3*glead3*g3;
gcompensated_fb = feedback(gcompensated_3,1);
pole(gcompensated_fb)
zero(gcompensated_fb)

%% PID Controller
s = tf('s');
num3 = (-s+0.2353)*0.17*(s+0.05)*(s+0.08);
den3 = (s+0.2353)*(s+0.1)*(s+0.03)*(s+0.004)*(s+5*10^-4);
gpade3 = num3/den3

%angle of deficiency is equal to 34.536 
z1=0.0004;
angle3=(deficiency*(pi/180)+angle(s3)-angle(s3+z1))*180/pi;
z2=((imag(s3)/tand(angle3))-real(s3));

GPID3s3=((s3+z1)*(s3+z2))/s3;
K3=(1/(abs(gpade3_D)*abs(GPID3s3)))

GPID3=K3*((s+z1)*(s+z2))/s;
Gcompensated_PID3=GPID3*g3;
Compensated_PID_fb3=feedback(Gcompensated_PID3,1); %%Closed loop feedback for
%the system

pole(Compensated_PID_fb3)
zero(Compensated_PID_fb3)

%% Verification
% Simulation to step input
figure(16)
hold on
step(gcompensated_fb)
grid on
hold off
stepinfo(gcompensated_fb)

stepinfo(Compensated_PID_fb3)
figure (17)
hold on
step(Compensated_PID_fb3)
grid on
hold off

%Comparison to second order
%_so means second order 
den_so = conv([1 -s3],[1 -conj(s3)]);
num_so = den_so(3);
g_so = tf(num_so,den_so);

figure(18)
step(gcompensated_fb,2000)
hold on
step(Compensated_PID_fb3,2000)
step(g_so,2000)
legend('Lag Lead Compensated', 'PID Compensated', 'Second Order')
grid on
hold off

% Phase and Gain Margins
[Gm_LL,Pm_LL,Wcg_LL,Wcp_LL] = margin(gcompensated_3)
[Gm_PID,Pm_PID,Wcg_PID,Wcp_PID] = margin(Gcompensated_PID3)



