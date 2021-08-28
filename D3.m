h = 1.5; %inches
l = 11.2; %inches
w = 2500 * (pi / 30); %radians per second
r = 10;
d = 0.8746;
radius = d / 2;
weight = 3;
mass = weight / 32.17; %lbf / ft/s^2
bore = 3.5;
Theta = 0:1:720; %degrees
Phi = asind(h .* sind(Theta) ./ l); %degrees
Area = pi * (bore / 2) ^ 2;
Vtdc = 1.023;
vert = [391.2467 589.6080];
horiz = [360 360];
I = pi / 64 * d .^ 4;
J = pi / 32 * d .^ 4;

%Material Properties
Sy = 90.65; %kpsi
Sut = 116.03; %kpsi

Point1 = 0;
Point2 = 45;
Point3 = 90;
Point4 = 135;
Point5 = 180;

A = -h * w^2 * (cosd(Theta) + (h * cosd(2 * Theta) ./ (l * cosd(Phi)) + (h .^ 3 * (sind(2 * Theta)) .^ 2) ./ (4 * l .^ 3 * (cosd(Phi)) .^ 3))); 
Acceleration = A ./ 12; %ft/s^2

Vtheta = Vtdc + ((h + l) - (h .* cosd(Theta) + l .* cosd(Phi))) .* pi ./ 4 * bore .^ 2;

Pcomp = (400 .* Vtdc .^ 1.3) ./ Vtheta(180:359) .^ 2;
Pexp = (600 .* Vtdc .^ 1.3) ./ Vtheta(361:540) .^ 2;
PcompZero = zeros(1,181) + Pcomp(1);
PexpZero = zeros(1,180) + Pexp(180);

Pressure = [PcompZero Pcomp Pexp PexpZero];

Force = Pressure .* Area;
Fcr = -(mass .* Acceleration + Force) ./ cosd(Phi);

Ox = -Fcr .* cosd(Phi);
Oy = Fcr .* sind(Phi);
T = h .* Fcr .* (cosd(Theta) .* sind(Phi) + sind(Theta) .* cosd(Phi));

Ax = -Ox;
Ay = -Oy;

Mres = 1 / 4 .* sqrt(Ax .^ 2 + Ay .^ 2); %in-lbf
Tres = -T; %in-lbf
Beta = atan2d(Ay, Ax);

c1 = abs(radius .* sind(Beta + Theta - Point1));
c2 = abs(radius .* sind(Beta + Theta - Point2));
c3 = abs(radius .* sind(Beta + Theta - Point3));
c4 = abs(radius .* sind(Beta + Theta - Point4));
c5 = abs(radius .* sind(Beta + Theta - Point5));

angleAdd = Theta + Beta;

Sigma1 = Mres .* c1 ./ I;
Sigma2 = Mres .* c2 ./ I;
Sigma3 = Mres .* c3 ./ I;
Sigma4 = Mres .* c4 ./ I;
Sigma5 = Mres .* c5 ./ I;

Tau1 = Tres .* radius ./ J;
Tau2 = Tres .* radius ./ J;
Tau3 = Tres .* radius ./ J;
Tau4 = Tres .* radius ./ J;
Tau5 = Tres .* radius ./ J;

% for i = 1:720
%     if angleAdd(i) > 90 && angleAdd(i) < 630 && (angleAdd(i) < 265 || angleAdd(i) > 450)
%         Sigma3(i) = -1 .* Sigma3(i);
%     end
% end
% 
% for i = 1:720
%     if angleAdd(i) > 90 && angleAdd(i) < 630 && (angleAdd(i) < 265 || angleAdd(i) > 450)
%         Tau3(i) = -1 .* Tau3(i);
%     end
% end

SigmaA = max(Sigma3) - min(Sigma3);
SigmaM = median(Sigma3);

TauA = max(Tau3) - min(Tau3);
TauM = median(Tau3);

SigmaPrimeA = (SigmaA ^ 2 + 3 * TauA ^ 2) ^ (1/2) / 1000; %kpsi
SigmaPrimeM = (SigmaM ^ 2 + 3 * TauM ^ 2) ^ (1/2) / 1000; %kpsi
SigmaPrimeMax = ((SigmaM + SigmaA) ^ 2 + 3 * (TauM + TauA) ^ 2) ^ (1/2) / 1000; %kpsi

ny = Sy / SigmaPrimeMax;

SePrime = 0.5 * Sut; %kpsi
% Ka = 39.9 * Sut ^ (-0.995); %As forged
Ka = 2.7 * Sut ^ (-0.265); %Machined
Kb = 0.879 * d ^ (-0.107);
Ke = 0.897;

Se = SePrime * Ka * Kb * Ke; %kpsi

modGoodman = ((SigmaPrimeA / Se) + (SigmaPrimeM / Sut)) ^ (-1);

figure(1);
plot(Theta, Phi);
grid on;
title('Angle Phi Vs. Crank Angle Theta');
xlabel('Crank Angle Theta (Degrees)');
ylabel('Angle Phi (Degrees)');

figure(2);
plot(Theta, A);
grid on;
title('Piston Acceleration Vs. Crank Angle Theta');
xlabel('Crank Angle Theta (Degrees)');
ylabel('Acceleration (in/s^2)');

figure(3)
plot(Theta, Vtheta);
grid on;
title('Volume vs Crank Angle Theta');
xlabel('Crank Angle Theta');
ylabel('Cylinder Volume in^3');

figure(4)
plot(Theta, Pressure);
grid on;
title('Cylinder Pressure vs Crank Angle Theta');
xlabel('Crank Angle Theta');
ylabel('Cylinder Pressure psi');

figure(5)
plot(Theta, Force);
grid on;
title('Force from Pressure vs Crank Angle Theta');
xlabel('Crank Angle Theta');
ylabel('Force lbf');

figure(6)
plot(Theta, Fcr);
grid on;
title('Force from Connecting Rod vs Crank Angle Theta');
xlabel('Crank Angle Theta');
ylabel('Force lbf');

figure(7)
plot(Theta, Ax);
grid on;
title('O_x vs Crank Angle Theta');
xlabel('Crank Angle Theta');
ylabel('Force lbf');

figure(8)
plot(Theta, Ay);
grid on;
title('O_y vs Crank Angle Theta');
xlabel('Crank Angle Theta');
ylabel('Force lbf');

figure(9)
plot(Theta, T);
grid on;
title('Torque vs Crank Angle Theta');
xlabel('Crank Angle Theta');
ylabel('Torque in-lbf');

figure(10)
hold on;
plot(Theta, Sigma3);
yline(median(Sigma3))
yline(max(Sigma3))
yline(min(Sigma3))
grid on;
title('Point 3 Normal Stress vs Crank Angle Theta');
xlabel('Crank Angle Theta');
ylabel('Normal Stress Psi');

figure(11)
hold on;
plot(Theta, Tau3);
yline(median(Tau3))
yline(max(Tau3))
yline(min(Tau3))
grid on;
title('Point 3 Torsional Shear Stress vs Crank Angle Theta');
xlabel('Crank Angle Theta');
ylabel('Torsional Shear Stress Psi');