clear; clc;
% This script is emulating the signal processing of MXR Distortion+ guitar pedal in digital


%sample rate and period
Fs = 44100;
Ts = 1/Fs;

% Op-amp stage

% component values of the DK substitution circuit
C1 = 10e-9;
R1 = Ts/(2*C1); 

C2 = 47e-9;
R2 = Ts/(2*C2);

R3 = 10e3; 
R4 = 1e6;
R5 = 4.7e3;
R8 = 1e6;

% "Distortion" RV-POT position from 0 to 1
% 1-Mega Ohm reverse-log pot
pot = 1; % from 0 to 1
k = 6; % sensitivity factor of the reverse-log pot
R6 = (exp(-k * pot) - exp(-k)) / (1 - exp(-k)) * 1000000; 
Rn = R5 + R6; % from 4.7k + 0  to  4.7k + 1Meg

% transfer function's G-values
Ga = 1 + (R3/R1);  
Gb = 1 + (R2/Rn);
Gh = 1 + (R4/(Gb*Rn)); 
Gx = (1/R8) + (1/(Ga*R1));

% transfer function's coefficients
b0 = Gh/(Ga*R1*Gx);
b1 = ((R3/(Ga*R1))-1)*(Gh/Gx);
b2 = (-R2*R4)/(Gb*Rn); 


% input signal: audio file
[input, Fs] = audioread('testaudio.wav'); 
N = length(input);

% output signal setting
y = zeros(N,1); 

% Initial state value for DK substitution circuit of capacitors
x1 = 0; x2 = 0;

% opamp-stage sample-by-sample processing
for n = 1:N
    Vin = input(n,1); % discrete input signal
    Vout = b0 * Vin + b1 * x1 + b2 * x2;

    % some values included in state-update equations
    Vx = (1/Gh)*Vout + ((R2*R4)/(Gb*Gh*Rn))*x2;
    VR1 = (Vin-Vx+(R3*x1))/Ga; 

    VRn = (Vx-(R2*x2))/Gb; 
    VR2 = (R2/Rn)*VRn + (R2*x2);
  
    % state-update equations
    x1 = (2/R1)*VR1 - x1;
    x2 = (2/R2)*VR2 - x2;

    % create a vector 'y' for discrete output signals
    y(n,1) = Vout; 
end


% clipping stage

% Germanium Diode parameters
Is = 10e-6; % saturation current 
Vt = 0.026; % thermal voltage
eta = 2;% emission coefficient

% clipping-stage component values w/ DK substitution
Rb = 10000;
Ca = 1e-9;
Ra = Ts/(2*Ca);

% "OUTPUT" POT position from 0 to 0.99
outputpot = 0.9; 

Re = 10000 * (1 - log10(1 + 9 * (1 - outputpot)));
Rd = 10000 - Re;



Gg = (1/Rd) + (1/Re);


Vd = 0; % initial guess of Vd
Vout2 = Vd/(Rd*Gg); % initial Vout2 val.

TOL = 1e-10; % a very small value close enough to zero
xa=0; % initial state value for DK substitution circuit of a capacitor

% clipping-stage sample-by-sample processing
for n = 1:N

    Vin2 = y(n,1); % discrete input signal
    
    fVd = 2*Is*sinh(Vd/(eta*Vt))+(Vd/Rb)+(-Vin2/Rb)+(Vd/Ra)-xa + (Vd/Rd)+(-Vout2/Rd);
   
    count=0;

    % a nested loop to find the Vd value for each input 
    while((abs(fVd) > TOL) && (count<5))  


        der = ((2*Is/(eta*Vt)) * cosh(Vd/(eta*Vt))) + (1/Ra) + (1/Rb) + (1/Rd);

        Vd = Vd - fVd/der; % N-R update equation

        fVd = 2*Is*sinh(Vd/(eta*Vt))+(Vd/Rb)+(-Vin2/Rb)+(Vd/Ra)-xa + (Vd/Rd)+(-Vout2/Rd); %find fVd with updated Vd val.
        count = count+1;

    end

    Vout2 = Vd/(Rd*Gg); 
    xa = ((2/Ra) * Vd) - xa; % state-update equation 

    % create a vector 'z' for the final discrete output signals
    z(n,1) = Vout2; 

end


audiowrite('input.wav', input, Fs);
[input, Fs] = audioread('input.wav');

audiowrite('output.wav', z, Fs); 
[output,Fs] = audioread('output.wav');

% audioplayer for output
output=audioplayer(output,Fs);
input=audioplayer(input,Fs);


