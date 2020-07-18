function r = wabpC(Araw, Offset,Scale, Fs)
% r = wabp(Araw,Offset,Scale, Fs);
% Entrada: Araw (125Hz sampled) Forma de onda según el formato wfdb-MIT  
%        Offset, Scale
% Salida: Instantes de tiempo de la onda ABP de entrada
% Por defecto Offset = 1600; Scale=20; Fs=125; 
%
% Gnu Public License Applies
% 
% James Sun Feb 09 2005 with some changes from Gari Clifford
% based upon wabp.c by Wei Zong (www.physionet.org)

% Si la señal no tiene una frecuencia de 125HZ entonces hay que
% muestrearla a dicha frecuencia.

if nargin<4
Fs=125;
end

if nargin < 3
Scale = 20;
end

if nargin<2
Offset = 1600;
end

if Scale==0
Scale = 20;
Offset = 1600;
end

% Si la señal no tiene una frecuencia de 125HZ entonces hay que
% muestrearla a dicha frecuencia.
if Fs~=125
Q=round(Fs)
P=round(125);
Araw = resample(Araw, P, Q);
end


ArawReal = (Araw+Offset)/Scale;

% LPF 
A = filter([1 0 0 0 0 -2 0 0 0 0 1],[1 -2 1],Araw)/24+30;
A = (A+Offset)/Scale;

A = A(4:end);  

% Función Slope-sum 
x = zeros(size(A));

dyneg = [A' 0] - [0 A'];
dyneg(find(dyneg>0)) = 0;

dypos = [A' 0] - [0 A'];
dypos(find(dypos<0)) = 0;
h = ones([16 1]);
ssf = conv(h,dypos);
ssf = [0 0 ssf]'; %'
%plot(ssf);


% Regla de decisión
avg0 = sum(ssf(1:1000))/1000;    
Threshold0 = 3*avg0;            % Umbral de decisión inicial 


lockout = 0;    
timer = 0;
BeatTime = 0;
z=zeros(1,100000);
z1=z;
counter=0;

for t= 50:length(ssf)-17
    lockout = lockout -1;
    timer = timer + 1;      
    
    if (lockout < 1) & (ssf(t) > avg0+5)      
        timer = 0;
        maxSSF = max(ssf(t:t+16));  % Se encuentra el máximo local de la función SSF
        minSSF = min(ssf(t-16:t));  % Se encuentra el mínimo local de la función SSF
        if maxSSF > (minSSF + 10)
            onset = 0.01*maxSSF ;  % Se toma como instante inicial cuando la función SSF excede 0.01*maxSSF
            
            tt = t-16:t;
            dssf = ssf(tt) - ssf(tt-1);
            BeatTime = max(find(dssf < onset))+t-17;
            counter = counter+1;

            if isempty(BeatTime)
                counter = counter-1;
            else
            z(counter) = BeatTime;
        end

            Threshold0 = Threshold0 + 0.1*(maxSSF - Threshold0);  % Se ajusta el umbral dínamicamente
            avg0 = Threshold0 / 3;        % Se ajusta el parámetro avg
            
            lockout = 32;   
        end
    end
    
    if timer > 312  % Se baja el umbral dinámico si no existe detección del pulso
        Threshold0 = Threshold0 -1;
        avg0 = Threshold0 / 3;
    end
end

r = z(find(z))'; 
r = r-2;         