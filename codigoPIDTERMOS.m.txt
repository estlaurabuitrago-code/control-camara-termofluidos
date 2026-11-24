clc;
clear all;
close all;

%%%%%%%------------ INICIO PID CONTINUO- DISCRETO ------------%%%%%%%%%%%
%%\ Parámetros del sistema
numerator = [4.738e-10];
denominator = [1 0.001363 2.153e-06 2.535e-10];
num=numerator %[output:37d5afa0]
den=denominator %[output:3679fbcb]
TF = tf(num, den) %[output:71433da6]
nump= cell2mat(TF.num) %[output:02523a69]
denp= cell2mat(TF.den) %[output:2e2be209]
%%\ Simulación de la respuesta al escalón
t = 0:10:50000; % Tiempo de simulación
[y, t] = step(TF, t); % Respuesta al escalón

% Graficar la respuesta
figure; %[output:30cc8f46]
plot(t, y, 'LineWidth', 2); %[output:30cc8f46]
grid on; %[output:30cc8f46]
title('Respuesta del sistema a una entrada escalón'); %[output:30cc8f46]
xlabel('Tiempo (s)'); %[output:30cc8f46]
ylabel('Salida (y(t))'); %[output:30cc8f46]
%%\ 
%%\ 

%%\ Diseño del controlador PID continuo
% Parámetros deseados
zita = 0.9 % Coeficiente de amortiguamiento deseado %[output:58ace067]
settling_time = 30000 %[output:4ef2414a]
ts = 2400% Tiempo de establecimiento deseado %[output:71e936eb]

% Cálculo de frecuencia
%wn = 4 / (ts * (zita - sqrt(zita^2 - 1)));
wn = 4/(ts*zita) %[output:28f455f4]

% Polinomio deseado (2do orden)
%         4.738e-10                               kd2*s^3 + Kd*s^2 + Kp*s + Ki
%----------------------------------------------* ----------------------
%1*s^3 + 0.001363*s^2 + 2.153e-06*s + 2.535e-10              s^1

syms s;
PD = [1 2*zita*wn wn^2] %[output:371c3460]
PD = conv(PD,[1 10*zita*wn]) %[output:2d094491]
PD = vpa(conv(PD,[1 10*zita*wn]),4) %[output:9ecc9272]
PD = vpa(PD(1)*s^4 + PD(2)*s^3 + PD(3)*s^2 + PD(4)*s^1 + PD(5),4) %[output:4247710a]
% Función de transferencia de la planta
NP = 4.738e-10 %[output:57829534]
DP = 1*s^3 + 0.001363*s^2 + 2.153e-06*s + 2.535e-10 %[output:64134517]

% Controlador PI
syms Kp Ki Kd Kd2;

NC = Kd2*s^3 + Kd*s^2 + Kp*s + Ki %[output:7193e880]
DC = s^1 %[output:03907c0a]

% Ecuación característica del sistema en lazo cerrado
PC = vpa(collect(expand(DC*DP + NC*NP))) %[output:3e0e0278]

% Coeficientes para igualar
eq1 = coeffs(PC, s, 'All') %[output:5c8c94a2]
eq2 = coeffs(PD, s, 'All') %[output:1a5892f2]
eqs = eq1 == eq2 %[output:6844c8f0]

% Solución para las ganancias del controlador
solutions = solve(eqs, [Kp, Ki,Kd, Kd2]);
Kp = double(vpa(solutions.Kp)) %[output:23229880]
Ki = double(vpa(solutions.Ki)) %[output:247e7b44]
Kd = double(vpa(solutions.Kd)) %[output:9cb67a9c]
Kd2 = double(vpa(solutions.Kd2)) %[output:91718d43]
syms z s
%%\ Discretización del controlador
Ts = 0.05; % Tiempo de muestreo
TFcont = (Kd2*s^3 + Kd*s^2 + Kp*s + Ki) / s %[output:71758da0]
% num  = [Kd  Kp  Ki];
% den  = [1 0];
% 
% TF = tf(num, den)
% Gz=c2d(TF,0.006) % discretiza la planta por retenedor de orden cero
% numdisc=cell2mat(Gz.num);
% dendisc=cell2mat(Gz.den);
% Método de Euler hacia atrás
seuler = (z - 1) / (Ts * z);
TFdis = subs(TFcont, s, seuler) %[output:82160c90]
TFdis = simplify(TFdis) %[output:16821bef]

% Obtener coeficientes discretos
[num, den] = numden(TFdis);
numdis = sym2poly(num);
dendis = sym2poly(den);

% Normalizar coeficientes
numdis = sym2poly(num) %[output:8b4947da]
dendis = sym2poly(den) %[output:7c42d01b]

% Función de transferencia discreta
TFdiscrete = tf(numdis, dendis, Ts) %[output:580d9789]

numdis=double(vpa(numdis/dendis(1),4)) %[output:4e096393]
dendis=double(vpa(dendis/dendis(1),4)) %[output:6bc408e5]

syms z
%%%%%%%------------ FINAL PID CONTINUO- DISCRETO ------------%%%%%%%%%%%
fcnum=vpa(2.9821 -8.9447*z^-1 + 8.9431*z^-2 -2.9805*z^-3,4) %[output:9c7803d9]
fcden=vpa(1 -1*z^-1 ,4) %[output:66755bc7]
%ECUACION EN DIFERENCIAS
syms e ep epp eppp epppp eppppp epppppp u up upp uppp upppp uppppp upppppp ea Ua
NUMdiff= fcnum;
NUMdiff=subs(NUMdiff,z^(-6),epppppp);
NUMdiff=subs(NUMdiff,z^(-5),eppppp);
NUMdiff=subs(NUMdiff,z^(-4),epppp);
NUMdiff=subs(NUMdiff,z^(-3),eppp);
NUMdiff=subs(NUMdiff,z^(-2),epp);
NUMdiff=subs(NUMdiff,z^(-1),ep);
NUMdiff=vpa(subs(vpa(NUMdiff,4),vpa(2.9821,4),2.9821*ea),4) %[output:37e07119]

DENdiff=fcden %[output:9ee831c5]
DENdiff=subs(DENdiff,z^(-6),upppppp);
DENdiff=subs(DENdiff,z^(-5),uppppp);
DENdiff=subs(DENdiff,z^(-4),upppp);
DENdiff=subs(DENdiff,z^(-3),uppp);
DENdiff=subs(DENdiff,z^(-2),upp);
DENdiff=subs(DENdiff,z^(-1),up);
DENdiff=vpa(subs(vpa(DENdiff,4),vpa(1,4), Ua),4) %[output:727a23c4]


ec_en_diferencias=vpa(vpa(NUMdiff,4)==vpa(DENdiff,4),4) %[output:180b0379]
%ECUACION PARA MICROCONTROLADOR
Ua_solve=vpa(solve(ec_en_diferencias,Ua),5) %[output:6479a8d3]
u=Ua_solve;



%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:37d5afa0]
%   data: {"dataType":"textualVariable","outputData":{"name":"num","value":"4.7380e-10"}}
%---
%[output:3679fbcb]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"den","rows":1,"type":"double","value":[["1.0000","0.0014","0.0000","0.0000"]]}}
%---
%[output:71433da6]
%   data: {"dataType":"textualVariable","outputData":{"name":"TF","value":"                   4.738e-10\n  --------------------------------------------\n  s^3 + 0.001363 s^2 + 2.153e-06 s + 2.535e-10\n \nContinuous-time transfer function.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 78 117 109 101 114 97 116 111 114 58 32 123 91 48 32 48 32 48 32 52 46 55 51 56 48 101 45 49 48 93 125 10 32 32 32 32 32 68 101 110 111 109 105 110 97 116 111 114 58 32 123 91 49 32 48 46 48 48 49 52 32 50 46 49 53 51 48 101 45 48 54 32 50 46 53 51 53 48 101 45 49 48 93 125 10 32 32 32 32 32 32 32 32 86 97 114 105 97 98 108 101 58 32 39 115 39 10 32 32 32 32 32 32 32 32 32 73 79 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 10 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>"}}
%---
%[output:02523a69]
%   data: {"dataType":"matrix","outputData":{"columns":4,"exponent":"-9","name":"nump","rows":1,"type":"double","value":[["0","0","0","0.4738"]]}}
%---
%[output:2e2be209]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"denp","rows":1,"type":"double","value":[["1.0000","0.0014","0.0000","0.0000"]]}}
%---
%[output:30cc8f46]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAX8AAAEfCAYAAACpsMY9AAAQAElEQVR4AeydB6AVNdbHzwVpCiggRbDQBSsWXBQUXV1hde0N0BUrNlwbigWs2PuKZRWkSLOAIoqI+CkqimCBteBKU5AioCigoILvm1\/ey31z5936bpu598DLTSY5SU7+yfyTSTKZKiX6TxFQBBQBi8Cf1qF2oSNQRfSfIqAIKAIWgZB1+M0u8ZtCgddHyT\/wVagFUATSRCAQ0X3bKwUCvWhKKvlHQ0X9FIG0EdCRatoQagJZRUDJP6vwauLFi0AmR6rakeSlHRU47Er+eWlVfs1U9fInApnsSPxZQl9qVeCwK\/n7stWpUoqAIlDYCOT\/sSLA5J9\/8Aq7cWrpFIFsI5Dneziv2ef\/sSIW+We71jOQfm7BmzZtmjzxxBMZ0Du5JMirZcuWgiHv5GKVS3322Wey9957m\/j9+vUrD0jTlcl0KWNlypZmEdKOnlfOSFt7PyUQeQ\/nvD042W\/YsEF69epl7pOuXbvKqlWr\/ARQVnWpQsEBIKu5BDxxyLNPnz4BL4V\/1Ke90e7uuece\/yiVgiYOZ6QgHV1UO5ByXILVHgqn5qrMnDlTMOVVoS43AoxMJ0yY4PbKjTsTDJMFTffcc0+ZO3euLFq0SO67775K5TBq1Kiib3M+rd7I+szRVbDaQ+HUnJn2WbBgQY6qWbNJGoHCGWAkXWQVVAQUgdwhEJX87WMY883WMPURTS3m6ayMtaPNnRG\/ZYuWZh6aJw0e+608c9PMJbvTN\/LOnLc3jDk50icuaaCrjYcbP8KsIR0b7raj631oxJwfcd3TPUxTkC5xbVo8GeDnNuiHnlYmke1NozT92OxPeDr5RdMH\/MHanS5udHPLu+XAxx3GNXHchvqgXpDD5hocucaAL\/LufHDjZw16kS\/y1lgZG+bNm3BksW062GCHv9uQNukQbg3X+Lvl4rm9+ZMOfvHi2DDkkPfmSRuiLREGbuBHHMqAH4byua\/x4xo5tyEuaRDuNtFk3fHcbrc+No1o8b3lsdc2DjqTrtUpVnsAfzAhHmlgcGNsvm4Z\/DHEwZ88vMadBnKff\/65VyR8bfUjTbexeYcF4zgMZod2NWsKNo1Y8d26WVnqDD2iZUE6Vs7aFlu3PH423Nq0K0P+rVu3Dsui7NFHH13hsZypD68iKEvFeWlq6dKlcsopp7iItEzCeWJav369WWChA7CZ4nfccccJSlq\/VO3K6C1SpldZZkuXLvHoXRZgrUhxswAMedlgay91yv+3v\/1NYjVAK4dNBXrTANN77rmX4ArGYu4OSCU\/dzzrBnfwpx6sn7XRbdq0afYyqk3jpG3QRrwC1DNpIOMN8+JPeDQ80Av90BMZt7Fh3rzJ8x\/\/+Idgu+XvcdYZyAM\/qpM0SZt08LOGa\/wJt37RbMoVq+zoRBgy0eJmwo\/yUSZ3Wly79Y51bxAHWdoU7niGtkybpq255eLFtxiCgzvOVVddldS94Y5DGhj86tSpI507dzZ8QR2RD\/7WcI2\/GwPqgLqwaSCLHH60UQlzAa1CDHfBgx840+LIuk28MrvlwpgtWRpOnXBv\/Gi6IYdBN+oYGa4xuNGbdLh2G2Rt+8YfN3643YZ6rGKBtAEkSMBOO+1kOgDmdp988kkTjCLMz3FBg5o9ezZOObBTJ\/nvf\/9r5oGtLGmQlhEQh\/VLHeaXPCdOnBghT8Btt91mQMedqiEv8kxF706dDoyrN3Patjzoc03\/a4zOF154odFz7NixeMuJJ55o\/MHqmmuuMX40rBkzZhh3xE9p2zJeNI7HH3\/cuGNhYgLLfmjMtvGiF\/lRJ5SZ\/O68806hYZSJJ21NmTLFyLp1sOkSULFeIutz8eLF8sUXXyAaFQvSot3Url1bxowZIxYjcdqFLccRRxxhCMHigQzlo111ctqXOP8q6uF4On9uvYnneJm\/L7\/80uRFOrQ35Ah45513DE6\/bNggTz\/9NF5CHuSFLDoZT+fHYuM4o\/5RNgyBxCM+6ZAefuACPrgza0obEmWibOTrLrtbb9oM94ZbFp1pN+jEfcz9jDuaoU3RtmhjlIvyufMjfdpmtLjkQV7Ic58gQzrcGxXbg4jFkPaArNvYMNacWrRokVLd0f7Qg\/TcZbA6idMWxfwrbduUCczqOh2NxZf4lAexVDH7rIwfbR2RvsWM9kE7IV30ASuMlSVf9Cccgxs\/3Mggi47UL37cQ3ALdRqPo6rcf\/\/9wiIekYjwxhtv4JSePXtKo0aNjJuKQCku7I3z\/fffy08\/\/YSXufEpABfIogwG8sTPa9x5Im\/TBmwallc+0XWu9HbrATbTp083pO8up\/spyi0fdpe2LXPJDcCNwMVFF10Urgc3JoRhuAHdREUDxh89Bg4ciDOiHoxHEj+ku3z5ciOJLujEBena8mFzjX+pKSWeUrcIbYG4XLtvCjpJ2gEGN+HxzIgRI4R0uMFsm4AgrrvuOqFhx2ofjEhtG2ZEiCz5gNEZZ5yBUxo3bizbbbedcS9btkx+\/fVXIW06I\/TD5hoBZG0aXMcz1BXxMbiRJZ2mTZvizKIpbUjudgNmYEem1Cl1ixvs0Q\/StDhtvfXW0qxZcjpCNBjSOueccwxuuMEWjHG7OxuuraFt2rbTu3dvU4+EpbrOSD4Y4mLAmDqjXNhc4x+t7sAB3iKceqU9WXnI02JGuDWxMWtmReLa4IVBKBFmqdw\/bkIHD+qAPKhX7gP8uIe5BvfpcTiqim2wJOBWgpG0nR\/CpqcSEbE3Dj3v7rvvTjRzw\/KYhRyGRw0TEOUH8Kkgd1D37t3Dl6k2CiJW1LtFeI4tU3qTTzxDpTCPFu0RK1Y8W1YwgbTcct5OBLICe2RoVHvttVe4jDZPiBMskEnWcBMccsghYXF3vfNoyY0TDozhQJeddt7JhELQNEDaAcaObkxgnB\/ygbAQ8aZB26JshC1YWHFzghcr5DAQMOXDXUp28W9cdKDM7vyIm4qhvJTbtrtU4qYmW9oBRy17aVDM5Lg\/qbOZMz+MKeMOsO0UP9oa5cOUpjETb6HuwM9clP3Qrr33ellQypa7LqNFJu9Ydee+d+AsuMumkUy7QLYcs9Ly4hcP5lQwA0fbAXnbPu2JvKzh\/raDbi8mDELdHaGNY20vR5k5fxuYjE3GKMBN9dBDD4lV2h3XEohXcbeM1x0PSK9s4uvSUZFbLht6Uz5uAgyER8W588yH293oks2fUQ6jRq\/8TGe+k4aZaE6YEcbgRwaHR3XudCALFtZ4OnP7J3THaBAL5lck\/4RpxRFAL\/SjDikrZY4jHjUIfIiPobxRhRJ5xihv7GgV23hY1hNkiRH9MNyfYdkMOZYtW26epjKUXFLJZKLuYmWUDGYemGMlVeofpX4ZzNExmftncIbvn9JczboIdY7xclQVeoMyuQiLxyEeqbzG\/eiI0vaxIhp5xJqjjcio7CIlIMviRLNypTcNr5+zcGV1oPxgxbyk9cuWTSUyPUZ+XgORVyZfRg2kFU1\/RrF0dPHS5TGTtoFe6OeWZdTOnDE3lNs\/lpsBxcwPZ5opNXRyG\/SMFS9Vf9p+3759zZMrcdEb\/c38ad06eCU0jAjBxwqC38JFi8zah\/VLyk71BohCJrHyufnmm836HeGMxikf5aS8+KViKN8ip3zuOsE9ffrb4WniVNJDNtWiEydu3Tnz9MikYzKO2eJFFdrzdGdKBg5Fz0zeP6SHgaNYXMeN8XJUFfcogEc0GgeCqY4guSlpBO5GZUfbpGcNRMCclL3Gds8Xeh9jkedJAzkMbtLFbU3W9Y5yo1GGdevXGxW4ISi\/uUjhx5aVMpKeO6oXf\/fjqR0xuOUz5WYakHrE0JHadL36WH+vzRMhj57EZxQNkSOTSGfi8RiLLPVLPePOpqGt2qc1yore6GHyjFLnxt\/1Q2dm55IpJ+UFP8gsieiulEqdtAN3uXGDRWmo55dMPF7RLiFJ1mEIg+xpZxAN18ka206RT7YdIJusqQxWcevOk7H73vniiy\/Frk8ixsibtonbmnxiRvujHca6f+A6u3blnWpjIMLUM\/pTFuqaNoU7GkdVYdRiR3XMhTEnhrDbnx7EPhrziEs4cXiUwFg\/\/N1goiTK4u82djUaP9K2i8w0Tgz+7gZnOwduNkaQpQUqbzJZ1zvBjWZvCEDnaQf9kzHM89vOllV54hMPbMEftzU0Cjs3D2HZThtMesU5m6QcJZtSpE2eNBjqEZtrK2HLxbW7Prh2Gxod8TG4bZibvJo1aybchDYsmm3Xfqhf6pmyIUf7Im3aIO0FP7dJUD1u0ZhuW1byJG90iCkcJQCSprwEUX8vpvBWuBvb6G2dVNM3kBz3Jym5d4xwHc\/s6awv0bkh4753KSf1gnHXO3K5NJF1d0f4Sc7q4L531q9fJ2wssGHcR9xP9tprVxYzphCTxQzswBCD2+pAe6JdcW3vH54U2IyDH4MN6hE39y0cQlk6\/aWTme7B3xqLEXKWo8ycP7tIaPSAxEq4JSTmL1HILoDhz4o9CTLC4TECN0SFHAbyRgH83Sv9XFuz3hkxkyby2FwT5l4VdxOjTR9AKTCyEt6aJWb3gVvv8\/v0MYuhNu109KbzIr44\/2go6EwFufWz\/u6yO+JiAccdzTACY7cGYWBGfNIHd\/y8hpV9ZPCPwOSD0kUoGgWNg3BrEhEj8tQT8m4d0IM88KeeqW\/ckab0yq2XxYL4bvypH9oXMdxkR1mRhUjIg7yQoZ6pb8ImTBiPl7CbAczMheunxOVOxUn69galrOSFH3nbdLyjK+uPTXlsh0wbprykQZkIx+DPTYw7lnG3pVh6xIqbyJ\/67dixoxFz1y\/1ZDydHwgmno6NGzUS20Yoj7ecYGjrzUku5b9Y7SFeQtQT+SITiVn5Ira77txt1MpTV7hJw23CmDkNK4xZi5aSCmakkSxmbt3IA70w4Lze6azgH\/f9E00eXkBXynHiSScK95K7Xdl03XKG\/Gnstgfh5uJxASESsqbTgZ0Ef8KtH1MdPE7Ya2ujLPOKKGD9rG3D3I3F+rnlyQedCLNx0enNN98UW+nWHxt59EPGTXhc4084chirt1O3XIYNeXn1Jh5bU8NCjgNSx9+rH3qBpS0bj9v0tE6UmH\/M0XsxZAoC440E2fBI6A2rU7eOoDdpeeMkcw3u6G137LjjkBd4uf287lh6IQcmPNmBF9cY8iNd3NaAKW7y8uIhTkePH2GSwX\/coM8\/\/7ygo03WtgGrH\/uv3dMEVs7aYG5lrR\/1b\/DcqXQHlB3N2\/BIu8Rs8fW2JdpteVv3ttTIFBJdgRs6ueXQmTZDedc7gzHuEXe4102ductkw0nXPXdt\/VOxSRt93HFse3D7ud2p1p1to+hr06Hs3E9gbf2k7HUsg5lDovgb9B1SQceUMXMGZu72RXro4MbM6kb6hLvNTjvtLN77J568+z7hnvO2K3ShHtGhCnNLGBqxzdQmjr81Y0aPMSNsK2NtKs7KWJuFPzK2MtFswE0kTxqkZeVGjx4ty1csl+22DmvEvQAAEABJREFU3U6qVq0qsz6cJaeeeqq8P+N9KSkpMfpRmVYem2vK49UBvRd7Fq7IizyjyZKWNeiODLLEsf62QgnHz14jG8+gC\/LWUBcYe024O747DBl0QBe3DNf4E44+7rBobm6m6W+XvrdAHGvIyy0fL11kbTxrx8LAK8u1zYfy2vjWxs+GY3Ntw9xxY+lHG6AtEMetkym3s\/CGPwbMSIM03dfkGctYWeQx4O1Ol+tYccXp2MT5R57kTXwMujKdia6LFi02L8dRBkdU3PmBA34Yd57Et\/KEocMiV3snDXeeXCMXz7jTt2mRrjcOfoRTHvKw4bjxIwwZ649N\/vhbw3U8eeJ49SFt4hCXdOw1staQL2EYwiF+sOIarBs1amxFzcGF+FueIF3SJx7+XIeFHYfpJBzb\/deocSMhXeStQQe3jHWTnpWxNnEbOU9eVsZtR5N3twdk3fqSpk0PHczIH6E0TU6iT3ltipx3znlS4vTCgwYNkkF33G7yPffccyT+6MqI6Y8ioAgoAhlHwJK+Q0sZTzubCQaG\/Hk0HTp0iOy9z94yetQoOe2006SHY0aNGi3td2svQ4YMrbDQk03gNG1FQBFQBEAgaKSPzhifk7\/tU0VWrFgpq9eskS6du0S8TFS\/fj3Zf\/+OsuaH1bJx40bKpCYTCJRDn4nUNI1iQEDLGCgEckb+zDEx58R8GfNQyaFU3qe2bdvGzJ3xUo477sZNG2X+11\/LNltvI1tttZU7SN3pIFAOfTqpaFy\/I6CdvN9rKGv6JSR\/v7eNt\/\/vLfngw5ly5JFHSv369bMGlCasCBQkAtrJF2S1JlOohOTv57bx3nvvyY033SQHdDxAzj333GTKWyrj9x6tVMsA\/6rqioAikAsE0qGyhOSfiwKkmgfbOl+dPFkuuuhCadWqlbAPnz270dLhZYkKplVL8xJYBf+W6q+YaBvQNpCdNhDvfZFo3JWMXzqD88CR\/5YtW+TJp56Uy\/91mXTteqg89dRT0rBhw7g4sdagZlGFg6UWLYzi59oLninMqJxMpaXpVK7OtA4qh1sm29vhhx9ONfjGRCH\/dB4kKpQrox5\/\/PGH3H333XL\/fffL2eecLQ8++GDEzp+MZlYMiaUzbCgGfLSMikABIxCF\/P3JCEz1PD1smAwfPlyu6tdPOOuiWrVqBVw1fimafwcDfkFI9VAEgohAFPL3ZzHmz58vTw8dKttuu60s+PprueGGG+T6668PG54I1q1b50\/lA61V+oMBzqgJNAQFoHxe66AA8PNfEdIflPmU\/CsW7KPZH8nq1avkhx\/WyPgJE2TcuLHGjB07TjCTX3tNNm3c5L86Uo2yhkDFVpK1rAo2YcUwqFWb\/qDMp+RfsWC9Tu\/lLFguLjMs3pS6Fy9eJJjpb78tHKIU1KpUvVNHoGIrST2NYo+hGBZvC\/Ap+eeyQnTskzu0NSdFQBHwCwJK\/mVH6vqlQlQPRUARUARygUAgyV\/H6rloGpqHIlBJBPQGrSRwuY2WL\/IPl7Iy7STo85SVKXMYMHUoAgkQyHv7CvoNmgDfQgnOO\/kXYzspxjIXyg2Tv3IkT+navvJXS0HKOe\/kHySw\/KFr8iTgD31Vi8wgkD6lL\/lxk3jNjAU\/iduMnb1SvObu17+RaOaSsfMkljn20U8lkekw6APpkILJDI7lqRT7naTkX94WknPlvcWkTwLJFVSl8oGAm5wtKbvJ2JKwm3S9JGsJtf6Vb4k1+z3yTVSiPeaxT8Vt3Ola992vL3bIv6Kxes1Y+FOFDuM9xy+RcZc1GTflymSdFN6dlBo5Kfmn2poKr8WkioDKJ4GAJTMvgUPellQtaUNqlqRxW2NJ2cpjWyK2xIvtJVmbdxJqZkSE\/DKSkJPIzvVrOr+x\/zKZl0hqZBlbK7+EpEZOSv5+qbei1iM4hYd8MJAuxpK5l8hjETjkTTzMkrWbBOImveAgkF1NwSJaB4DfyxfvI9iZ0yA1ssxcvv5IScnfH\/WgWvgEAciH0Trk7CZ2yJzROTaGUTjGkjkkjiF+skVJRTbZNOPJQZxe06XVduI2PTs2kWimf7cW4jWP9mwvscwkh6gTmTkDDpRoBpKP5o\/u8cqnYakhoOSfGl4qXQAIQLqW4CFwDCN3S+5Mt+DnJnbi5KLoEJw1lpTdZGwJ2E26XpK1xPnjA4eJNR9f2jw60V6yj7zsMu503e7+3Zo75B9p3Hp53Z1bbyeJjC1nMnYusC+2PPxG\/sWGv5Y3iwhA2BhG8ZB5NIInDMMUTCZVsYTmJfD+zgjakqolbTdZWze2JWUrj92\/jITdZOslWZt3JsujaRUeAkr+hVenRVkiSB4DkVuiZ3oGwzX+TMvEAoe4scLc\/pZYIXUIuL+HzCFtRtvYGC+B93fIm3gYS9o2TXc+6lYEso2Akn+2Edb0s4IAZA2hQ+yM6I997FOzlZFr\/CH6yuzlsEQMOWMYbdsRupfUCfOSOfGzUuBiSjQwZa1MC\/NP4ZT8\/VMXqkkcBLxk7x3RE+6NHmsvBwSNgdwxkDgE7yZ3\/DCE2xG6N\/2cXQebY3IGU+4zitXCcq9JZXJU8q8Mahon6whA5izKsuOGkb2X7K0CkLh1R7MJh8AxkDkkz3QMhmsMYRB8tPi+8PMVx2hP5Is2kQEllPwzAKImkRkEIHymbJi6geyPeexTGTt7hdkLHysH4tgwS\/QQOiRvR\/JcY3xP8rYgvrZ91RP5GqncK5dax+wz8k9N+dyDqzlmGgHIG8L3ju5tPoRbt9uG6DEQOsTuJnr8fD2SdxdE3YpAxhBIrWP2GfmnpnzGMNOEcooAhA7h2xE+Ngu0iZSwZM+onmkbDMQP2SeKq+GKgCIQiUDa5J+jsXqk1noVSAQgfebwmdKB8OkA4hXEkj0E7x7Z66g+Hmo+CFNS8EElJFYhbfLXsXpikDMtEaR7C8KfNG+DOd4X0uet2Xh4QPj9u7UQO7qH+HVkHw8xH4YpKfiwUiqqlDb5V0zSfz5BIstk0AvCvQXp21H+zdPWxF209RJ+\/27NzdEAyWChMkWEgBY1owgUBfkHgSwzWqt5TAzSZ0on0SgfwmdEb0f4Svh5rDTNuigRKAryT7ZmC+0JIdlyZ0IO0rc7dmLN5VvCZyrHLtbq\/H0m0Nc0FIHUEcgq+QeNTPUJIfUG5Cb9WDt2mtbdSvo78\/gc1QvxM+JPPSd3DHUrAopAughklfyVTNOtHv\/Gh\/Tt9E4s0mekD+lP6r2jQ\/7NM\/whDv9iEwzNgjY0CwaqQdIyq+QfJCBU1+QQgPTtQm686R1G+Ezt9HcWb5NLWaVyi4AOzXKLt\/9yCzr5+w\/RAtXIkj6nZ8barslI3y7g6tROgTYELVbBIJAX8tcHzmC1H4jfkj5ur\/Zu0tcFXC86XGuLBwU1\/kIgL+SvD5z+agSxtIHo7RQPbq8cpG+nd5T0vei4rwuvxUd2Z5FX7pLnxK2ZVAqBvJB\/pTTVSDlFgOOU7WjfmzGk379bC2FOX6d3vOgUx3VkdxZ5lQgB7SoSIZSb8ECRvzaa7DcKRviM9o957FPB7c0R0mfLZn9dyPVCo9dJIpBaV5FkoiqWMgKBIn9tNCnXb0oRIPt4o30WcyF9Rv6xE85nFx1bKw1RBBSBSAQCRf6RqutVJhFgtM+RDHQA7nQh+v5lUzzJzetrF+3GT92KgF8RUPL3a83kSC\/IHuKPtn0T4n+0R3vzglaO1Cm+bPRBqfjq3CclLlTy9wm8\/lYD4o81zZPaaN\/f5fS1dvqg5LPqKZ7eWMnfZ00vV+rwdm68aZ7+uqCbq6rQfHyFQPH0xnkl\/z\/++EOuv\/566devX1LV\/8svv8iZZ54pLVu2dJkW0sK5TjaNpDIqcCGmeTiXx1tMpnl0J48XFb0OLAJ+VdwnDxd5I\/+NGzfKPffcI88++1zSVbTBIf9vv10ie++1l\/To0aPM9JSejrtjx45Jp1OsgkzzQPzR5vf7ly3q0gEUKz5abkUgJwj45OGiSk4K68qkpKRE5s2bJ7169ZKhQ4dKScmfrtD4zpUrVshPa3+UPhdcIHfccUeEOe200+JHLvJQiJ+pHi\/xQ\/b9HeLvr9M8Rd5CtPjFhkDOyX\/16tVy4YUXyoIFC+Tqq6+WA\/7yl6QxX758uUioiuy4446i\/5JHAOLvO3aeFC\/xJ4+VSioCxYJA5sg\/yXmsUCgkHTp0kPHjx8s\/\/\/lPqeJcJwv2woULpX79ejJ79mw5sFMnadWqtRx11FEyY8b7zhNEkgokm1mByEH87Ojxnrm\/c\/2awvy+Hs9QIBWtxVAEUkQgc+QfSi7nhg0bysMPPyxt27ZNLkKZ1ObNm+Wrr76Sb7\/9Vl6aOFEuu+IKufHGgcIi8HnnniNTpkwpk1TLImCJH9v6YVvix+ZajSKgCBQfApkj\/yxj9+uvG2Xjxl+lS5cuMnrUKOnhzPH37t1bxk8YL+3atzcdClNKCdQommAInxE\/trvQED4HsmG7\/dWtCCgCxYWAD8g\/uUeGunXrOAvET8vIkSOlTp064VravsH20r1bd1myZIl88803YX+3w24NXbx4sRSDmTF3vhz179kVDmY7pn1tefH0JjnH4Lvvvst5nsVQz6mUUesgP\/e+5R5sNyf5wZ0U+bMfn1H1nDlz5NVXXxVsrrds2ZKBMiQ\/V09+6OLNdOdddpZNmzbJzz\/\/7A0y14sWLRJMixYtpNBN1W13kEsmrZHl6zabstsf5vZHnN8xL+Vngb7Qcfd7+bQOXPd+DnkA3rHG3ot+sWOS\/8aNG+Wll16SE044QXbbbTf5y1\/+IieeeKJceumlxua6Xbt2csopp8i0adMkGiknU8jkxv0in332mey3377y0EMPV0j2i8+\/kO3q1ZOmTZtWCCsmD6Z42NWD7S53z45NhI+uuP3UrQgoAsWNQBVv8SH9ESNGyMEHHyxXXnml2ZK5V9lLVddee6088sgjpgM49dRTnUXbXQ0p9+nTRw488EAhHvG9aca9TpL9mzdvbnb3TJ36uixbtiycJI+zr015Tfbacy8zqg0H+NGR\/ENOytpD+BC\/d1ePEn\/KUGoERSDjCGTx1q+0rhHkP3fuXPnHP\/4h9913rxx55JHhKR62Zd5xxx0CyR999NFyxRVXyF133eWEv2LIn\/BDDz1U7rzzThN\/5syZSStUUhKd\/XmaaNmypTzxxBMmLeb5zzvvPFm6dKmceWZvGTNmjAwfMVxOdRZ+1\/28Ti6\/\/DKpVauWkfXtT\/SiZkRdXuDyEn+XVtvpiD8uuhqoCOQGgSze+pUuQJj8V61aJdddd50cfvjhZt88ZN++fXupWrVq3MSrVasm++yzj9Nh3CezZs2SY489Vm699da4cSICQ8n3id27d5ehTz8tNWpUl4EDB8rdTge0m6Pjc88\/Z94diEi3iC6iHdnAbp6XL9mniFDQoioCikAqCITJv54zZz5u3Dhz0FrdunWTT8PF3cS77LLLhHSSSaB27doyZvQY03F45Y844h3X0XsAABAASURBVAizSHvhhReGg0KhkHQ+6CCZPHmy8MLXvHlfyZAhQ\/w\/3RMuQeYdjPijvbnLds7M56YpKgLBQ8BFUcFTPosah8mfETzknXJeUZ5nKpVOyhlnNEIgE1vy4ybxns5pRvwXxx7x640QyKpWpdNAIApFpZFa4UQNk3+0Iv32228yffp0ubpfPzMdxPQO8\/DYrAnccMMN8u6771Z6p0+0PNUvOQQgfl7icktD\/Hx5C9vt73brjeBGQ92KQPEiUE7+riHhj2vXmhMz999\/fzn77LNl\/IQJwgslwMR+YY5a4GC2sWPHSu\/evWW\/\/faTu+++W9Y68ZBRk10EIH529mC7c+rfrYUk951ddyx1KwKKgEGgyH7Kyd8ZErJN88knn5Suhxxi3qRlB8+wYcOEl7p4UeHTTz+Vd955x+zwYc79ww8\/lMGDB5sFX+QOceIRf9PGjUUGY26Lyzy\/d2cPxM+2ztxqorkpAopAUBEIk\/+6devk\/PPPF8h7wIABAtGzp79r164SbQ4\/FAoJh7Rxqib7+z\/66CMh3tNPPy1nOU8LQQXE73pD\/N4FXki\/f7fmfldd9VMEFIEcIOCaxImbW5j8OTrh5JNPNnP4fBgl1T3z7Nwh3ltvvWW+sBU3Vw2sFAJLYizwBuLt3WRbZKWQyXQkTU8RCC4CziROUsqHyZ+tnscff3zaL0rRaZBOUrmrUNIIQPzRFng5kz\/pRPIpmGyLzKeOmrciUIgIxBh4hcnfW2amgVjE\/fzzz71BEdeE9+3bVzjoLSKggC9iYJnVEjPVQwfgzqS\/s8Abb2ePW1bd2UUgH20iuyXS1AsGgRgDr5jkzymZvEy1cuXKuBh88snHZgE4lW\/xxk3QX4FRtYmBZVTZTHgyz49xp8U8P8btp+78IZDrNpG\/kmrOhYJAmPzZ6XPN1VcLO3YwnOC5fPly88Yv14cc0rUszNqHmMPf7rjjTtl+++2d6aKtCwUTX5WD0X60F7kCMc\/vKyRVGUUggAhk8ZEyTP7M1Z\/Zu7fQCXBSJsTPIvCaNWuE6+++W+qxvxNkmjVrJldddWXEB1YCCLFvVWY\/v1s5pnl4kQu\/LLYLklejCKSAQAG3xhRQyLhoFh8pw+SP0nvssYfMnj3bnKnDyZw77bST2frJHv9ohr3+b775phx0UGeiq8kwAkz1ePfz9+y4Q\/hFriy2iwyXRJPzNQIZ4W1tjb6u4yjKRZC\/O7xBgwbywgsvmKkdt7+6c4MA0z0s8rpz44jm\/rqf3w2JujOBgPJ2JlDMfRppdtph8md3D3v0mfahFBzlzEtcNWrU4DJpwxe93nvvvaTlVTA6An3Hzov4Bi\/TPYN7to8urL45QkCzUQR8hECanXaY\/Jnf\/\/e\/\/y2cmT9p0qSUD2uD9InHUcz333+\/jxAKnirRpnv667bO4FWkaqwI+BiBMPnzkhdfxzr99NPNDh8Oa+PUzo8\/\/lg2bNgQtQj4M8q\/+uqrZd999zWffezWrZv5ylbUCOqZEAGme7y7e5ju6dmxScK4lRdI8\/mx8hlrzAwjoDWZYUALOLkw+VNGdvyc36ePTJkyRTizZ8L48eYD7XzDt327dmb+\/5BD2OJ5iPDxdvzP\/OeZwoifzzu+\/fbbpuMgHdIrYJO1onnn+XMz3ZPm82PW0NCEU0VAazISMe0MI\/FwX0WQPwE0HrZv3nXXXTL7o4\/MqZ10BA2231544YttnytXrjB7+\/Ef\/Ohg4VA35IlHGmoqh8CMBT8JUz7u2D077iB0AG4\/dSsCikByCMBnyUkWn1QF8ndDwGFthuAHDxamd+bPn2+2gWJzzXHOhCPnjqfu1BFguifaqL+\/7u5JHUyNoQjkAoGA5xGT\/Dmr57777jNHO7MYnO1yFvvj2YyFP4l7Tz+jffsyV7axz3T6OavLnGWUaYTykZ6OgfOBej7zTHR7xCT\/kpISM5d\/0kknSadOncyXvZYtW5a1shRz02TU713k7dxqu\/DLXFkDPUsJ56wuc5ZRloDKabKJqCCnymhmOUAg0e0Rk\/wbNWokL730ktx2221Sv149GTr0abPge\/jhh8vo0aMTfLJRG1oqdRttukfP7kkFwXzLav6KQPAQiEn+FIXtn2z9fH3qVJk58wO57rrr8JaBAwfKAQccIGeeeaZMmzYtyjsBifoc0X9lCDDq9y7y9u\/WoixULUXATwjooM5PtZGuLnHJ3504b\/vymUfO8uHbvbcNGiQ\/\/\/yzXHDBBbLnnntKv3795Ov5XwvTRe546o6PAG\/yuiWY6++Z1T397tzUXcwIpE7lOqgrpPaSNPnbQvNiF4e\/vfbaZPnf\/0rJnn39dAp\/7\/53Oeuss4rqwy4WlzI7JYutne5FXiIHdZEX3dUECwGl8mDVV6a1TYr8IXw+7MI0zz4d9hG+3PXJx5\/IMcf8Q8aPH2\/2+bPXf8iQIfLll1+aaaFMK1po6THd453r503ezq23829RUx8qZrwsPlAh42XSBBWBfCAQk\/w5q+e1114z8\/r77LOPXNK3r3zwwQfS5eAuMmLECPnkk0\/k3nvvFcI4BA5z2GGHSdeuXWXevHn5KEug8oy2tbO\/3+f6fTBU9IEKgWpnqmweEAhIljHJf+3atcJbuzNmzJDddttN7i\/b8z9s2DCz66datWpRi1izZk3Zddddo4apZykC0Ub9Qd7aWVoq\/VUEChGBwn3WjEn+kPjll18ufNRl4sSJwmcdk3mTd5CzEPzkk08WYivIWJkY9dMB2ARZ5E131F+4TdSi5Cdb0fZTbWRXl8J91gyTP2\/xMtVjgaxbt66ccMIJwi4f65es7U4n2TjFIgfpe+f6GfXTAaSDQeE20XRQyVbcZNHOVv6ariKQPgJh8v\/hhx\/kuOOOk8cffzzmEc6JsmNhmFE\/nUYi2byH52nwlo1Rf96xVAUUAUUgcAiEyZ8R\/oABA8ynGw866CBznAMneCbat0\/44sWL5dZbbxXijRo1yhzr7Hsk8jB4y9ao3\/dYq4KKgCLgOwTC5B8KhQx5v\/LKq2Yr59ixY4Wz+\/lIC1s8Wfzl3P5XX33VnPnDNf6c6X\/4Xw+X559\/Xs4991x55ZVXTDq+K2l+FTK5B2nUn6cHI4OT\/igCikD2EQiTv82qVq2a0qdPH7Ot85577jFz\/uz4YTrnsssuk0svvVSwucZ\/hx12kHvuvcfI489agU1L7XIEGPWPnbWi3MNxZWKu30kmK395eDDKSjkCnaj2wIGuPr8rX4H8rcLs7Dn55JNl6tSp5sWtd999VyD8Rx55xNhc80IX4cghb+OqXRGBpT9ukiVrN4UDWOBNd4dPODF1FCYC2gMXVr36rDQxyd+tZ40aNYSvdB1xxBFy9NFHCzbX+Lvl1B0bAXb4MPq3En4e9Vsds2LraDYrsGqiikCqCCRF\/qkmqvKRCHjP8GHU36vjDpFCxXKlo9liqemycmpvXwaE7ywl\/xxUCaN+dzY716sZ2A+1uMuh7soiUEzxtLf3a20Hlvx5kez66683R0n7FVz0YqrHe3KnzvWDjBpFQBHIJwKBJP+NGzcKO5Gee+65fGKXVN7eUX\/OTu7Up+2k6keFFIFiRSBQ5M8LZZwY2qvX6TJ06FD5c8ufvq43Rv0JvtKVPf31aTt72KadsvbMaUOoCaSNQKDIf\/Wq1XLhRRfJggXz5eqrr5YDOv0lbQCymQAvdbnTZ6HX1+f1u5VVdxYR0J45i+Bq0kkikJD8Odp5ypQp5nhn5ti9hs839ujRQ1atWpVklpUXC1UJSYe9O5gPyPzzn\/+UqqGE6lc+szRjMur3TvnoXH+aoGr0DCKgTx8ZBDO9pPIUOy57zpkzR4488ki5+OKLzYtd48aNE6+ZMGGCcChcKJT90QznDz388EPStm1bA1f05hvd10TI4Y95qevHyJe6enZskl0N\/FH07JZRU88QAtm\/XzOkqCaTJQRikv\/mzZuFD7f8\/vvv5pA3zvRp3bqNXHDBBfLmm2\/KnXfeKS1atDAvfE16+WVzDESWdIyTbDS2q9ioW7ZsKRgOoMuVuXli5NfM9m68lWQ9728WZz+PxannwQGBWS97JfQqJp20DlJvt5loH\/CONXGILC9BMcmf6Z7PP\/\/cnOl\/2mmnSbt27WSPPXaT+fPnS9OmTQU\/Pue4aNEiGe+M\/vOifUWej6oGOmLorHJhqm67g3y8rHzUj1J9DmtjOstc5O+3PHbcccfkyt68RXJyzqDDb2XMjD7ZK3\/SdaDYZrQNwjvWwAN+MjHJn501jP47dOggoVDImLZtdzUjy3Xr1pky0KC6d+8ur7\/+uvzyyy\/GT39EvHP9OdveGXTwQ0EvgOqvCAQHgZjkHwqFZKuttjLz+bY4rVu3lp9++klWr15tvaR9+\/by7bff5of8S\/zHFiz05m17Z7hW1KEI5AiBaDOvOcpas0kPgZjkz9HMPMrOcRZ9eZuWbFhw\/fPPP2XhwoVcCk8HTA2Zizz8hMR\/LW\/Gwp\/KkXDUq+T2zvI01JU2Ak41pJ2GJhADAf+Nv2Iomo53YbagmOTPiZ2HHXaYTJ48Wc4++2xhqscuXNxyyy3mJav7779f+HJXmzZtZNttt00H3UDETaYJRJzZ79wYPYv1ADcf1ahTDT7SRlUJHgKF2YJikj8VdOqpp5oPu1SpUkUwderUkfPOO09+++03uf322+Wxxx6TmjVrmC9\/0VkQJ5cmGTLOpD6JmkC00zv7d2ueSRU0LUVAEfAgkGse8GSfvcsspxyX\/KtVqyb9+\/eXIUOGyDbbbGNU+fvf\/262ej7xxBNm7\/\/UqW8Ii8ImMIc\/fDxmzJgxct999+Uw1\/hZeRd6ObM\/fgwNVQQUgUoh4GL8RIOySqVfBJHikr8tf\/Xq1c1uH3vdqFEj8\/IXH3WpV6+e9S5qm4Ve7+mdRXtmf1G3BC18ThBQxk8b5jD5c1LmW2+9JbzMFdVMfjUibPKr5dfEI37a2gQ4Ae8On4oLva6hSoDLqarnAwHNUxHIPAJh8l+\/fr3cfPPNcumll8qlfR2D7TYev77hsL4mHvEzr14wUmTUP3b2ighl+3drEXEtokMVKfB\/2r0HsIKLuNLC5N+gQQN54YUXZObMmTLzQ8c49rvvvitnnHGGeePt4Ycflvfee09MuBPG0Q977723bL99Q7P4S\/wAVn1Y5XTaQF7O8Qlrrg6\/IKDdu19qIgU9irjSwuRftWpVcz4P8\/nWsM1zwYIF8tJLL8kxxxxjjnWwYV27dpWRI0eaM3Oef\/4Fs+c\/Bch9J5pOG8jhQq\/vcFOFFAFFIJgIhMnfq\/6PP\/4oL774ohx++OHCFk9vONf40yl8\/vlnwllA+BWbYcpHF3qDX+vpPPkFv\/TFVwKtb5GY5M+5PpzXY9\/ujdU8mOsPxQosAn\/vQm8y5\/how\/NfwyjmNuy\/2si+RoGubwNP+iwSk\/x5Y5c3d8ePHy\/Lli0z2Xl\/OPKU8N332EOKccsno37vQm\/PA3bwwlThOvgNr0KR1EPRTki1AAAQAElEQVQRUARyikD6LBKT\/Hljt3fv3rJ06VI57rjjZPDgwcICMB0B9q233mr8OeiN4x84BC6nZfdBZrrQ64NKUBUUAUWgUgjEJH9S69KliyF9iP2BBx4QOoODDz7Y2MOGDzcLxLz9m483fNEv30YXevNdA8Wev5ZfEag8AnHJPxQKmS91zZgxQziznw6gT58+gj319ddl6tSpwnbPymcf3JhM+aS90Jv+tF1wAfSr5lonfq0Z1SvDCMQlf5sX20CZ\/z\/++OPl2muvFWyu8bcyxWZ7F3orvtGbBCLpT9slkYmKpISA1klKcKlwjhDIwqAkTP4cz8AxDRjcGNzJHPWAHPI5gsEX2cxYsDZCj4pv9EYEZ\/tC01cEig+BLBCib0HMwqAkTP5s2eR4B8z69etkvfu4h\/BRDmXHPriPenDcpXHW+xa3TCsW7ejmnh2bZDobTS9QCBQTE\/mkYrJAiD4pWU7UCJM\/xzNwvAOmQYPtxV7b4xys\/cHM0qMf7DV2aZwGOVHYD5mM8Zzjs3O9mn5QS3XIKwLKRHmFXzMXSRGDMPkzf9+wYUOzgwc3hmt7nIO1GzdqJNZt7MaNwnFSzDuQ4iz0euf7dconkFWpSisCRY1AmPyLGoUUCs\/efrd4pRZ63QmoWxFQBBSBPCAQJv9Vq1YJh7XZ7\/QmtFu0NIe6IUc84udB\/5xneffriyPy1K91RcChF75AoFwJXYkox0JdkQiEyT8UCskOO+wgO+64Y3Jmp3I54oVChT\/nyZSP3dvPiB8o9WtdoKDGrwgU\/l3pV+T9r1eY\/JnfHzdunLzzzjspG+IRn+IW8khjxsKfKKIxdAR0AJ1bb2eu9SeLCBRyo8oibJq0IhAPgTD5xxNKFMbJn1u2bDFihTzSGDsr8mtdPTsmPsTNgJK\/nxRy9jHDFnKjSqGGVFQRSAqBJG\/luORfUlIin376qdx0003mYy6HHHKIeA3n\/3Tv3l1++OGHpPQKqhAjfTvlY8vA8c3WHVS7vJ0owwa1DlVvRSACgSRv5bjk\/+abb0qvnj3lmWeekS+++EK+++47WeoYbGs4859F35o1a0bkX2gXM1xTPpStUKZ8TDsp7wEomhpFQBEoJARilCUm+f\/222\/CXH6jxo3l5Zdflo8+\/tgc4nbJxRfL\/\/73P3njjTfMoW\/NmzeXe+65R+rWrRsji8LwLugpH9MDFEY9aSkUgWgI6PimIioxyf\/nn3+W+fPnm+mePfbYQ+rXqyft2rWTuf\/9r\/z+++\/SqlUrufvuu02KI0aMMHah\/hTqlE+h1peWSxHwIqDjGy8iEvszjlZ09913t05p3bq1LPn2W2GqB0++3sU3fjniYcOGDXgVpNEXuwqyWouwUFpkRaAcgZgjfz7gss0228i3Dtlb8V122UX4ctf3339vvcw7AStWrJBff\/017FdoDu+LXYW1y6fQakvLowgoAskgEJP8mcNnIfftt9+WtWtLjy9u2rSp0CmwA4jEN2\/eLB9++KHxC4VCeBWc0SmfgqtSLZAikHUEgrDGEJP8Ifnjjz9e5syZI7169ZIff\/xRWrRoIXvuuafcdtttcv7558sZZ5whzz33nOzprAkwBZRdRPMDZ4FO+WS3qjR1RaDIEQjCUDgm+VN3zOffe++9Zq6\/evXqUqtWLenXr5\/svPPOwjbQWbNmya677ipXX3ONGf0TJ3smP3B69\/brWT7Zq2FNWRFQBHKHQJV4WYVCIbPb55FHHpHatWsbUXb+TJkyxUz3MOUzadIkadasmQkrxB\/vF7v0LJ9CrGUtkyJQfAhUqUyRq1WrZs7w5zwfzv2vTBpBiOOd7y+UF7uCgL3qqAgoAmUIZGnGO+rInwXeZ599Vp566inhZS9U4OyekSNHyl5772WOcmbP\/y233CLr1q0juCBNhfl+H36xK0vtoiDrUwulCAQSgSzNeEeQP2f5TJgwwZzfc91115o3fC25P\/300wLZ\/7LhF\/MlL3YD8XJX3759zfd+UwM1GJTl3eLZuXW9JIqZ27JlqV0kUU4VKRwEtCTFiEAE+bOz55Zbb5EaNWrItddeJyOfeUbq168vS5cuFUb9TPcMHjxYeKmL+f5rr71WZsyYIS+++GKK2AWDspas3RRRruQOcgtG2SIKpheKgCJQdAiEyZ9R\/5gxY0ScgSvTPX369JFmTZsKc\/qc5bN8+XLp1KmTHHrooQakUCgkp512muy1117m\/H87PWQCC+BnyY+bBGOLovP9Fgm1FYFkEHCIJBmxQpAJaFHD5M\/0zoIFC+Tggw8W9vK764TRPZ1Dp04Hmu2eNmzbbbc1Z\/xwBhBnAVn\/QrArnOLpw\/n+LOCsSSoCGUKgiJ6AA1rUMPkzcudFLo5m5gUv2wIgdaaDmPLZY4\/drXda9uLFi6VHjx7Spk0badeuvfCUwRHRiRLlTKEzzzzTLDjz9nHLFi1L3Y7N+weJ4qcS7j3FM7n5\/lRyUFlFQBFQBKIhkJtHiTD5h0KhqC9qrVq1SpYtWyYc7cBJnm5VOd5h06ZNss02W0eN65a17q+\/\/toQ\/4oVK2XAgAHSt+8l8sknH8t5550nTC1ZuWg25M9ZQ3vvvbdJo0fPHmF7\/44do0WptF+F+f6W+rnGSoOpERUBRSAFBHLzKBEmf3bvtGjeQhZ\/s9js3rGafvDBB7JmzRozt8++fuuPzYFuc+fOdUbfrZI6z5+pI9YVeIoYPnyY9O7d2yH\/vvLUU0OEtDgqgnRjGWTYhsqTwh133CFu08NZf\/DGK\/F6JHnNXD\/Gipv5\/jZK\/hYPtRUBRSD4CITJnx0+h\/31MPnsv\/+VBx98UDii+fPPPxe2c1arXk34VKN7Omjjxo3y2GOPmdF6t27dkhr5Q9yzPpwlBx54oPARGAvfbrvtJvvv31E4LoJ8rb\/X5skgFAqZk0S9YdGuK9t\/BmF\/f7Tyqp8ioAgoAskiECZ\/Ihx77LHSuXMXGT58uBnpc838\/EknnmS+2oUM5MyXuw477DDhRbDOnTuHdwARHs+sXLlSVq76Xlq1biWhUDk10\/G0a99Oli79Lu5LYwsXLpR69evJ7Nmzzc6j1q1ayVFHHSXvv\/++8FQRL+9Uwiqc55PU\/v5UclBZRSAICKiOhYxABPnXqVPHjOYfeOAB6dKli+yzzz7y0EMPmZe7mKoBCN70ZSrohx9+kJNPPlnY9088whIZ4m7Zsll22XmXCqJ1atcxxE+6FQIdD9YX5n01T5Z8u0QmTpwoV1xxhQy8caD5sMw555wjU6a85khl5m\/GgrURCSW3vz8iil4oAoqAIuBrBCLIH005uZOjnHmpa\/z48cLo3xI\/4awNDB06VDjTnycArvFPxvz2++9S8meJuNOz8UK8YFD+MGC9wzYfi9n06ybTKY0aNcq8Y9C791kywdGxffv28tDD\/5bVq1eH5d0OszOoZUvhKSaRmTF3vnhH\/k2rrk0qbqK0sxburNNkLe3Fi9MuOzu5\/KxfMeimdZB+O65MO7Hcg+3mJD+4K5B\/IqVCoZB569ee8plI3h1eo3p1M93zxx9\/uL2N2yzOlphfc+39oZMZ+vRQ86ax+0mjwfbbC2sOS5cskW+++cYbzVwvWrRIMC1atDDfJIhnyzYNTRz7w2JvPHlfhDkL9VnWIyFu8fLfcccd04ofL20NS9ymwcjPddA8ifuSMgTGuO5HeMcayyl+sVMm\/3QU521hFo3ZrulNh7WEunW3lQYNGniDwtdbtmyRaB0Hn5fctGmj8E5CWLiSDu+oX8\/vrySQGi2LCMQeJGUx06wlHeeBP2t5ZjXhbBYog1WfU\/Jv0qSJNG7cWFi4dS\/Q8oIZ8\/nNmjWNuWX0s88+k3333VcefvjhCvX2xRdfSL169c27CBUCU\/SoMN+vi70pIqji2Ucgm+ySfe01hzQQyGDVp0z+aajtEHQ9OeCAA4QFY\/cUzZdffikfzf7I7OCJNZ3UvHlzadO6jbz++uvmpTOrB3OZr732mjmSgsdC619Z2\/ty1856rENloYweL4Mjl+gZqK8ioAgkg0BOyT8UCpkdQoz0zzjjdPMOAbuFzj\/\/fHNM9CknnxLWedq0acIiyX+eeML4Mc9\/7nnnCieM9j7zTBk7ZowMHzHCLPxyLtHll18ece6QiRTxk5h1eLELY6Mx39+59Xb2Uu1MIJDBkUsm1NE0FIFiRSDr5O+lXD4DOXr0aNlpp51l0KBBZqso0znDhg2TZjs2i6wHJ7LzF\/bjRTO+K1C9Rg0ZMHCg3HXnncJOH94M7tChQ1guuiMx6+jLXdGRU19FIBIBvSoEBLJO\/tEot23btuZDMZwG+tVXX8mTTz5Z4a3dI444QhYtXiQXXnhhGOdQKCQHHXSQTJ482awbEHfIkCFmJ0lYKA1HhcVene9PA02NWhAIuEdfBVEgLYRFIOvkbzNK2c5Do6uw2NtKp3xSrjeNUFgIRBu9FVYJi7Y0\/iX\/PDQ678i\/aFtF8gVXSUVAEQgoAj4k\/zwM+Z3Kcy\/0Opeii72goEYRUAQKFQEfkn8ehvxO7epirwNCNv\/y06dns0SatiIQaAQyRv6BRsFR3jvl09m3i70BZdH89OlOzeqfIuATBHx26yr5l7WL4Cz2KouWVZlaikCwEPDLrVvWCSn5lzUf78i\/zLvcKgOs3ENdgUNA6zDLVabJBwKBsk5Iyd+praQWe8sAc8T1L6gIaB0GteZU7ywgoOTvgKqLvQ4I+qcIKAJFhYCSv1Pd3ikf\/y72OsoG40+1zBYCOnWVLWSLLl0lf6fKg7PY6yhb6H9KbvFruFimroqiHeS3kEr+zq3mPcbZ8dK\/fCFQLOSWL3yDkm9RtIP8FjLr5O\/3tsZiL8bqqW\/2WiTUVgQyg0B+x7eZKUPhpVIiRU\/+uthbeM1aS+QvBPI1vtVOJ147CCn562Kv6D9FIMsI5Cf5fHU6+Slt6rkW\/cjfu9jLtE\/qMGoMRSD\/COhIN\/91ECQNip78vZWl3+z1IqLXQUFAR7pBqSl\/6FnU5M9Cr3vah1F\/5wTf7NXRVVoNVyMrAoqATxAoavL3LvYmUyc6ukoGJZVRBBQBvyNQ1OTvHvVTUTrlAwpqFAGLgD7nWiQK0c4b+fsBTO9ib2ffnuHvB7RUh+JDQJ9zC7nOi5r8vRXbpdV2Xi8RHfxUxCQPPtmvhuznkAfYNEtFICYCRUv+SS\/26uAnZuPJZUD2qyH7OeQSr2DkpVrmE4GiJf\/KLPbms6I07yJCQB9Ciqiy81fUoiV\/XezNX6PTnBMgoA8hCQDS4EwgULTkr4u9lW8+WRiYVl4ZjakIKAKVQqBoyd97jHPUxd5KQVr4kXRgWvh1rCUsfASKkvxZ7MXY6k3mzV4rG0hbh+qBrLbESmvFJsZIJWIh4Dvyj6VoJv29i7325a5c3kq5zEt0qJ7J5uOjtLRifVQZgVOlKMnfu9hrX+7K5a2Uy7wC1ypVYUVAEcg6AkVJ\/t7FXp3vz3o70wwUgRQQyKVoTp\/Bc1mwdFLQywAAEABJREFUhHkVHfkz1+8d+e9Uv2ZCoFRAEVAEChGB4n0Gj0n+hdofVpjvd4ifBd9CbNZaJkVAEVAEYiEQk\/yD0B9WpoPyjvp1yidW0\/CFvyqhCCgCWUIgJvlnKb+MJluZDso739+z4w4Z1UkTUwQUAUUgCAgEmvxTBVjn+1NFTOUVAUWgUBEIDPlnogKKd76\/MhNkmUBc01AEFAG\/IhA48l+8eLH07NFDWrdpI+3atZML+lwg3333XVL43v364gi54pnyqcwEWQRUeqEIKAIFhkCgyP\/rr782xL9ixQoZOGCA9O17iXz0ycdy7rnnyfLly+NWTbQpn3QXe3U8HRdyDVQEMoyAJpdJBAJD\/iUlJTJmzBjZqlo1GT5iuPTu3dsh\/0tlyJCnZOX3K+W5556Li8uMhT9FhLO9s3PrKF\/uipCKf6Hj6fj4aKgioAj4F4HAkP\/atWvlww9nyUEHHSi77NI8jOhu7XeTjvvvL7NnzZYNGzaE\/d0ORv3FO+XjRkLdySKgT3XJImXlFDGLRFDswJD\/ypUrZdWqldKyZSsJhcrH3DVq1JB2u7aTpUuXyLp166LiPnb2SqEDsIGM+vt3K+9ArL\/agUEg64qWt7CsZ1UgGShiQavI3JJ\/KoMDj+yWLVtks2N22WWXChjXrlNbfnaI\/4cffogIg\/TXHXm3eEf9nVulN90TkYleJETg8MMPTyijAtlFQOsgu\/gGMfXckn8qgwOP7O+\/\/y4lW0qkmjPnXwFoI2t+wkEQ\/yVj58mfW28f9sNR5dc18toNRztPEC3VtMwNBuDeMkd5aT7R67So66BFdExy3VaoAz+Z3JJ\/GiWvXr26hKqE5I8\/\/ohIhYuQ85QQiuR+vCsYpnsm9vubLFq0SI1ioG2gWNrAYv\/c7xVIKY8egSH\/qlWrStWtqsqSb7+tANe6Deulbp060qBBg3BYz45N5NGe7c01pM\/1yxfvI6ns8HH6FBM\/mD+xtY8dEsySqtaKgCKQOgJVUo+SnxhNmjSRJo0by4KFC4Vtn1aL3377Tb766n\/StFkzqVu3rvU2NoQ\/Z8CBgqEjoBMwAUn+xH+Y8DuFxtY+dkiSwKhYDAT83iZiqB14by1AZRAIDPnXq1dPDuh4gHzwwQfyzbffhMv65ZdfykezZ0unTp2kdu3aYX\/rSJXwbbzEtlJoYoyKTULbRLHVeJDLGxjyD4VCcvIppwgj\/TNOP0NGDB8hgwc\/Iuedf740btRITnHCglwRwdBdR7bBqCfVUhFIjEBgyJ+i7LHHHjJ69GjZeaedZNDtg2Two4\/JfvvuK08PGybNnGkfZNRkEwHfj2yzVHjt9LIErCabRwQCRf7g1LZtWxk7bpzMnz9fvpo3T5588knZcccdCQqsUWrxe9UVa6fn93pR\/SqPQIkEjvyTKSwnf\/bo0UPatGktnPzZp0+fpE\/+TCb9TMsUGrWwIP\/+++\/Lcccd59RBG2ndqpUcddRRgh9hmcZP00uMwLvvvmvuhWnTpiUWVomMIbB69Wq59tprZc899zTvFXXp0kVefvll4aXVjGVSqYRChUf+9uRPTv4cMGCg9O3bVz75+BM577zEJ3+K\/ssIAlOmTJFzzjnHpDVo0G1y+x13GPc555wrrzlh5qKyP\/qYlDJyy5Z9J7feeqv8\/tvvKcfVCJVHANxP79VL6HAvuKCPPPDAA9KiRQu54oorZZwze1H5lDMTM7sj\/5LMKJlsKowqOfmzarVqMnz48LKTP\/vKU0OfEjqDRCd\/JptPdLkcFza6Enn3Xb9+vQwZMkQ6dOggo0aNktNO6+GY02T0qNHSvn17GeqEIVNpRQvtManSQCQXkZciH3roYVmzZo1Uq14tuUgqlREEnn\/+BVn700\/mfujb91I5\/vjj5cmnnpSDD+4iEyZMkLTug4QaJuaj7JJ\/jm9UTv6cNWuWHHTggdK8efMwPObkz44dhbBYJ3+GhSvtyHFhK61ndiPSyUI0PN7WqVMnnFm9+vWkY8f95QeHhDZu3Bj2V0f2EGAwNOzpYfLpp5\/KoNtvl5o1a2YvM1+mnD+lfvzxR5k69XUz8OnQoUNYkVo1a8lwZ2A6fvx4cd8fYYGMORLzUXbJP2MFSS4hTv78\/vuVZo45FCovvDn5s92u8t3SpTFP\/kwuB5VKhAAL8tOnTzfTbW7ZjZs2ClNytbbZRrbaait3kLqzhMDcuXNlxMgRctVVV5odclnKRpONgsCyZctk2bIVsvvuu8vIESOlU6cDDS\/5ae2roMifRZTNm7fIztFO\/qxdxxC\/9+TPKPWmXllA4O233paZMz+UbkceKfXr189CDpqkGwEWGq+\/\/nrp1q2bdO\/+d3eQunOAAE+\/Gzasl9sGDTJTPRdddIHcettt8tvvvznrYefKtDffzIEW8bMoKPJnQYtH3agnf0pISuJjoaFZQuC9996TG2+8UQ44oKOce+65WcolbrJFFcg8\/4MPPijbb7+9s7h4hYRC5U\/BRQVEqoXNIEEwEC0p+dNMtbG427v3WdLLWfwdN3actGzZXP7z+ONZnvNPXPiCIv\/qNaqbhk7jr1j0DNZsxcTVJwoCJSUlMnnyZLnoooukVetWcv\/992dxnlPr11bBpFcmybRpb8iVV16ZRbxtbgVkZ7yPDEnXQw6JeA+pYcOGctRRR8uib76Rpc40dD7RKyjyr1q1qplP\/jbKyZ8s9G5bt27EyZ\/pAq90ExtBRj5PPfmUXHbZZdK1a1fBTcOPHSPdkFC6CRREfNr5C8+\/IGvW\/CAnnniiM8rkLPsW5p0Ldpfwzsvee+8tn332WUGU16+F2HbbbZ1Rfy1p3LhxBRVZ89qyeXPe9\/pXqaBZgD2aNGliwF4Y9eTPr8zJn96TP9MprtJNdPR48rr77rvl3vvuk7PPPluYgsjuzoboehSjLzt6\/vWvf8kjjzziMoOl\/zXXSM1aNeWCCy6Qe++9V3bcKdhvxfu9biH9evW2Ey8XoTdrARxCSQfBdb5MQZF\/PU7+POAAeZ+TP53HKgsqJ3\/Onv2Rs+Ie\/eRPK6d2+ggw1TNs2NNmO9vVV\/eT6667LvrX19LPSlOIggCjyk6dOsnRRx8dYQ7q3FmqbVVN9ttvXznSWXSvt129KLELzCuPj+ZNmzaV\/ffvKG+\/\/baw60rK\/n333Xfy2muvSds2baRho4ZlvvmxCor8Q6GQnHzyyaUnf55xhowYMVwGDx4s559\/vjTSkz9z0sI4c2no0GFS13nsZWvnDTfcIOw6seYe54lg3bp1OdFFM4mGQBE9r8YsajRcMutHJ3z22WeZRDldAB4aO2aM9O7d2+w6vOzyy4U9\/0YgTz8FRf5gyMmfY0ePlp04+XPQIEP+++67rzAa1ZM\/QSi7ZvZHs2X16lXCltoJ4yeY19jZ7WDNq84C8KZNm7KrhKYeHYE8kmF0hQrbt0OHDsKpAtiQ\/8CbbjLbnIcNG27egM936QuO\/AG0Tdu2hnTmz18gX331lRTCyZ+UKwjm9F6nm2\/jLub7sFG+ncoLYDyFBaEshaQjB4vNnTNXjjjiiEIqlu\/Lwlk+HHcCDy2YP1+ef\/55YYDqB8ULkvz9AKzqEAgEVElFICYCeVwyiKlTJgOKkPwLvUoz2Tw0LUWgeBEo9FmyIiT\/Qq\/S4r1ZteSKgCKQPAJFSP7RwVFfRUARUASKCYEMk79Oqfin8Whd+KcuVBNFwH8IZJj8dUrFP1Ucpy60X\/BPNakmqSOQtfabuipBjpFh8g8yFEWke5x+oYhQ0KIGFQFtvxmpOSX\/jMCoiSgCioAiUI4A3+3lQMNVq1aVe\/rMpeTvswpRdXyBgCqhCFQaAQj\/tttuq3T8XEVU8s8V0pqPPPHEE2VHDHPMcHzDhy9mzZolHD\/cr18\/n6CXm8lmjl4+88wzDV6pFJxjmg8\/\/HCZM2dOKtFUNsMITJgwoew7CiFff0BKyT\/DFa\/JxUZg1113lR49eoTNSSedJJzxz\/G2xx9\/fNgfmcMOO0z4PkPs1PIRkv3JZk5FHTVqlKxZs0bAJ5VScmwAUw0coc25\/qnEVdnMIEAH\/OmnnwqHuYlD\/dlvMZXXW8k\/AXYanDkEIPQ77rhDrLnpppukVatWwlHc1157bdifcE5i3W+\/\/cxxuPfdd1\/mlPB5SpyK+swzzxjyoGNMRd1QKGQ6UNKYOnVqKlHTlM3NE1GaSiaMzlHLf\/vb34RpG68wnSlTOZyRRJvlWOy33npL6KzdsiNGjJBTTjlFGNC4\/f3oVvKPUyuF0aTjFFCDfIUARDJu3LPCB4cOPvjgSunWpk0b6dy5s9CBrF27tlJppB7Jz+PbxKUBdxZoGYD8\/vvvToTIO3\/jxo3mk5ic0MmR8QxOatWqaT5P+ua0Nx350j\/SwBWUw\/OU\/KmtGCbYTTpGoXzrHXnDoSaP0NHm\/JctWybcqIzCWrZoKR07dhS+E8DojHgY3L16nS6sF8yePduMxiDGdu3ayaWXXiqrV68WPqzBZw3xIwz38uXLiW4MN3PLli3l5ZdfFu+oj1EipGEEy374gtkLL7wgPOEwOiRNRoGsXXhly6JEWOjzxhtTzUeH+Pi6O5BvIzCVgK7otM8++wgYUA63XCgUkn\/84x8yb95X8t5777mD1B0FATrIm2++2RA5ay2lIpF3\/qxZs+Sdd96RW265RXhaPe2002TUqNGyT4cO8vgTjwvxaG+c2Ml5\/aVp+P83u+Rf8X72PyK+1LDAgYy81+LWAF9FOvGkE+XVV181H+55+N8PyyGHHCJPDRkiF198sbkRyxMokY8\/\/lj69u0rO++8s9xz771GdvLkyQLRc6PWqrW13HXXXXLUUUfJm2++KRABI73yNEQGDRokEydOlAv6XGA6gS1btpg0OarXykH8EEP\/\/v2F70YwOhwwYICZQmDxlpF4og6Ahdrvv\/9eunTpIqFQOSiff\/65nH766cJ0znXXXicPP\/yw\/PWvf5Xx48ebckBgVg9sOp2GDbc3hJUoT+QLzYAjHTykbMtm62fSpFcipmpuv\/12GTV6lPz9738XPnFZKh95v\/E1rh122EEOPPDA0mDnl8+SHnPMMbJg4QL55ptvZPHixTJz5kw57rjjzKYG2tfSpUuFaSQGMU4U3\/2Vkn9kWTOnZHn7zVyaRZlSgQNp2l+sMpb7b9y4SR599FHZ5NgjR440ozBuwPvvv18eeOAB+fDDD2Xo0KERLWTJkiVyySWXCOsGJziLyny\/dq+99pI5c+cawn\/ooQeFxWbSgADmzpljngjciUCgEH3fS\/tKr169hNE9aYwc+YxwgyM7ZcoUefbZZ4VR4bBhw4xN5wJBt2\/fXh597FFD3sjGMp988onZJcInAN0ykM9vv\/0m\/\/73v6X3Wb2FMlOeCy+80HQuCxYscIsL34Zt0aKFfPHFF+LtGCIEC\/RixYoV8sYbb5g1JDpyOmvqb9y4cbJy5YqIUtevX18efOBBYZEc3MQ0N\/Nj5H755RfzHd4dmu5g6sZ4lqFefJkAAAgJSURBVP3w3ZA\/fv\/DtBeeQhmYLOI7Fo7hGyJ8UAo9CCuL4iurlPzLy+or5VQZRYAdExaFBQvmG4I\/stuR0sF55Lb+2Iz+d9ttN5k+fbr8\/PPPeBnTpEkTOfTQQyUUKm3k3OBMydSpXdt82CQUKvXns3vtnCmhTQ7Jer80BqGX5WfSZIG695m9Zdmy78zoGmKeOPElIS8IuVq1akaOH6Zv6ARWf79a3n33XbyiGksyLBQSxytEHpA5REZYKBSSq666SmbMmGGmvfCzZpttthE+mLPSeYpYuXKlFNu\/7t27C538iy++KLfecqs8\/vjjhtzB67zzzwu3BXC5\/vrr5dhjjy3fWWYGIoSUGjr+zZs3S5PGTQRcS31Lf2vWrCnVqleTb7\/9ttQjYL+l5B8wpVVdHyLguWlS1zBxAkyJ8Ci\/ZfMWYeqGqZ9SM9mQ8J9\/\/mlGwm7yh9S5Sb360Ak0bNgwwhvZ0q4gwrtCR0Noq9atzEiQKZl169Y5j\/3fmE+HMpIk3G1MZ1O3jjMPP8\/tHeGGYCgbU0Zbb711RBidV10nPtNIdEL\/+te\/hPUI5pkjBF0XdESb\/\/hDmO5weefAmbges61EKBQy6x4Q+4svvWieClmoPfvssyXk\/E8lf9oUg4EqVapIKBSqGDVGcVn0ZSBCJ1wxkj98qvhDDdUi8AhEuS9SK1PiBOz0xksvvWQWbZnXvbTvpY67r1x22WXC3OqGDb9EjPx3bNZMvGSKXtzMkD1ut\/Hey4zEvZ1EuXzIkKsdHZJmeVhqLkb2P\/74Y9RI7N9\/\/vkXwvv+X3nlFTPXz6LvLc4i5IYNGyrEq127juDvXRCuIJhxj8T1GC1LL+7RZFL1g7Tp\/EKhUKW3XlKnDB7oBKjnVHXIt3w8XJX8K1k7Gi33COyyyy5CY2bnjZ1bXbR4kYTdzlzr3LlzxD3HinwFTaN6VpAyHrEItJQM\/jSj\/1AoJHQk+JlIlfipUaOGRHtqsEk1b97cTGWwmPn666+brYcNGjSQESNGCPPZVs7aGzasN4QXu+Oykv6wK9dlRNcdkmathXWgU089Va644goZPny4wclOm0WPaX3LG0goFJKqW20lK79fKUzNieuf7Vxoly5vXznj4ark76uqKlJl4rVQKQ9kLryWM88KAXKDu9Fiqoc3Ynv06CGxRtBh+fIkw17xHPaJwy2zePFiM7JmKyr78llgpRNisdEth3vhwoVmFxILv1xHM3Qe7CBZtmyZ\/Prrr2ERCIcnHBalKVfVqlWF3TzsYHrqqadM58MW0XCEMgdz\/aw9YMq8isZi8X3AgIFywgknyI033mi2cdIBsLDLhgBv26kITHkDYZ6\/datWsmL5CqcOI5+w5n\/9ten0vQv0FdPzp0\/eyb+8j\/UnQKpVDhCI2wjKAyFP5rx5e\/Wjjz4KK8bNPHHiRKFTYNGWBdlwYAYcjCIhZZsURy8w4ma\/Pe8YMGpni98PP\/wg7PRhqsEr26B+A7OF0\/p7bUiGtQGeNEjfhlv\/efPmmcVd64\/NGgF50WlwbQ1p8L5C48aNhbl\/618sNtsyjzqquzDnX6tWLbOYyzsSPXv1FHbghELl5B4dk\/I2R\/ihhx4qdOoffPA+l8aA\/aRJk0xHTDswngH7yTv5J6qGgOGp6mYcgfIWUrt2bbnmmmukevXqcs4555iXbiZNmmRe5GIqiJuQYyFCofI4sdRJLFEek+2ibPGE8DE8Ycz7ap7ZbcMCLZLsMGFX0LPPPStnOwuLbPu0siwK97mgj7Rt2xbRmIY3c3\/66acKu0d4UYz3FCAz5vhZ5B48eLCw8Avx9+zZMyJNnhDYgbL77rubozMiAovgggHCAw88aJ6KbHF5AgI7tvNav9h2ZOvgmJG\/\/OUvpjPhXQ7q9owzzpBP58wROhXqIHZa\/g2pEtnH+VdR1SwFBAqqUiMLw43Ngu\/RRx9t9tuz0MuTAIQ8ZswY84JVMkhFpho\/BtsG2b3B3npe+GIHxzMjnzEvBtmYllxuuvEm4SkBor7zrjvNyHv06NGGJEKhSFKxca0NWUPy\/\/d\/\/xfxIhIdzHBnzpoXu3jHgGmg\/\/znP2arKidIejsV3gbmKYQ3fUOh+HnavNUW1wSjRPxj0PHII4+YdzzAm7rlnRO2kNLpRwgH6KKKNo0A1VayquahUlMhU1sMbioIe\/r06WZfuvU3q7rOBQu3c+fONS9oOZfhP8iQt3LZ3cM8OzbX7sVNmzbp4w5HdhyQeIU8HX\/26JMf+TqX4T92Cw0cONDsJuItW17jZ7onLCBinHQAvM3LgV\/M83817yvz4tcBBxwgoVDiSmGKBsLmTVHvPD5l5s1eyuouM\/4m87IftozyNERHsu+++5b5qpUMAhdceKF5T4TO3SvPug7rB+BP3U6d+ro5xiMUSlyv3rTydu25SfM+7ZM3IDTjjCKQ0Vsgo4lltJhZTSwUCpm3d\/\/4\/XehA6lMZgucxeUPPvjAHAcR1OmIypRb4ySBgOe+UvJPAjMVCRoCniFOgNRnJ0\/vs84SjiJIdY8+C9\/PPfustG7dWjhyOEDFVlXzgICSf4ZB1+T8gIBniOMHlZLUISQhYTGRba2sbSQZzYixsMx0Fovi3qkuI5Dpn+D2sZlGIpDpKfkHstqyq7Te06X4ssjL\/Dp2qU8OfkNidqmMHDlS2LmUSo6sVXAyKYviqcSrtKyja6XjasS8I6Dkn\/cq8J8Cek\/7r05UoyAhEAxdlfyDUU+qpSKgCCgCGUXg\/wEAAP\/\/Mrlw\/wAAAAZJREFUAwDySW3HMBVB9AAAAABJRU5ErkJggg==","height":191,"width":255}}
%---
%[output:58ace067]
%   data: {"dataType":"textualVariable","outputData":{"name":"zita","value":"0.9000"}}
%---
%[output:4ef2414a]
%   data: {"dataType":"textualVariable","outputData":{"name":"settling_time","value":"30000"}}
%---
%[output:71e936eb]
%   data: {"dataType":"textualVariable","outputData":{"name":"ts","value":"2400"}}
%---
%[output:28f455f4]
%   data: {"dataType":"textualVariable","outputData":{"name":"wn","value":"0.0019"}}
%---
%[output:371c3460]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"PD","rows":1,"type":"double","value":[["1.0000","0.0033","0.0000"]]}}
%---
%[output:2d094491]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"PD","rows":1,"type":"double","value":[["1.0000","0.0200","0.0001","0.0000"]]}}
%---
%[output:9ecc9272]
%   data: {"dataType":"symbolic","outputData":{"name":"PD","value":"$$ \\left(\\begin{array}{ccccc}\n1.0 & 0.03667 & 0.0003923 & \\text{1.04e-6} & \\text{9.526e-10}\n\\end{array}\\right) $$"}}
%---
%[output:4247710a]
%   data: {"dataType":"symbolic","outputData":{"name":"PD","value":"$$ 1.0\\,s^4 +0.03667\\,s^3 +0.0003923\\,s^2 +\\text{1.04e-6}\\,s+\\text{9.526e-10} $$"}}
%---
%[output:57829534]
%   data: {"dataType":"textualVariable","outputData":{"name":"NP","value":"4.7380e-10"}}
%---
%[output:64134517]
%   data: {"dataType":"symbolic","outputData":{"name":"DP","value":"$$ s^3 +\\frac{3142864021558265\\,s^2 }{2305843009213693952}+\\frac{5083627518809173\\,s}{2361183241434822606848}+\\frac{612925390544617}{2417851639229258349412352} $$"}}
%---
%[output:7193e880]
%   data: {"dataType":"symbolic","outputData":{"name":"NC","value":"$$ {\\textrm{Kd}}_2 \\,s^3 +\\mathrm{Kd}\\,s^2 +\\mathrm{Kp}\\,s+\\mathrm{Ki} $$"}}
%---
%[output:03907c0a]
%   data: {"dataType":"symbolic","outputData":{"name":"DC","value":"$$s$$"}}
%---
%[output:3e0e0278]
%   data: {"dataType":"symbolic","outputData":{"name":"PC","value":"$$ s^4 +{\\left(0.00000000047379999999999995617945672987019\\,{\\textrm{Kd}}_2 +0.0013630000000000000622002449546244\\right)}\\,s^3 +{\\left(0.00000000047379999999999995617945672987019\\,\\mathrm{Kd}+0.0000021529999999999999692765293574093\\right)}\\,s^2 +{\\left(0.00000000047379999999999995617945672987019\\,\\mathrm{Kp}+0.00000000025350000000000000348407182282091\\right)}\\,s+0.00000000047379999999999995617945672987019\\,\\mathrm{Ki} $$"}}
%---
%[output:5c8c94a2]
%   data: {"dataType":"symbolic","outputData":{"name":"eq1","value":"$$ \\left(\\begin{array}{ccccc}\n1 & 0.00000000047379999999999995617945672987019\\,{\\textrm{Kd}}_2 +0.0013630000000000000622002449546244 & 0.00000000047379999999999995617945672987019\\,\\mathrm{Kd}+0.0000021529999999999999692765293574093 & 0.00000000047379999999999995617945672987019\\,\\mathrm{Kp}+0.00000000025350000000000000348407182282091 & 0.00000000047379999999999995617945672987019\\,\\mathrm{Ki}\n\\end{array}\\right) $$"}}
%---
%[output:1a5892f2]
%   data: {"dataType":"symbolic","outputData":{"name":"eq2","value":"$$ \\left(\\begin{array}{ccccc}\n1.0 & 0.036666666666633318527601659297943 & 0.00039231824417029770302178803831339 & 0.0000010402377686319352889654510363471 & 0.00000000095259868922379351864615504297262\n\\end{array}\\right) $$"}}
%---
%[output:6844c8f0]
%   data: {"dataType":"symbolic","outputData":{"name":"eqs","value":"$$ \\left(\\begin{array}{ccccc}\n1.0=1.0 & 0.00000000047379999999999995617945672987019\\,{\\textrm{Kd}}_2 +0.0013630000000000000622002449546244=0.036666666666633318527601659297943 & 0.00000000047379999999999995617945672987019\\,\\mathrm{Kd}+0.0000021529999999999999692765293574093=0.00039231824417029770302178803831339 & 0.00000000047379999999999995617945672987019\\,\\mathrm{Kp}+0.00000000025350000000000000348407182282091=0.0000010402377686319352889654510363471 & 0.00000000047379999999999995617945672987019\\,\\mathrm{Ki}=0.00000000095259868922379351864615504297262\n\\end{array}\\right) $$"}}
%---
%[output:23229880]
%   data: {"dataType":"textualVariable","outputData":{"name":"Kp","value":"2.1950e+03"}}
%---
%[output:247e7b44]
%   data: {"dataType":"textualVariable","outputData":{"name":"Ki","value":"2.0106"}}
%---
%[output:9cb67a9c]
%   data: {"dataType":"textualVariable","outputData":{"name":"Kd","value":"8.2348e+05"}}
%---
%[output:91718d43]
%   data: {"dataType":"textualVariable","outputData":{"name":"Kd2","value":"7.4512e+07"}}
%---
%[output:71758da0]
%   data: {"dataType":"symbolic","outputData":{"name":"TFcont","value":"$$ \\frac{\\frac{5000398828688115\\,s^3 }{67108864}+\\frac{3536823477727487\\,s^2 }{4294967296}+\\frac{2413412401994369\\,s}{1099511627776}+\\frac{4527356586979633}{2251799813685248}}{s} $$"}}
%---
%[output:82160c90]
%   data: {"dataType":"symbolic","outputData":{"name":"TFdis","value":"$$ \\frac{z\\,{\\left(\\frac{12067062009971845\\,{\\left(z-1\\right)}}{274877906944\\,z}+\\frac{88420586943187175\\,{{\\left(z-1\\right)}}^2 }{268435456\\,z^2 }+\\frac{625049853586014375\\,{{\\left(z-1\\right)}}^3 }{1048576\\,z^3 }+\\frac{4527356586979633}{2251799813685248}\\right)}}{20\\,{\\left(z-1\\right)}} $$"}}
%---
%[output:16821bef]
%   data: {"dataType":"symbolic","outputData":{"name":"TFdis","value":"$$ \\frac{1343026164261655689563476273\\,z^3 -4028336569421644712069079040\\,z^2 +4027594744925276410898022400\\,z-1342284339760760031805440000}{45035996273704960\\,z^2 \\,{\\left(z-1\\right)}} $$"}}
%---
%[output:8b4947da]
%   data: {"dataType":"matrix","outputData":{"columns":4,"exponent":"27","name":"numdis","rows":1,"type":"double","value":[["1.3430","-4.0283","4.0276","-1.3423"]]}}
%---
%[output:7c42d01b]
%   data: {"dataType":"matrix","outputData":{"columns":4,"exponent":"16","name":"dendis","rows":1,"type":"double","value":[["4.5036","-4.5036","0","0"]]}}
%---
%[output:580d9789]
%   data: {"dataType":"textualVariable","outputData":{"name":"TFdiscrete","value":"  1.343e27 z^3 - 4.028e27 z^2 + 4.028e27 z - 1.342e27\n  ---------------------------------------------------\n              4.504e16 z^3 - 4.504e16 z^2\n \nSample time: 0.05 seconds\nDiscrete-time transfer function.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 78 117 109 101 114 97 116 111 114 58 32 123 91 49 46 51 52 51 48 101 43 50 55 32 45 52 46 48 50 56 51 101 43 50 55 32 52 46 48 50 55 54 101 43 50 55 32 45 49 46 51 52 50 51 101 43 50 55 93 125 10 32 32 32 32 32 68 101 110 111 109 105 110 97 116 111 114 58 32 123 91 52 46 53 48 51 54 101 43 49 54 32 45 52 46 53 48 51 54 101 43 49 54 32 48 32 48 93 125 10 32 32 32 32 32 32 32 32 86 97 114 105 97 98 108 101 58 32 39 122 39 10 32 32 32 32 32 32 32 32 32 73 79 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 46 48 53 48 48 10 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>"}}
%---
%[output:4e096393]
%   data: {"dataType":"matrix","outputData":{"columns":4,"exponent":"10","name":"numdis","rows":1,"type":"double","value":[["2.9821","-8.9447","8.9431","-2.9805"]]}}
%---
%[output:6bc408e5]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"dendis","rows":1,"type":"double","value":[["1","-1","0","0"]]}}
%---
%[output:9c7803d9]
%   data: {"dataType":"symbolic","outputData":{"name":"fcnum","value":"$$ \\frac{8.943}{z^2 }-\\frac{8.945}{z}-\\frac{2.98}{z^3 }+2.982 $$"}}
%---
%[output:66755bc7]
%   data: {"dataType":"symbolic","outputData":{"name":"fcden","value":"$$ 1.0-\\frac{1.0}{z} $$"}}
%---
%[output:37e07119]
%   data: {"dataType":"symbolic","outputData":{"name":"NUMdiff","value":"$$ 2.982\\,\\mathrm{ea}-8.945\\,\\mathrm{ep}+8.943\\,\\mathrm{epp}-2.98\\,\\mathrm{eppp} $$"}}
%---
%[output:9ee831c5]
%   data: {"dataType":"symbolic","outputData":{"name":"DENdiff","value":"$$ 1.0-\\frac{1.0}{z} $$"}}
%---
%[output:727a23c4]
%   data: {"dataType":"symbolic","outputData":{"name":"DENdiff","value":"$$ \\mathrm{Ua}-1.0\\,\\mathrm{up} $$"}}
%---
%[output:180b0379]
%   data: {"dataType":"symbolic","outputData":{"name":"ec_en_diferencias","value":"$$ 2.982\\,\\mathrm{ea}-8.945\\,\\mathrm{ep}+8.943\\,\\mathrm{epp}-2.98\\,\\mathrm{eppp}=\\mathrm{Ua}-1.0\\,\\mathrm{up} $$"}}
%---
%[output:6479a8d3]
%   data: {"dataType":"symbolic","outputData":{"name":"Ua_solve","value":"$$ 2.9821\\,\\mathrm{ea}-8.9447\\,\\mathrm{ep}+8.9431\\,\\mathrm{epp}-2.9805\\,\\mathrm{eppp}+\\mathrm{up} $$"}}
%---
