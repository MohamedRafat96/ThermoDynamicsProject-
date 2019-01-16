% The Value of compression isentropic efficiency
Eff_Isentropic  = (60:1:100)/100;
To=300; 
Tsource=2000;
%% Graph Thermal and Second Law Efficiency 
%To Calculate T4
% T4c = T4 = T2

for x = 1: length(Eff_Isentropic)
    T4c(x) = (110.62+(300*Eff_Isentropic(x)))/Eff_Isentropic(x) ;
end

%The Constant Value of Temperature  
T1c=300;      T3c=T1c;    %  T1 = T3;
T5c=1600;     T7c=T5c;    %  T5 = T7;
T6c=1168.95;  T8c=T6c;    %  T6 = T8;
for x = 1: length(Eff_Isentropic)
    %Law of Thermal Efficiency
    Eff_Thermal(x) = (T5c-T6c+T7c-T8c-T4c(x)+T1c-T4c(x)+T3c)/(T5c-T4c(x)+T7c-T6c);
end
figure(1)
plot(Eff_Isentropic,Eff_Thermal);
xlabel('isentropic efficiency');
ylabel('thermal efficiency');

for x = 1:length(Eff_Isentropic)
    % Law of Second Law of Thermodynamics Efficiency
    %T2c(x) = T4c(x);
    Eff_SecondLaw(x) = (T5c-T6c+T7c-T8c-T4c(x)+T1c-T4c(x)+T3c)/((1-(To/Tsource))*(T5c-T4c(x)+T7c-T6c));
end
figure(2)
plot(Eff_Isentropic,Eff_SecondLaw);
xlabel('isentropic efficiency');
ylabel('second law efficiency');

%% Graph Thermal and Second Law Efficiency With Maximum Temperature
Tmax = (1000:50:1500)+273; % T5 = T7 -> Maximum Temperature (degrees Kelven)

T6M=1168.95; T8M=T6M; % T6 = T8 
T4M=447.49;  T2M=T4M; % T4 = T2
T1M=300;     T3M=T1M; % T1 = T3

for x = 1:length(Tmax)
    %Thremal efficiency With Maximum Temprature 
    Eff_Thermal_Tmax(x) = (Tmax(x)-T6M+Tmax(x)-T8M-T2M+T1M-T4M+T3M)/(Tmax(x)-T4M+Tmax(x)-T6M);
end
figure(3)
plot(Tmax-273,Eff_Thermal_Tmax);
xlabel('Tmax (degree celsius)');
ylabel('thermal efficiency');

for x = 1:length(Tmax)
    % Law of Second Low Efficiency With Maximum Temprature 
    Eff_SecondLaw_Tmax(x) = (Tmax(x)-T6M+Tmax(x)-T8M-T2M+T1M-T4M+T3M)/((1-(To/Tsource))*(Tmax(x)-T4M+Tmax(x)-T6M));
end
figure(4)
plot(Tmax-273,Eff_SecondLaw_Tmax);
xlabel('Tmax (degree celsius)');
ylabel('Second Low efficiency');

%% Graph Thermal and Second Law Efficiency With Minimum Temperature
Tmin = (0:1:50)+273; % T1 = T3 -> Minimum Temperature (degrees Kelven)

T6N=1168.95;  T8N=T6N; % T6 = T8
T4N=447.49;   T2N=T4N; % T4 = T2
T5N=1600;     T7N=T5N; % T5 = T7

for x = 1:length(Tmin)
    % Law of Thremal Efficiency With Minimum Temprature 
    Eff_Thermal_Tmin(x) = (T5N-T6M+T7N-T8N-T2N+Tmin(x)-T4N+Tmin(x))/(T5N-T4N+T7N-T6N);
end
figure(5)
plot(Tmin-273,Eff_Thermal_Tmin);
xlabel('Tmin (degree celsius)')
ylabel('thermal efficiency')

for x = 1:length(Tmin)
    % Law of Second Low Efficiency With Minimum Temprature 
    Eff_SecondLaw_Tmin(x) = (T5N-T6N+T7N-T8N-T2N+Tmin(x)-T4N+Tmin(x))/((1-(To/Tsource))*(T5N-T4N-T6N+T7N));
end
figure(6)
plot(Tmin-273,Eff_SecondLaw_Tmin);
xlabel('Tmin (degree celsius)')
ylabel('Second Law efficiency')

%% Brayton cycle with four stages
N(1)= 2;
N(2)= 3;
N(3)= 4;

Rp = 9;
for x = 1:length(N)
    Rp_Stage(x) = Rp^(1/N(x));
end
Tmin_NN = 300 ; % 300 Kelven degree -- Tmin = T1
Tmax_NN = 1600; % 1600 Kelven degree -- Tmax = T5
P_hight = 9*(10)^5; % The Hight Pressure = 9 bar 
P_low = (10)^5;     % The Low Pressure = 1 bar 
y = 1.4; % The specific Heat Ratio for Air 
Eff_ise_Compressor  = 75;
Eff_ise_Turbine = 100;

%Brayton cycle with N = 2
T1_N2 = Tmin_NN;
T3_N2 = T1_N2;
T5_N2 = Tmax_NN;
T7_N2 = T5_N2;
T8_N2 = 1168.95;
T4_N2 = 447.49;

%Process (1)--->(2s) In Comperssor 
T2s_N2 = Tmin_NN*(Rp_Stage(1))^((y-1)/y);
T2a_N2 = (((T2s_N2-Tmin_NN))/(Eff_ise_Compressor/100))+Tmin_NN;

%Process (5)--->(6s) In Turbine 
T6s_N2 = Tmax_NN/((Rp_Stage(1))^((y-1)/y));
T6a_N2 = T6s_N2; % The Efficiency of Turbine = 100%

Eff_Thermal_N2 = (T5_N2-T6a_N2+T7_N2-T8_N2-T2a_N2+T1_N2-T4_N2+T3_N2)/(T5_N2-T4_N2+T7_N2-T6a_N2);
Eff_SecondLaw_N2 = (T5_N2-T6a_N2+T7_N2-T8_N2-T2a_N2+T1_N2-T4_N2+T3_N2)/((1-(To/Tsource))*(T5_N2-T4_N2+T7_N2-T6a_N2));

%Brayton cycle with N = 3

T1_N3 = Tmin_NN ; 
T3_N3 = T1_N3 ;
T5_N3 = T3_N3;
T7_N3 = Tmax_NN;
T9_N3 = T7_N3;
T11_N3 = T9_N3;

%Process (1)--->(2)  
T2s_N3 = Tmin_NN*(Rp_Stage(2))^((y-1)/y);
T2a_N3 = (((T2s_N3-Tmin_NN))/(Eff_ise_Compressor/100))+Tmin_NN;

%Process (3)--->(4s)
T4s_N3 = Tmin_NN*(Rp_Stage(2))^((y-1)/y);

%Process (3)--->(4)
T4a_N3 = (((T2s_N3-Tmin_NN))/(Eff_ise_Compressor/100))+Tmin_NN;
T6_N3 = T4a_N3;

%Process (7)--->(8s)
T8s_N3 = ((1/Rp_Stage(2))^((y-1)/y))*Tmax_NN;
T8_N3 = T8s_N3;

%Process (9)--->(10)
T10_N3 = T8_N3;
T12_N3 = T10_N3;

Eff_Thermal_N3 = (T7_N3-T8_N3+T9_N3-T10_N3+T11_N3-T12_N3-T6_N3+T5_N3-T4a_N3+T3_N3-T2a_N3+T1_N3)/(T7_N3-T6_N3+T9_N3-T8_N3+T11_N3-T10_N3);
Eff_SecondLaw_N3 = (T7_N3-T8_N3+T9_N3-T10_N3+T11_N3-T12_N3-T6_N3+T5_N3-T4a_N3+T3_N3-T2a_N3+T1_N3)/((1-(To/Tsource))*(T7_N3-T6_N3+T9_N3-T8_N3+T11_N3-T10_N3));


%Brayton cycle with N = 4
T1_N4 = Tmin_NN ; 
T3_N4 = T1_N4 ;
T5_N4 = T3_N4;
T7_N4 = T5_N4;
T9_N4 = Tmax_NN ; 
T11_N4 = T9_N4 ;
T13_N4 = T11_N4;
T15_N4 = T13_N4;

%Process (1)--->(2) 
T2s_N4 = Tmin_NN*(Rp_Stage(3))^((y-1)/y);
T2_N4 = (((T2s_N4-Tmin_NN))/(Eff_ise_Compressor/100))+Tmin_NN;

%Process (3)--->(4)
T4_N4 = (((T2s_N4-Tmin_NN))/(Eff_ise_Compressor/100))+Tmin_NN;

%Process (5)--->(6)
T6_N4 = T4_N4;

%Process (7)--->(8)
T8_N4 = T6_N4;

%Process (9)--->(10)
T10_N4 = Tmax_NN/(Rp_Stage(3))^((y-1)/y);

%Process (11)--->(12)
T12_N4 = T10_N4;

%Process (13)--->(14)
T14_N4 = T12_N4;

%Process (15)--->(16)
T16_N4 = T14_N4;

Eff_Thermal_N4 = ((T9_N4-T10_N4)+(T11_N4-T12_N4)+(T13_N4-T14_N4)+(T15_N4-T16_N4)-(T8_N4-T7_N4)-(T4_N4-T3_N4)-(T2_N4-T1_N4)-(T6_N4-T5_N4))/((T9_N4-T8_N4)+(T11_N4-T10_N4)+(T13_N4-T12_N4)+(T15_N4-T14_N4));
                    
                    
Eff_SecondLow_N4 = ((T9_N4-T10_N4)+(T11_N4-T12_N4)+(T13_N4-T14_N4)+(T15_N4-T16_N4)-(T8_N4-T7_N4)-(T4_N4-T3_N4)-(T2_N4-T1_N4)-(T6_N4-T5_N4))/((1-(To/Tsource))*((T9_N4-T8_N4)+(T11_N4-T10_N4)+(T13_N4-T12_N4)+(T15_N4-T14_N4)));



Eff_Thermal_Stages(1) = Eff_Thermal_N2;
Eff_Thermal_Stages(2) = Eff_Thermal_N3;
Eff_Thermal_Stages(3) = Eff_Thermal_N4;

figure(7)
plot(N,Eff_Thermal_Stages);

Eff_SecondLow_Stages(1) = Eff_SecondLaw_N2;
Eff_SecondLow_Stages(2) = Eff_SecondLaw_N3;
Eff_SecondLow_Stages(3) = Eff_SecondLow_N4;

figure(8)
plot(N,Eff_SecondLow_Stages);





