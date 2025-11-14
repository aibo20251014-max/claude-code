clc
clear
% 固定电流密度，改变压力
j = 0.35;             % 固定电流密度为0.35 A/cm²
P = 1:0.2:10;         % 压力从1到10 bar变化

t = 80;               % 温度，摄氏度
T = t + 273.15;       % 温度，开尔文
T_ref = 298.15;       % 参考温度，K
F = 96485;            % 法拉第常数，C/mol
R = 8.314;            % 常用气体常数，J/(mol*K)
M_H2O = 18;           % 水的摩尔质量，g/mol
M_O2 = 32;            % 氧气的摩尔质量，g/mol
M_H2 = 2;             % 氢气的摩尔质量，g/mol

P_an = P;             % 阳极压力，bar
P_cat = P;            % 阴极压力，bar
l_anc = 0.125;        % 阳极通道-催化剂距离,cm
l_catc = 0.125;       % 阴极通道-催化剂距离,cm

% O2,H2压强
epslion_an = 0.3;     % 阳极孔隙率
epslion_cat = 0.3;    % 阴极孔隙率
tau_m = 2.18;         % 膜的迂曲度
tau_an = tau_m;       % 阳极的迂曲度
tau_cat = tau_m;      % 阴极的迂曲度

% D_eff求解
r = 4.5e-6;           % 电解槽平均孔隙半径,cm
sigma_O2 = 3.467e-8;
sigma_H2 = 2.827e-8;
sigma_H2O = 2.461e-8;
epslion_O2 = 106.7;
epslion_H2 = 59.7;
epslion_H2O = 809.1;
sigma_O2_H2O = (sigma_O2 + sigma_H2O) / 2;
sigma_H2_H2O = (sigma_H2 + sigma_H2O) / 2;
epslion_O2_H2O = sqrt(epslion_O2 * epslion_H2O);
epslion_H2_H2O = sqrt(epslion_H2 * epslion_H2O);
tau_O2_H2O = T / epslion_O2_H2O;
tau_H2_H2O = T / epslion_H2_H2O;
Tau_D_O2 = 1.06/(tau_O2_H2O)^0.156 + 0.193/exp(0.476*tau_O2_H2O) + 1.036/exp(1.53*tau_O2_H2O) + 1.765/(3.894*tau_O2_H2O);
Tau_D_H2 = 1.06/(tau_H2_H2O)^0.156 + 0.193/exp(0.476*tau_H2_H2O) + 1.036/exp(1.53*tau_H2_H2O) + 1.765/(3.894*tau_H2_H2O);

% 初始化数组
V_oc = zeros(size(P));
V_act = zeros(size(P));
V_con = zeros(size(P));
V_ohm = zeros(size(P));
theta_array = zeros(size(P));

% 对每个压力值进行计算
for i = 1:length(P)
    P_current = P(i);

    D_eff_O2 = 0.00133*((1/M_O2+1/M_H2O)^0.5)*T^1.5/(P_current*(sigma_O2_H2O)^2*Tau_D_O2);
    D_eff_H2 = 0.00133*((1/M_H2+1/M_H2O)^0.5)*T^1.5/(P_current*(sigma_H2_H2O)^2*Tau_D_H2);
    D_eff_H2OK = 4/3*r*sqrt(8*R*T/(pi*M_H2O));
    D_eff_an = 1/(epslion_an/tau_an*((1/D_eff_O2)+(1/D_eff_H2OK)));
    D_eff_cat = 1/(epslion_cat/tau_cat*((1/D_eff_H2)+(1/D_eff_H2OK)));
    X_H2O_an = ((epslion_an/tau_an)*exp(R*T*l_anc*j/(2*F*P_current*D_eff_an)));
    X_H2O_cat = ((epslion_cat/tau_cat)*exp(R*T*l_catc*j/(2*F*P_current*D_eff_cat)));
    X_O2 = 1 - X_H2O_an;
    X_H2 = 1 - X_H2O_cat;

    %% 开路过电位（可逆电位）
    m = 8;  % 摩尔浓度
    wt = 35; % 溶液浓度
    P_H2O_KOH_sat = 10^(-0.01508*m - 0.0016788*m^2 + 2.25887e-5*m^3 + ...
                    (1 - 0.0012062*m + 5.6024e-4*m^2 - 7.8228e-6*m^3) * ...
                    (35.4462 - 3343.93/T - 10.9*log10(T) + 0.0041645*T));
    P_O2 = ((1/X_H2O_an) - 1) * P_H2O_KOH_sat;
    P_H2 = ((1/X_H2O_cat) - 1) * P_H2O_KOH_sat;
    V_std0 = 1.229;
    delta_S_0 = 0.9e-3;
    alpha_H20_KOH = 0.82;

    m_calc = wt*(183.1211 - 0.56845*(T+273.15) + 984.5679*exp(wt/115.96277))/5610.5;
    P_H2O_0 = exp(37.043 - 6275.7/T - 3.4159*log(T));
    P_H2O = exp(0.016214 - 0.13802*m_calc + 0.19330*sqrt(m_calc) + 1.0239*log(P_H2O_0));

    V_oc(i) = V_std0 + (T - T_ref)*delta_S_0 + (R*T/(2*F))*log((P_current - P_H2O)^(3/2)/(P_H2O/P_H2O_0));

    %% 活化过电位
    n = 2;
    j_an = j;
    j_cat = j;
    j_lim = 30;
    gamma_a = 1.25;
    gamma_c = 1.05;
    j_refa_0 = 1.34535e-5;
    j_refc_0 = 1.8456e-3;
    deltaG_an = 41500;
    deltaG_cat = 23450;
    theta = (-97.25 + 182*T/T_ref - 84*(T/T_ref)^2)*(j/j_lim)^0.3 * P_current/(P_current - P_H2O_KOH_sat);
    theta_array(i) = theta;  % 保存theta值
    j0_an = gamma_a*j_refa_0*exp(-deltaG_an/R*(1/T - 1/T_ref));
    j0_cat = gamma_c*j_refc_0*exp(-deltaG_cat/R*(1/T - 1/T_ref));
    alpha_a = 0.0675 + 0.00095*T;
    alpha_c = 0.1175 + 0.00095*T;
    b_an = R*T/(n*F*alpha_a);
    b_cat = R*T/(n*F*alpha_c);

    if j_an > 0 && j0_an > 0
        V_act_an = b_an*log(j_an/j0_an) + b_an*log(1/(1-theta));
    else
        V_act_an = 0;
    end

    if j_cat > 0 && j0_cat > 0
        V_act_cat = b_cat*log(j_cat/j0_cat) + b_cat*log(1/(1-theta));
    else
        V_act_cat = 0;
    end

    V_act(i) = V_act_an + V_act_cat;

    %% 浓差极化过电位
    t_an = 0.2;
    t_cat = 0.2;
    n_O2 = j/(4*F);
    n_H2 = j/(2*F);
    C_O2el_an = (P_current*X_O2)/(R*T) + (t_an*n_O2)/D_eff_an;
    C_H2el_cat = (P_current*X_H2)/(R*T) + (t_cat*n_H2)/D_eff_cat;

    % 标况下的
    P_ref = 1;
    t_ref = 85 + 273.15;
    D_eff_O20 = 0.00133*((1/M_O2+1/M_H2O)^0.5)*t_ref^1.5/(P_ref*(sigma_O2_H2O)^2*Tau_D_O2);
    D_eff_H20 = 0.00133*((1/M_H2+1/M_H2O)^0.5)*t_ref^1.5/(P_ref*(sigma_H2_H2O)^2*Tau_D_H2);
    D_eff_H2OK0 = 4/3*r*sqrt(8*R*t_ref/(pi*M_H2O));
    D_eff_an0 = 1/(epslion_an/tau_an*((1/D_eff_O20)+(1/D_eff_H2OK0)));
    D_eff_cat0 = 1/(epslion_cat/tau_cat*((1/D_eff_H20)+(1/D_eff_H2OK0)));
    X_H2O_an0 = ((epslion_an/tau_an)*exp(R*t_ref*l_anc*j/(2*F*P_ref*D_eff_an0)));
    X_H2O_cat0 = ((epslion_cat/tau_cat)*exp(R*t_ref*l_catc*j/(2*F*P_ref*D_eff_cat0)));
    X_O20 = 1 - X_H2O_an0;
    X_H20 = 1 - X_H2O_cat0;
    C_O20_an = (P_ref*X_O20)/(R*t_ref) + (t_an*n_O2)/D_eff_an0;
    C_H20_cat = (P_ref*X_H20)/(R*t_ref) + (t_cat*n_H2)/D_eff_cat0;

    V_con(i) = R*T/(4*F)*log(C_O2el_an/C_O20_an) + R*T/(2*F)*log(C_H2el_cat/C_H20_cat);

    %% 欧姆过电位
    A = 15000;
    A_e = A;
    A_m = A;
    rou_an = 6.4e-6;
    rou_cat = 6.4e-6;
    k_an = 0.00586;
    k_cat = 0.00586;
    d_an = 0.125;
    d_cat = 0.125;
    t_m = 0.05;
    R_an = (rou_an/(1-epslion_an)^1.5)*(t_an/A_e)*(1 + k_an*(T-T_ref));
    R_cat = (rou_cat/(1-epslion_cat)^1.5)*(t_cat/A_e)*(1 + k_cat*(T-T_ref));
    R_e = R_an + R_cat;

    sigma_KOH_free = -2.04*m - 0.0028*m^2 + 0.005332*m*T + 207.2*m/T + 0.001043*m^3 - 0.0000003*m^2*T^2;
    rou_el = 1/sigma_KOH_free;
    R_el_free = rou_el*(d_an/A_e + d_cat/A_e);
    R_el_bubble = R_el_free*(1/((1-2/3*theta)^1.5) - 1);
    R_el = R_el_free + R_el_bubble;

    omega_m = 0.85;
    epslion_m = 0.42;
    R_m = rou_el*tau_m^2*t_m/(omega_m*epslion_m*A_m);

    I = j*A_e;
    V_ohm(i) = I*(R_m + R_e + R_el);
end

% 绘制5张图
figure('Position', [100, 100, 1400, 900]);

% 图1：压力对可逆电位的影响
subplot(3, 2, 1);
plot(P, V_oc, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('压力 (bar)', 'FontSize', 12);
ylabel('可逆电位 V_{oc} (V)', 'FontSize', 12);
title('压力对可逆电位的影响', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 11);

% 图2：压力对活化过电位的影响
subplot(3, 2, 2);
plot(P, V_act, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('压力 (bar)', 'FontSize', 12);
ylabel('活化过电位 V_{act} (V)', 'FontSize', 12);
title('压力对活化过电位的影响', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 11);

% 图3：压力对浓差极化过电位的影响
subplot(3, 2, 3);
plot(P, V_con, 'g-d', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('压力 (bar)', 'FontSize', 12);
ylabel('浓差极化过电位 V_{con} (V)', 'FontSize', 12);
title('压力对浓差极化过电位的影响', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 11);

% 图4：压力对欧姆过电位的影响
subplot(3, 2, 4);
plot(P, V_ohm, 'm-^', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('压力 (bar)', 'FontSize', 12);
ylabel('欧姆过电位 V_{ohm} (V)', 'FontSize', 12);
title('压力对欧姆过电位的影响', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 11);

% 图5：压力对气泡覆盖率的影响
subplot(3, 2, 5);
plot(P, theta_array, 'k-p', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('压力 (bar)', 'FontSize', 12);
ylabel('气泡覆盖率 \theta', 'FontSize', 12);
title('压力对气泡覆盖率的影响', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 11);

% 添加总标题
sgtitle('电流密度固定为0.35 A/cm² 时压力对各过电位和气泡覆盖率的影响', 'FontSize', 16, 'FontWeight', 'bold');

% 输出一些关键数据
fprintf('电流密度: %.2f A/cm²\n', j);
fprintf('温度: %.0f °C\n', t);
fprintf('压力范围: %.1f - %.1f bar\n', min(P), max(P));
fprintf('\n在压力为 %.0f bar 时:\n', P(1));
fprintf('  可逆电位: %.4f V\n', V_oc(1));
fprintf('  活化过电位: %.4f V\n', V_act(1));
fprintf('  浓差极化过电位: %.4f V\n', V_con(1));
fprintf('  欧姆过电位: %.4f V\n', V_ohm(1));
fprintf('  气泡覆盖率: %.4f\n', theta_array(1));
fprintf('  总电压: %.4f V\n', V_oc(1)+V_act(1)+V_con(1)+V_ohm(1));
fprintf('\n在压力为 %.0f bar 时:\n', P(end));
fprintf('  可逆电位: %.4f V\n', V_oc(end));
fprintf('  活化过电位: %.4f V\n', V_act(end));
fprintf('  浓差极化过电位: %.4f V\n', V_con(end));
fprintf('  欧姆过电位: %.4f V\n', V_ohm(end));
fprintf('  气泡覆盖率: %.4f\n', theta_array(end));
fprintf('  总电压: %.4f V\n', V_oc(end)+V_act(end)+V_con(end)+V_ohm(end));
