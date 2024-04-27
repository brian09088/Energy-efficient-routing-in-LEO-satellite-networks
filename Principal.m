clear all;
clc;

%% 產生模擬時間表：日期、總模擬時間
%%%Segundos, Snapshot %%%
%模擬 = 模擬('2018-06-01 12:30:00',12200,100);
             %Simulate(日期, TempoTotalSimulacao, TempoIntervaloSimulacao);
             %TotalSimulation時間：1圈=6100秒； 2 圈 =
             %12200秒； 5圈=30000秒；
             %SimulationIntervalTime：快照間隔（場景），100
             % 秒；
tic

Simulacao = Simular('2020-07-20 13:00:00',90000,100); 

Simulacao.Lb = 0.3;     % 0.0045 路線上至少 100 J；
Simulacao.Lc = 0;       % 0.001
Simulacao.TipoAdaptacao = 0; % 不考慮日食

Simulacao.Pesos = [0.35 0.35 0.30];       % 權重
Simulacao.PesosLimiar = [0.15 0.70 0.15];   % 權重Threshold
Simulacao.constanteA = 0.8;

Simulacao.TotalFontes = 1000;       % 要模擬的總源 - 要模擬的終端數量
%Simulacao.CBRs = [1.5 1 0.5];      
Simulacao.CBRs = [1];               % 要模擬的 % CBR 速率（資料速率）：以 Mbps 為單位
%Simulacao.Metricas = {'TP';'LASER';'PROPOSTA'};   
Simulacao.Metricas = {'PROPOSTA'};   % 指標
Simulacao.PotenciaCP = 117000;       % 電池電量 => 117 KJ
Simulacao.PotenciaTX = 7;            % 傳輸功耗 => 瓦特 = J/s
Simulacao.PotenciaRX = 3;            % 接收功耗 => 瓦 = J/s
Simulacao.PotenciaON = 4;            % 正常運作時的功耗 => 瓦特 = J/s
Simulacao.PotenciaCG = 20;           % 能量捕獲功率：20 瓦 = 20 J/s
Simulacao.ISLCap = 10;               % 連結容量（Mbps)
%TotalSimulacoes = 30;               % 信賴區間的總模擬百分比
TotalSimulacoes = 1;                 % 模擬總數 = 1


%% II - 儲存模擬資料 - 場景
save(strcat('Resultados/',num2str(Simulacao.TotalFontes),'_fontes','/Simulacao.mat'),'Simulacao','-v7.3');

%% III - 執行路由
for i=1:TotalSimulacoes
  Resultados = Simulacao.Roteamento(i); % Roteamento
  save(strcat('Resultados/',num2str(Simulacao.TotalFontes),'_fontes','/Resultados_',num2str(i),'.mat'),'Resultados','-v7.3');
end

%% IV - 合併數據以產生圖形
DadosConsolidados = Consolidar(strcat('Resultados/',num2str(Simulacao.TotalFontes),'_fontes/'),TotalSimulacoes);

time = toc ;
save(strcat('Resultados/',num2str(Simulacao.TotalFontes),'_fontes','/DadosConsolidados','.mat'),'DadosConsolidados','time','-v7.3');




