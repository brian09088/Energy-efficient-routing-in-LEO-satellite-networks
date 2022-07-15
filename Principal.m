clear all;
clc;

%% I - GERA CN�RIO DA SIMUL�AO: Data, Tempo Total da Simula��o em
%%%Segundos, Snapshot %%%
%Simulacao = Simular('2018-06-01 12:30:00',12200,100);
            %Simular(Data, TempoTotalSimulacao, TempoIntervaloSimulacao);
            % TempoTotalSimulacao:  1 volta = 6100 segundos; 2 voltas =
            % 12200 segundos; 5 voltas = 30000 segundos;
            % TempoIntervaloSimulacao: intervaldo do snapshot (cen�rios), 100
            % segundos;
tic

Simulacao = Simular('2020-07-20 13:00:00',90000,100); 

Simulacao.Lb = 0.3;  % 0.0045 Minimo de 100 J na Rota;
Simulacao.Lc = 0;  % 0.001
Simulacao.TipoAdaptacao = 0; % 2 - sem considerar o ecliplse

Simulacao.Pesos = [0.35 0.35 0.30];
Simulacao.PesosLimiar = [0.15 0.70 0.15];
Simulacao.constanteA = 0.8;

Simulacao.TotalFontes = 1000;       % Total de fontes para simular - n�meros de terminais para simula��o
%Simulacao.CBRs = [1.5 1 0.5];      % Taxas CBR (Taxa de dados)para Simular: em Mbps - Megabits por segundo
Simulacao.CBRs = [1];         % Taxas CBR (Taxa de dados)para Simular: em Mbps - Megabits por segundo
%Simulacao.Metricas = {'TP';'LASER';'PROPOSTA'};   % M�tricas
Simulacao.Metricas = {'PROPOSTA'};   % M�tricas
Simulacao.PotenciaCP = 117000;       % Capacidade de pot�ncia da bateria => 117 KJ
Simulacao.PotenciaTX = 7;            % Consumo de pot�ncia na transmiss�o  => watt = J/s
Simulacao.PotenciaRX = 3;            % Consumo de pot�ncia na recep��o  => watt = J/s
Simulacao.PotenciaON = 4;            % Consumo de pot�ncia na opera��o normal => watt = J/s
Simulacao.PotenciaCG = 20;           % Pot�ncia de capta��o de energia: 20 watt = 20 J/s
Simulacao.ISLCap = 10;              % Capacidade dos links (Mbps - Megabits por segundo)
%TotalSimulacoes = 30;               % Total de simula��es para intervalo de confian�a
TotalSimulacoes = 1;


%% II - SALVA DADOS DA SIMULA��O - CEN�RIO
save(strcat('Resultados/',num2str(Simulacao.TotalFontes),'_fontes','/Simulacao.mat'),'Simulacao','-v7.3');

%% III - REALIZA O ROTEAMENTO
for i=1:TotalSimulacoes
  Resultados = Simulacao.Roteamento(i); % Roteamento
  save(strcat('Resultados/',num2str(Simulacao.TotalFontes),'_fontes','/Resultados_',num2str(i),'.mat'),'Resultados','-v7.3');
end

%% IV - CONSOLIDA DADOS PARA GERAR GR�FICOS 
DadosConsolidados = Consolidar(strcat('Resultados/',num2str(Simulacao.TotalFontes),'_fontes/'),TotalSimulacoes);

time = toc ;
save(strcat('Resultados/',num2str(Simulacao.TotalFontes),'_fontes','/DadosConsolidados','.mat'),'DadosConsolidados','time','-v7.3');




