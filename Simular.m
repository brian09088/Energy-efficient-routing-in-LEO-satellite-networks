classdef Simular

     properties
        Cenario;            % 模擬場景數據
        TotalFontes;        % 模擬終端數量
        CBRs;               % 
        Metricas;           % 模擬指標
        PotenciaCP;         % 額定運轉功耗
        PotenciaTX;         % 傳輸功耗
        PotenciaRX;         % 接收期間的功耗百分比
        PotenciaON;         % 正常運作時的功耗百分比
        PotenciaCG;         % 能量捕獲功率
        ISLCap;             % 鏈路容量 
        Resultados;         % 結果變數
        Lb;
        Lc;
        TipoAdaptacao;
        Pesos;
        PesosLimiar;
        constanteA;
    end
    
    methods
        
        % 構造函數
        function obj = Simular(Data,TempoTotalSimulacao,TempoIntervaloSimulacao)
          % 產生模擬場景
          obj.Cenario = Cenario(Data,TempoTotalSimulacao,TempoIntervaloSimulacao);
        end

        function obj = set.Cenario(obj,val)
          obj.Cenario = val;
        end
        
        function obj = set.TotalFontes(obj,val)
          obj.TotalFontes = val;
        end     
        
        function obj = set.CBRs(obj,val)
          obj.CBRs = val;
        end  
        
        function obj = set.Metricas(obj,val)
          obj.Metricas = val;
        end  
        
        function obj = set.PotenciaCP(obj,val)
          obj.PotenciaCP = val;
        end   
        
        function obj = set.PotenciaON(obj,val)
          obj.PotenciaON = val;
        end 
        
        function obj = set.PotenciaTX(obj,val)
          obj.PotenciaTX = val;
        end        
        
        function obj = set.PotenciaRX(obj,val)
          obj.PotenciaRX = val;
        end
        
        function obj = set.PotenciaCG(obj,val)
          obj.PotenciaCG = val;
        end  
        
        function obj = set.ISLCap(obj,val)
          obj.ISLCap = val;
        end        
        
        function obj = set.Resultados(obj,val)
          obj.Resultados = val;
        end
        
        function obj = set.Lc(obj,val)
          obj.Lc = val;
        end
        
        function obj = set.Lb(obj,val)
          obj.Lb = val;
        end
        
        function Resultados = Roteamento(obj,loopId,Interacoes)
                  
            if nargin == 2
                Interacoes = size(obj.Cenario.Propagacao,2);
            end   
            Resultados = {};
            Edges = {};
            % ---------------- 及時模擬 -------------------------------------
            LoopId = 1;
            
            %fprintf(' Teste -- Simula誽o %d (Roteamento %d/%d)\n',loopId,LoopId,Interacoes);  
            
            while LoopId <= Interacoes 
                fprintf('Simula誽o %d (Roteamento %d/%d)\n',loopId,LoopId,Interacoes);       
                % ---------------------- 生成來源/目的地 -----------------   
                % 存取：結果（終端 ID）.來源
                % 對 - 選擇原產地大陸
                %CkSrc = randsrc(obj.TotalFontes,1,[1:max(obj.Cenario.Zonas.Continentes(:)); ...
                %obj.Cenario.Zonas.ContinentesHosts./sum(obj.Cenario.Zonas.ContinentesHosts(:))]);
                CkSrc = randsrc(obj.TotalFontes,1,1:max(obj.Cenario.Zonas.Continentes(:)));
                for fonte_id=1: obj.TotalFontes
                    % 原產地
                    Pares(fonte_id).ContinenteOrigem = CkSrc(fonte_id); 
                    % 有流量需求的來源區域的 id
                    IdCkSrc = find(obj.Cenario.Zonas.Continentes...
                        .*(obj.Cenario.Zonas.DensidadeHosts>0) == CkSrc(fonte_id));
                    % 主區域選擇
                    %ZkSrc = randsrc(1,1,[IdCkSrc'; obj.Cenario.Zonas.DensidadeHosts(IdCkSrc)'...
                    %   ./sum(obj.Cenario.Zonas.DensidadeHosts(IdCkSrc))]);  
                    ZkSrc = randsrc(1,1,IdCkSrc');  
                    
                    Pares(fonte_id).ZonaOrigem = ZkSrc; 

                    % 起源大陸的衛星
                    SatSrc = obj.Cenario.Propagacao(LoopId).Zonas.Satelites(ZkSrc)';
                    Pares(fonte_id).SateliteOrigem = SatSrc;
                    % 目的地大陸
                    CkDst = randsrc(1,1,[1:max(obj.Cenario.Zonas.Continentes(:));...
                        obj.Cenario.Propagacao(LoopId).Zonas.Fluxo(CkSrc(fonte_id),:)]);                  
                    Pares(fonte_id).ContinenteDestino = CkDst;      
                    % 目的地區域/衛星選擇 - 避免選擇甚至
                    while(1)
                        IdCkDst = find(obj.Cenario.Zonas.Continentes...
                            .*(obj.Cenario.Zonas.DensidadeHosts>0) == CkDst);
                         %ZkDst = randsrc(1,1,[IdCkDst'; obj.Cenario.Zonas.DensidadeHosts(IdCkDst)'...
                         %   ./sum(obj.Cenario.Zonas.DensidadeHosts(IdCkDst))]);  
                        ZkDst = randsrc(1,1,IdCkDst'); 
                        Pares(fonte_id).ZonaDestino = ZkSrc;
                        SatDst = obj.Cenario.Propagacao(LoopId).Zonas.Satelites(ZkDst)';
                        Pares(fonte_id).SateliteDestino = SatDst;
                        % 避免為出發地和目的地選擇相同的衛星
                        if SatSrc ~= SatDst
                            break;
                        end       
                    end
                end
                Resultados(LoopId).Fontes = Pares; 
                % ----------------------- 活動連結 ------------------------------
                % 每個連結的傳播時間矩陣 - 毫秒
                MatDelay = ((obj.Cenario.Propagacao(LoopId).Enlaces...
                ./299792458).*1000000);   %%  velocidade da luz no v塶uo = 299.792.458 m/s
                %% Tempo = Enlace (km)/ V (m/s) = Enlace*10^3 m / 299.792.458 m/ 10^3 milissegundos 
                %% Tempo = Enlace (km)/ V (m/s) = Enlace*10^3 * 10^3 milissegundos / 299.792.458
                %% Tempo = Enlace (km)/ V (m/s) = Enlace*1000000 milissegundos / 299.792.458
               
                % 網路鄰接矩陣
                MatRede  = (MatDelay>0); % 二進位矩陣 ={ 1, for MatDelay>0; 否則為 0} 
          
                % --------------------- 位元率 --------------------------------
                for cbr_id = 1: size(obj.CBRs,2)
                    cbr = obj.CBRs(cbr_id);
                    % --------------------- 指標 ----------------------------------
                    for metrica_id=1:size(obj.Metricas,1)
                        Metrica = num2str(cell2mat(obj.Metricas(metrica_id)));   
                        % --------------- 能量矩陣 ---------------------------
                         ConsumoEnergia = zeros(obj.Cenario.TotalSatelites,1);
                         CicloVida = zeros(obj.Cenario.TotalSatelites,1);
                         GD = zeros(obj.Cenario.TotalSatelites,1);
                         DOD = zeros(obj.Cenario.TotalSatelites,1);
                         DOD_ant = zeros(obj.Cenario.TotalSatelites,1);
                         
                        if LoopId == 1 % ######### 初始條件  
                            EnergiaInicio = ones(obj.Cenario.TotalSatelites,1)*obj.PotenciaCP;
                            DOD = zeros(obj.Cenario.TotalSatelites,1);
                            %ConsumoEnergia = zeros(obj.Cenario.TotalSatelites,1);
                        else % ######### 非初始條件
                            
                            % 時間間隔開始時的能量百分比等於上一個間隔結束時的剩餘能量
                           
                            EnergiaInicio = Resultados(LoopId-1).(Metrica).EnergiaFinal(:,cbr_id);
                            DOD_ant =  Resultados(LoopId-1).(Metrica).DOD(:,cbr_id); %DOD da Iteracao anterior
                            %ConsumoEnergia = Resultados(LoopId-1).(Metrica).ConsumoEnergia(:,cbr_id);
                        end 
                        Resultados(LoopId).(Metrica).EnergiaInicio(:,cbr_id) = EnergiaInicio;
                        
                        Resultados(LoopId).(Metrica).ConsumoEnergia(:,cbr_id) = ConsumoEnergia;
                        
                        % 能量容量矩陣已扣除能量 
                        % 衛星標稱運行百分比
                        EnergiaON = ones(obj.Cenario.TotalSatelites,1)...
                            .*obj.PotenciaON.*obj.Cenario.TempoIntervaloSimulacao;
                        %%% 衛星當前剩餘能量 %%%%%
                        Energia =  EnergiaInicio  -  EnergiaON;
                        % 更新能源消耗
                        ConsumoEnergia = EnergiaON; 
                        
                        % --- 容量、需求和需求滿足矩陣 ------
                        % 鏈路容量
                        Capacidade = MatRede.*obj.ISLCap; % Mbps
                        % 每個鏈路傳輸的資料百分比
                        MatTrafego = MatRede.*0;
                        % 出發地/目的地需求百分比
                        MatDemanda = MatRede.*0;       
                        % 按需服務出發地/目的地
                        MatDemandaAtendida = MatRede.*0;
                       
                        % ----  掃描終端對（來源和目標） --------  
                        for fonte_id=1: obj.TotalFontes   
                            Origem  = Resultados(LoopId).Fontes(fonte_id).SateliteOrigem;
                            Destino = Resultados(LoopId).Fontes(fonte_id).SateliteDestino;
                            Demanda = cbr; % CBR, Mbps
                            % --------------------- 指標的應用 -------------
                            switch Metrica    
                                case 'TP'    
                                    G = digraph(MatDelay);
                                    % 傳播時間                    
                                    Tmin = min(G.Edges.Weight); % 烏托邦  
                                    Tmax = max(G.Edges.Weight); % 最低點  
                                    % 歸一化傳播時間
                                    TPNorm = (G.Edges.Weight - Tmin)/(Tmax - Tmin);
                                    G.Edges.Weight = TPNorm;           
                                case 'LASER'        
                                    G = digraph(MatDelay);
                                    % 傳播時間                   
                                    Tmin = min(G.Edges.Weight); % 烏托邦 
                                    Tmax = max(G.Edges.Weight); % 最低點  
                                    % 歸一化傳播時間
                                    TPNorm = (G.Edges.Weight - Tmin)/(Tmax - Tmin);
                                    [Si,Sj] = findedge(G);
                                    Ei = [obj.Cenario.Propagacao(LoopId).Satelites(Si).EclipseStatus]';
                                    Ej = [obj.Cenario.Propagacao(LoopId).Satelites(Sj).EclipseStatus]'; 
                                    % 電池電量標準化
                                    Bi = Energia(Si)./obj.PotenciaCP;
                                    Bj = Energia(Sj)./obj.PotenciaCP;
                                    Dij = Ei./Bi + Ej./Bj;
                                    Dij(isnan(Dij)) = 0; 
                                    Dmin = min(Dij); % 烏托邦  
                                    Dmax = max(Dij); % 最低點 
                                    DijNorm = (Dij-Dmin)/(Dmax-Dmin);
                                    G.Edges.Weight = 0.5.*TPNorm+0.5.*DijNorm;   
                                    
                                case 'PROPOSTA' 
                                    G = digraph(MatDelay);
                                    % 傳播時間                    
                                    Tmin = min(G.Edges.Weight); % utopia   
                                    Tmax = max(G.Edges.Weight); % nadir     
                                    % 歸一化傳播時間
                                    TPNorm = (G.Edges.Weight - Tmin)/(Tmax - Tmin);
                                    [Si,Sj] = findedge(G);
                                    Ei = [obj.Cenario.Propagacao(LoopId).Satelites(Si).EclipseStatus]';
                                    Ej = [obj.Cenario.Propagacao(LoopId).Satelites(Sj).EclipseStatus]'; 
                                    Bi = Energia(Si)./obj.PotenciaCP;
                                    Bj = Energia(Sj)./obj.PotenciaCP;
                                    % 新模型，僅考慮
                                    Peso = ones(size(G.Edges.Weight,1),3);
                                    Peso(:,1) = obj.Pesos(1);
                                    Peso(:,2) = obj.Pesos(2);
                                    Peso(:,3) = obj.Pesos(3);
                                    
                                    xij = zeros(size(G.Edges.Weight,1),1);
                                    pos = zeros(size(G.Edges.Weight,1),1);
                                    if ( obj.Lb > 0)
                                        BijBmax = [Bi(:,:) Bj(:,:)];
                                        xij = any (BijBmax < obj.Lb, 2);
                                        xij = double(xij(:));
                                        pos = find( xij == 1); % 要移除的衛星的位置百分比
                                    end
                                      
                                    if (obj.TipoAdaptacao == 2)
                                       Dij =  -(Bi + Bj); 
                                    else
                                       Dij = Ei./Bi + Ej./Bj; 
                                    end
                                    
                                    Cij = zeros(size(G.Edges.Weight,1),1);
                                    CijResidual = zeros(size(G.Edges.Weight,1),1);
                                    for s_ij = 1:size(Cij,1)
                                        CijResidual(s_ij,1) = Capacidade(Si(s_ij),Sj(s_ij)) / obj.ISLCap ;
                                        Cij(s_ij,1) = 1 - CijResidual(s_ij,1); 
                                        %Cij(s_ij,1) = 1-(Capacidade(Si(s_ij),Sj(s_ij))/obj.ISLCap);
                                    end 
                                    
                                    % 從連結中刪除 - 衛星和
                                    % 電池連接未滿足
                                    % 電池容量和連結容量
                                    
                                    Zij = zeros(size(G.Edges.Weight,1),1);
                                    posZij = zeros(size(G.Edges.Weight,1),1);
                                    if ( obj.Lc > 0)
                                      Zij = any (CijResidual(:,:) < obj.Lc, 2); 
                                      Zij = double(Zij(:));
                                      posZij = find(  Zij == 1);
%                                       pos = union(pos, posZij)
                                    end
                                    
                                    if (sum (pos) > 0 && sum (posZij) > 0 )
                                         pos = union(pos, posZij);
                                    elseif ( sum (pos) == 0 && sum (posZij) > 0)
                                        pos = posZij;
                                    end
                                    
                                    if (sum (pos) > 0 )
                                         Cij(pos, :) = NaN; % 分配 NaN 以從標準化中刪除
                                         Dij(pos, :) = NaN; 
                                         Peso(pos,1) = obj.PesosLimiar(1);
                                         Peso(pos,2) = obj.PesosLimiar(2);
                                         Peso(pos,3) = obj.PesosLimiar(3);
                                    else
                                         Dij(isnan(Dij)) = 0;
                                    end
                                    
                                    Dmin = min(min(Dij,[],'omitnan')); % utopia
                                    Dmax = max(max(Dij,[],'omitnan')); % nadir 
                                    DijNorm = (Dij-Dmin)/(Dmax-Dmin);

                                    Cmin = 0;
                                    Cmax = 1;
                                    CijNorm = (Cij - Cmin)./(Cmax-Cmin);
                                    
                                    if (sum (pos) > 0 )
                                      CijNorm(pos, :) = inf; % 刪除沒有能力的連結
                                      DijNorm(pos, :) = inf; % 無需電池即可用衛星去除邊緣
                                      TPNorm(pos, :) = inf;  % 無需電池即可用衛星去除邊緣  
                                      %Se for NaN, Not a Number, atribui
                                      %infinito..
%                                       CijNorm(isnan(CijNorm)) = inf; 
%                                       DijNorm(isnan(DijNorm)) = inf ; 
%                                       TPNorm(isnan(TPNorm)) = inf; 
                                    end
                                      % 如果不是數字則禁用路線
                                      
                                     CijNorm(isnan(CijNorm)) = 0; 
                                     DijNorm(isnan(DijNorm)) = 0 ; 
                                     TPNorm(isnan(TPNorm)) = 0;  
                                    
                                    %%%% 建議指標：w1、w2、w3 %%%%%%%
                                    %G.Edges.Weight = 0.35.*TPNorm+0.35.*DijNorm+0.3.*CijNorm;
                                    TDC = [TPNorm DijNorm CijNorm];
                                    G.Edges.Weight = diag( Peso(:,:)  * TDC(:,:)');
                                    %Resultados(LoopId).(Metrica).Fontes.demanda(fonte_id,cbr_id) = Demanda;
                             end
                            
                             DemandaAtendida = Demanda; % Mbps
                             % 時間間隔內出發地/目的地需求的百分比矩陣
                              MatDemanda(Origem,Destino) = MatDemanda(Origem,Destino) ...
                                + Demanda*obj.Cenario.TempoIntervaloSimulacao;  
                              
                              [rota, custo] = shortestpath(G,Origem,Destino); %%%%% 該費用未使用 %%%%%%
                                % 路徑上衛星的能量
                                EnergiaRota =  Energia(rota);  
                                % 最小路徑能量
                                
                                MinEnergiaRota = min(EnergiaRota);
                                %%%% 能量 Tx、Rx 是衛星 k-l 之間的每個連結
                                % 傳輸所需能量
                                EnergiaTx = ((DemandaAtendida/obj.ISLCap)*obj.PotenciaTX)*obj.Cenario.TempoIntervaloSimulacao;
                                % 接收所需能量百分比
                                EnergiaRx = ((DemandaAtendida/obj.ISLCap)*obj.PotenciaRX)*obj.Cenario.TempoIntervaloSimulacao;
                                % 傳輸和接收所需的總能量
                                EnergiaTxRx = EnergiaTx + EnergiaRx;

                                % 檢查是否有足夠的功率來傳輸
                                if MinEnergiaRota < EnergiaTxRx
                                   DemandaAtendida = 0;
                                end  
                                
                            CapRota = [];
                            % 路由跳數總計百分比
                            saltos = size(rota,2)-1;
                            % 找到一條有容量的路線
                            for st=1:saltos        
                                % 跳數連結容量
                                CapEnlace = Capacidade(rota(st),rota(st+1)); % Mbps
                                CapRota = cat(1,CapRota,CapEnlace);
                                if CapEnlace < DemandaAtendida % 如果仍有連結容量
                                   % 限制容量需求
                                   DemandaAtendida = CapEnlace;
                                end
                            end 
                            
                            % 連結傳播時間（以毫秒為單位）
                            Delay = 0;
                            % 如果有資料要在該路線上行駛，則執行路由
                            if DemandaAtendida > 0
                                
                                for st=1:saltos        
                                    % 更新能量 - 傳輸
                                    Energia(rota(st)) = Energia(rota(st)) - EnergiaTx;
                                    % 更新能量 - 接收
                                    Energia(rota(st+1)) = Energia(rota(st+1)) - EnergiaRx; 
                                    % 更新鏈路容量
                                    Capacidade(rota(st),rota(st+1))= Capacidade(rota(st),rota(st+1))...
                                        - DemandaAtendida;
                                    % 連結流量矩陣
                                    MatTrafego(rota(st),rota(st+1))...
                                        = MatTrafego(rota(st),rota(st+1))+DemandaAtendida...
                                        *obj.Cenario.TempoIntervaloSimulacao; % (Mb) = 需求 (MB/s)* 模擬間隔時間 (s)；
                                    % 路由延遲（傳播時間）
                                    Delay = Delay+MatDelay(rota(st),rota(st+1));   
                                    
                                    %%%%%#### 添加能量傳輸
                                    %%%%% 接待花費
                                    ConsumoEnergia(rota(st)) = ConsumoEnergia(rota(st)) + EnergiaTx; % 能量消耗轉移
                                    ConsumoEnergia(rota(st+1)) =  ConsumoEnergia(rota(st+1)) + EnergiaRx; % 接收功耗
                                    %%%%%%%
                                end
                            else
                                % 沒有傳送資料    
                                saltos = 0;          
                            end
                             % 時間間隔內滿足的需求矩陣百分比
                             MatDemandaAtendida(Origem,Destino) = MatDemandaAtendida(Origem,Destino)...
                                + DemandaAtendida*obj.Cenario.TempoIntervaloSimulacao; 
                             % 來源摘要
                             Resultados(LoopId).(Metrica).Fontes.rota(fonte_id,cbr_id) = {rota};
                             Resultados(LoopId).(Metrica).Fontes.demanda(fonte_id,cbr_id) = Demanda;
                             Resultados(LoopId).(Metrica).Fontes.capacidadesrota(fonte_id,cbr_id) = {CapRota};
                             Resultados(LoopId).(Metrica).Fontes.enlacessaturados(fonte_id,cbr_id) = size(find(CapRota==0),1);
                             Resultados(LoopId).(Metrica).Fontes.demandaatendida(fonte_id,cbr_id)...
                                 = DemandaAtendida;
                             Resultados(LoopId).(Metrica).Fontes.saltos(fonte_id,cbr_id) = saltos;
                             %if Delay ==0
                                % Delay = inf;
                             %end
                             Resultados(LoopId).(Metrica).Fontes.tempopropagacao(fonte_id,cbr_id) = Delay;   
                             Resultados(LoopId).(Metrica).Fontes.eclipserota(fonte_id,cbr_id)...
                                 = size(find([obj.Cenario.Propagacao(LoopId).Satelites(rota).EclipseStatus] ==1),2); 
                            
                            %%%%% 每個衛星連結的能耗 %%%% 
%                           Consumo =  Resultados(LoopId).(Metrica).ConsumoEnergia(Origem,cbr_id);
%                           Consumo =  EnergiaON(Origem) + TotalTxRx; %+ Consumo;
%                           Resultados(LoopId).(Metrica).ConsumoEnergia(Origem,cbr_id) = Consumo;                             
                        end
                      
                        % 需求
                        %Resultados(LoopId).(Metrica).Demanda(:,cbr_id) = {MatDemanda};
                        %Resultados(LoopId).(Metrica).DemandaAtendida(:,cbr_id) = {MatDemandaAtendida};     
                        % 間隔結束時的連結容量
                        Resultados(LoopId).(Metrica).Trafego(:,cbr_id) = {MatTrafego};  % MBytes           
                        % 間隔結束時的連結容量
                        Resultados(LoopId).(Metrica).Capacidade(:,cbr_id) = {Capacidade};
                        % 間隔能量捕獲
                        CaptacaoEnergetica = obj.PotenciaCG*max(0,obj.Cenario.TempoIntervaloSimulacao-[obj.Cenario.Propagacao(LoopId).Satelites.EclipseTempoIntervalo]');           
                        %%%%% 能量捕獲 ...
                        %%% 添加以計算能量
                        %CaptacaoEnergetica = Resultados(LoopId).(Metrica).CaptacaoEnergetica(:,cbr_id) + CaptacaoEnergetica
                        %%%% 原始計算 %%%
                        Resultados(LoopId).(Metrica).CaptacaoEnergetica(:,cbr_id) = CaptacaoEnergetica;
                        % 更新能量矩陣
                        
                        % 剩餘能量 = 能量 
                          A = obj.constanteA;
                        % Energia(:) =  Energia(:) - Energia(:)*0.25; % teste diminui 25%
                          Dt2 =  (obj.PotenciaCP - Energia(:) ) ./ obj.PotenciaCP;
                          DOD(:,:) =  Dt2(:,:);
                          Dt1 = DOD_ant(:,:);
                         
                         for k=1:obj.Cenario.TotalSatelites
                            if( Dt2(k) > Dt1(k))
                                CicloVida(k,:) = Dt2(k) * 10^(A*(Dt2(k)-1)) - Dt1(k) * 10^(A*(Dt1(k)-1));
                            else
                                CicloVida(k,:) = 0;
                            end
                            GD(k)= Dt2(k) * 10^(A*(Dt2(k)-1));
                         end
                        
                        Energia = min(obj.PotenciaCP, Energia + CaptacaoEnergetica);
                        Energia = max(0, Energia); % 添加以避免能量為負值
                       
                        Resultados(LoopId).(Metrica).EnergiaFinal(:,cbr_id)  = Energia;
                        Resultados(LoopId).(Metrica).ConsumoEnergia(:,cbr_id) = ConsumoEnergia;
                        
                        Resultados(LoopId).(Metrica).DOD(:,cbr_id) = DOD;
                        Resultados(LoopId).(Metrica).CicloVida(:,cbr_id) = CicloVida;
                        Resultados(LoopId).(Metrica).GD(:,cbr_id)  = GD;
                       
                    end 
                end   
                LoopId = LoopId+1;
            end                       
        end   
    end   
end