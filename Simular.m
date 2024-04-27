classdef Simular

     properties
        Cenario;            % ���������ƾ�
        TotalFontes;        % �����׺ݼƶq
        CBRs;               % 
        Metricas;           % ��������
        PotenciaCP;         % �B�w�B��\��
        PotenciaTX;         % �ǿ�\��
        PotenciaRX;         % �����������\�Ӧʤ���
        PotenciaON;         % ���`�B�@�ɪ��\�Ӧʤ���
        PotenciaCG;         % ��q����\�v
        ISLCap;             % ����e�q 
        Resultados;         % ���G�ܼ�
        Lb;
        Lc;
        TipoAdaptacao;
        Pesos;
        PesosLimiar;
        constanteA;
    end
    
    methods
        
        % �c�y���
        function obj = Simular(Data,TempoTotalSimulacao,TempoIntervaloSimulacao)
          % ���ͼ�������
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
            % ---------------- �ήɼ��� -------------------------------------
            LoopId = 1;
            
            %fprintf(' Teste -- Simula��o %d (Roteamento %d/%d)\n',loopId,LoopId,Interacoes);  
            
            while LoopId <= Interacoes 
                fprintf('Simula��o %d (Roteamento %d/%d)\n',loopId,LoopId,Interacoes);       
                % ---------------------- �ͦ��ӷ�/�ت��a -----------------   
                % �s���G���G�]�׺� ID�^.�ӷ�
                % �� - ��ܭ첣�a�j��
                %CkSrc = randsrc(obj.TotalFontes,1,[1:max(obj.Cenario.Zonas.Continentes(:)); ...
                %obj.Cenario.Zonas.ContinentesHosts./sum(obj.Cenario.Zonas.ContinentesHosts(:))]);
                CkSrc = randsrc(obj.TotalFontes,1,1:max(obj.Cenario.Zonas.Continentes(:)));
                for fonte_id=1: obj.TotalFontes
                    % �첣�a
                    Pares(fonte_id).ContinenteOrigem = CkSrc(fonte_id); 
                    % ���y�q�ݨD���ӷ��ϰ쪺 id
                    IdCkSrc = find(obj.Cenario.Zonas.Continentes...
                        .*(obj.Cenario.Zonas.DensidadeHosts>0) == CkSrc(fonte_id));
                    % �D�ϰ���
                    %ZkSrc = randsrc(1,1,[IdCkSrc'; obj.Cenario.Zonas.DensidadeHosts(IdCkSrc)'...
                    %   ./sum(obj.Cenario.Zonas.DensidadeHosts(IdCkSrc))]);  
                    ZkSrc = randsrc(1,1,IdCkSrc');  
                    
                    Pares(fonte_id).ZonaOrigem = ZkSrc; 

                    % �_���j�����ìP
                    SatSrc = obj.Cenario.Propagacao(LoopId).Zonas.Satelites(ZkSrc)';
                    Pares(fonte_id).SateliteOrigem = SatSrc;
                    % �ت��a�j��
                    CkDst = randsrc(1,1,[1:max(obj.Cenario.Zonas.Continentes(:));...
                        obj.Cenario.Propagacao(LoopId).Zonas.Fluxo(CkSrc(fonte_id),:)]);                  
                    Pares(fonte_id).ContinenteDestino = CkDst;      
                    % �ت��a�ϰ�/�ìP��� - �קK��ܬƦ�
                    while(1)
                        IdCkDst = find(obj.Cenario.Zonas.Continentes...
                            .*(obj.Cenario.Zonas.DensidadeHosts>0) == CkDst);
                         %ZkDst = randsrc(1,1,[IdCkDst'; obj.Cenario.Zonas.DensidadeHosts(IdCkDst)'...
                         %   ./sum(obj.Cenario.Zonas.DensidadeHosts(IdCkDst))]);  
                        ZkDst = randsrc(1,1,IdCkDst'); 
                        Pares(fonte_id).ZonaDestino = ZkSrc;
                        SatDst = obj.Cenario.Propagacao(LoopId).Zonas.Satelites(ZkDst)';
                        Pares(fonte_id).SateliteDestino = SatDst;
                        % �קK���X�o�a�M�ت��a��ܬۦP���ìP
                        if SatSrc ~= SatDst
                            break;
                        end       
                    end
                end
                Resultados(LoopId).Fontes = Pares; 
                % ----------------------- ���ʳs�� ------------------------------
                % �C�ӳs�����Ǽ��ɶ��x�} - �@��
                MatDelay = ((obj.Cenario.Propagacao(LoopId).Enlaces...
                ./299792458).*1000000);   %%  velocidade da luz no v�cuo = 299.792.458 m/s
                %% Tempo = Enlace (km)/ V (m/s) = Enlace*10^3 m / 299.792.458 m/ 10^3 milissegundos 
                %% Tempo = Enlace (km)/ V (m/s) = Enlace*10^3 * 10^3 milissegundos / 299.792.458
                %% Tempo = Enlace (km)/ V (m/s) = Enlace*1000000 milissegundos / 299.792.458
               
                % �����F���x�}
                MatRede  = (MatDelay>0); % �G�i��x�} ={ 1, for MatDelay>0; �_�h�� 0} 
          
                % --------------------- �줸�v --------------------------------
                for cbr_id = 1: size(obj.CBRs,2)
                    cbr = obj.CBRs(cbr_id);
                    % --------------------- ���� ----------------------------------
                    for metrica_id=1:size(obj.Metricas,1)
                        Metrica = num2str(cell2mat(obj.Metricas(metrica_id)));   
                        % --------------- ��q�x�} ---------------------------
                         ConsumoEnergia = zeros(obj.Cenario.TotalSatelites,1);
                         CicloVida = zeros(obj.Cenario.TotalSatelites,1);
                         GD = zeros(obj.Cenario.TotalSatelites,1);
                         DOD = zeros(obj.Cenario.TotalSatelites,1);
                         DOD_ant = zeros(obj.Cenario.TotalSatelites,1);
                         
                        if LoopId == 1 % ######### ��l����  
                            EnergiaInicio = ones(obj.Cenario.TotalSatelites,1)*obj.PotenciaCP;
                            DOD = zeros(obj.Cenario.TotalSatelites,1);
                            %ConsumoEnergia = zeros(obj.Cenario.TotalSatelites,1);
                        else % ######### �D��l����
                            
                            % �ɶ����j�}�l�ɪ���q�ʤ��񵥩�W�@�Ӷ��j�����ɪ��Ѿl��q
                           
                            EnergiaInicio = Resultados(LoopId-1).(Metrica).EnergiaFinal(:,cbr_id);
                            DOD_ant =  Resultados(LoopId-1).(Metrica).DOD(:,cbr_id); %DOD da Iteracao anterior
                            %ConsumoEnergia = Resultados(LoopId-1).(Metrica).ConsumoEnergia(:,cbr_id);
                        end 
                        Resultados(LoopId).(Metrica).EnergiaInicio(:,cbr_id) = EnergiaInicio;
                        
                        Resultados(LoopId).(Metrica).ConsumoEnergia(:,cbr_id) = ConsumoEnergia;
                        
                        % ��q�e�q�x�}�w������q 
                        % �ìP�кٹB��ʤ���
                        EnergiaON = ones(obj.Cenario.TotalSatelites,1)...
                            .*obj.PotenciaON.*obj.Cenario.TempoIntervaloSimulacao;
                        %%% �ìP��e�Ѿl��q %%%%%
                        Energia =  EnergiaInicio  -  EnergiaON;
                        % ��s�෽����
                        ConsumoEnergia = EnergiaON; 
                        
                        % --- �e�q�B�ݨD�M�ݨD�����x�} ------
                        % ����e�q
                        Capacidade = MatRede.*obj.ISLCap; % Mbps
                        % �C������ǿ骺��Ʀʤ���
                        MatTrafego = MatRede.*0;
                        % �X�o�a/�ت��a�ݨD�ʤ���
                        MatDemanda = MatRede.*0;       
                        % ���ݪA�ȥX�o�a/�ت��a
                        MatDemandaAtendida = MatRede.*0;
                       
                        % ----  ���y�׺ݹ�]�ӷ��M�ؼС^ --------  
                        for fonte_id=1: obj.TotalFontes   
                            Origem  = Resultados(LoopId).Fontes(fonte_id).SateliteOrigem;
                            Destino = Resultados(LoopId).Fontes(fonte_id).SateliteDestino;
                            Demanda = cbr; % CBR, Mbps
                            % --------------------- ���Ъ����� -------------
                            switch Metrica    
                                case 'TP'    
                                    G = digraph(MatDelay);
                                    % �Ǽ��ɶ�                    
                                    Tmin = min(G.Edges.Weight); % �Q����  
                                    Tmax = max(G.Edges.Weight); % �̧C�I  
                                    % �k�@�ƶǼ��ɶ�
                                    TPNorm = (G.Edges.Weight - Tmin)/(Tmax - Tmin);
                                    G.Edges.Weight = TPNorm;           
                                case 'LASER'        
                                    G = digraph(MatDelay);
                                    % �Ǽ��ɶ�                   
                                    Tmin = min(G.Edges.Weight); % �Q���� 
                                    Tmax = max(G.Edges.Weight); % �̧C�I  
                                    % �k�@�ƶǼ��ɶ�
                                    TPNorm = (G.Edges.Weight - Tmin)/(Tmax - Tmin);
                                    [Si,Sj] = findedge(G);
                                    Ei = [obj.Cenario.Propagacao(LoopId).Satelites(Si).EclipseStatus]';
                                    Ej = [obj.Cenario.Propagacao(LoopId).Satelites(Sj).EclipseStatus]'; 
                                    % �q���q�q�зǤ�
                                    Bi = Energia(Si)./obj.PotenciaCP;
                                    Bj = Energia(Sj)./obj.PotenciaCP;
                                    Dij = Ei./Bi + Ej./Bj;
                                    Dij(isnan(Dij)) = 0; 
                                    Dmin = min(Dij); % �Q����  
                                    Dmax = max(Dij); % �̧C�I 
                                    DijNorm = (Dij-Dmin)/(Dmax-Dmin);
                                    G.Edges.Weight = 0.5.*TPNorm+0.5.*DijNorm;   
                                    
                                case 'PROPOSTA' 
                                    G = digraph(MatDelay);
                                    % �Ǽ��ɶ�                    
                                    Tmin = min(G.Edges.Weight); % utopia   
                                    Tmax = max(G.Edges.Weight); % nadir     
                                    % �k�@�ƶǼ��ɶ�
                                    TPNorm = (G.Edges.Weight - Tmin)/(Tmax - Tmin);
                                    [Si,Sj] = findedge(G);
                                    Ei = [obj.Cenario.Propagacao(LoopId).Satelites(Si).EclipseStatus]';
                                    Ej = [obj.Cenario.Propagacao(LoopId).Satelites(Sj).EclipseStatus]'; 
                                    Bi = Energia(Si)./obj.PotenciaCP;
                                    Bj = Energia(Sj)./obj.PotenciaCP;
                                    % �s�ҫ��A�ȦҼ{
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
                                        pos = find( xij == 1); % �n�������ìP����m�ʤ���
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
                                    
                                    % �q�s�����R�� - �ìP�M
                                    % �q���s��������
                                    % �q���e�q�M�s���e�q
                                    
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
                                         Cij(pos, :) = NaN; % ���t NaN �H�q�зǤƤ��R��
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
                                      CijNorm(pos, :) = inf; % �R���S����O���s��
                                      DijNorm(pos, :) = inf; % �L�ݹq���Y�i�νìP�h����t
                                      TPNorm(pos, :) = inf;  % �L�ݹq���Y�i�νìP�h����t  
                                      %Se for NaN, Not a Number, atribui
                                      %infinito..
%                                       CijNorm(isnan(CijNorm)) = inf; 
%                                       DijNorm(isnan(DijNorm)) = inf ; 
%                                       TPNorm(isnan(TPNorm)) = inf; 
                                    end
                                      % �p�G���O�Ʀr�h�T�θ��u
                                      
                                     CijNorm(isnan(CijNorm)) = 0; 
                                     DijNorm(isnan(DijNorm)) = 0 ; 
                                     TPNorm(isnan(TPNorm)) = 0;  
                                    
                                    %%%% ��ĳ���СGw1�Bw2�Bw3 %%%%%%%
                                    %G.Edges.Weight = 0.35.*TPNorm+0.35.*DijNorm+0.3.*CijNorm;
                                    TDC = [TPNorm DijNorm CijNorm];
                                    G.Edges.Weight = diag( Peso(:,:)  * TDC(:,:)');
                                    %Resultados(LoopId).(Metrica).Fontes.demanda(fonte_id,cbr_id) = Demanda;
                             end
                            
                             DemandaAtendida = Demanda; % Mbps
                             % �ɶ����j���X�o�a/�ت��a�ݨD���ʤ���x�}
                              MatDemanda(Origem,Destino) = MatDemanda(Origem,Destino) ...
                                + Demanda*obj.Cenario.TempoIntervaloSimulacao;  
                              
                              [rota, custo] = shortestpath(G,Origem,Destino); %%%%% �ӶO�Υ��ϥ� %%%%%%
                                % ���|�W�ìP����q
                                EnergiaRota =  Energia(rota);  
                                % �̤p���|��q
                                
                                MinEnergiaRota = min(EnergiaRota);
                                %%%% ��q Tx�BRx �O�ìP k-l �������C�ӳs��
                                % �ǿ�һݯ�q
                                EnergiaTx = ((DemandaAtendida/obj.ISLCap)*obj.PotenciaTX)*obj.Cenario.TempoIntervaloSimulacao;
                                % �����һݯ�q�ʤ���
                                EnergiaRx = ((DemandaAtendida/obj.ISLCap)*obj.PotenciaRX)*obj.Cenario.TempoIntervaloSimulacao;
                                % �ǿ�M�����һݪ��`��q
                                EnergiaTxRx = EnergiaTx + EnergiaRx;

                                % �ˬd�O�_���������\�v�Ӷǿ�
                                if MinEnergiaRota < EnergiaTxRx
                                   DemandaAtendida = 0;
                                end  
                                
                            CapRota = [];
                            % ���Ѹ����`�p�ʤ���
                            saltos = size(rota,2)-1;
                            % ���@�����e�q�����u
                            for st=1:saltos        
                                % ���Ƴs���e�q
                                CapEnlace = Capacidade(rota(st),rota(st+1)); % Mbps
                                CapRota = cat(1,CapRota,CapEnlace);
                                if CapEnlace < DemandaAtendida % �p�G�����s���e�q
                                   % ����e�q�ݨD
                                   DemandaAtendida = CapEnlace;
                                end
                            end 
                            
                            % �s���Ǽ��ɶ��]�H�@�����^
                            Delay = 0;
                            % �p�G����ƭn�b�Ӹ��u�W��p�A�h�������
                            if DemandaAtendida > 0
                                
                                for st=1:saltos        
                                    % ��s��q - �ǿ�
                                    Energia(rota(st)) = Energia(rota(st)) - EnergiaTx;
                                    % ��s��q - ����
                                    Energia(rota(st+1)) = Energia(rota(st+1)) - EnergiaRx; 
                                    % ��s����e�q
                                    Capacidade(rota(st),rota(st+1))= Capacidade(rota(st),rota(st+1))...
                                        - DemandaAtendida;
                                    % �s���y�q�x�}
                                    MatTrafego(rota(st),rota(st+1))...
                                        = MatTrafego(rota(st),rota(st+1))+DemandaAtendida...
                                        *obj.Cenario.TempoIntervaloSimulacao; % (Mb) = �ݨD (MB/s)* �������j�ɶ� (s)�F
                                    % ���ѩ���]�Ǽ��ɶ��^
                                    Delay = Delay+MatDelay(rota(st),rota(st+1));   
                                    
                                    %%%%%#### �K�[��q�ǿ�
                                    %%%%% ���ݪ�O
                                    ConsumoEnergia(rota(st)) = ConsumoEnergia(rota(st)) + EnergiaTx; % ��q�����ಾ
                                    ConsumoEnergia(rota(st+1)) =  ConsumoEnergia(rota(st+1)) + EnergiaRx; % �����\��
                                    %%%%%%%
                                end
                            else
                                % �S���ǰe���    
                                saltos = 0;          
                            end
                             % �ɶ����j���������ݨD�x�}�ʤ���
                             MatDemandaAtendida(Origem,Destino) = MatDemandaAtendida(Origem,Destino)...
                                + DemandaAtendida*obj.Cenario.TempoIntervaloSimulacao; 
                             % �ӷ��K�n
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
                            
                            %%%%% �C�ӽìP�s������� %%%% 
%                           Consumo =  Resultados(LoopId).(Metrica).ConsumoEnergia(Origem,cbr_id);
%                           Consumo =  EnergiaON(Origem) + TotalTxRx; %+ Consumo;
%                           Resultados(LoopId).(Metrica).ConsumoEnergia(Origem,cbr_id) = Consumo;                             
                        end
                      
                        % �ݨD
                        %Resultados(LoopId).(Metrica).Demanda(:,cbr_id) = {MatDemanda};
                        %Resultados(LoopId).(Metrica).DemandaAtendida(:,cbr_id) = {MatDemandaAtendida};     
                        % ���j�����ɪ��s���e�q
                        Resultados(LoopId).(Metrica).Trafego(:,cbr_id) = {MatTrafego};  % MBytes           
                        % ���j�����ɪ��s���e�q
                        Resultados(LoopId).(Metrica).Capacidade(:,cbr_id) = {Capacidade};
                        % ���j��q����
                        CaptacaoEnergetica = obj.PotenciaCG*max(0,obj.Cenario.TempoIntervaloSimulacao-[obj.Cenario.Propagacao(LoopId).Satelites.EclipseTempoIntervalo]');           
                        %%%%% ��q���� ...
                        %%% �K�[�H�p���q
                        %CaptacaoEnergetica = Resultados(LoopId).(Metrica).CaptacaoEnergetica(:,cbr_id) + CaptacaoEnergetica
                        %%%% ��l�p�� %%%
                        Resultados(LoopId).(Metrica).CaptacaoEnergetica(:,cbr_id) = CaptacaoEnergetica;
                        % ��s��q�x�}
                        
                        % �Ѿl��q = ��q 
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
                        Energia = max(0, Energia); % �K�[�H�קK��q���t��
                       
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