classdef Consolidar
    properties
        Diretorio;
        TotalSimulacoes;
        Simulacao;
        Estimadores;
    end
    
    methods   
        

        function obj = Consolidar(Diretorio,TotalSimulacoes)
          obj.Diretorio = Diretorio;
          obj.TotalSimulacoes = TotalSimulacoes;
          load(strcat(Diretorio,'Simulacao.mat'));
          obj.Simulacao = Simulacao;
          obj.Estimadores.Amostras = obj.CalcAmostras();
        end
        
        function obj = set.Estimadores(obj,val)
          obj.Estimadores = val;
        end       

        function Estimadores =  CalcAmostras(obj)    
                   
            for simulacao_id=1:obj.TotalSimulacoes              
                % 交通
                Resultados ={};
                fprintf('Consolidar Amostras (%d/%d)\n',simulacao_id,obj.TotalSimulacoes);
                load(strcat(obj.Diretorio,'Resultados_',num2str(simulacao_id),'.mat'));  
                for metrica_id=1:size(obj.Simulacao.Metricas,1)  
                    Metrica = num2str(cell2mat(obj.Simulacao.Metricas(metrica_id)));   
                    for loopId =1:size(Resultados,2)  
                        % 需求滿足百分比
                        DemandaAtendida.(Metrica)(loopId,:) ...
                            = mean(Resultados(loopId).(Metrica)...
                            .Fontes.demandaatendida)./obj.Simulacao.CBRs;   
                         % 被封鎖的字體總數百分比
                         FontesBloqueadas.(Metrica)(loopId,:)  ...
                             = sum(Resultados(loopId).(Metrica).Fontes.demandaatendida  == 0);
                         % 平均延遲百分比
                         Delay.(Metrica)(loopId,:) = sum(Resultados(loopId)...
                             .(Metrica).Fontes.tempopropagacao)./sum((Resultados(loopId)...
                             .(Metrica).Fontes.tempopropagacao>0));          
                         % 平均跳躍
                         Saltos.(Metrica)(loopId,:) = sum(Resultados(loopId)...
                             .(Metrica).Fontes.saltos)./sum((Resultados(loopId)....
                             .(Metrica).Fontes.saltos>0));                
                          % 飽和連結總數百分比
                          EnlacesSaturados.(Metrica)(loopId,:)  = sum(Resultados(loopId)...
                              .(Metrica).Fontes.enlacessaturados)./size(find(obj.Simulacao...
                              .Cenario.Propagacao(loopId).Enlaces >0),1);     
                          % 日蝕路線
                          EclipseRota.(Metrica)(loopId,:)  = sum(Resultados(loopId)...
                               .(Metrica).Fontes.eclipserota)...
                               ./sum(Resultados(loopId).(Metrica).Fontes.eclipserota >0) ;                           
                    end
                    Estimadores.(Metrica).Delay(simulacao_id,:) = mean(Delay.(Metrica));
                    Estimadores.(Metrica).DemandaAtendida(simulacao_id,:) = mean(DemandaAtendida.(Metrica));
                    Estimadores.(Metrica).Saltos(simulacao_id,:) = mean(Saltos.(Metrica));
                    Estimadores.(Metrica).FontesBloqueadas(simulacao_id,:) = mean(FontesBloqueadas.(Metrica));
                    Estimadores.(Metrica).EnlacesSaturados(simulacao_id,:) = mean(EnlacesSaturados.(Metrica));
                    Estimadores.(Metrica).EclipseRota(simulacao_id,:) = mean(EclipseRota.(Metrica));    
                end  
            end        
        end         
    end
end