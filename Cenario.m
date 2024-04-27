classdef Cenario
	properties (Constant)
        mu = 398600.5;              % 地球引力常數
        rp=6371;                    % 地球半徑 - km 
        rs = 695500;                % 太陽半徑 - km
	end
     
	properties  ( Access = public )
        TempoTotalSimulacao;        % 總模擬時間 seconds
        TempoIntervaloSimulacao;	% 模擬間隔 seconds      
        TotalSatelites = 66;        % 星座中衛星總數
        TotalPlanos = 6;            % 軌道平面總數
        TotalSatPlanos = 11;        % 每個平面的衛星總數
        BordaPolar = 75;            % 極座標緯度
        Zonas;                      % 將土地劃分為區域
        Constelacao;                % 星座      
        Propagacao;                 % 傳播
        Eclipse;                    % 日蝕
    end
    
    methods
        
        %% Construtor
        function obj = Cenario(Data,TempoTotalSimulacao,TempoIntervaloSimulacao)
            
            obj.Constelacao	= obj.InfoTLE();
            obj.Zonas =  Zonas();
            obj.TempoTotalSimulacao = TempoTotalSimulacao; 
            obj.TempoIntervaloSimulacao = TempoIntervaloSimulacao;	
            obj.Propagacao = obj.Propagar(datevec(Data));  
            obj.Eclipse = obj.DadosEclipse();
        end

        %% 從 TLE 檔案中提取數據
        function DadosOrbitais = InfoTLE(obj)        
                fd = fopen('TLE.txt', 'rb');
                A0 = fgetl(fd);
                A1 = fgetl(fd);
                A2 = fgetl(fd);
                DadosOrbitais = {};     % 軌道資料
                cont=1;                 % 記數
                while ischar(A2)          
                    % 偏心率
                    DadosOrbitais(cont).epoca = obj.Epoca2Data(str2num(A1(19:32)));	 
                    % 傾角 inclination
                    DadosOrbitais(cont).e = str2double(A2(27:33))/ (1e7);
                    % 升交點黃經
                    DadosOrbitais(cont).i = str2double(A2(9:16));  
                    % 近地點角
                    DadosOrbitais(cont).raan = str2double(A2(18:25)); 
                    % 平均異常
                    DadosOrbitais(cont).w = str2double(A2(35:42));  
                    % 平均移動百分比
                    DadosOrbitais(cont).M = str2double(A2(44:51));    
                    % 軌道資料續
                    n = str2double(A2(53:63));  
                    DadosOrbitais(cont).n = n;   
                    % 半長軸
                    DadosOrbitais(cont).a = (obj.mu/(n*2*pi/(24*3600))^2)^(1/3);	
                    cont = cont + 1; 
                    A0 = fgetl(fd);
                    A1 = fgetl(fd);
                    A2 = fgetl(fd); 
                end
                fclose(fd);               
        end
        
        %% 傳回 TLE 紀元的日期/時間向量
        function vetor_data = Epoca2Data(obj, tle_epoca )
            ymd=floor(tle_epoca);
            yr=fix(ymd/1000);
            dofyr=mod(ymd,1000);
            if (yr < 57)
                 year=  yr+ 2000;
            else
                 year=  yr+ 1900;
            end; 
            decidy=round((tle_epoca-ymd)*10^8)/10^8;
            temp=decidy*24;
            hh=fix(temp);
            temp=(temp-hh)*60;
            mm=fix(temp);
            temp=(temp-mm)*60;
            ss=floor(temp);
            nd = eomday(year,1:12);
            temp = cumsum(nd);
            month=find(temp>=dofyr, 1 );
            temp=temp(month)-dofyr;
            date=nd(month)-temp;
            vetor_data=[year,month,date,hh,mm,ss];
        end
 
        %% 傳回 ECI 中的位置向量，從軌道元素開始
        function rECI = Oe2Eci(obj,SatNum,Data)  
                DeltaT = etime(Data,obj.Constelacao(SatNum).epoca);
                j2 = 1.0826359*10^-3;                     % Pertuation J2（地球變平）           
                a = obj.Constelacao(SatNum).a;            % 半長軸       
                e = obj.Constelacao(SatNum).e;            % 偏心率       
                i = obj.Constelacao(SatNum).i;            % 傾角 inclination
                raan = obj.Constelacao(SatNum).raan;      % 升交點黃經    
                w = obj.Constelacao(SatNum).w;            % 近地點角
                M = obj.Constelacao(SatNum).M;            % 平均異常百分比          
                p = a*(1-e^2);                            %       
                T = 2*pi*sqrt(a^3/obj.mu);                % 週期,seconds      
                Tp = (M/360*T);                 % 透過近地點的時間
                ti = Tp+DeltaT;                 % 考慮近地點的模擬初始時間
                j2_raan = -(3/2)*(j2)*((obj.rp/p)^2)*sqrt(obj.mu/(a^3))*cosd(i);	   % 升交點經度的擾動     
                j2_w = (3/4)*(j2)*((obj.rp/p)^2)*sqrt((obj.mu/a^3))*(5*cosd(i)^2-1) ;  % 近地點角擾動
                M = mod(360*ti/T,360);          % 平均異常更新
                E = obj.CalcKeplerEq(M,e);    	% 偏心異常更新             
                v = mod(atan2((sind(E)*(1-e^2)^.5),(cosd(E)-e)),2*pi)*180/pi;  % 真實異常更新
                RAAN = ((raan*pi/180)+j2_raan*DeltaT)*180/pi;                  % 升交點經度更新
                w = ((w*pi/180)+j2_w*DeltaT)*180/pi;  % 近地點參數更新     
                
                Quaternio = quatnormalize(angle2quat(RAAN*pi/180,i*pi/180,v*pi/180+w*pi/180,'ZXZ'));
                % 計算四元數旋轉矩陣
                Rot_m = quat2rotm(Quaternio) ;  % 四元數旋轉矩陣
                % 軌道半徑
                r = a*((1-(e*cosd(E))));
                % ECI 中軌道向量「r」的計算
                rECI = Rot_m*[1;0;0]*r;    
        end   
        
        %% 傳回 ECEF 中的位置向量
        function [rECEF] = Eci2Ecef(obj,rECI,Data)
            JD = juliandate(Data);
            Tj = (JD - 2451545.0)/36525.0;  
            GMST = mod(280.46061837 + 360.98564736629*(JD - 2451545) + 0.000387933*Tj^2- Tj^3/38710000,360);
            R3 = [cosd(GMST),sind(GMST),0;-sind(GMST),cosd(GMST),0;0,0,1];
            rECEF = R3*rECI; 
        end
        
        %% 回傳偏心異常 - E
        function [E] = CalcKeplerEq(obj,M,e)
        tol = 10^-8;
        M = M*pi/180; % 弧度
        Etemp = M;
        ratio = 1;
        while abs(ratio) > tol
            f_E = Etemp - e*sin(Etemp) - M;
            f_Eprime = 1 - e*cos(Etemp);
            ratio = f_E/f_Eprime;
            if abs(ratio) > tol
                Etemp = Etemp - ratio;
            else
                E = Etemp;
            end
        end
        E = mod(E,2*pi)*180/pi; % Graus
        end       
        
        %% 根據 ECEF 中的位置向量計算緯度和經度
        function [Lat,Long] = LatLong(obj,rECEF)
            r_delta = norm(rECEF(1:2));
            sinA = rECEF(2)/r_delta;
            cosA = rECEF(1)/r_delta;
            Long = atan2(sinA,cosA);
            if Long < -pi
                Long = Long + 2*pi;
            end
            Lat = asin(rECEF(3)/norm(rECEF));    
            Lat = Lat*180/pi;
            Long = Long*180/pi;
        end
        
        %% 計算本影中心與本影與半影末端的距離
        function [Eu Kp] = CalcEuKp(obj,vetor_rs,Sps)
            % Xu = 本影錐頂點與地心之間的距離 - Km
            Xu = (obj.rp*2*Sps)/(obj.rs*2-obj.rp*2); 
            AlphaU = asin((2*obj.rp)/(2*Xu)); % 本影角
            Eu = (Xu-norm(vetor_rs))*tan(AlphaU);
            % Xp = 半影錐頂點與地球中心之間的距離 - Km
            Xp = (2*obj.rp*Sps)/(2*obj.rs+2*obj.rp);
            AlphaP = asin((2*obj.rp)/(2*Xp)); % 半影角 
            Kp = (Xp+norm(vetor_rs))*tan(AlphaP);
        end 
        
        %% 計算日食發生
        function [EclipseStatus] = CalcEclipse(obj,Data,R)
            JD = juliandate(Data);
            % 從 2000 年起的儒略世紀
            T = (JD - 2451545.0)/36525.0;  
            % 太陽的幾何平均經度
            L0 = mod(280.46645 + 36000.76983*T + 0.0003032*T^2,360);
            % 平均太陽異常
            M = mod(357.52910 + 35999.05030*T - 0.0001559*T^2 - 0.00000048*T^3,360);
            % 地球軌道偏心率
            e = mod(0.016708617 - 0.000042037*T - 0.0000001236*T^2,360);
            % 太陽中心方程，相對於其平均距平
            C = mod((1.914600-0.004817*T-0.000014*T^2)*sind(M)+(0.019993-0.000101*T)*sind(2*M)+0.000290*sind(3*M),360);
            % 太陽的黃道經度（lambda）
            Ls = mod((L0 + C),360)*pi/180;  
            % 太陽的真實異常百分比
            f=mod(M+C,360);
            % 地球 <-> 太陽距離
            RS = ((1.000001018*(1 - e^2))/(1 + e*cosd(f)))*149597870.70;
            % 黃道傾角
            Ep = mod(23 + 26/60 + 21.448/3600 - 46.8150/3600*T - (0.00059/3600)*T^2+ (0.001813/3600)*T^3,360)*pi/180;
            % Unit 向量太陽
            vetor_solar = [cos(Ls),sin(Ls)*cos(Ep),sin(Ls)*sin(Ep)];   
            vetor_rs = dot(R,vetor_solar)*vetor_solar;
            vetor_projecao = R'-vetor_rs; 
            [Eu Kp] = obj.CalcEuKp(vetor_rs,RS);
                % 日食狀態：0 - sem 日食，1 = 半影，2 本影 
                if(dot(R,vetor_solar) < 0);   
                    NORM = norm(vetor_projecao);
                    if(NORM > Kp) 
                        EclipseStatus = 0;
                    end    
                    if(Eu<NORM && NORM <Kp)  
                        EclipseStatus = 1; % Penumbra                   
                    end   
                    if(NORM < Eu)            
                        EclipseStatus = 1; % Umbra
                    end 
                else
                    EclipseStatus = 0; 
                end   
        end      

        %% 衛星隨時間傳播
        function Propagacao = Propagar(obj,Data,Interacoes)
            if nargin == 2
                Interacoes = ceil(obj.TempoTotalSimulacao/obj.TempoIntervaloSimulacao) ;
            end
            LoopId = 1;
            while LoopId <= Interacoes  
                fprintf('場景 - 第 %d 次交互，共 %d\n',LoopId,Interacoes);  
                Propagacao(LoopId).Data = Data;
                Satelites = [];
               
                for SatNum=1:obj.TotalSatelites   
                    rECI  = obj.Oe2Eci(SatNum,Data);
                    Satelites(SatNum).rECI = rECI;
                    %## ECEF 中的位置向量
                    rECEF = obj.Eci2Ecef(rECI,Data);   
                    %## 確定緯度和經度
                    [Lat,Long] = obj.LatLong(rECEF);    
                    Satelites(SatNum).Latitude = Lat;          
                    Satelites(SatNum).Longitude = Long;
                    %## 日食發生
                    EclipseStatus =  obj.CalcEclipse(Data,rECI);
                    % 日食狀況：1 - 日食，0 - 太陽
                    Satelites(SatNum).EclipseStatus = EclipseStatus;     
                    % 日食中的曝光時間百分比
                    TempoEclipseIntervalo = 0;   
                    if EclipseStatus ~= 0 % 衛星在 日蝕 啟動
                        % 檢查間隔內日食中剩餘的時間
                        for i=1:obj.TempoIntervaloSimulacao 
                            % 日期增量
                            DataIntervalo = datevec(addtodate(datenum(Data), i, 'second'));
                            % 計算新日期的軌道半徑
                            R = obj.Oe2Eci(SatNum,DataIntervalo);
                            % 計算狀態
                            Eclipse = obj.CalcEclipse(DataIntervalo,R);  
                            if Eclipse ~= EclipseStatus % 已更改
                                TempoEclipseIntervalo = TempoEclipseIntervalo+1;
                                break
                            else
                                TempoEclipseIntervalo = TempoEclipseIntervalo+1;
                            end
                        end
                    else % 衛星從有太陽情形下開始      
                        % 檢查時間間隔內留在陽光下的時間
                        for i=1:obj.TempoIntervaloSimulacao 
                            % 日期增加
                            DataIntervalo = datevec(addtodate(datenum(Data), i, 'second'));
                            % 計算新日期的軌道半徑
                            R = obj.Oe2Eci(SatNum,DataIntervalo);
                            % 計算狀態
                            Eclipse = obj.CalcEclipse(DataIntervalo,R);  
                            if Eclipse ~= EclipseStatus % 已更改                         
                                 TempoEclipseIntervalo = obj.TempoIntervaloSimulacao-i; 
                                 Satelites(SatNum).EclipseStatus = Eclipse; 
                                break
                            end
                        end
                    end
                     Satelites(SatNum).EclipseTempoIntervalo =  TempoEclipseIntervalo;   
                    if LoopId == 1 % 初始條件
                        Satelites(SatNum).EclipseTempoTotal	= TempoEclipseIntervalo;        
                    else
                        if(Satelites(SatNum).EclipseStatus ~=0) 
                             % 在日蝕中
                             % 累計迭代的日食時間
                             % 上一個與當前時間
                            Satelites(SatNum).EclipseTempoTotal = Propagacao(LoopId-1).Satelites(SatNum).EclipseTempoTotal+TempoEclipseIntervalo;
                            %fprintf(" SatNum %d Satelites(SatNum).EclipseTempoTotal %d \n", SatNum, Satelites(SatNum).EclipseTempoTotal);  
                        else
                            % 更新累積的日蝕時間
                            Satelites(SatNum).EclipseTempoTotal = Propagacao(LoopId-1).Satelites(SatNum).EclipseTempoTotal; % Estava em zero(0);
                        end  
                     end     
                end
                Propagacao(LoopId).Satelites = Satelites;
                % 活躍連結
                q = obj.TotalSatPlanos; % 每個平面上的衛星數量
                Adj = zeros(obj.TotalSatelites); % 鄰接矩陣分配
                for k=1:obj.TotalPlanos  % 掃描所有平面
                    for j=1:q     % 掃描單一平面上的衛星            
                        SatNum = (k*q)-(q-j); % 衛星 1 的編號順序... t    
                        satAnt   = SatNum-1;  % 軌道前的衛星
                        satPost  = SatNum+1;  % 軌道後的衛星
                        planAnt  = SatNum-q;  % 先前平面上的衛星百分比
                        planPost = SatNum+q;  % 後平面上的衛星
                        %%%%% 添加 GSO - GeoStationary Satllite 的連結 %%%%%% 
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if(j==1) satAnt   = k*q; end            % 軌道前第一顆衛星
                        if(j==q) satPost  = k*q-(q-1); end      % 軌道後最後一顆衛星       
                        if(k==1) planAnt  = ((obj.TotalPlanos-1)*q)+j; end    % 前一條軌道的第一顆衛星
                        if(k == obj.TotalPlanos) planPost = j; end            % 前一條軌道的最後一顆衛星            
                        % ----  同一方案的 ISL - 固定成本           
                        La = sqrt(2)*norm(Satelites(SatNum).rECI)*sqrt(1-cosd(360/q)); % 連結成本
                        Adj(SatNum,satAnt)  =  La; % 來自同一軌道的前一顆衛星
                        Adj(SatNum,satPost) =  La; % 同一軌道的後來衛星          
                        if abs(Satelites(SatNum).Latitude) < obj.BordaPolar % 衛星位於“Simulacao.BordaPolar”，不存在 ISL 平面內                 
                            Le = sqrt(2)*norm(Satelites(SatNum).rECI)*sqrt(1-cosd(360/(2*obj.TotalPlanos)))*cosd(Satelites(SatNum).Latitude); % 維持平面之間的ISL
                            if  abs(Satelites(planAnt).Latitude) < obj.BordaPolar % 緯度 坐在前面                 
                               if(k>1) % seam 
                                     Adj(SatNum,planAnt) = Le;  % 衛星給出先前的軌道
                               end
                            end            
                            if abs(Satelites(planPost).Latitude) < obj.BordaPolar % 緯度 坐在後面
                                if(k<obj.TotalPlanos) % 接縫
                                   Adj(SatNum,planPost) =  Le;  % 衛星給出後軌道                                                                              
                                end
                            end 
                       end             
                    end
                end          
                
                Propagacao(LoopId).Enlaces = Adj;
                
                %% 各地區和當地時間的衛星百分比
                tam = size(obj.Zonas.LatLong);
                Sats = zeros(tam);
                Horas = zeros(tam);
                SatLatLong = [[Propagacao(LoopId).Satelites.Latitude]', [Propagacao(LoopId).Satelites.Longitude]'];
                diff = (360/tam(2))/2; % 考慮經絡的度數差異
                for i=1: tam(1)
                    for j=1:tam(2)           
                        % 計算最近的衛星
                        latlon=obj.Zonas.LatLong(i,j);
                        [k, dist] = dsearchn([latlon{1}(1),latlon{1}(2)],SatLatLong);   
                        [dist, SatNum] = min(dist);
                        Sats(i,j) = SatNum;   
                        lon = latlon{1}(2); %longitude       
                        t = datenum(Data);
                        zd = timezone(lon+diff,'degrees');
                        t = addtodate(t, -zd, 'hour');
                        hora = hour(t);
                        if hora ==0 
                         hora = 24;
                        end
                        Horas(i,j) = hora;            
                    end
                end
                Propagacao(LoopId).Zonas.Satelites = Sats; % 按地區分列的衛星覆蓋率             
                % 依小時概況的需求流量 500TB/天
                % 原程式碼中並沒有解釋500000000000000這個值，我們研究了這個值，就是每天500TB = 5x10^14 bytes/day。
                % 發送封包的平均大小為 210 bytes
                [ti,tj] = size(obj.Zonas.MatrizTrafego);
                Fk = zeros(max(obj.Zonas.Continentes(:))); % 各大洲流量
                Tij = zeros(ti,tj);
                for i=1:ti
                   for j=1:tj       
                       Tij(i,j) = (obj.Zonas.MatrizTrafego(i,j)/...
                           sum(obj.Zonas.MatrizTrafego(:)))*(500000000000000/3600)...
                           *((obj.Zonas.AtividadeHora(Horas(obj.Zonas.ZonaId(i))))/100);         
                       % 原產地大陸
                       CkOrg = obj.Zonas.Continentes(obj.Zonas.ZonaId(i));
                       % 目的地大陸 
                       CkDst = obj.Zonas.Continentes(obj.Zonas.ZonaId(j));   
                       Fk(CkOrg,CkDst)  = Fk(CkOrg,CkDst)+Tij(i,j); 
                   end        
                end              
                % 各大洲之間的流量
                tk = max(obj.Zonas.Continentes(:));
                for k=1:tk
                    Propagacao(LoopId).Zonas.Fluxo(k,:) = Fk(k,:)/sum(Fk(k,:)); 
                end 
                % 根據模擬間隔更新模擬日期/時間
                Data = datevec(addtodate(datenum(Data), obj.TempoIntervaloSimulacao, 'second'));                  
                % 增量循環控制
                LoopId = LoopId+1;       
            end
            %% 計算總日蝕
            
            
        end
        
        function Eclipse = DadosEclipse(obj)
            tam = size(obj.Propagacao,2); % 場景數量
            Saida=zeros(obj.TotalSatelites,1);
            Entrada=zeros(obj.TotalSatelites,1);
            Tempo=zeros(obj.TotalSatelites,1);
            for s=1:obj.TotalSatelites
                   ciclo = 1; % 代表公轉(繞地球轉)
                   for i=2:tam  
                       % 日食結束
                       % 累積日食時間
                       if(obj.Propagacao(i).Satelites(s).EclipseStatus  == 0 && obj.Propagacao(i-1).Satelites(s).EclipseStatus ~=0)  
                           Saida(s,ciclo) = i-1;
                           EclipseTempo = obj.Propagacao(i-1).Satelites(s).EclipseTempoTotal ;  
                           Entrada(s,ciclo) = Saida(s,ciclo)-round(EclipseTempo/obj.TempoIntervaloSimulacao); 
                           Tempo(s,ciclo) = ceil(EclipseTempo/60);   
                           ciclo = ciclo+1;
                       end
                   end
            end
            Eclipse.Entrada = Entrada;
            Eclipse.Saida = Saida;
            Eclipse.Tempo = Tempo;
            Eclipse.TempoMaximo = max([Tempo]')';         
        end
    end    
end