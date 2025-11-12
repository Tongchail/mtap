function  [H2Osat]  =  fluidsat(var)

MW  = [60.0855, 79.88, 101.96, 71.85, 40.3, 56.08, 61.98, 94.2, 18.02]; 

P    = min(var.P,0.5)*1e9/1e5;
T    = var.T+273.15;                   % 温度从 °C 转为 K。
X    = var.X./MW ./ sum(var.X./MW,2);  % 把输入的质量分数 var.X 转成摩尔分数
X    = X./sum(X(:,1:end-1),2);         % 再次归一化：只对“无水氧化物”部分（第 1–8 列）做和为 1 的归一化。H2O不参与归一化
ln_fH2O = log(P);                      % 取 fH₂O 的自然对数。这里近似 fH₂O ≈ P (bar)，即把水的逸度当作压力（理想近似）

a = 2565;
b = [-1.997,-0.9275,2.736];
c = 1.171;
d = -14.21;

XH2O   = exp((a./T + sum(b.*X(:,[3,4,7]),2).*(P./T) + c.*ln_fH2O + d)/2); % Moore et al., 1998 公式计算溶解于熔体中的 H₂O 的摩尔分数 XH2O（注意：这里是“溶解的水”这一成分在熔体中的 mole fraction）。

H2Osat = XH2O.*MW(end) ./ sum([X(:,1:end-1), XH2O].*MW,2);  % 把 XH2O（水的摩尔分数）转换成质量分数（wt fraction），输出的 H2Osat 就是饱和水含量（wt fraction）；如果要 wt%，乘以 100。

end