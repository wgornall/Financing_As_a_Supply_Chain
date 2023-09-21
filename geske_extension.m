
baseEconomy = struct('rf', 0.025*5,'rfPerYear', 0.025, 'tau', 0.25, 'alphab', .2, 'T',5,'ForceSpread',-Inf,...
    'bout',0,'bout2',0,'depi',0,'fracown',.5,'b2prob',.5,...
    'CBaselIRB',-Inf,'CBaselS',-Inf,  'CFirmLevU',-Inf, 'CFirmLevL',-Inf, 'CBankLevU',-Inf,  'CBankLevL' ,-Inf, 'Assets', []  );

baseAssetWeights = [    1.008512262742020   1.041302859831466   0.221789624087612  ];

baseAssets = struct();
baseAssets.corporateLoans  = struct('weight',  baseAssetWeights(1),'AssetSpecificSpread',-Inf,'rho',.2,'sigma',.4,'tau',.25,'alpha',0.1,'DefaultProtected',0,'FloatingLeverage',true,'maxload',1,'Repayment',.3,'riskWeight',1,'corrFormula','c');
baseAssets.mortgage  = struct('weight', baseAssetWeights(2),'AssetSpecificSpread',-Inf,'rho',.2,'sigma',.25,'tau',0,'alpha',0.25,'DefaultProtected',0.5,'FloatingLeverage',false,'maxload',-Inf,'Repayment',.8*exp( -baseEconomy.rf ),'riskWeight',.5,'corrFormula','r');
baseAssets.cash  = struct('weight', baseAssetWeights(3),'AssetSpecificSpread',-Inf,'rho',.5,'sigma',.5,'tau',0,'alpha',0,'DefaultProtected',1,'FloatingLeverage',false,'maxload',-Inf,'Repayment',1,'riskWeight',0,'corrFormula','0');

baseBankH = setfield(baseEconomy,  'Assets',[baseAssets.mortgage]);

baseBankH = setfield(FSC_H_findCS(baseBankH),'val',0)


%% Basic Setup

Periods = 1;
Weights = ones(1,Periods)/Periods*baseAssetWeights(2);
Repayments = .55540;
nnn = 1001;
BankLoanIRate = 0.0021;


MedianCumulativeReturns = linspace(-10,10,nnn)';

DiscountFactor = exp(-baseEconomy.rfPerYear*baseBankH.T*(1:Periods)/Periods);

OneA = ones(nnn,1);

LoanRepayment = @(PeriodOfPortfolio) .80*DiscountFactor(PeriodOfPortfolio);
VolAtPeriod = @(PeriodOfPortfolio) baseAssets.mortgage.sigma* sqrt(baseBankH.T*PeriodOfPortfolio/Periods);

BankPortfolioCashFlow = @(y,PeriodOfPortfolio) (baseAssets.mortgage.DefaultProtected* LoanRepayment(PeriodOfPortfolio) + ...
    (1- baseAssets.mortgage.DefaultProtected) *(...
    LoanRepayment(PeriodOfPortfolio) .* normcdf( ( - log( LoanRepayment(PeriodOfPortfolio))  - 0.5 .* VolAtPeriod(PeriodOfPortfolio)^2 + sqrt(baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio) .* y/PeriodOfPortfolio^.5  ) ./ ( sqrt(1-baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio)) ) ...
    +  (1-baseAssets.mortgage.alpha).*(1-baseAssets.mortgage.tau).*exp( sqrt(baseAssets.mortgage.rho) .* y/PeriodOfPortfolio^.5  *VolAtPeriod(PeriodOfPortfolio)   - 0.5 .* baseAssets.mortgage.rho .* VolAtPeriod(PeriodOfPortfolio)^2 ) ...
    .* normcdf( ( + log(LoanRepayment(PeriodOfPortfolio)) - (0.5  - baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio)^2 - sqrt(baseAssets.mortgage.rho)  .* y/PeriodOfPortfolio^.5 *VolAtPeriod(PeriodOfPortfolio) ) ./ ( sqrt(1-baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio) ) )));



PortfolioCashFlows = [];
for iter0 = 1:size(MedianCumulativeReturns,1)
    for iter1 = 1:Periods
        PortfolioCashFlows(iter0,iter1) =  BankPortfolioCashFlow(MedianCumulativeReturns(iter0), iter1)*Weights(iter1);
    end
end

LiquidityDefault = PortfolioCashFlows<OneA*Repayments;

TransitionPDF = eps+diff(normcdf([-Inf ;MedianCumulativeReturns(1:end-1)/2+MedianCumulativeReturns(2:end)/2 ;Inf]));

TransitionMatrix = fliplr(flipud(spdiags(fliplr((ones(nnn,1)*TransitionPDF')'+eps))));

TransitionMatrix(:,1) = sum(TransitionMatrix(:,1:1+(nnn-1)/2),2);
TransitionMatrix(:,2:1+(nnn-1)/2) = [];

TransitionMatrix(:,end) = sum(TransitionMatrix(:,end-(nnn-1)/2:end),2);
TransitionMatrix(:,end-(nnn-1)/2:end-1) = [];

StatProb = TransitionPDF;
for iter1 = 2:Periods
    StatProb(:,iter1) = StatProb(:,iter1-1)'*TransitionMatrix;
end

InDefault = LiquidityDefault+eps;
for iter1 = 2:Periods
    InDefault(:,iter1) = max(LiquidityDefault(:,iter1)', InDefault(:,iter1-1)'*TransitionMatrix);
end

CashFlowToDebt = (1-InDefault).*(OneA*Repayments)+InDefault.*PortfolioCashFlows*(1-baseEconomy.tau)*(1-baseEconomy.alphab);

DebtValue = DiscountFactor.*sum(StatProb.*CashFlowToDebt);

BankInterest = Repayments-DebtValue;

PortfolioCosts = sum(StatProb.*PortfolioCashFlows).*DiscountFactor.*exp(-BankLoanIRate*(1:Periods)*baseEconomy.T/Periods);

EquityValue = DiscountFactor*sum((1-InDefault).*StatProb.*max(0,PortfolioCashFlows-OneA*Repayments - baseEconomy.tau*max(0,PortfolioCashFlows-OneA*BankInterest-OneA*PortfolioCosts)))';

BankValue = sum(DebtValue)+EquityValue;
Leverage = sum(DebtValue)/BankValue;

[BankValue/sum(sum(PortfolioCashFlows.*StatProb))*100 Leverage]




%% Basic Setup - Symetrical

Periods = 1;
Weights = ones(1,Periods)/Periods*baseAssetWeights(2);
Repayments = .55190;
nnn = 1001;
BankLoanIRate = 0.0021;


MedianCumulativeReturns = linspace(-10,10,nnn)';

DiscountFactor = exp(-baseEconomy.rfPerYear*baseBankH.T*(1:Periods)/Periods);

OneA = ones(nnn,1);

LoanRepayment = @(PeriodOfPortfolio) .80*DiscountFactor(PeriodOfPortfolio);
VolAtPeriod = @(PeriodOfPortfolio) baseAssets.mortgage.sigma* sqrt(baseBankH.T*PeriodOfPortfolio/Periods);

BankPortfolioCashFlow = @(y,PeriodOfPortfolio) (baseAssets.mortgage.DefaultProtected* LoanRepayment(PeriodOfPortfolio) + ...
    (1- baseAssets.mortgage.DefaultProtected) *(...
    LoanRepayment(PeriodOfPortfolio) .* normcdf( ( - log( LoanRepayment(PeriodOfPortfolio))  - 0.5 .* VolAtPeriod(PeriodOfPortfolio)^2 + sqrt(baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio) .* y/PeriodOfPortfolio^.5  ) ./ ( sqrt(1-baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio)) ) ...
    +  (1-baseAssets.mortgage.alpha).*(1-baseAssets.mortgage.tau).*exp( sqrt(baseAssets.mortgage.rho) .* y/PeriodOfPortfolio^.5  *VolAtPeriod(PeriodOfPortfolio)   - 0.5 .* baseAssets.mortgage.rho .* VolAtPeriod(PeriodOfPortfolio)^2 ) ...
    .* normcdf( ( + log(LoanRepayment(PeriodOfPortfolio)) - (0.5  - baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio)^2 - sqrt(baseAssets.mortgage.rho)  .* y/PeriodOfPortfolio^.5 *VolAtPeriod(PeriodOfPortfolio) ) ./ ( sqrt(1-baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio) ) )));



PortfolioCashFlows = [];
for iter0 = 1:size(MedianCumulativeReturns,1)
    for iter1 = 1:Periods
        PortfolioCashFlows(iter0,iter1) =  BankPortfolioCashFlow(MedianCumulativeReturns(iter0), iter1)*Weights(iter1);
    end
end

LiquidityDefault = PortfolioCashFlows<OneA*Repayments;

TransitionPDF = diff([0 ;normcdf(MedianCumulativeReturns(1:end-1)/2+MedianCumulativeReturns(2:end)/2 ); 1]);

TransitionMatrix = fliplr(flipud(spdiags(fliplr((ones(nnn,1)*TransitionPDF')'+eps))));

TransitionMatrix(:,1) = sum(TransitionMatrix(:,1:1+(nnn-1)/2),2);
TransitionMatrix(:,2:1+(nnn-1)/2) = [];

TransitionMatrix(:,end) = sum(TransitionMatrix(:,end-(nnn-1)/2:end),2);
TransitionMatrix(:,end-(nnn-1)/2:end-1) = [];

StatProb = TransitionPDF;
for iter1 = 2:Periods
    StatProb(:,iter1) = StatProb(:,iter1-1)'*TransitionMatrix;
end

InDefault = LiquidityDefault+eps;
for iter1 = 2:Periods
    InDefault(:,iter1) = max(LiquidityDefault(:,iter1)', InDefault(:,iter1-1)'*TransitionMatrix);
end

CashFlowToDebt = (1-InDefault).*(OneA*Repayments)+InDefault.*PortfolioCashFlows*(1-baseEconomy.tau)*(1-baseEconomy.alphab);

DebtValue = DiscountFactor.*sum(StatProb.*CashFlowToDebt);

BankInterest = Repayments-DebtValue;

PortfolioCosts = sum(StatProb.*PortfolioCashFlows).*DiscountFactor.*exp(-BankLoanIRate*(1:Periods)*baseEconomy.T/Periods);

EquityValue = DiscountFactor*sum((1-InDefault).*StatProb.*max(0,PortfolioCashFlows-OneA*Repayments - baseEconomy.tau*(PortfolioCashFlows-OneA*BankInterest-OneA*PortfolioCosts)))';

BankValue = sum(DebtValue)+EquityValue;
Leverage = sum(DebtValue)/BankValue;

[BankValue/sum(sum(PortfolioCashFlows.*StatProb))*100 Leverage]



%% Liquidity Default - N Maturities at T/N, 2T/N, ... T - symetrical

Periods = 1;
Weights = ones(1,Periods)/Periods*baseAssetWeights(2);
Repayments = 0.48;[.182 .185 .185];
nnn = 1001;


MedianCumulativeReturns = linspace(-10,10,nnn)';

DiscountFactor = exp(-baseEconomy.rfPerYear*baseBankH.T*(1:Periods)/Periods);

OneA = ones(nnn,1);

LoanRepayment = @(PeriodOfPortfolio) .80*DiscountFactor(PeriodOfPortfolio);
VolAtPeriod = @(PeriodOfPortfolio) baseAssets.mortgage.sigma* sqrt(baseBankH.T*PeriodOfPortfolio/Periods);

BankPortfolioCashFlow = @(y,PeriodOfPortfolio) (baseAssets.mortgage.DefaultProtected* LoanRepayment(PeriodOfPortfolio) + ...
    (1- baseAssets.mortgage.DefaultProtected) *(...
    LoanRepayment(PeriodOfPortfolio) .* normcdf( ( - log( LoanRepayment(PeriodOfPortfolio))  - 0.5 .* VolAtPeriod(PeriodOfPortfolio)^2 + sqrt(baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio) .* y/PeriodOfPortfolio^.5  ) ./ ( sqrt(1-baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio)) ) ...
    +  (1-baseAssets.mortgage.alpha).*(1-baseAssets.mortgage.tau).*exp( sqrt(baseAssets.mortgage.rho) .* y/PeriodOfPortfolio^.5  *VolAtPeriod(PeriodOfPortfolio)   - 0.5 .* baseAssets.mortgage.rho .* VolAtPeriod(PeriodOfPortfolio)^2 ) ...
    .* normcdf( ( + log(LoanRepayment(PeriodOfPortfolio)) - (0.5  - baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio)^2 - sqrt(baseAssets.mortgage.rho)  .* y/PeriodOfPortfolio^.5 *VolAtPeriod(PeriodOfPortfolio) ) ./ ( sqrt(1-baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio) ) )));

PortfolioCashFlows = [];
for iter0 = 1:size(MedianCumulativeReturns,1)
    for iter1 = 1:Periods
        PortfolioCashFlows(iter0,iter1) =  BankPortfolioCashFlow(MedianCumulativeReturns(iter0), iter1)*Weights(iter1);
    end
end

LiquidityDefault = PortfolioCashFlows-OneA*Repayments - baseEconomy.tau*(PortfolioCashFlows-OneA*taxShieldOLD)<0;

TransitionPDF = diff([0 ;normcdf(MedianCumulativeReturns(1:end-1)/2+MedianCumulativeReturns(2:end)/2 ); 1]);

TransitionMatrix = fliplr(flipud(spdiags(fliplr((ones(nnn,1)*TransitionPDF')'+eps))));
TransitionMatrix = TransitionMatrix -eps;

TransitionMatrix(:,1) = sum(TransitionMatrix(:,1:1+(nnn-1)/2),2);
TransitionMatrix(:,2:1+(nnn-1)/2) = [];

TransitionMatrix(:,end) = sum(TransitionMatrix(:,end-(nnn-1)/2:end),2);
TransitionMatrix(:,end-(nnn-1)/2:end-1) = [];

StatProb = TransitionPDF;
for iter1 = 2:Periods
    StatProb(:,iter1) = StatProb(:,iter1-1)'*TransitionMatrix;
end

InDefault = LiquidityDefault+eps;
for iter1 = 2:Periods
    InDefault(:,iter1) = max(LiquidityDefault(:,iter1)', InDefault(:,iter1-1)'*TransitionMatrix);
end

CashFlowToDebt = (1-InDefault).*(OneA*Repayments)+InDefault.*PortfolioCashFlows*(1-baseEconomy.tau)*(1-baseEconomy.alphab);

DebtValue = DiscountFactor.*sum(StatProb.*CashFlowToDebt);

BankInterest = Repayments-DebtValue;

EquityValue = DiscountFactor*sum((1-InDefault).*StatProb.*max(0,PortfolioCashFlows-OneA*Repayments - baseEconomy.tau*(PortfolioCashFlows-OneA*taxShieldOLD)))';

BankValue = sum(DebtValue)+EquityValue;
Leverage = sum(DebtValue)/BankValue;

[BankValue Leverage]

taxShieldOLD =  ((BankValue+BankInterest))*ones(1,Periods)/Periods; %FIXXX %baseEconomy.tau*(OneA*BankInterest+OneA*PortfolioCosts);


%% Solvency Default - N Maturities at T/N, 2T/N, ... T - symetrical

% abc = [];

for Periods = 1:1;
    baseBankH.T = 10*Periods/mean(1:Periods);
    
    Weights = DiscountFactor.^-(1:Periods)/Periods;
    nnn = 501;
    LiquidityDefault = false;
    
    firstTry = [Weights*.72 0.14 ];
    
    DiscountFactor = exp(-baseEconomy.rfPerYear*baseBankH.T/Periods);
    
    MedianCumulativeReturns = linspace(-6*Periods^.5,6*Periods^.5,nnn)';
    
    TransitionPDF = diff([0 ;normcdf(MedianCumulativeReturns(1:end-1)/2+MedianCumulativeReturns(2:end)/2 ); 1]);
    
    TransitionMatrix = fliplr(flipud(spdiags(fliplr((ones(nnn,1)*TransitionPDF')'+eps))));
    TransitionMatrix = max(0,TransitionMatrix -eps);
    
    TransitionMatrix(:,1) = sum(TransitionMatrix(:,1:1+(nnn-1)/2),2);
    TransitionMatrix(:,2:1+(nnn-1)/2) = [];
    
    TransitionMatrix(:,end) = sum(TransitionMatrix(:,end-(nnn-1)/2:end),2);
    TransitionMatrix(:,end-(nnn-1)/2:end-1) = [];
    
    LoanRepayment = @(PeriodOfPortfolio) .80*DiscountFactor^PeriodOfPortfolio;
    VolAtPeriod = @(PeriodOfPortfolio) baseAssets.mortgage.sigma* sqrt(baseBankH.T*PeriodOfPortfolio/Periods);
    
    BankPortfolioCashFlow = @(y,PeriodOfPortfolio) (baseAssets.mortgage.DefaultProtected* LoanRepayment(PeriodOfPortfolio) + ...
        (1- baseAssets.mortgage.DefaultProtected) *(...
        LoanRepayment(PeriodOfPortfolio) .* normcdf( ( - log( LoanRepayment(PeriodOfPortfolio))  - 0.5 .* VolAtPeriod(PeriodOfPortfolio)^2 + sqrt(baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio) .* y/PeriodOfPortfolio^.5  ) ./ ( sqrt(1-baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio)) ) ...
        +  (1-baseAssets.mortgage.alpha).*(1-baseAssets.mortgage.tau).*exp( sqrt(baseAssets.mortgage.rho) .* y/PeriodOfPortfolio^.5  *VolAtPeriod(PeriodOfPortfolio)   - 0.5 .* baseAssets.mortgage.rho .* VolAtPeriod(PeriodOfPortfolio)^2 ) ...
        .* normcdf( ( + log(LoanRepayment(PeriodOfPortfolio)) - (0.5  - baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio)^2 - sqrt(baseAssets.mortgage.rho)  .* y/PeriodOfPortfolio^.5 *VolAtPeriod(PeriodOfPortfolio) ) ./ ( sqrt(1-baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio) ) )));
    
    PortfolioCashFlows = [];
    for iter0 = 1:size(MedianCumulativeReturns,1)
        for iter1 = 1:Periods
            PortfolioCashFlows(iter0,iter1) =  BankPortfolioCashFlow(MedianCumulativeReturns(iter0), iter1)*Weights(iter1);
        end
    end
    
    x = fmincon(@(x) -geske_HelperFunction(baseEconomy,PortfolioCashFlows,x(1:end-1),x(end),Periods,1,TransitionMatrix,LiquidityDefault,false) ,firstTry,...
        [-eye(Periods+1);eye(Periods+1)],[zeros(Periods+1,1);[median(PortfolioCashFlows),1]'],[],[],[],[],@(x) geske_HelperFunction(baseEconomy,PortfolioCashFlows,x(1:end-1),x(end),Periods,1,TransitionMatrix,LiquidityDefault,true) );
    
    [~,b ] = geske_HelperFunction(baseEconomy,PortfolioCashFlows,x(1:end-1),x(end),Periods,1,TransitionMatrix,LiquidityDefault,false);
    b.Periods = Periods;
    b.nnn = nnn;
    b.AnnualDefault = 1-(1- b.BankDefProb)^(1/baseBankH.T);
    b.Years = mat2str(round(baseBankH.T*(1:Periods)/Periods,1));
    
    abc = [abc b];
    
    struct2table(abc)
end