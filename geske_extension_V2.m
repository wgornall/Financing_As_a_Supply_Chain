
baseEconomy = struct('rf', 0.025*5,'rfPerYear', 0.025, 'tau', 0.25, 'alphab', .2, 'T',5,'ForceSpread',-Inf,...
    'bout',0,'bout2',0,'depi',0,'fracown',0,'b2prob',0,...
    'CBaselIRB',-Inf,'CBaselS',-Inf,  'CFirmLevU',-Inf, 'CFirmLevL',-Inf, 'CBankLevU',-Inf,  'CBankLevL' ,-Inf, 'Assets', []  );

baseAssetWeights = [    1.008512262742020   1.041302859831466   0.221789624087612  ];

baseAssets = struct();
baseAssets.corporateLoans  = struct('weight',  baseAssetWeights(1),'AssetSpecificSpread',-Inf,'rho',.2,'sigma',.4,'tau',.25,'alpha',0.1,'DefaultProtected',0,'FloatingLeverage',true,'maxload',1,'Repayment',.3,'riskWeight',1,'corrFormula','c');
baseAssets.mortgage  = struct('weight', baseAssetWeights(2),'AssetSpecificSpread',-Inf,'rho',.2,'sigma',.25,'tau',0,'alpha',0.25,'DefaultProtected',0.5,'FloatingLeverage',false,'maxload',-Inf,'Repayment',.8*exp( -baseEconomy.rf ),'riskWeight',.5,'corrFormula','r');
baseAssets.cash  = struct('weight', baseAssetWeights(3),'AssetSpecificSpread',-Inf,'rho',.5,'sigma',.5,'tau',0,'alpha',0,'DefaultProtected',1,'FloatingLeverage',false,'maxload',-Inf,'Repayment',1,'riskWeight',0,'corrFormula','0');

baseBankH = setfield(baseEconomy,  'Assets',[baseAssets.mortgage]);

baseBankH = setfield(FSC_H_findCS(baseBankH),'val',0)


 abc = [];

%% Maturities at T/N, 2T/N, ... T - symetrical


for Periods = 1:9;
    baseBankH.T = 5*Periods/mean(1:Periods);
    
    Weights = DiscountFactor.^-(1:Periods)/Periods;
    nnn = 4001;
    LiquidityDefault = false;
    
    firstTry = [Weights*.73 0.14 ];
    
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
        .* normcdf( ( + log(LoanRepayment(PeriodOfPortfolio)) - (0.5  - baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio)^2 - sqrt(baseAssets.mortgage.rho)  .* y/PeriodOfPortfolio^.5 *VolAtPeriod(PeriodOfPortfolio) ) ./ ( sqrt(1-baseAssets.mortgage.rho) .* VolAtPeriod(PeriodOfPortfolio) ) ))...
        )/DiscountFactor^PeriodOfPortfolio;
    
   
    
    PortfolioCashFlows = [];
    for iter0 = 1:size(MedianCumulativeReturns,1)
        for iter1 = 1:Periods
            PortfolioCashFlows(iter0,iter1) =  BankPortfolioCashFlow(MedianCumulativeReturns(iter0), iter1)*Weights(iter1);
        end
    end
    
    x = fmincon(@(x) -geske_HelperFunction(baseEconomy,PortfolioCashFlows,x(1:end-1),x(end),Periods,DiscountFactor,TransitionMatrix,LiquidityDefault,false) ,firstTry,...
        [-eye(Periods+1);eye(Periods+1)],[zeros(Periods+1,1);[median(PortfolioCashFlows),1]'],[],[],[],[],@(x) geske_HelperFunction(baseEconomy,PortfolioCashFlows,x(1:end-1),x(end),Periods,DiscountFactor,TransitionMatrix,LiquidityDefault,true) );
    
    [~,b ] = geske_HelperFunction(baseEconomy,PortfolioCashFlows,x(1:end-1),x(end),Periods,DiscountFactor,TransitionMatrix,LiquidityDefault,false);
    b.Periods = Periods;
    b.nnn = nnn;
    b.AnnualDefault = 1-(1- b.BankDefProb)^(1/baseBankH.T);
    b.Years = mat2str(round(baseBankH.T*(1:Periods)/Periods,1));
    
    abc = [abc b];
    
    struct2table(abc)
end