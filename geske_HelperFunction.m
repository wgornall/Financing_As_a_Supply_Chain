function [ output, output2 ] = Untitled( baseEconomy , PortfolioCashFlows, Repayments, taxShieldTotal,Periods,DiscountFactor,TransitionMatrix,LiquidityDefault,const)

taxShield = taxShieldTotal*ones(1,Periods)/Periods;

EOPDebtCFIfD    = (1-baseEconomy.alphab)*PortfolioCashFlows;
EOPDebtCFIfND   = ones(size(PortfolioCashFlows,1),1)*Repayments;
EOPEquityCFIfND = PortfolioCashFlows- ones(size(PortfolioCashFlows,1),1)*Repayments-max(0,baseEconomy.tau*PortfolioCashFlows-ones(size(PortfolioCashFlows,1),1)*taxShield);

EOPDebtVIfBOPD    = zeros(size(PortfolioCashFlows)+[0,1]);
EOPDebtVIfBOPND   = zeros(size(PortfolioCashFlows)+[0,1]);
EOPEquityVIfBOPND = zeros(size(PortfolioCashFlows)+[0,1]);

EOPTriggerDefault      = zeros(size(PortfolioCashFlows)+[0,1]);

EventualDefaultProb = zeros(size(PortfolioCashFlows)+[0,1]);

for iter1 = Periods:-1:1
    %%
    temp = EOPEquityCFIfND(:,iter1)  +  DiscountFactor * TransitionMatrix*EOPEquityVIfBOPND(:,iter1+1);
        
    EOPTriggerDefault(:,iter1) = 0<(0>temp)+ LiquidityDefault.*(0>EOPEquityCFIfND(:,iter1));

    try
        firstSolvent = find(EOPTriggerDefault(:,iter1)==0,1);

        EOPTriggerDefault(firstSolvent,iter1) = interp1(temp(firstSolvent-1:firstSolvent),0:1,0);
    end
    
    EventualDefaultProb(:,iter1) = max(EOPTriggerDefault(:,iter1), TransitionMatrix *EventualDefaultProb(:,iter1+1) );
    
    EOPDebtVIfBOPD(:,iter1) = EOPDebtCFIfD(:,iter1)  +  DiscountFactor * TransitionMatrix * EOPDebtVIfBOPD(:,iter1+1) ;
    
    EOPEquityVIfBOPND(:,iter1) = (1-EOPTriggerDefault(:,iter1)).*temp;
        
    EOPDebtVIfBOPND(:,iter1) = EOPTriggerDefault(:,iter1).*EOPDebtVIfBOPD(:,iter1)...
        +(1-EOPTriggerDefault(:,iter1)).*(EOPDebtCFIfND(:,iter1)  +  DiscountFactor * TransitionMatrix * EOPDebtVIfBOPND(:,iter1+1));
    
    
end

Time0EquityValue = DiscountFactor * TransitionMatrix * EOPEquityVIfBOPND(:,1);
Time0DebtValue  = DiscountFactor * TransitionMatrix * EOPDebtVIfBOPND(:,1);
Time0DefaultProb = TransitionMatrix * EventualDefaultProb(:,1);

BankEquityValue =Time0EquityValue(end/2+.5);
BankDebtValue = Time0DebtValue(end/2+.5);
BankDefProb = Time0DefaultProb(end/2+.5);

BankValue = BankEquityValue+BankDebtValue;

Leverage = BankDebtValue/BankValue;

BankInterest = sum(Repayments)-BankDebtValue;

actualTaxShieldTotal =  baseEconomy.tau*(BankValue+BankInterest);


if const

    output = taxShieldTotal-actualTaxShieldTotal ;
    output2 = 0;     
else
    output = BankValue + (taxShieldTotal-actualTaxShieldTotal);
    
    output2=[];
    output2.BankValue = BankValue;
    output2.Leverage = Leverage;
    output2.BankEquityValue= BankEquityValue;
    output2.BankDebtValue= BankDebtValue;
    output2.Repayments = Repayments;
    output2.BankDefProb = BankDefProb;

end
end

