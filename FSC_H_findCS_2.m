% Financing as a Supply Chain: The Capital Structure of Banks and Firms - William Gornall and Ilya A. Strebulaev
% Supporting Code
%
% Author: William Gornall
% email: wrgornall@gmail.com
% 2012; Last revision: Feb 24 2014


function res=FSC_H_findCS2(inputSet, vvv)
%finddebtlevelsCOMBINED - Finds the optimal capital structure for banks and firms
%
% Syntax:  [output1,output2] = finddebtlevels(ModelParameters,Assets)
% Syntax:  [output1,output2] = finddebtlevels(ModelParameters,Assets,'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2,...)
%
% Inputs:
%    ModelParameters - Structure containing initial parameters
%    Assets - Structure containing bank's asset base
%
%    'PropertyName',PropertyValue - Optional parameters that change the
%       model behavior:
%
%         SkipOptimization - Outputs model results for a firm with no debt
%
%         Tol - Optimization tolerance
%
%         StartingPoint - Starting capital structure parameter values in
%         terms of bank  equity, bank debt
%
%         OutputEfficiencyCosts - Governs wether efficiency costs are
%           calculated and included in model output
%
%         ReplaceIffBetter - ModelParameters includes a
%
% Outputs:
%    res - Structure containing input parameters and capital structure
%    outputs
%

%------------- BEGIN CODE --------------

lastX = [];
val = [];

LoanPrices = [];
FEval = [];
BEval = [];
BDval = [];
BankFailurePoint = [];

DefaultOptions = struct('StartingPoint',NaN,'FixRfd',-1,'BONDS',false,'GetSP',0,'ForceFirmBorrowing',0,'EquityPenalty',0,'DebtBenefit',0,'Tol',1e-15,'OutputEfficiencyCosts',0,'ReplaceIfBetter',false,'SkipSolve',false,'OutputOnlyVal',false,'BoundsOff',false);
Options = setstructfields(DefaultOptions, vvv);


if(isfield(inputSet, 'BONDS'))
    if(inputSet.BONDS >0)
        Options.BONDS = true;
    end
end



OptimizationTol = Options.Tol;

InputParams = inputSet;
Assets = inputSet.Assets;



for itera = 1:size(Assets(:),1)
    Assets(itera).sigma =  Assets(itera).sigma *  InputParams.T^.5;
end

FloatingLeverage =  cumsum([Assets.FloatingLeverage]).*[Assets.FloatingLeverage];

if InputParams.CFirmLevU>-Inf
    Assets([Assets.FloatingLeverage]).Repayment = min(InputParams.CFirmLevU*.8,Assets([Assets.FloatingLeverage]).Repayment);
end

BankRepayment = sum([Assets.weight].*[Assets.Repayment])*.5;

BankFailurePoint = NaN;

fastAccessParameterTable = [[Assets.sigma]; [Assets.rho]; [Assets.tau];[Assets.weight];[Assets.alpha];[Assets.DefaultProtected];[Assets.Repayment];[Assets.maxload]];


if isnan(Options.StartingPoint)
    TrialVector = [BankRepayment fastAccessParameterTable(end-1,[Assets.FloatingLeverage])  0 ];
else
    TrialVector = Options.StartingPoint;
end

lb = min(1e-10,abs(TrialVector));
lb(end) = -1;
ub = 1+0*TrialVector;
if numel(FloatingLeverage)>0
    ub(2:end-1) = fastAccessParameterTable(end,fastAccessParameterTable(end,:)>-Inf);
end
ub(1) = sum(max(fastAccessParameterTable(end-1:end,:)).*[Assets.weight]);

if Options.BoundsOff
    ub = max(ub,1000);
end

spread = 0;

if InputParams.CFirmLevL ~= -Inf
    ub(2) = Inf;
end
if InputParams.CBankLevL ~= -Inf
    ub(1) = Inf;
end


%if ~(InputParams.ForceSpread>-Inf)
if ~Options.SkipSolve
    try
        CapitalStructResults = fmincon(@(x) -objective(x),TrialVector,[],[],[],[],lb,ub,@(x) constraint(x),optimset('MaxFunEvals',5000,'TolFun',eps,'TolX',eps, 'TolCon', OptimizationTol,'Display','off','Algorithm','active-set'));
    catch
        try
            sp = fmincon(@(x) sum(max(constraint(x),0)),TrialVector,[],[],[],[],lb,ub,[],optimset('MaxFunEvals',5000,'TolFun',eps,'TolX',eps, 'TolCon', OptimizationTol,'Display','off','Algorithm','active-set'));
            CapitalStructResults = fmincon(@(x) -objective(x),sp,[],[],[],[],lb,ub,@(x) constraint(x),optimset('MaxFunEvals',5000,'TolFun',eps,'TolX',eps, 'TolCon', OptimizationTol,'Display','off','Algorithm','active-set'));
        catch
            sp = fminsearchcon(@(x) sum(max(constraint(x),0)),TrialVector,lb,ub,[],[],[],optimset('MaxFunEvals',5000,'TolFun',eps,'TolX',eps, 'TolCon', OptimizationTol));
            CapitalStructResults = fminsearchcon(@(x) -objective(x),sp,lb,ub,[],[],@(x) constraint(x),optimset('MaxFunEvals',5000,'TolFun',eps,'TolX',eps, 'TolCon', OptimizationTol));
        end
    end
else
    CapitalStructResults = fmincon(@(x) -objective([Options.StartingPoint x]),0,[],[],[],[],lb(end),ub(end),@(x) constraint([Options.StartingPoint x]),optimset('MaxFunEvals',5000,'TolFun',eps,'TolX',eps, 'TolCon', OptimizationTol,'Display','off','Algorithm','active-set'));
    CapitalStructResults = [Options.StartingPoint(1:end-1) CapitalStructResults];
end
%
% else
%     CapitalStructResults = fmincon(@(x) -objective([x InputParams.ForceSpread] ),TrialVector(1:end-1),[],[],[],[],lb(1:end-1),ub(1:end-1),@(x) constraint([x InputParams.ForceSpread] ),optimset('MaxFunEvals',5000,'TolFun',eps,'TolX',eps, 'TolCon', OptimizationTol,'Display','off','Algorithm','active-set'));
%     CapitalStructResults = [CapitalStructResults InputParams.ForceSpread];
% end


if Options.OutputOnlyVal
    res = val;
elseif ~Options.ReplaceIfBetter || ~isfield(inputSet,'val') || val > inputSet.val
    res = append(CapitalStructResults);
else
    res = inputSet;
end


    function ret = normcdffast(x)
        ret =  0.5 * erfc(-x./ sqrt(2));
    end

    function a = WT ( y )
        
        AssetValue = ones(size(Assets,2),size(y,2));
        
        for iter1 = 1:size(Assets,2)
            
            s=fastAccessParameterTable(1,iter1); %sigma
            r=fastAccessParameterTable(2,iter1); % rho
            t=fastAccessParameterTable(3,iter1); % tau
            w=fastAccessParameterTable(4,iter1); % weight
            a=fastAccessParameterTable(5,iter1); % alpha
            p=fastAccessParameterTable(6,iter1); % default protected
            d=fastAccessParameterTable(7,iter1); % default point
            q=fastAccessParameterTable(8,iter1); % repayment
            
            AssetValueTemp = q .* normcdffast( ( - log( d)  - 0.5 .* s^2 + sqrt(r) .* s .* y ) ./ ( sqrt(1-r) .* s ) ) ...
                +  (1-a).*(1-t).*exp( sqrt(r) .* s .* y  - 0.5 .* r .* s^2 ) ...
                .* normcdffast( ( + log( d) - (0.5  - r) .* s^2   - sqrt(r) .* s .* y ) ./ ( sqrt(1-r) .* s ) );
            AssetValueTemp = AssetValueTemp * (1-p) + q * p;
            AssetValue(iter1,:) = w*AssetValueTemp;
        end
        a = sum(AssetValue,1);
    end

    function a = WTBONDS ( y )
        
        AssetValue = ones(size(Assets,2),size(y,2));
        
        for iter1 = 1:size(Assets,2)
            
            
            s=fastAccessParameterTable(1,iter1); %sigma
            r=fastAccessParameterTable(2,iter1); % rho
            t=fastAccessParameterTable(3,iter1); % tau
            w=fastAccessParameterTable(4,iter1); % weight
            a=fastAccessParameterTable(5,iter1); % alpha
            p=fastAccessParameterTable(6,iter1); % default protected
            d=fastAccessParameterTable(7,iter1); % default point
            q=fastAccessParameterTable(8,iter1); % repayment
            
            if (iter1>1)||(size(Assets,2)<3)
                
                AssetValueTemp = q .* normcdffast( ( - log( d)  - 0.5 .* s^2 + sqrt(r) .* s .* y ) ./ ( sqrt(1-r) .* s ) ) ...
                    +  (1-a).*(1-t).*exp( sqrt(r) .* s .* y  - 0.5 .* r .* s^2 ) ...
                    .* normcdffast( ( + log( d) - (0.5  - r) .* s^2   - sqrt(r) .* s .* y ) ./ ( sqrt(1-r) .* s ) );
                
                AssetValue(iter1,:) = w*(AssetValueTemp * (1-p) + q * p);
                
            else
                
                q = q*(1-InputParams.BONDS);
                
                d = min(d,  q/((1-a).*(1-t)));
                
                fac = min(1/(1-InputParams.BONDS), 1/((1-a).*(1-t)) );
                
                AssetValueTemp = q .* normcdffast( ( - log( d)  - 0.5 .* s^2 + sqrt(r) .* s .* y ) ./ ( sqrt(1-r) .* s ) ) ...
                    +  (1-a).*(1-t)*fac.*exp( sqrt(r) .* s .* y  - 0.5 .* r .* s^2 ) ...
                    .* normcdffast( - log(fac) + (   + log( d) - (0.5  - r) .* s^2   - sqrt(r) .* s .* y ) ./ ( sqrt(1-r) .* s ) );
                
                AssetValue(iter1,:) = w*(AssetValueTemp * (1-p) + q * p);
            end
            
        end
        a = sum(AssetValue,1);
    end


    function yfail = CalcBankFailurePoint( offset)
        if WT(-10)*offset>BankRepayment
            yfail = -10;
        elseif WT(10)*offset < BankRepayment
            yfail = 10;
        else
            yfail = fzero(@(x)WT(x)*offset-BankRepayment,1);
        end
    end

    function a = BEpo( W  )
        
        a = (W - BankRepayment -  InputParams.tau * max(0, W - sum(LoanPrices) - (BankRepayment -BDval)  ) ).*(W>=BankRepayment) + ...
            (1-InputParams.fracown)*InputParams.b2prob*(W - BankRepayment +W*InputParams.bout2  ).*(W>= BankRepayment-InputParams.bout2).*(W<BankRepayment) ;
        
        
        
    end

    function a =  BDpo ( W  )
        
        gteeRAW = 0;
        for iter2 = 1:numel(Assets)
            AiQ = Assets(iter2);
            
            if (AiQ.AssetSpecificSpread == -Inf)
                spreadToUse = spread;
            else
                spreadToUse = AiQ.AssetSpecificSpread ;
            end
            
            gteeRAW =  gteeRAW + LoanPrices(iter2)*exp( spreadToUse  * InputParams.T);
            
        end
        
        
        
        
        
        
        gteebase = min(BankRepayment,gteeRAW);
        a =   BankRepayment +  (W<BankRepayment).*( 1 -  InputParams.b2prob* (W>BankRepayment-InputParams.bout2) ).* ((1-InputParams.bout)*max((1-InputParams.alphab).* W  , gteebase* InputParams.depi ) +InputParams.bout*gteebase - BankRepayment );
        
        
    end

    function Fval=objective(trialCS)
        
        
        
        if ~isequal(trialCS,lastX)
            lastX=trialCS;
            
            BankRepayment = trialCS(1);
            spread = trialCS(end);
            
            for iter2 = 1:numel(Assets)
                if ( FloatingLeverage(iter2) >0)
                    Assets(iter2).Repayment = trialCS(1+FloatingLeverage(iter2));
                end
                
                Assets(iter2).DefaultPoint = max(1e-1000, Assets(iter2).Repayment);
                
                if Assets(iter2).tau>0
                    AiQ = Assets(iter2);
                    
                    if (AiQ.AssetSpecificSpread == -Inf)
                        spreadToUse = spread;
                    else
                        spreadToUse = AiQ.AssetSpecificSpread ;
                    end
                    
                    for iterj = 1:5
                        AiQ.DefaultPoint = max(1e-100, AiQ.Repayment) +  AiQ.tau/(1-AiQ.tau) * exp(-spreadToUse  * InputParams.T - InputParams.rf) * ( AiQ.DefaultProtected + (1-AiQ.DefaultProtected)* (...
                            + AiQ.Repayment * (1-normcdf( ( + AiQ.sigma^2/2 + log(AiQ.DefaultPoint)   )/AiQ.sigma)) ...
                            + (1- AiQ.tau)*(1- AiQ.alpha)*normcdf( ( - AiQ.sigma^2/2 + log(AiQ.DefaultPoint)   )/AiQ.sigma)  ...
                            ));
                    end
                    Assets(iter2) = AiQ;
                end
                
            end
            
            fastAccessParameterTable(7,:) = [Assets.DefaultPoint];
            fastAccessParameterTable(8,:) = [Assets.Repayment];
            
            
            LoanPrices = ones(numel(Assets),1);
            for iter2 = 1:numel(Assets)
                
                AiQ = Assets(iter2);
                
                if (AiQ.AssetSpecificSpread == -Inf)
                    spreadToUse = spread;
                else
                    spreadToUse = AiQ.AssetSpecificSpread ;
                end
                
                % CORR WG 5/15
                LoanPrices(iter2,:) = AiQ.weight*exp(-spreadToUse  * InputParams.T - InputParams.rf) * ( AiQ.Repayment *AiQ.DefaultProtected + (1-AiQ.DefaultProtected)* (...
                    + AiQ.Repayment * (1-normcdf( ( + AiQ.sigma^2/2 + log(AiQ.DefaultPoint)   )/AiQ.sigma)) ...
                    + (1- AiQ.tau)*(1- AiQ.alpha)*normcdf( ( - AiQ.sigma^2/2 + log(AiQ.DefaultPoint)   )/AiQ.sigma)  ...
                    ));
                
                
            end
            
            
            AssetValue = zeros(numel(Assets),1);
            for iter1 = 1:numel(Assets)
                if(FloatingLeverage(iter1)>0)
                    AiQ = Assets(iter1);
                    AssetValue(iter1) = (1-AiQ.tau)*AiQ.weight*exp(-InputParams.rf)*(normcdf( (AiQ.sigma^2/2 - log( AiQ.DefaultPoint ))/AiQ.sigma) ...
                        -AiQ.DefaultPoint* normcdf( ( -AiQ.sigma^2/2 - log(AiQ.DefaultPoint ))/AiQ.sigma));
                    
                end
            end
            
            FEval = sum(AssetValue);
            
            
            
            BankFailurePoint = CalcBankFailurePoint(1);
            
            warning('off','MATLAB:integral:MinStepSize');
            
            if ~Options.BONDS
                BDval =  exp(InputParams.T*Options.DebtBenefit)*exp(-InputParams.rf)* (BankRepayment*normcdf(-BankFailurePoint-.01) + integral(  @(y)  BDpo( WT( y) ) .* normpdf(y),-8,BankFailurePoint+.01,'RelTol',OptimizationTol,'AbsTol',OptimizationTol ,'Waypoints',BankFailurePoint ));
                
                BEval =   exp(-InputParams.T*Options.EquityPenalty)* exp(-InputParams.rf)* integral(  @(y)  BEpo( WT( y) ) .* normpdf(y),-8,8,'RelTol',OptimizationTol,'AbsTol',OptimizationTol ,'Waypoints',BankFailurePoint );
                
                val =   FEval + BEval + BDval;
            else
                BDval =  exp(InputParams.T*Options.DebtBenefit)*exp(-InputParams.rf)* (integral(  @(y)  BDpo( WTBONDS( y) ) .* normpdf(y),-8,8,'RelTol',OptimizationTol,'AbsTol',OptimizationTol ,'Waypoints',BankFailurePoint ));
                
                BEval =   exp(-InputParams.T*Options.EquityPenalty)* exp(-InputParams.rf)* integral(  @(y)  BEpo( WTBONDS( y) ) .* normpdf(y),-8,8,'RelTol',OptimizationTol,'AbsTol',OptimizationTol ,'Waypoints',BankFailurePoint );
                
                BONDSval =   exp(-InputParams.rf)* integral(  @(y)  ( WT( y) -  WTBONDS( y) ) .* normpdf(y),-8,8,'RelTol',OptimizationTol,'AbsTol',OptimizationTol ,'Waypoints',BankFailurePoint );
                
                val =   FEval + BEval + BDval + BONDSval;
            end
            
            
            
        end
        
        Fval=val;
        
    end

    function [Cvali , Cvale]=constraint(trialCS)
        
        if ~isequal(trialCS,lastX)
            objective(trialCS);
        end
        
        CRBasel = [];  CRBaselS = []; CRFirmLevU = []; CRFirmLevL = []; CRBankLevU = []; CRBankLevL = [];
        
        if InputParams.CBankLevU ~= -Inf
            CRBankLevU = BDval/(BDval+BEval) - InputParams.CBankLevU;
        end
        if InputParams.CBankLevL ~= -Inf
            CRBankLevL = -( BDval/(BDval+BEval) - InputParams.CBankLevL);
        end
        
        if InputParams.CFirmLevU ~= -Inf
            CRFirmLevU = 1-FEval/(FEval+LoanPrices(1)) - InputParams.CFirmLevU;
        end
        if InputParams.CFirmLevL ~= -Inf
            CRFirmLevL = -(1- FEval/(FEval+LoanPrices(1)) - InputParams.CFirmLevL) ;
        end
        if InputParams.CBaselS ~= -Inf
            BaselReq = InputParams.CBaselS*LoanPrices.*.7;%[Assets.riskWeight]';
            CRBaselS = sum(BaselReq)-BEval;
        end
        
        
        if InputParams.ForceSpread == -Inf
            BankProfit =BEval+BDval-sum(LoanPrices);
        else
            BankProfit =BEval+BDval-sum(LoanPrices)-InputParams.ForceSpread*(numel(Assets)==3);
        end
        
        if InputParams.CBaselIRB ~= -Inf
            
            
            BaselReq = ones(numel(Assets),1);
            %%
            for iter1 = 1:numel(Assets)
                
                %%
                AiQ = Assets(iter1);
                
                if AiQ.corrFormula == '0'
                    BaselReq(iter1,:)  = 0;
                else
                    
                    LGD = 1-(1- AiQ.tau)*(1- AiQ.alpha).*max(1e-10,normcdf( ( - AiQ.sigma^2/2 + log(AiQ.DefaultPoint))/AiQ.sigma))/max(1e-10,normcdf( ( + AiQ.sigma^2/2 + log(AiQ.DefaultPoint)   )/AiQ.sigma)/AiQ.Repayment);
                    
                    PD = (1-AiQ.DefaultProtected)*normcdf( (AiQ.sigma^2/2 + log(AiQ.DefaultPoint) )/AiQ.sigma);
                    
                    if AiQ.corrFormula == 'r'
                        LGD = 0.35;
                        R = 0.15;
                        PD = PD/30;
                        
                    elseif AiQ.corrFormula == 'c'
                        R= 0.12 * (1 - exp( -50 * PD)) / (1- exp(-50) ) + .24 * (1 - (1 - exp( -50 * PD)) / (1- exp(-50) ));
                        LGD = 0.45;
                        PD = PD/InputParams.T;
                    end
                    
                    PD = max(PD,1e-4);
                    
                    b= (0.11852 - 0.05478 * log(PD))^2;
                    BaselReq(iter1,:) = LoanPrices(iter1)*LGD * (normcdf(  norminv(PD)/sqrt(1-R) + norminv(0.999) * sqrt(R)/sqrt(1-R) ) -  LGD*PD )* (1+(InputParams.T-2.5) * b) / (1-1.5* b);
                end
            end
            
            BaselReq(isnan(BaselReq)) = LoanPrices(isnan(BaselReq));
            
            CRBasel = InputParams.CBaselIRB*sum(BaselReq)-BEval;
        end
        
        JJ = [Assets.weight];
        
        xtrc = trialCS(1)-sum(trialCS(2:end).*JJ(FloatingLeverage>0))-sum((1-[Assets.FloatingLeverage]).*[Assets.Repayment].*[Assets.weight]);
        if Options.BoundsOff
            xtrc = [];
        end
        Cvali = [   xtrc CRBasel   CRBaselS  CRFirmLevU CRFirmLevL CRBankLevU  CRBankLevL ];
        Cvale=BankProfit;
        
    end


    function res = append(trialCS)
        
        if  ~isequal(trialCS,lastX) %update stored quantities if x has changed
            objective(trialCS);
        end
        
        if ~Options.OutputOnlyVal
            
            if InputParams.bout2<=0
                BDEF    =  normcdf(BankFailurePoint);
            else
                BDEF    = normcdf(CalcBankFailurePoint(1+InputParams.bout2)) * InputParams.b2prob + normcdf(BankFailurePoint)*(1-InputParams.b2prob);
            end
            
            CapitalStructResults = trialCS;
            
            BankProfit = (BEval+BDval-sum(LoanPrices))/sum(LoanPrices);
            BankMLeverage = BDval/(BEval+BDval);
            BankBankLeverage = BDval/sum(LoanPrices);
            FirmLeverage = (LoanPrices(1)+(val - FEval - BEval - BDval))/(LoanPrices(1)+FEval+(val - FEval - BEval - BDval));
            
            FDEFprob = normcdf( ( + Assets(1).sigma^2/2 + log(Assets(1).DefaultPoint)   )/Assets(1).sigma);
            
            [aee bee] = constraint(CapitalStructResults);
            if max([aee abs(bee)])>OptimizationTol*10
                val = 0;
            end
            
            res = setstructfields(InputParams,struct('StartingPoint',CapitalStructResults,'CSVals',[BDval,BEval,FEval],'BankProfit',BankProfit,'LoanPrices',LoanPrices,'val',val,'BankLev',BankMLeverage,'BankLevBook',BankBankLeverage,'FirmLev',FirmLeverage,'AnnBDEF',(1-(1-BDEF)^(1/InputParams.T)),'AnnFDEF',(1-(1-FDEFprob)^(1/InputParams.T))));
            
        end
    end


end