% Financing as a Supply Chain: The Capital Structure of Banks and Firms - William Gornall and Ilya A. Strebulaev
% Supporting Code
%
% Author: William Gornall
% email: wrgornall@gmail.com
% 2012; Last revision: Feb 28 2014


%% Model Setup

savpath = 'imagesnew';

baseEconomy = struct('rf', 0.025*5, 'tau', 0.25, 'alphab', .2, 'T',5,'ForceSpread',-Inf,...
    'bout',0,'bout2',0,'depi',0,'fracown',.5,'b2prob',.5,...
    'CBaselIRB',-Inf,'CBaselS',-Inf,  'CFirmLevU',-Inf, 'CFirmLevL',-Inf, 'CBankLevU',-Inf,  'CBankLevL' ,-Inf, 'Assets', []  );

baseAssetWeights = [    1.008512262742020   1.041302859831466   0.221789624087612  ];

baseAssets = struct();
baseAssets.corporateLoans  = struct('weight',  baseAssetWeights(1),'AssetSpecificSpread',-Inf,'rho',.2,'sigma',.4,'tau',.25,'alpha',0.1,'DefaultProtected',0,'FloatingLeverage',true,'maxload',1,'Repayment',.3,'riskWeight',1,'corrFormula','c');
baseAssets.mortgage  = struct('weight', baseAssetWeights(2),'AssetSpecificSpread',-Inf,'rho',.2,'sigma',.25,'tau',0,'alpha',0.25,'DefaultProtected',0.5,'FloatingLeverage',false,'maxload',-Inf,'Repayment',.8*exp( -baseEconomy.rf ),'riskWeight',.5,'corrFormula','r');
baseAssets.cash  = struct('weight', baseAssetWeights(3),'AssetSpecificSpread',-Inf,'rho',.5,'sigma',.5,'tau',0,'alpha',0,'DefaultProtected',1,'FloatingLeverage',false,'maxload',-Inf,'Repayment',1,'riskWeight',0,'corrFormula','0');

baseBank = setfield(baseEconomy,  'Assets',[baseAssets.corporateLoans,  baseAssets.mortgage, baseAssets.cash ]);
baseBankC = setfield(baseEconomy,  'Assets',[baseAssets.corporateLoans]);
baseBankH = setfield(baseEconomy,  'Assets',[baseAssets.mortgage]);

baseBank =  setfield(FSC_H_findCS(baseBank),'val',0);
baseBankC = setfield(FSC_H_findCS(baseBankC),'val',0);
baseBankH = setfield(FSC_H_findCS(baseBankH),'val',0);

bankAssetMix = baseBank.LoanPrices/sum(baseBank.LoanPrices);
 % baseAssetWeights = baseAssetWeights./bankAssetMix'.*[.2 .6 .2]


outputCharts = false;

%% Figure 1: Impact of Seniority and Diversification on Distribution of Returns

if true
    
    trials = 10000;
    pointsToTry = 0.5/trials:1/trials:1;
    
    tempAsset = baseAssets.corporateLoans;
    
    tempBank = baseBankC;
    tempBank.tau=0;
    tempBank.alphab=0;
    
    tempBank = FSC_H_findCS(setfield(tempBank,'CFirmLevU',0.25));
    tempAsset.DPS  = tempBank.StartingPoint(2) + tempAsset.tau/(1-tempAsset.tau)*tempBank.LoanPrices;
    tempAsset.Repayment = tempBank.StartingPoint(2)/(1-tempAsset.tau);
    
    tempAsset.T = baseEconomy.T;
    tempAsset.rf = baseEconomy.rf;
    tempAsset.sigma = tempAsset.sigma*tempAsset.T^.5;
    
    loanPayoff = @(z) (z>=tempAsset.DPS)*tempAsset.Repayment + z.*(z<tempAsset.DPS)* (1 - tempAsset.alpha);
    
    payoffMatrix = [
        exp(tempAsset.sigma*norminv(pointsToTry)-tempAsset.sigma^2/2);
        exp(sqrt(tempAsset.rho)*tempAsset.sigma*norminv(pointsToTry)-tempAsset.rho*tempAsset.sigma^2/2  );
        loanPayoff( exp(tempAsset.sigma*norminv(pointsToTry)-tempAsset.sigma^2/2 ) );
        integral(@(x) loanPayoff( exp(sqrt(tempAsset.rho)*tempAsset.sigma*norminv(pointsToTry)+sqrt(1-tempAsset.rho)*tempAsset.sigma*norminv(x)-tempAsset.sigma^2/2 )),0,1,'ArrayValued',true,'AbsTol',0.00001);
        ];
    payoffMatrixRaw = payoffMatrix;
    
    mean(payoffMatrix');
    
    payoffMatrix(3,floor(trials:-1/tempAsset.DefaultProtected:1)) = tempAsset.Repayment;
    payoffMatrix(4,:) = payoffMatrix(4,:)*(1-tempAsset.DefaultProtected) + tempAsset.Repayment*tempAsset.DefaultProtected;
    
    payoffMatrix = payoffMatrix./(mean(payoffMatrix,2)*ones(1,size(payoffMatrix,2)));
    payoffMatrix = log(payoffMatrix);%/tempAsset.T^.5;
    
    xAxis = -2:.0005:2;
    
    plot1 = ksdensity(payoffMatrix(1,:),xAxis);
    plot2 = ksdensity(payoffMatrix(2,:),xAxis);
    plot3 = ksdensity(payoffMatrix(4,:),xAxis);
    
    
    %%
    close(figure(2));           figure(2)
    set(gcf, 'PaperPositionMode', 'auto');    set(gcf, 'Position', [0 0 800 400]);        set(gcf, 'color', [1 1 1]);
    
    plot(xAxis,plot1,'k:',xAxis,plot2,'k--',xAxis,plot3,'k','LineWidth',2);
    axis([-1.5,1.5, 0,13]);  set(gca, 'FontSize', 14);  xlabel('Five Year Log Return');         ylabel('Probability Density');
    set(gca,'XTick',-2:.5:2)
    
    legend('Single Firm','Pool of Firms','Pool of Loans','Orientation','Horizontal','Location','SouthOutside');      legend('boxoff');
    
    arrowp1 = [.268 .515];      arrowp2 = [.410 .63];    arrowp3 = [.42 .77];
    
    
    arrowd1 = .045 * [0 -2];    arrowd2 = .09 * [1.2 -2];    arrowd3 = .073 * [1.2 0];
    
    annotation('textarrow', [ arrowp1(1) arrowp1(1)+arrowd1(1)] , [ arrowp1(2) arrowp1(2)+arrowd1(2)] , 'String' , ['Single Firm, 5-year vol ' num2str(100*std(payoffMatrix(1,:)),2) '%'],'Fontsize',12);
    annotation('textarrow', [ arrowp2(1) arrowp2(1)+arrowd2(1)] , [ arrowp2(2) arrowp2(2)+arrowd2(2)] , 'String' , ['Pool of Firms, 5-year vol ' num2str(100*std(payoffMatrix(2,:)),2) '%'],'Fontsize',12);
    annotation('textarrow', [ arrowp3(1) arrowp3(1)+arrowd3(1)] , [ arrowp3(2) arrowp3(2)+arrowd3(2)] , 'String' , ['Pool of Loans, 5-year vol ' num2str(100*std(payoffMatrix(4,:)),2) '%'],'Fontsize',12);
   
    if outputCharts
        export_fig([savpath '\figure_cdfplots_paper.pdf'])
    end
end


%%

% Bank with 1\% default rate.
inp = baseBankC;
inp.CBankLevU= .9+.0001;
inp.CBankLevL=.8925;
inp.CFirmLevL;
inp = FSC_H_processPoints(inp);        
inp.AnnBDEF


%% Bank with lognormally distributed cash flow
inp2 = baseBankC;
inp2.CFirmLevL=1;
inp2.Assets.rho=.99999;
inp2.Assets.sigma=.057/5^.5;
inp2.CBankLevL=.8925;
inp2 = FSC_H_findCS(inp2);
inp2.AnnBDEF


%% Simple Leverage Charts
% Figure 2: Optimal Firm Leverage for Given Bank Leverage
% Figure 3: Optimal Bank Leverage for Given Firm Leverage

if true
    for runtype = 'bf'
        %%
        rn = ['SimpleLeverageChartsRunOn' runtype];
        ri.(rn).time = now();
        ri.(rn).date = date();
        
        nn = 51;
        
        pvrange =  [0.01/nn (1:nn-1)/nn 1-0.01/nn];
        
        rr.(rn) = struct( baseBankC );
        
        if runtype == 'b'
            ri.(rn).var = 'BankLev';
            rr.(rn) = [arrayfun( @(x)  setfield(baseBankC,'CBankLevU',x), pvrange(baseBankC.BankLev>pvrange)');  arrayfun( @(x)  setfield(baseBankC,'CBankLevL',x), pvrange(baseBankC.BankLev<pvrange)')] ;
        else
            ri.(rn).var = 'FirmLev';
            rr.(rn) =  [arrayfun( @(x)  setfield( setfield(baseBankC,'CFirmLevL',x-.05),'CFirmLevU',x), pvrange(baseBankC.FirmLev>pvrange)' );  arrayfun( @(x)  setfield(baseBankC,'CFirmLevL',x), pvrange(baseBankC.FirmLev<pvrange)' )] ;
        end
        
        ri.(rn).name = [ri.(rn).var '_basic_paper.pdf'];
        rr.(rn) = FSC_H_processPoints(rr.(rn),[1 11]);
        
        if runtype == 'b'
            rr.(rn)(1).FirmLev = 0;  rr.(rn)(1).BankLev = 0;
            rr.(rn)(end).FirmLev = 0; rr.(rn)(end).BankLev = 1;
        else
            rr.(rn)(1).FirmLev = 0;  rr.(rn)(1).BankLev = 1;
        end
        
        %%
        
        xAxis = [rr.(rn)(:).( ri.(rn).var)];
        xAxis(end)=1;
        
        plot1 = [rr.(rn)(:).FirmLev];
        plot2 = [rr.(rn)(:).BankLev];
        
        if runtype == 'b'
            xAxis(1:end-1) = xAxis(1:end-1)*.998;
            plot1(end) = 0;
        end
        
        close(figure(2));           figure(2)
        set(gcf, 'PaperPositionMode', 'auto');    set(gcf, 'Position', [0 0 800 400]);        set(gcf, 'color', [1 1 1]);
        plot(xAxis,plot1,'k:',xAxis,plot2,'k','LineWidth',2);
        set(gca,'Position',[  0.1300    0.2225    0.7750    0.7025])
        
        xlabh = legend('Loan Leverage','Bank Leverage','Location','SouthOutside','Orientation','horizontal');      legend('boxoff');
        
        set(xlabh,'Position',get(xlabh,'Position') + [-.125 -.25 0.25 0.25]);
        xlabh2 = get(gca,'XLabel');
        set(xlabh2,'Position',get(xlabh2,'Position') - [0 0.125 0 ]);
         
        xlabel('Leverage'); 
         
        if runtype == 'b'
            xlabel('Bank Leverage');
            hold on 
         %    plot(baseBankC.BankLev*[1,1],[0,1],'-')
            hold off
        else
            xlabel('Loan Leverage');
            hold on 
       %      plot(baseBankC.FirmLev*[1,1],[0,1],'-')
            hold off
        end
        
        axis([0,1, 0,1]);  set(gca, 'FontSize', 14);         ylabel('Leverage');
        
        
        
        set(gca,'XTick',[0,1]);     set(gca,'YTick',[0,1]);
        export_fig([savpath '\figure_' ri.(rn).name])
        
    end
    
end





%% Figure 4: Impact of Borrower Leverage on Bank Default Rates

if true
    
    rn = 'FirmLeverageAndBankDefault';
    %    rr = rmfield(rr,rn)
    ri.(rn).time = now();
    ri.(rn).date = date();
    ri.(rn).name = ['FirmLeverageAndBankDefault_paper.pdf'];
    
    nn = 21;
    paramrange =  [0.01/nn (1:nn-1)/nn 1-0.01/nn];
    
    
    pl = [0.85 0.9 0.95 ];
    
    close(figure(2));           figure(2)
    set(gcf, 'PaperPositionMode', 'auto');    set(gcf, 'Position', [0 0 800 400]);        set(gcf, 'color', [1 1 1]);
    
    for iter1 = 1:3
        inp = setfield(setfield(baseBankC, 'CBankLevU',pl(iter1)+.001), 'CBankLevL',pl(iter1)-.001);
        rr.(rn)(iter1,:)  = [    arrayfun( @(x)  setfield( inp ,'CFirmLevU',x), paramrange )];
        rr.(rn)(iter1+3,:) = [   arrayfun( @(x)  setfield( inp ,'CFirmLevL',x), paramrange )  ];
    end
    
    rr.(rn) = FSC_H_processPoints(   rr.(rn),[]);
    
    rr.(rn)= reshape(rr.(rn),[6 numel(paramrange)]);
    
    %%
    for iter1 = 1:3
        rr.(rn)(iter1+6,:) = [ [ rr.(rn)(iter1,[ true   diff([rr.(rn)(iter1+3,:).FirmLev])<1e-6      ] )] [ rr.(rn)(iter1+3,~[ true  diff([rr.(rn)(iter1+3,:).FirmLev])<1e-6  ])] ];
    end
    
    xAxis = 0:.001:1;
    
    plot1 = interp1(paramrange,[rr.(rn)(7,:).AnnBDEF]',xAxis,'PCHIP');
    plot2 = interp1(paramrange,[rr.(rn)(8,:).AnnBDEF]',xAxis,'PCHIP');
    plotc = interp1(paramrange,[rr.(rn)(9,:).AnnBDEF]',xAxis,'PCHIP');
    
    plot1 = [0 cumsum(max(0,diff(plot1))) ];
    plot2 = [0 cumsum(max(0,diff(plot2))) ];
    plotc = [0 cumsum(max(0,diff(plotc))) ];
    
    
    %%
    
    h1=     subplot(1,2,1);
    
    plot(0,0,'w',xAxis,plotc,'k:',xAxis,plot2,'k--',xAxis,plot1,'k','LineWidth',2)
    axis([0 1 0 .4])
    set(gca, 'FontSize', 14);    xlabel('Firm Leverage');   ylabel('Annual Bank Default Probability');
    
    set(gca,'XTick',0:.2:1)
    
    
    
    %%
    rn2 = [rn '2'];
    % rr = rmfield(rr,rn2)
    for iter1 = 1:3
        
        inp = setfield(setfield(baseBankH, 'CBankLevU',pl(iter1)+0.001), 'CBankLevL',pl(iter1)-0.001);
        
        inputruns = [    arrayfun( @(x)  setfield( baseAssets.mortgage ,'Repayment',x / exp( baseEconomy.rf)) , paramrange )];
        
        rr.(rn2)(iter1,:) = [ arrayfun( @(x)  setfield(inp,'Assets',inputruns(x)), 1:numel(inputruns))];
    end
    
    rr.(rn2) = FSC_H_processPoints(   rr.(rn2),[1]);
    
    rr.(rn2)=    reshape(rr.(rn2),[3 numel(paramrange)]);
    %%
    xAxis = 0:.001:1;
    plot1 = interp1(paramrange,[rr.(rn2)(1,:).AnnBDEF]',xAxis,'PCHIP');
    plot2 = interp1(paramrange,[rr.(rn2)(2,:).AnnBDEF]',xAxis,'PCHIP');
    plotc = interp1(paramrange,[rr.(rn2)(3,:).AnnBDEF]',xAxis,'PCHIP');
    %%
    
    h2=     subplot(1,2,2);
    
    
    plot(xAxis,plotc,'k:',xAxis,plot2,'k--',xAxis,plot1,'k','LineWidth',2)
    axis([0 1 0 0.4])
    set(gca, 'FontSize', 14);    xlabel('Mortgage LTV');   ylabel('Annual Bank Default Probability');
    
    set(gca,'XTick',0:.2:1)
    
    
    
    hL = legend('  95% Bank Leverage','  90% Bank Leverage','  85% Bank Leverage','Location','North','Orientation','horizontal');
    legend('boxoff');
    
    
    newPosition = [0.4 0.05 0.2 0.05];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits);
    
    abc = get(h1,'Position');
    abc(2) =abc(2)+abc(4)/5;
    abc(4) = abc(4)*4/5;
    set(h1,'Position', abc,'Units', newUnits);
    abc = get(h2,'Position');
    abc(2) =abc(2)+abc(4)/5;
    abc(4) = abc(4)*4/5;
    set(h2,'Position', abc,'Units', newUnits);
    
    
    export_fig([savpath '\figure_altmap.pdf'])
end




%%  Main Results Charts
% Figure 5: Impact of Systematic Risk on Leverage and Default Rates
% Figure 6: Impact of Asset Volatility on Leverage and Default Rates
% Figure 7: Impact of an Insured Deposit Base Proportional to Liabilities on Leverage and Default Rates
% Figure 8: Impact of Debt Guarantee Expectations on Leverage and Default Rates
% Figure 9: Impact of Equity Injection Expectations on Leverage and Default Rates
% Figure 10: Impact of Bank Leverage Limits on Leverage and Default Rates

if true
     
    for paramset = [...
%             struct('var','rho','lb',0,'ub',1,'propername',  '\rho' )
%             struct('var','sigma','lb',0,'ub',1.5,'propername',  '\sigma'   )
%             struct('var','depi','lb',0,'ub',1,'propername', 'Ins. Dep. as Fraction of Bank Liabilities '  )
%            struct('var','bout','lb',0,'ub',1,'propername', 'Probability of Debt Guarantee'  )
%            struct('var','bout2','lb',0,'ub',0.2,'propername', 'Equity Injection Size'  )
%            struct('var','CBaselS','lb',0,'ub',0.5,'propername', 'Equity Requirement h (%)'  )
%            struct('var','depih','lb',0,'ub',1,'propername', 'Ins. Dep. as Fraction of Bank Liabilities '  )
%            struct('var','tau','lb',0,'ub',1,'propername',  '\tau'   )
%             struct('var','alpha','lb',0,'ub',1,'propername',  '\alpha_f'   )
%             struct('var','alphab','lb',0,'ub',1,'propername',  '\alpha_b'   )  
           struct('var','BONDS','lb',0,'ub',1,'propername',  'Bonds as a % of Firm Debt'   )     
        %    struct('var','ForceSpread','lb',0,'ub',.05,'propername',  'Bank Profit \delta'   )  
            ]';
        %%
        
        rn = ['LeveragePlotsRunOn' paramset.var  ];
        ri.(rn).time = now();
        ri.(rn).date = date();
        ri.(rn).var = paramset.var;
        ri.(rn).name = [paramset.var '_paper.pdf'];
        
        nprm = baseBank;
        
        nn = 41;
        paramrange =  [0.01/nn (1:nn-1)/nn 1-0.01/nn];
        
        pvrange = paramset.lb + (paramset.ub - paramset.lb) * paramrange;
        
        if strcmp( paramset.var,'depih')
            for iter1 = 1: numel(nprm.Assets)
                nprm.Assets(iter1).sigma = nprm.Assets(iter1).sigma *1.5;
                nprm.Assets(iter1).rho = nprm.Assets(iter1).rho *1.5;
            end
            ri.(rn).var =  'depi';
            pvrange = 0.805:0.01:1;
            paramset.var = 'depi';
        end
        
        if strcmp( paramset.var,'depi')
            pvrange = [0 0.8+.2*paramrange];
        end
        
        
        breakpt = -Inf;
        
        rr.(rn) = arrayfun( @(x)  setfield(setstructfields(nprm,paramset),ri.(rn).var ,x), pvrange');
        
        for iter2 = 1: numel(pvrange);
            for iter1 = 1: numel(rr.(rn)(iter2).Assets)
                rr.(rn)(iter2).Assets(iter1).( ri.(rn).var  )=rr.(rn)(iter2).( ri.(rn).var  );
            end
        end
        
         rr.(rn) = FSC_H_processPoints(rr.(rn),[1 11]);
        
        %%
        
        bankres =  rr.(rn);
        %%
        xAxis = [bankres.(ri.(rn).var)];
        axisbd = [bankres(1).lb,bankres(1).ub, 0,1];
        
        plot1 = [bankres.FirmLev];
        y1b = [bankres.BankLev];
        y2a = [bankres.AnnFDEF];
        y2b = [bankres.AnnBDEF];
        y2b(~(y2b>-Inf)) = 1;
        
        range1 =y2b>Inf;
        
        if strcmp(bankres(1).var,'depi')
            range1 =y2b>.05;
        elseif strcmp(bankres(1).var,'bout')
            range1 =y2b>.05;
        end
        
        close(figure(2));           figure(2)
        set(gcf, 'PaperPositionMode', 'auto');    set(gcf, 'Position', [0 0 800 400]);        set(gcf, 'color', [1 1 1]);
        
        subplot(1,2,1)
        
        plot(xAxis(range1),plot1(range1),'k:',xAxis(range1),y1b(range1),'k','LineWidth',2); hold on
        plot(xAxis(~range1),plot1(~range1),'k:',xAxis(~range1),y1b(~range1),'k','LineWidth',2);
        legend('Firm Leverage','Bank Leverage','Location','SouthOutside');          legend('boxoff');
        axis(axisbd);  set(gca, 'FontSize', 14);  xlabel(bankres(1).propername);         ylabel('Leverage');
        
        
        if( strcmp(bankres(1).var,'CBaselS') )
            set(gca,'XTick',[0 0.04 0.08 0.13 0.2:0.1:1])
            set(gca,'XTickLabel',[0 0.04  0.08 0.13 0.2:0.1:1]*100)
            line([.04,0.04],[0,1],'Color',[0.5,0.5,0.5])
            line([.13,0.13],[0,1],'Color',[0.5,0.5,0.5])
        end
        
        subplot(1,2,2)
        
        plot(xAxis(range1),y2a(range1),'k:',xAxis(range1),y2b(range1),'k','LineWidth',2); hold on
        plot(xAxis(~range1),y2a(~range1),'k:',xAxis(~range1),y2b(~range1),'k','LineWidth',2);
        legend('Firm Default Probability','Bank Default Probability','Location','SouthOutside');            legend('boxoff');
        axis(axisbd);  set(gca, 'FontSize', 14);  xlabel(bankres(1).propername);         ylabel('Annual Default Probability');
        
        if( strcmp(bankres(1).var,'CBaselS') )
            set(gca,'XTick',[0  0.08 0.25:0.25:1])
            set(gca,'XTickLabel',[0  0.08 0.25:0.25:1]*100)
            line([.04,0.04],[0,1],'Color',[0.5,0.5,0.5])
            line([.13,0.13],[0,1],'Color',[0.5,0.5,0.5])
        end
        
        export_fig([savpath '\figure_'  ri.(rn).name  ])
    end
end


%% Figure 5: Impact of Systematic Risk on Leverage and Default Rates

if true
    %%
    
    close(figure(2));           figure(2)
    set(gcf, 'PaperPositionMode', 'auto');    set(gcf, 'Position', [0 0 800 400]);        set(gcf, 'color', [1 1 1]);
    
    subplot(1,2,1)
    bankres =  rr.LeveragePlotsRunOnrho;
    
    xaxis = [bankres.(bankres(1).var)];
    axisbd = [bankres(1).lb,bankres(1).ub, 0,1];
    
    y1a = [bankres.FirmLev];
    y1b = [bankres.BankLev];
    
    plot(xaxis,y1a,'k:',xaxis,y1b,'k','LineWidth',2);
    legend('Firm Leverage','Bank Leverage','Location','SouthOutside');          legend('boxoff');
    axis(axisbd);  set(gca, 'FontSize', 14);  xlabel(bankres(1).propername);         ylabel('Leverage');
    
    
    
    subplot(1,2,2)
    bankres =  rr.LeveragePlotsRunOnsigma;
    
    xaxis = [bankres.(bankres(1).var)];
    axisbd = [bankres(1).lb,bankres(1).ub, 0,1];
    
    y1a = [bankres.FirmLev];
    y1b = [bankres.BankLev];
    
    plot(xaxis,y1a,'k:',xaxis,y1b,'k','LineWidth',2);
    legend('Firm Leverage','Bank Leverage','Location','SouthOutside');          legend('boxoff');
    axis(axisbd);  set(gca, 'FontSize', 14);  xlabel(bankres(1).propername);         ylabel('Leverage');
    
    export_fig([savpath '\figure_sigmarho1.pdf'])
    export_fig([savpath '\figure_sigmarho2.pdf'])
    
end



%% Figure 11: Impact of Capital Regulation on Mortgage Interest Rates

if true
    
    rn = 'CostOfRegulationRun';
    ri.(rn).time = now();
    ri.(rn).date = date();
    ri.(rn).name = [rn '_paper.pdf'];
    
    pvrange = [-Inf 0 0.1:.02:0.5];
    rr.(rn) = arrayfun( @(x)  setfield(baseBankH,'CBaselS',x), pvrange);
    rr.(rn) =    FSC_H_processPoints( rr.(rn),1)   ;
    
    spreads = [rr.(rn).LoanPrices];
    spreads = (spreads(1,:)/baseBankH.Assets(1).weight/baseBankH.Assets(1).Repayment).^(-1/ baseBank.T)-1;
    
    %%
    xAxis = pvrange;
    plot1 = 100*spreads;
    
    close(figure(2));           figure(2)
    set(gca,'XTick',[0 0.04 0.08 0.13 0.199999:0.1:0.5])
    set(gca,'XTickLabel',[0 0.04 0.08 0.13 0.2:0.1:0.5]*100)
    
    
    set(gcf, 'PaperPositionMode', 'auto');    set(gcf, 'Position', [0 0 800 400]);        set(gcf, 'color', [1 1 1]);
    plot(xAxis,plot1,'k','LineWidth',2);
    axis([0,0.5, 4,5]);  set(gca, 'FontSize', 14);  ylabel('Mortgage Interst Rate (%)');         xlabel('Equity Requirement h (%)');
    
    line([.04,0.04],[0,10],'Color',[0.5,0.5,0.5])
    line([.13,0.13],[0,10],'Color',[0.5,0.5,0.5])
    
    set(gca,'XTick',[ 0 .04 .08 .13 .2-0.0000000000000001:.1:1])
    set(gca,'XTickLabel',[ 0 .04 .08 .13 .2:.1:1]*100)
    
    export_fig([savpath '\figure_costofreg_paper.pdf'])
    
end




%% Figure 12: Bank Gambling and Deposit Insurance or Debt Guarantees

if true
    
    for gambleTrial = { 
            struct('var','depi','ub',0,'lb',.3,'mode','CBaselS','rnn','Deposit Insurance Level','rnr','Bank Leverage Limit','bank',baseBankH)
            struct('var','bout','ub',0,'lb',.3,'mode','CBaselS','rnn','Debt Guarantee Probability','rnr','Bank Leverage Limit','bank',baseBankH)
            }'
        gambleTrial = gambleTrial{1};
        
        rn = ['mHorseRace' gambleTrial.var];
        ri.(rn).time = now();
        ri.(rn).date = date();
        ri.(rn).name = ['horserace' gambleTrial.var '_paper.pdf'];
        
        rr.(rn).gambleTrial = gambleTrial;
        
        trials = 41;
        rr.(rn).range = gambleTrial.ub - (gambleTrial.ub-gambleTrial.lb)*[0:1/trials:1];
        
        zeroRhoBank = gambleTrial.bank;
        
        oneRhoBank = gambleTrial.bank;
        for ii1 = 1:numel(oneRhoBank.Assets)
            oneRhoBank.Assets(ii1).rho = 1-1e-5;
        end
        
        %%
        zeroRhoBankSP = [.2, 0.0];
        
        for ii1 = 1:numel(rr.(rn).range)
            
            tempZeroRhoBank = setfield(zeroRhoBank,gambleTrial.mode,rr.(rn).range(ii1));
            rr.(rn).baseline(ii1,1) = FSC_H_findCS( tempZeroRhoBank,'StartingPoint',zeroRhoBankSP );
            
            oneRhoBankSP = rr.(rn).baseline(ii1,1).StartingPoint;
            tempOneRhoBank = setfield(oneRhoBank,gambleTrial.mode,rr.(rn).range(ii1));
            
            %%
            trialAdjustment = -0.1;
            pastSuccess = 1-trialAdjustment;
            %%
            for ii2 = 2:100
                if abs(trialAdjustment)<1e-10
                    break;
                end
                
                tempOneRhoBank.(gambleTrial.var) = pastSuccess + trialAdjustment;
                tempZeroRhoBank.(gambleTrial.var) = pastSuccess + trialAdjustment;
                
                rr.(rn).baseline(ii1,ii2) = FSC_H_findCS( tempZeroRhoBank,'StartingPoint',zeroRhoBankSP ,'SkipBankDef',true);
                rr.(rn).gamble(ii1,ii2) = FSC_H_findCS(tempOneRhoBank,'StartingPoint',oneRhoBankSP,'SkipBankDef',true,'BoundsOff',true);
                
                if ii1>1
                    rr.(rn).gamble(ii1,ii2) = FSC_H_findCS(rr.(rn).gamble(ii1,ii2),'StartingPoint',rr.(rn).gamble(ii1-1,1).StartingPoint,'SkipBankDef',true,'ReplaceIfBetter',true,'BoundsOff',true);
                end
                
                %%
                if  [[rr.(rn).gamble(ii1,ii2).val]-rr.(rn).baseline(ii1,ii2).val]>0
                    oneRhoBankSP = rr.(rn).gamble(ii1,ii2).StartingPoint;
                    
                    pastSuccess = pastSuccess + trialAdjustment;
                    trialAdjustment = trialAdjustment*1.5;
                else
                    trialAdjustment = trialAdjustment/3.5;
                end
                
                
            end
            rr.(rn).gamble(ii1,1) = rr.(rn).gamble(ii1,ii2-1);
            rr.(rn).baseline(ii1,1) = rr.(rn).baseline(ii1,ii2-1);
            
        end
        
       
        %%
        
        plot1 =  [ rr.(rn).gamble(:,1).((rr.(rn).gambleTrial.var))]' ;
        xAxis =  rr.(rn).range(1:numel(plot1));
        
        close(figure(2));           figure(2)
        set(gcf, 'PaperPositionMode', 'auto');    set(gcf, 'Position', [0 0 800 400]);        set(gcf, 'color', [1 1 1]);
        
        rectangle('Position',[.04,0,.13-.04,1],'FaceColor',[0.75,0.75,0.75],'LineStyle','none');hold on;
        set(gca,'XTick',[0 0.04 0.08 0.13 0.199999:0.1:0.5])
        set(gca,'XTickLabel',[0 0.04 0.08 0.13 0.2:0.1:0.5]*100)
        line([.08,0.08],[0,1],'Color',[0 0 0]);
        
        h = area(xAxis,[plot1 1-plot1]);
        
        set(h(2),'FaceColor',[0.75,0.75,0.75]);
        set(h(1),'FaceColor',[1,1,1]);
        legend([rr.(rn).gambleTrial.rnn ' Above Which Bank Gambles'],'Location','SouthOutside');
        legend('boxoff');
        ybound=0;
        if strcmp(rr.(rn).gambleTrial.var,'depi')
            ybound=.7;
        end
        axis([0,.26, ybound,1]);  set(gca, 'FontSize', 14);        ylabel(rr.(rn).gambleTrial.rnn);
        
        xlabel('Equity Requirement h (%)');    
        
        set(gca,'XTick',[0 0.04 0.08 0.13 0.199999:0.05:0.5])
        set(gca,'XTickLabel',[0 0.04 0.08 0.13 0.2:0.05:0.5]*100)
        
        
        line([.04,0.04],[0,1],'Color',[0.5,0.5,0.5])
        line([.13,0.13],[0,1],'Color',[0.5,0.5,0.5])
        
        text(.13,.5+plot1(floor(12))/2,'Gambling Strategy','Fontsize',18,'HorizontalAlignment','Center')
        text(.13,ybound/2+ plot1(10)/2,'Safe Strategy','Fontsize',18,'HorizontalAlignment','Center')
        
        export_fig([savpath '\figure_' ri.(rn).name])
        
    end
end



%% Table 1: Impact of Seniority and Diversification on Return Moments

if true
    rn = 'ReturnMo';

    trials = 10000;
    pointsToTry = 0.5/trials:1/trials:1;
    
    res = zeros(10,4);
        
    tempAsset = baseAssets.corporateLoans;
    
    tempBank = baseBankC;
    tempBank.tau=0;
    tempBank.alphab=0;
    
    levs = [.25 .15 .35];
    for iter1 = 1:3
        tempBank = FSC_H_findCS(setfield(tempBank,'CFirmLevU',levs(iter1)));
        trialDefaultPoints(iter1) = tempBank.StartingPoint(2) + tempAsset.tau/(1-tempAsset.tau)*tempBank.LoanPrices;
        trialRepayments(iter1) = tempBank.StartingPoint(2)/(1-tempAsset.tau);
    end
    
    tempAsset.Repayment = trialRepayments(1);
    tempAsset.DPS = trialDefaultPoints(1);
    
    trialAssetsC = [ arrayfun( @(x)  setfield(tempAsset,'rho',x), tempAsset.rho*[1 0.5 2])  arrayfun( @(x) setfield(setfield(tempAsset,'DPS',trialDefaultPoints(x)),'Repayment',trialRepayments(x)), 2:3 )];
    
    
    tempAsset = baseAssets.mortgage;
    trialRepayments = [0.8 0.6 1]/exp(baseEconomy.rf);
    trialDefaultPoints = trialRepayments;
    
    tempAsset.Repayment = trialRepayments(1);
    tempAsset.DPS = trialDefaultPoints(1);
    
    trialAssetsH = [ arrayfun( @(x)  setfield(tempAsset,'rho',x), tempAsset.rho*[1 0.5 2])  arrayfun( @(x) setfield(setfield(tempAsset,'DPS',trialDefaultPoints(x)),'Repayment',trialRepayments(x)), 2:3 )];
    
    
    trialAssets = [trialAssetsC trialAssetsH];
    
    for iter1 = 1:10
        
        tempAsset = trialAssets(iter1);
        tempAsset.T = baseEconomy.T;
        tempAsset.rf = baseEconomy.rf;
        tempAsset.sigma = tempAsset.sigma*tempAsset.T^.5;
        
        loanPayoff = @(z) (z>=tempAsset.DPS)*tempAsset.Repayment + z.*(z<tempAsset.DPS)* (1 - tempAsset.alpha);
        
        payoffMatrix = [
            exp(tempAsset.sigma*norminv(pointsToTry)-tempAsset.sigma^2/2);
            exp(sqrt(tempAsset.rho)*tempAsset.sigma*norminv(pointsToTry)-tempAsset.rho*tempAsset.sigma^2/2  );
            loanPayoff( exp(tempAsset.sigma*norminv(pointsToTry)-tempAsset.sigma^2/2 ) );
            integral(@(x) loanPayoff( exp(sqrt(tempAsset.rho)*tempAsset.sigma*norminv(pointsToTry)+sqrt(1-tempAsset.rho)*tempAsset.sigma*norminv(x)-tempAsset.sigma^2/2 )),0,1,'ArrayValued',true,'AbsTol',0.00001)
            ];
        
        mean(payoffMatrix');
        
        payoffMatrix(3,floor(trials:-1/tempAsset.DefaultProtected:1)) = tempAsset.Repayment;
        payoffMatrix(4,:) = payoffMatrix(4,:)*(1-tempAsset.DefaultProtected) + tempAsset.Repayment*tempAsset.DefaultProtected;
        
        payoffMatrix = payoffMatrix./(mean(payoffMatrix,2)*ones(1,size(payoffMatrix,2)));
        payoffMatrix = log(payoffMatrix)/tempAsset.T^.5;
        
        res(iter1,:) = std(payoffMatrix')';
    end
    
    disp(['Table: ' rn])
    res
    
    
    rr.(rn) = res;
    
end



%% Table 2:  Bank Leverage with Varying Portfolios

if true
    rn = 'BankLevAsPortfoliosVary';

    rr.(rn)= [baseBank, ...
        arrayfun( @(x)  setfield(baseBankC,'CFirmLevU',x ), [ .15 .25 ]),...
        arrayfun( @(x)  setfield(baseBankC,'CFirmLevL',x ), [ .35 .55 .75]),...
        arrayfun( @(x)  setfield(baseBankH,'Assets','Repayment',x) , [0.6 0.7 0.8 .9 1.0]*exp( -baseEconomy.rf )  ) ];
    
    rr.(rn) = [rr.(rn), ...
        arrayfun( @(x)  setfield(rr.(rn)(x),'T',2.5), [ 1:numel(rr.(rn)) ])];
    
    rr.(rn) = FSC_H_processPoints(rr.(rn),[]);
    
    
    disp(['Table: ' rn])
    [rr.(rn).BankLev; rr.(rn).AnnBDEF]'
    
end




%%  TABLE 3: Capital Structure of Banks and Firms

if true
    rn = 'CapStructTable';
    
    rr.(rn).rows = [baseBank
        setfield(baseBank,'Assets',[ arrayfun( @(x)  setfield(baseBank.Assets(x),'rho',0.1) , [1:3] )])
        setfield(baseBank,'Assets',[ arrayfun( @(x)  setfield(baseBank.Assets(x),'rho',0.4) , [1:3] )])
        setfield(baseBank,'Assets',[ arrayfun( @(x)  setfield(baseBank.Assets(x),'sigma',0.1) , [1:3] )])
        setfield(baseBank,'Assets',[ arrayfun( @(x)  setfield(baseBank.Assets(x),'sigma',0.5) , [1:3] )])
        setfield(         setfield(baseBank,'tau',0.1),'Assets',[ arrayfun( @(x)  setfield(baseBank.Assets(x),'tau',0.1) , [1:3] )])
        setfield(         setfield(baseBank,'tau',0.35),'Assets',[ arrayfun( @(x)  setfield(baseBank.Assets(x),'tau',0.35) , [1:3] )])
        setfield(baseBank,'rf',0.01*5)
        setfield(baseBank,'rf',0.05*5)
        setfield(baseBank,'T',1)
        setfield(baseBank,'T',2.5)
        setfield(baseBank,'Assets',[ arrayfun( @(x)  setfield(baseBank.Assets(x),'alpha',0.05) , [1:3] )])
        setfield(baseBank,'Assets',[ arrayfun( @(x)  setfield(baseBank.Assets(x),'alpha',0.2) , [1:3] )])
        setfield(baseBank,'alphab',0.1)
        setfield(baseBank,'alphab',0.4)
                setfield(baseBank,'depi',0.85)
        setfield(baseBank,'depi',0.9)
        setfield(baseBank,'depi',0.95)];
    
    rr.(rn).rows2 = rr.(rn).rows;
    for iter00 = 1:numel(rr.(rn).rows)
        rr.(rn).rows2(iter00).tau = 0;
        rr.(rn).rows2(iter00).alphab = 0;
        rr.(rn).rows2(iter00).Assets = rr.(rn).rows2(iter00).Assets(1);                     
    end
    
    rr.(rn).rows = FSC_H_processPoints(rr.(rn).rows,[]);
    rr.(rn).rows2 = FSC_H_processPoints(rr.(rn).rows2,[]);
    
    
    disp(['Table: ' rn])
    [rr.(rn).rows.FirmLev; rr.(rn).rows.AnnFDEF;  rr.(rn).rows.BankLev; rr.(rn).rows.AnnBDEF; rr.(rn).rows2.FirmLev; rr.(rn).rows2.AnnFDEF]'*100
    
    
end



%%  TABLE 4: Capital Structure of Banks and Firms Under Extensions

if true
    rn = 'CapStrExt';
    
    rr.(rn) = [baseBank
        setfield(baseBank,'depi',0.85)
        setfield(baseBank,'depi',0.90)
        setfield(baseBank,'depi',0.95)
        setfield(baseBank,'depi',0.98)
        setfield(baseBank,'bout',0.25)
        setfield(baseBank,'bout',0.50)
        setfield(baseBank,'bout',0.75)
        setfield(baseBank,'bout',0.95)
        setfield(baseBank,'bout2',0.05)
        setfield(baseBank,'bout2',0.10)
        setfield(baseBank,'bout2',0.20)
        setfield(baseBank,'bout2',0.40)
        ];
    
    rr.(rn) = [rr.(rn)'
        [ arrayfun( @(x)  setfield(rr.(rn)(x),'CBaselS',0.08) , 1:numel(rr.(rn)) )]
        [ arrayfun( @(x)  setfield(rr.(rn)(x),'CBaselIRB',1) , 1:numel(rr.(rn)) )]
        [ arrayfun( @(x)  setfield(rr.(rn)(x),'CBaselS',0.16) , 1:numel(rr.(rn)) )]
        [ arrayfun( @(x)  setfield(rr.(rn)(x),'CBaselIRB',2) , 1:numel(rr.(rn)) )]
        ]';
    
    rr.(rn) = FSC_H_processPoints(rr.(rn),1);
    
    
    disp(['Table: ' rn])
    [rr.(rn).BankLev; rr.(rn).AnnBDEF]'
end



%%  IN TEXT NUMBERS
if true
    
    disp(['In text numbers'])
    %Need to run depi chart, depih chart(hidden), bout chart, cost of reg, gamble charts,
    
    defaultUnderReg = reshape([rr.CapStrExt.AnnBDEF],size(rr.CapStrExt));
    tableddefaultUnderReg =  max(defaultUnderReg( 1==[0 1 1 1 0 1 1 1 0 1 1 1 0],:));
    
    itn = [];
    
    itn{1} = baseBank.BankLev;
    itn{2} = baseBank.FirmLev;
    
    itn{3} = rr.('ReturnMo')(1,4);
    
    safePoints =  [rr.LeveragePlotsRunOndepi.AnnBDEF]<0.05;
    unsafePoints =  [rr.LeveragePlotsRunOndepi.AnnBDEF]>=0.05;
    valuesForPoints = [rr.LeveragePlotsRunOndepi.val];
    valuesForPoints2 = [rr.LeveragePlotsRunOndepi.depi];
    (valuesForPoints(cumsum(unsafePoints)==1) - max(valuesForPoints(safePoints)))/diff(valuesForPoints(abs(cumsum(unsafePoints)-1.5)==0.5));
    itn{4} = valuesForPoints2(cumsum(unsafePoints)==1) - (valuesForPoints(cumsum(unsafePoints)==1) - max(valuesForPoints(safePoints)))/diff(valuesForPoints(abs(cumsum(unsafePoints)-1.5)==0.5))*diff(valuesForPoints2(abs(cumsum(unsafePoints)-1.5)==0.5));
    
    
    safePoints =  [rr.LeveragePlotsRunOnbout.AnnBDEF]<0.05;
    unsafePoints =  [rr.LeveragePlotsRunOnbout.AnnBDEF]>=0.05;
    valuesForPoints = [rr.LeveragePlotsRunOnbout.val];
    valuesForPoints2 = [rr.LeveragePlotsRunOnbout.bout];
    (valuesForPoints(cumsum(unsafePoints)==1) - max(valuesForPoints(safePoints)))/diff(valuesForPoints(abs(cumsum(unsafePoints)-1.5)==0.5));
    itn{5} = valuesForPoints2(cumsum(unsafePoints)==1) - (valuesForPoints(cumsum(unsafePoints)==1) - max(valuesForPoints(safePoints)))/diff(valuesForPoints(abs(cumsum(unsafePoints)-1.5)==0.5))*diff(valuesForPoints2(abs(cumsum(unsafePoints)-1.5)==0.5));
    
    itn{6} = max(1-tableddefaultUnderReg(4:5)./tableddefaultUnderReg(2:3) );
    
    spreads = [rr.CostOfRegulationRun.LoanPrices];
    spreads = (spreads(1,:)/baseBankH.Assets(1).weight/baseBankH.Assets(1).Repayment).^(-1/ baseBank.T)-1;
    spreads = diff(spreads)./diff([rr.CostOfRegulationRun.CBaselS])*100;
    itn{7} = spreads(end);
    
    itn{8} = interp1( [rr.mHorseRacedepi.gamble(:,1).CBaselS],[rr.mHorseRacedepi.gamble(:,1).depi],0.08);
    
    itn{9} = 'See Table 1';
    
    
    safePoints =  [rr.LeveragePlotsRunOndepih.AnnBDEF]<0.05;
    unsafePoints =  [rr.LeveragePlotsRunOndepih.AnnBDEF]>=0.05;
    valuesForPoints = [rr.LeveragePlotsRunOndepih.val];
    valuesForPoints2 = [rr.LeveragePlotsRunOndepih.depi];
    (valuesForPoints(cumsum(unsafePoints)==1) - max(valuesForPoints(safePoints)))/diff(valuesForPoints(abs(cumsum(unsafePoints)-1.5)==0.5));
    itn{10} = valuesForPoints2(cumsum(unsafePoints)==1) - (valuesForPoints(cumsum(unsafePoints)==1) - max(valuesForPoints(safePoints)))/diff(valuesForPoints(abs(cumsum(unsafePoints)-1.5)==0.5))*diff(valuesForPoints2(abs(cumsum(unsafePoints)-1.5)==0.5));
    
    itn{11} =  rr.LeveragePlotsRunOnbout2(end).BankLev;
    
    itn{12} = interp1([rr.FirmLeverageAndBankDefault(7,:).FirmLev]', [rr.FirmLeverageAndBankDefault(7,:).AnnBDEF]',0.30);
    itn{13} = interp1([rr.FirmLeverageAndBankDefault(7,:).FirmLev]', [rr.FirmLeverageAndBankDefault(7,:).AnnBDEF]',0.60);
    itn{14} = itn{13}/itn{12};
    
    jjj = [rr.FirmLeverageAndBankDefault2(1,:).Assets];
    jjj = [jjj.Repayment]'*exp(baseEconomy.rf);
    
    itn{15} =  interp1(jjj, [rr.FirmLeverageAndBankDefault2(1,:).AnnBDEF]',0.8);
    itn{16} =  interp1(jjj, [rr.FirmLeverageAndBankDefault2(2,:).AnnBDEF]',1,'linear','extrap');
    
    
    jjj = baseBank;
    jjj.Assets(1).sigma =  jjj.Assets(1).sigma *1.5;
    jjj.Assets(2).sigma =  jjj.Assets(2).sigma *1.5;
    jjj = FSC_H_findCS(jjj,'SkipSolve',true,'StartingPoint',baseBank.StartingPoint);
    
    itn{17} = jjj.AnnBDEF;
    itn{18} = baseBank.AnnBDEF;
    
    jjj = baseBank;
    jjj.Assets(1).rho =  0.4;
    jjj.Assets(2).rho =  0.4;
    jjj = FSC_H_findCS(jjj,'SkipSolve',true,'StartingPoint',baseBank.StartingPoint);
    
    itn{19} = jjj.AnnBDEF;
    
    jjj = baseBank;
    
    jjj.Assets(1).rho =  0.4;
    jjj.Assets(2).rho =  0.4;
    jjj.CBankLevL =  jjj.BankLev;
    jjj.CFirmLevL =  jjj.FirmLev;
    
    jjj = FSC_H_findCS(jjj);
    
    
    
    itn{20} = tableddefaultUnderReg(3);
    itn{21} = tableddefaultUnderReg(5);
    
    itn{22} = rr.mHorseRacedepi.gamble(1).depi;
    itn{23} = rr.mHorseRacebout.gamble(1).bout;
    
    itn{24} = interp1([[rr.mHorseRacedepi.gamble(:,1).CBaselS]'], [rr.mHorseRacedepi.gamble(:,1).depi]',0.08);
    itn{25} = interp1([[rr.mHorseRacebout.gamble(:,1).CBaselS]'], [rr.mHorseRacebout.gamble(:,1).bout]',0.08);
    
    
    hhh = [rr.CostOfRegulationRun.StartingPoint];
   hh2 = [rr.CostOfRegulationRun.CBaselS];
   
    
    itn{26} = interp1( hh2(2:end), hhh(4:2:end),0.16)-interp1( hh2(2:end), hhh(4:2:end),0.08);
    
    itn{27} = FSC_H_findCS(setfield(baseBank,'tau',0.3));
   
    
    itn{28} = (baseBank.CSVals(1)*exp(0.025*5)-baseBank.StartingPoint(1)*(1-baseBank.AnnBDEF)^5)/(1-(1-baseBank.AnnBDEF)^5)/sum(baseBank.CSVals(1:2));
    
    itn{29} = baseBank.StartingPoint(end);
    itn{30} = interp1([[rr.LeveragePlotsRunOnForceSpread.ForceSpread]'], [rr.LeveragePlotsRunOnForceSpread.BankLev]',0.05);
    
    itn{31} = [rr.('SimpleLeverageChartsRunOnf')(end).BankLev];
    itn{32} = rr.LeveragePlotsRunOnForceSpread(end).BankLev;
    
    itn{33} = interp1([rr.SimpleLeverageChartsRunOnf.BankLev],[rr.SimpleLeverageChartsRunOnf.FirmLev],.5542);
    
    for i=1:numel(itn)
        disp([num2str(i) ':']);
        disp(itn{i});
    end
    
end

%% APPENDIX

%% Alternative Debt Benefits.

rn = 'AppAltDB';

tempBank = baseBankH;
tempBank.tau = 0;


;
results = [   FSC_H_findCS(setfield(baseBank,'tau',0.15))
     FSC_H_findCS(setfield(baseBank,'tau',0.25))
      FSC_H_findCS(setfield(baseBank,'tau',0.35))
     FSC_H_findCS(tempBank,'DebtBenefit',0.0025)
 FSC_H_findCS(tempBank,'DebtBenefit',0.005)
 FSC_H_findCS(tempBank,'DebtBenefit',0.01)
 FSC_H_findCS(tempBank,'EquityPenalty',0.05)
 FSC_H_findCS(tempBank,'EquityPenalty',0.1)
 FSC_H_findCS(tempBank,'EquityPenalty',0.15)
 ];


    disp(['Table: ' rn])

[[results.BankLev]' [ results.AnnBDEF]']


%
% ps2 = struct();
%
% ps2.asset =
% ps2.bank = [...
%     struct('var','tau','values',[0.1,0.35],'lb',0,'ub',1,'assetlevel',false,'propername', '\tau'  ),...
%     struct('var','alphab','values',[0.05,0.2],'lb',0,'ub',1,'assetlevel',false,'propername','\alpha_B'  ),...
%     struct('var','CBankLevU','values',[0.75,0.84,0.87,0.92,.96],'lb',0,'ub',1,'assetlevel',false,'propername', 'Maximum Bank Leverage'  ),...
%     struct('var','ForceSpread','values',[0,0.05],'lb',0,'ub',0.10,'assetlevel',false,'propername','\theta'  ),...
%     struct('var','rf','values',[0.025,0.1],'lb',0,'ub',.2,'assetlevel',false,'propername', 'r_f'  )
%     struct('var','spread','values',[0.01,0.02,0.03],'lb',0,'ub',.05,'assetlevel',true,'propername',  'Bank Bargaining Power'  )
%     struct('var','alpha','values',[0.1,0.5],'lb',0,'ub',1,'assetlevel',true,'propername',  '\alpha_F' ),...
%     ];
% ps2.moralhazard = [...
%     ];
%
%     abc =  [ps2.asset ps2.bank ps2.moralhazard];
%
%     for paramset = [...
%             struct('var','rho','values',[0.1,0.4],'lb',0,'ub',1,'assetlevel',true,'propername',  '\rho'  ),...
%             struct('var','sigma','values',[0.2*ps.modelParams.T^.5,.8],'lb',0,'ub',1.5,'assetlevel',true,'propername',  '\sigma'   ),...
%             struct('var','depi','values',[1,0.95,0.9,.85,.8],'lb',0,'ub',1,'assetlevel',false,'propername', 'Ins. Dep. as Fraction of Bank Liabilities '  ),...
%             struct('var','bout','values',[0.25:0.25:0.75],'lb',0,'ub',1,'assetlevel',false,'propername', 'Probability of Debt Guarantee'  ),...
%             struct('var','bout2','values',[0.01,0.02,0.04],'lb',0,'ub',0.15,'assetlevel',false,'propername', 'Equity Injection Size'  ),...
%             struct('var','CBaselS','values',[1],'lb',0,'ub',1,'assetlevel',false,'propername', 'Equity Capital Requirement'  ),...
%             ];

 
