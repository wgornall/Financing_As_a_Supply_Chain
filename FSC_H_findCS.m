% Financing as a Supply Chain: The Capital Structure of Banks and Firms - William Gornall and Ilya A. Strebulaev
% Supporting Code
%
% Author: William Gornall
% email: wrgornall@gmail.com
% 2012; Last revision: Feb 24 2014


function res=FSC_H_findCS(inputSet, varargin)
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



if numel(inputSet.Assets ) == 1
    bankToTest = inputSet;
else
    splitBank = inputSet;
    splitBank.Assets = splitBank.Assets(2:end);
	splitBank.Assets(1).AssetSpecificSpread = -Inf;
    splitBank = FSC_H_findCS(splitBank, struct(varargin{:}));
    splitBank.Assets(1).AssetSpecificSpread = splitBank.StartingPoint(end);

    bankToTest = splitBank;
    bankToTest.Assets = [inputSet.Assets(1) splitBank.Assets];
end

    bankToTest.Assets(1).AssetSpecificSpread = -Inf;
    
    inputArgs = struct(varargin{:});
    
    if numel(inputSet.Assets ) < 3 
        inputArgs.StartingPoint=[];
        inputArgs = rmfield(inputArgs,'StartingPoint' );
        inputArgs.SkipSolve=[];
        inputArgs = rmfield(inputArgs,'SkipSolve' );
    end
    
    
     
    res = FSC_H_findCS_2(bankToTest, inputArgs);
end