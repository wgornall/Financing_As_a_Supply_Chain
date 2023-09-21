% Financing as a Supply Chain: The Capital Structure of Banks and Firms - William Gornall and Ilya A. Strebulaev
% Supporting Code
%
% Author: William Gornall
% email: wrgornall@gmail.com
% 2012; Last revision: Feb 24 2014


function [ outp ] = FSC_H_processPoints(inp,stren)
len3 = numel(inp);
outp=inp;

sp = floor(len3/2);
tocheck = [sp:len3 sp-1:-1:1];
checkorder = tocheck*0;
sen = 1e-6;

iter0 = 1;

h = waitbar(1); 
verb0z = @(p,x) waitbar(p,h,x);

while(iter0<=numel(tocheck))
    
    eiq = tocheck(iter0);
    fff = checkorder(iter0);
    
    
    if not(or(eiq > len3, eiq<1))
        try
            if(fff == 0)
                trial = FSC_H_findCS(outp(eiq),'BoundsOff',true);
            elseif( abs(fff) == 1)
                trial = FSC_H_findCS(outp(eiq),'StartingPoint',outp(min( max(eiq-fff,1),len3)).StartingPoint,'BoundsOff',true);
            elseif( abs(fff) == 11)
                trial = FSC_H_findCS(outp(eiq),'StartingPoint',(outp(min( max(eiq-1,1),len3)).StartingPoint + outp(min( max(eiq+1,1),len3)).StartingPoint)/2,'BoundsOff',true);
            else
                trial.val = -Inf;
            end
            
            if trial.val > outp(eiq).val + sen;
                verb0z(iter0/length(tocheck),[num2str(eiq) ' moved to ' num2str(trial.val) ' from ' num2str(outp(eiq).val) ' on round ' num2str(fff) ]);
                outp(eiq) = trial;
                tocheck = [tocheck (eiq+1+stren*0) (eiq-1+stren*0)];
                checkorder = [checkorder stren -stren];
            else
                verb0z(iter0/length(tocheck),[num2str(eiq) ' stayed the same at ' num2str(outp(eiq).val) ' on round ' num2str(fff) ]);
            end
        catch
            verb0z(iter0/length(tocheck),'error')
        end
    end
    
    iter0 = iter0+1;
end

verb0z(0,'done...')
end

