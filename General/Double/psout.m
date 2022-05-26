function stop = psout(optimValues,state,psoname)
stop=false;
switch state
    
    case 'iter'
        
    if mod(optimValues.iteration,200)==0 %save swarm every 200 iterations
        nit=optimValues.iteration/200; %swarm number
        temp=['swarm',num2str(nit),'=optimValues.swarm;']; %extract swarm
        eval(temp);
        if nit==1 %first time, create output mat-file
        clear nit temp state optimValues
        save(psoname)
        else %later, append into existing mat-file 
        clear nit temp state optimValues
        save(psoname,'-append')
        end
    end
    case 'done'
        swarmfinal=optimValues.swarm; %save final swarm when algorithm converged
        
        if optimValues.iteration<200 %if first time writing (convergence in fewer than 200 iterations)
            clearvars -except stop psoname swarmfinal
            save(psoname)
        else %otherwise append final swarm
            clearvars -except stop psoname swarmfinal
            save(psoname,'-append')
        end
end
end