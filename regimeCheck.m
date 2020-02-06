function [regime,ind] = regimeCheck(vetor,samp_freq,we,n_ciclos,steady_time)

% Calculando o intervalo de tempo que demoram os ciclos desejados
ciclo = 1/we;
intervalo_t = ciclo*n_ciclos;

% Quantas amostras são requisitadas nesse momento
N = samp_freq*intervalo_t;

% Calcular o rms desse numero de ciclos fatiando o vetor
fatias = floor(length(vetor)/N);

rms_fatia_ant=0;
steady = 0;

for i=0:fatias-1
    ind_ini = N*i+1;
    ind_fin = N*(i+1);
    rms_fatia_atual = rms(vetor(ind_ini:ind_fin));
    
    % Calculo da diferença percentual do rms de uma fatia de n_ciclos pra
    % outra
    ddif = abs((rms_fatia_atual - rms_fatia_ant)/rms_fatia_atual);
    
    if ddif < 0.01
        steady = steady+1;
        %disp(i)
    
    elseif ddif > 0.01 && steady > 1 % Caso tenha algum erro no meio, resete a contagem
        steady = 0;
        %disp('reset')
    end
    
    % Caso o número de periodos com regime permanente for, no total, igual
    % ao tempo pedido, quebre o laço
    if steady_time/intervalo_t == steady
        regime = 'steady';
        ind = ind_fin;
        break
    else
        regime = 'trasient';
        ind = 0;
    end
    
    
    rms_fatia_ant = rms_fatia_atual;
end



