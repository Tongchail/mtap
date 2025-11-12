function state = checkOversat(var, cal, step, iter, time, T, state)
% state: struct('done',logical,'log',[step iter time T_C maxΔ])

    if nargin < 7 || isempty(state)
        state = struct('done',false,'log',[]);
    end
    if state.done, return, end

    if isstruct(var) && isfield(var,'H2Om') && isstruct(cal) && isfield(cal,'H2Osat')
        is_over = any(var.H2Om(:) > cal.H2Osat(:));
        if is_over
            excess = var.H2Om(:) - cal.H2Osat(:);
            maxOvershoot = max(excess(excess>0));
            T_mean_C = mean(T(:)) - 273.15; 

            fprintf('[FIRST oversat] step=%d iter=%d time=%.6g T=%.2f°C maxΔ=%.3e\n', ...
                step, iter, time, T_mean_C, maxOvershoot);

            state.log  = [step, iter, time, T_mean_C, maxOvershoot];
            state.done = true;
        end
    end
end