function overall_exit_status = check_termination(soc_pct,v_cell,params)

% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

exit_status = 0; % No abnormal condition has been reached. Simulation valid.

if soc_pct < params.CutoffSOC
    fprintf('Cell SOC is below lower cutoff\n');
    exit_status = [exit_status;1];
elseif soc_pct > params.CutoverSOC
    fprintf('Cell SOC is above upper cutoff\n');
    exit_status = [exit_status;2];
end

if v_cell < params.CutoffVoltage
    fprintf('Cell Voltage is below lower cutoff\n');
    exit_status = [exit_status;3];
elseif v_cell > params.CutoverVoltage
    fprintf('Cell Voltage is above upper cutoff\n');
    exit_status = [exit_status;4];
end

overall_exit_status = sum(exit_status); % is 0 if no cutoffs were hit

end
