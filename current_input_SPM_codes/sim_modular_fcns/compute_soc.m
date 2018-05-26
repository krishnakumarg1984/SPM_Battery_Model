function soc_pct = compute_soc(cs_avg_neg, param)

% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

soc_num = (cs_avg_neg/param.cs_max_n) - param.theta_min_neg;
soc_denom = param.theta_max_neg - param.theta_min_neg;

soc = soc_num/soc_denom; % as a fraction between 0 and 1
soc_pct = 100*soc; % in percentage [%]

end
