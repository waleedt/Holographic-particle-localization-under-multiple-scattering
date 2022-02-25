function [fid_cost,s_est,SF,E_rb_all] = evaluateFidCost_BTR(s,Holo2D,forwardObj)
[s_est,SF,E_rb_all,~] = forwardObj.A(s);
resid_2D = s_est - Holo2D;
fid_cost = 0.5*norm(resid_2D(:), 2)^2;
end