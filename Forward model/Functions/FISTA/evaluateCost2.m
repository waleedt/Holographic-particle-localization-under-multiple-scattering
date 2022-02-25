function [cost,fid_cost,reg_cost] = evaluateCost2(resid_2D,xhat,tau,tvCost)

fid_cost = 0.5*norm(resid_2D(:), 2)^2;
reg_cost = tau*tvCost(xhat);
cost = fid_cost + reg_cost;
end