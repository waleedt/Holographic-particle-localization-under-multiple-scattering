function [total_cost,reg_cost] = evaluateCost3(fid_cost,xhat,tau,tvCost)

reg_cost = tau*tvCost(xhat);
total_cost = fid_cost + reg_cost;
end