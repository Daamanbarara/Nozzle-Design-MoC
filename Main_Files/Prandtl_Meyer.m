%% Function for the computation of the nu using Prandtl Meyer's function
function [nu]=Prandtl_Meyer(M,gamma)
nu=(sqrt((gamma+1)/(gamma-1))*(atand(sqrt(((gamma-1)/(gamma+1))*((M^2)-1)))))-...
    atand(sqrt((M^2)-1));
end