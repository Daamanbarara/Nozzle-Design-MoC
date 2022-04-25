%% Function for the computation of the nu using Prandtl Meyer's function
function [M_cal]=Inverse_Prandtl_Meyer(nu,gamma)
M=linspace(1,5,10000);
for i=1:length(M)
    nu_cal(i)=(sqrt((gamma+1)/(gamma-1))*(atand(sqrt(((gamma-1)/(gamma+1))*((M(i)^2)-1)))))-...
    atand(sqrt((M(i)^2)-1));
end
% [j]=find(round(nu_cal,3) == nu)
A = repmat(nu_cal,[1 length(nu)]);
[minValue,closestIndex] = min(abs(A-nu'));
closestValue = nu_cal(closestIndex); 
M_cal=M(closestIndex);

end