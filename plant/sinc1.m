function y = sinc1(x)
% 
%              sin(x)
% sinc1(x) := --------
%                x
% 
    if abs(x) < 1e-2
        y = 1 - x^2/6 + x^4/120; % Taylor series
        % y = (1 - (7/60)*x^2) / (1 + (1/20)*x^2); % Pade approximant
    else
        y = sin(x)/x;
    end
end