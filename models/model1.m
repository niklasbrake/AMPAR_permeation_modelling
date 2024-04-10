function G = model1(Params,Ca,V)
% MODEL1 computes the steady-state block probability for model 1.
%   G = model1(Params,Ca,V) returns the normalized conductance (G), given
%   model parmaeters (Params), and a vector of calcium concentrations (Ca)
%   and voltages (V), in mM and volts, respectively.
%
%   Parameter order should be: k1*, k1r*, Cai, k2r*, d1.

    % Get function that returns steady-state condutance of channel,
    % as function of calcium (Ca) and voltage (V).
    N = solve_system(Params);
    % Loop through each Ca and V value
    G = nan(length(Ca),length(V));
    for i = 1:length(Ca)
        for j = 1:length(V)
            G(i,j) = N(Ca(i),V(j));
        end
    end
end

function N = solve_system(Params)
% A = solve_system(Params) returns the steady-state probability of calcium
% not blocking the pore for model 1, given parameters Params. Parameters
% and equations defined in the manuscript,
% Supplemental Information - Description of calcium binding models.
    F = 96485; % C/M
    R = 8.31; % J/MK
    T = 293; %K
    z = 2; % # electrons

    k1x = Params(1);
    k1rx = Params(2);
    k2rx = Params(4);
    k2x = k1x*k2rx/k1rx;
    Cai = 10.^(-10*Params(3));
    d1 = Params(5);

    k1 = @(V) k1x * exp(-z*d1*F*V/(2*R*T));
    k1r = @(V) k1rx * exp(z*d1*F*V/(2*R*T));
    k2 = @(V) k2x * exp(z*(1-d1)*F*V/(2*R*T));
    k2r = @(V) k2rx * exp(-z*(1-d1)*F*V/(2*R*T));


    a = @(Cao,V) k1(V)*Cao + k2(V)*Cai;
    b = @(Cao,V) k1r(V) + k2r(V);
    N = @(Cao,V) b(Cao,V)/(a(Cao,V)+b(Cao,V));
end