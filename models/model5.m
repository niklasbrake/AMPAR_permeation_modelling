function [G,E] = model5(Params,Ca,V)
% MODEL5 computes the steady-state block probability for model 5.
%   [G,E] = model5(Params,Ca,V) returns the normalized conductance (G),
%   and the steady-state occupancy of each pore state (E), given model
%   parmaeters (Params), and a vector of calcium concentrations (Ca) and
%   voltages (V), in mM and volts, respectively.
%
%   Parameter order should be: k1*, k1r*, k2*, k2r*, k3*, k4r*, k4*, d1, d3.

    % Compute transition matrix for the parameters,as function of
    % calcium (Ca) and voltage (V).
    A = markov(Params);

    % Loop through each Ca and V value and compute the null
    % space of the transition matrix to get steady-state values
    E = nan(4,length(Ca),length(V));
    for i = 1:length(Ca)
        for j = 1:length(V)
            try
                E(:,i,j) = null(A(Ca(i),V(j)));
            catch
                E(:,i,j) = 0;
            end
        end
    end
    % Define condutcance as 1 minus probability of B1 or B2 occupancy,
    % i.e., either state 2, 3 or 4.
    E = E./sum(E);
    G = 1-squeeze(sum(E(2:4,:,:),1));
end

function A = markov(Params)
% A = markov(Params) returns the transition matrix for the model given
% parameters Params. Parameters and equations defined in the manuscript,
% Supplemental Information - Description of calcium binding models.

    F = 96485; % C/M
    R = 8.31; % J/MK
    T = 293; %K
    z = 2; % # electrons
    Cai = 1e-4; % 100 nM

    k1x = Params(1);
    k1rx = Params(2);
    k2x = Params(3);
    k2rx = Params(4);
    k3x = Params(5);
    k3rx = k2rx*k1rx*k3x/(k2x*k1x);
    k4x = Params(6);
    k4rx = Params(7);
    d1 = Params(8);
    d3 = Params(9);
    d2 = abs(d1-d3);
    d4 = max(0,1-max(d1,d3));

    k1 = @(V) k1x * exp(-z*d1*F*V/(2*R*T));
    k1r = @(V) k1rx * exp(z*d1*F*V/(2*R*T));
    k2 = @(V) k2x * exp(-z*d2*F*V/(2*R*T));
    k2r = @(V) k2rx * exp(z*d2*F*V/(2*R*T));
    k3 = @(V) k3x*exp(-z*d3*F*V/(2*R*T));
    k4 = @(V) k4rx*exp(z*d4*F*V/(2*R*T));
    k3r = @(V) k3rx*exp(z*d3*F*V/(2*R*T)) + k4x*exp(-z*d4*F*V/(2*R*T));

    A = @(Ca,V) [   -Ca*k3(V)-Cai*k4(V)-Ca*k1(V),  k1r(V),    0,  k3r(V); ...
                    Ca*k1(V),  -k1r(V)-k2(V)-Ca*k3(V)-Cai*k4(V),  k3r(V),    k2r(V); ...
                    0,  Ca*k3(V)+Cai*k4(V),  -k3r(V)-k1r(V),   Ca*k1(V); ...
                    Ca*k3(V)+Cai*k4(V),  k2(V), k1r(V),    -k2r(V)-Ca*k1(V)-k3r(V)];
end