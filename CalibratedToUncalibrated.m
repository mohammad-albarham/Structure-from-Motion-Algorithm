function P_un = CalibratedToUncalibrated(P_ca, K)

N = length(P_ca);

P_un = cell(1,N);
for i=1:N
    P_un{i} = K * P_ca{i};
end

end