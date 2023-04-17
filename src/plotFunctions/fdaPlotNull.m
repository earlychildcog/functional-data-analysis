function fdaPlotNull(T, Tnulls, tvec, perms)

figure
hold on
for p = perms
    plot(tvec, abs(Tnulls(p,:)),'--', 'LineWidth', 1)
end

plot(tvec, abs(T'),'r-', 'LineWidth', 3)


title(sprintf('original data vs %d permutations', length(perms)))
xlabel('time (s)')
ylabel('T-statistic')

maxT = max(abs(T));
yline(maxT, 'LineStyle',':','LineWidth',2)
xlim([min(tvec) max(tvec)])


legend([{'permuted data'} repmat({''}, [1 length(perms)-1]) {'original data' ''}])