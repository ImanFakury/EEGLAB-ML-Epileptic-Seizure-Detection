entropyValues = cell(1, length(PSD));

for i = 1:length(PSD)
    probabilities = PSD{i} / sum(PSD{i});
    entropyValues{i} = -sum(probabilities .* log2(probabilities));
    disp(['Shannon Entropy of Signal ' num2str(i) ': ' num2str(entropyValues{i})]);
end
figure;
bar(cell2mat(entropyValues));
title('Shannon Entropy of Each Signal');
xlabel('Signal Number');
ylabel('Entropy (bits)');
