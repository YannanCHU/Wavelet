clc;
clear;
close all;

thresholds = 0:0.01:1;
PSNRs = zeros(size(thresholds));

parfor i = 1:length(thresholds)
    [~, PSNRi] = Main_SR(thresholds(i));
    PSNRs(i) = PSNRi;
end

[PSNR_optimal, threshold_optimal_index] = max(PSNRs);
threshold_optimal = thresholds(threshold_optimal_index);

%%
figure(1);
plot(thresholds, PSNRs, 'LineWidth', 2);
title("The PSNR against threshold value");
ylabel("PSNR value");   xlabel("Thresholds");

[SR_image, ~] = Main_SR(threshold_optimal);
figure();
image(SR_image);
title("The Super-resolved image with largest PSNR" + ...
        sprintf("\n(PSNR: %.2f dB; Threshold: ", PSNR_optimal) + " " + sprintf("%.2f", threshold_optimal) + ")");
