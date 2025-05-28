#!/bin/bash
#SBATCH --job-name=gwas_postprocess
#SBATCH --output=logs/postprocess.out
#SBATCH --error=logs/postprocess.err
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --dependency=singleton
#SBATCH --exclude=node50,node51,node52,node53,node54,node55,node56,node57,node58,node59,node60,node61,node62,node63,node64,node65,node66,node67,node68,node69

model="model1"

echo "==> Concatenating files..."
for i in {1..22}; do
  tail -n +2 mega_scores_chr"$i"_kinship_09_23_2024_5PCs_"$model".txt
done > manhattan_plot_input_sorted.txt

echo "==> Calculating beta and SE..."
awk -v OFS="\t" '
  BEGIN { FS=OFS }
  {
    beta = ($10 + 0 != 0) ? $9 / $10 : "NA"
    se   = ($10 + 0 != 0) ? 1 / sqrt($10) : "NA"
    split($2, loc, ":")
    print loc[1], loc[2], $0, beta, se
  }
' manhattan_plot_input_sorted.txt > temp_full.txt

echo "==> Filtering by allele frequency..."
awk -F'\t' '$10 >= 0.01 && $10 < 0.99' temp_full.txt > temp_filtered.txt

# Add header
echo -e "CHR\tBP\tCHR_full\tSNP\tcM\tPOS\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL\tBETA\tSE" > manhattan_plot_input_sorted_"$model".txt
cat temp_filtered.txt >> manhattan_plot_input_sorted_"$model".txt

echo "==> Creating input for plotting..."
awk -F'\t' '$14 != "NA" { print $4, $1, $2, $13 }' manhattan_plot_input_sorted_"$model".txt | tr ' ' '\t' > manhattan_plot_input_ready.txt

echo "==> Loading R and generating plots..."
module load R/4.2.2

# === Inline R for Manhattan and QQ plot generation ===
Rscript - <<'EOF'
#################################
#   By Basilio Cieza Huaman
#################################

# Manhattan Plot
library(CMplot)
data <- read.table("manhattan_plot_input_ready.txt", header = TRUE, sep = "\t")
CMplot(data,
       plot.type = 'm',
       cex = 1,
       band = 1,
       ylim = c(0, 15),
       col = c("grey30", "grey60"),
       threshold = c(5e-8),
       threshold.col = c("red"),
       threshold.lty = c(5),
       threshold.lwd = c(2),
       amplify = FALSE,
       LOG10 = TRUE)

# QQ Plot with Lambda
library(qqman)
p_values <- data$PVAL
observed_chisq <- qchisq(1 - p_values, df = 1)
lambda <- median(observed_chisq) / 0.4549
print(paste("Lambda:", round(lambda, 3)))

jpeg("qqplot_with_lambda.jpg", width = 800, height = 800, res = 300)
qq(p_values, main = "QQ Plot of P-values")
text(x = 0.3, y = max(-log10(p_values)) - 0.5,
     labels = paste0("Lambda = ", round(lambda, 3)),
     pos = 4, col = "blue", cex = 0.4)
dev.off()
EOF

echo "==> Cleaning up temporary files..."
rm temp_*.txt manhattan_plot_input_sorted.txt

echo "Post-processing completed successfully."
