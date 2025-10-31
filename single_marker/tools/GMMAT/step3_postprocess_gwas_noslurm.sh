#!/bin/bash
module load R/4.2.2
SCRIPT_DIR=$(dirname "$(realpath "$0")")

model="model1"

MANHATTAN_INPUT="$OUTPUT_DIR/manhattan_input.txt"
SORTED_OUTPUT="$OUTPUT_DIR/manhattan_plot_input_sorted_$model.txt"
SUM_STATS="$OUTPUT_DIR/sum_stats.txt"

awk '{print $2,$1,$4,$11}' "$OUTPUT_DIR/mega_scores_chr1_kinship_09_23_2024_5PCs_$model.txt" | sed 's/ /\t/g' > temp1
sed -i '1d' temp1
echo -e "SNP\tCHR\tBP\tPVAL" > "$MANHATTAN_INPUT"
cat temp1 >> "$MANHATTAN_INPUT"
rm temp1
mv "$OUTPUT_DIR/mega_scores_chr1_kinship_09_23_2024_5PCs_$model.txt" "$OUTPUT_DIR/mega_scores_chr1_kinship_09_23_2024_5PCs_sorted_$model.txt"


awk -v OFS="\t" '{
  beta = ($10 + 0 !=0)? $9/$10 : "NA"
  se   = ($10 + 0 !=0)? 1/sqrt($10) : "NA"
  split($2,loc,":")
  print loc[1], loc[2], $0, beta, se
}' "$OUTPUT_DIR/mega_scores_chr1_kinship_09_23_2024_5PCs_sorted_$model.txt" > "$OUTPUT_DIR/temp_full.txt"

awk -F'\t' '$10 >= 0.01 && $10 < 0.99' "$OUTPUT_DIR/temp_full.txt" > "$OUTPUT_DIR/temp_filtered.txt"

echo -e "CHR\tBP\tCHR_full\tSNP\tcM\tBP\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL\tBETA\tSE" > "$SORTED_OUTPUT"
cat "$OUTPUT_DIR/temp_filtered.txt" >> "$SORTED_OUTPUT"


export OUTPUT_DIR

module load R/4.2.2

echo "--------------"
echo " here start the R scripts "
echo "--------------"
pwd
echo "$SCRIPT_DIR/create_manhattan.R"
ls

Rscript "$SCRIPT_DIR/create_manhattan.R"  "$OUTPUT_DIR/manhattan_input.txt"
Rscript "$SCRIPT_DIR/create_qq.plot.R"  "$OUTPUT_DIR/manhattan_input.txt"
Rscript "$SCRIPT_DIR/create_circular_manhattan.R"  "$OUTPUT_DIR/manhattan_input.txt"
Rscript "$SCRIPT_DIR/create_density_plot.R"  "$OUTPUT_DIR/manhattan_input.txt"

rm "$OUTPUT_DIR"/temp_*.txt
mv "$SORTED_OUTPUT" "$SUM_STATS"
echo "Post-processing completed in $OUTPUT_DIR"

