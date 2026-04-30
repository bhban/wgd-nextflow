mkdir -p rooted_trees rooting_summaries

for tree in tree_*.nwk; do
    base=$(basename "$tree" .nwk)

    python root_gene_tree.py \
      --tree "$tree" \
      --genomes-tsv genomes.tsv \
      --output-tree "rooted_trees/${base}.rooted.nwk" \
      --summary-tsv "rooting_summaries/${base}.summary.tsv"
done

head -n 1 rooting_summaries/$(ls rooting_summaries | head -n 1) > all_rooting_summaries.tsv

for summary in rooting_summaries/*.summary.tsv; do
    tail -n +2 "$summary" >> all_rooting_summaries.tsv
done

cat all_rooting_summaries.tsv

cut -f3 all_rooting_summaries.tsv \
  | grep -o 'tier=[0-9]*' \
  | sort \
  | uniq -c
