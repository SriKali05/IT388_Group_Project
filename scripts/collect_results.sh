#!/bin/bash
#collect_results.sh - combine all the run outputs in results/raw/ into CSV files for Excel.
#
#   bash scripts/collect_results.sh [raw_dir] [out_dir]       (defaults: results/raw  results)
#
# Writes:
#   summary_avg.csv   one row per configuration (averaged over the repeats) - includes speedup and efficiency.
#    
#   summary.csv       one row per individual run (all repeats), for checking the spread.
#
# Raw files must be named <version>_r<ranks>_t<threads>_n<steps>_s<searches>_rep<k>.out
# (pfold_seq.sb and pfold_par.sb do this). (manually) Runs from the project folder, after the jobs finish.
#
# speedup    = mean total time of the sequential runs (same steps and searches) / mean total time
# efficiency = speedup / cores
# pct_moves, pct_energy = share of the available core time (total time x cores) spent making
#   moves / computing energy; the rest is everything else, including waiting and idle cores.

# AI Usage ACK: In line with the project's AI use scale, Claude, a generative AI tool, was used at a level 3
# to brainstorm effecient ways to process and calculate our experiment run results
# AI helped explain what grep, awk, and sed are and how/why they could be used to organize our output

cd "$(dirname "$0")/.." || exit 1
RAW=${1:-results/raw}
OUTDIR=${2:-results}
RUNS="$OUTDIR/summary.csv"
AVG="$OUTDIR/summary_avg.csv"
mkdir -p "$OUTDIR"

echo "version,ranks,threads,cores,proteins,searches,steps,rep,total_s,moves_s,energy_s,comm_s,imbalance,verified,moves_per_sec,file" > "$RUNS"

n=0
for f in "$RAW"/*.out; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .out)
    [[ $b =~ ^([a-z]+)_r([0-9]+)_t([0-9]+)_n([0-9]+)_s([0-9]+)_rep([0-9]+)$ ]] || continue
    ver=${BASH_REMATCH[1]}; r=${BASH_REMATCH[2]}; t=${BASH_REMATCH[3]}
    steps=${BASH_REMATCH[4]}; searches=${BASH_REMATCH[5]}; rep=${BASH_REMATCH[6]}

    total=$(sed -n 's/^Total time: \([0-9.]*\) sec.*/\1/p' "$f" | head -1)
    [ -n "$total" ] || { echo "skipping $b (no Total time line: run unfinished?)" >&2; continue; }
    proteins=$(sed -n 's/^[A-Za-z]*: \([0-9]*\) proteins.*/\1/p' "$f" | head -1)
    moves=$(sed -n 's/^Time in moves: \([0-9.]*\) sec, in energy: \([0-9.]*\) sec.*/\1/p' "$f" | head -1)
    energy=$(sed -n 's/^Time in moves: \([0-9.]*\) sec, in energy: \([0-9.]*\) sec.*/\2/p' "$f" | head -1)
    comm=$(sed -n 's/^Communication time on rank 0: \([0-9.]*\) sec.*/\1/p' "$f" | head -1)
    imb=$(sed -n 's/^Load imbalance (slowest \/ average): \([0-9.]*\).*/\1/p' "$f" | head -1)
    grep -q "Verification: PASSED" "$f" && ok=1 || ok=0
    mps=$(awk -v p="${proteins:-0}" -v s="$searches" -v n="$steps" -v t="$total" 'BEGIN{printf "%.0f", p*s*n/t}')

    echo "$ver,$r,$t,$((r * t)),$proteins,$searches,$steps,$rep,$total,$moves,$energy,$comm,$imb,$ok,$mps,$b.out" >> "$RUNS"
    n=$((n + 1))
done

# one row per configuration, averaged over the repeats
echo "version,ranks,threads,cores,steps,searches,reps,mean_total_s,min_total_s,max_total_s,speedup,efficiency,mean_moves_per_sec,mean_imbalance,mean_comm_s,pct_moves,pct_energy,all_verified" > "$AVG"
awk -F, 'NR > 1 {
    key = $1 "," $2 "," $3 "," $4 "," $7 "," $6
    if (!(key in cnt)) { order[++k] = key; mn[key] = $9; mx[key] = $9; okall[key] = 1 }
    cnt[key]++; tot[key] += $9
    if ($9 < mn[key]) mn[key] = $9
    if ($9 > mx[key]) mx[key] = $9
    mv[key] += $10; en[key] += $11
    if ($12 != "") { cm[key] += $12; ncm[key]++ }
    if ($13 != "") { im[key] += $13; nim[key]++ }
    mps[key] += $15
    if ($14 != 1) okall[key] = 0
    if ($1 == "seq") { bs[$7 "," $6] += $9; bn[$7 "," $6]++ }
}
END {
    for (i = 1; i <= k; i++) {
        key = order[i]; split(key, a, ",")
        mean = tot[key] / cnt[key]
        bk = a[5] "," a[6]
        if (bn[bk] > 0) { sp = (bs[bk] / bn[bk]) / mean; ef = sp / a[4]; sps = sprintf("%.3f", sp); efs = sprintf("%.3f", ef) }
        else { sps = ""; efs = "" }
        coretime = mean * a[4] * cnt[key]
        pm = (coretime > 0) ? sprintf("%.1f", 100 * mv[key] / coretime) : ""
        pe = (coretime > 0) ? sprintf("%.1f", 100 * en[key] / coretime) : ""
        cms = (ncm[key] > 0) ? sprintf("%.4f", cm[key] / ncm[key]) : ""
        ims = (nim[key] > 0) ? sprintf("%.2f", im[key] / nim[key]) : ""
        o = (a[1] == "seq") ? 1 : (a[1] == "openmp") ? 2 : (a[1] == "mpi") ? 3 : 4
        printf "%d,%s,%s,%s,%s,%s,%s,%d,%.4f,%.4f,%.4f,%s,%s,%.0f,%s,%s,%s,%s,%d\n", o, a[1], a[2], a[3], a[4], a[5], a[6], cnt[key], mean, mn[key], mx[key], sps, efs, mps[key] / cnt[key], ims, cms, pm, pe, okall[key]
    }
}' "$RUNS" | sort -t, -k1,1n -k6,6n -k5,5n -k3,3n | cut -d, -f2- >> "$AVG"

echo "wrote $n runs to $RUNS"
echo "wrote $(($(wc -l < "$AVG") - 1)) configurations to $AVG"