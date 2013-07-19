for i in */*.str_summary.txt; do
    echo $i
    bsub -q short -W 4:00 "./metrics_msitools.py $i > ${i/.str_summary.txt/}.metrics.txt"
    bsub -q short -W 4:00 "./metrics_msitools.py --min-mapq 60 --min-units 5 --min-supp 20 --max-ref-diff 80 $i > ${i/.str_summary.txt/}.minmapq60_minunits5_minsupp20_maxrefdiff80.metrics.txt"
    bsub -q short -W 4:00 "./candidate_loci.py --min-mapq 60 --min-units 5 --min-supp 20 --max-ref-diff 80 $i > ${i/.str_summary.txt/}.minmapq60_minunits5_minsupp20_maxrefdiff80.candidates.txt"
done
