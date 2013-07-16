for i in ../sputnik/*.sputnik.txt.gz; do
    x=$(basename $i .rmdup.nomapq0.sputnik.txt.gz)
    bsub -q short -W 12:00 -o $x/$x.log -R 'rusage[mem=32000]' "gunzip -c $i | ./msitools.pl --input /dev/stdin --outprefix $x/$x.flank10 --resource_path `pwd` --flank_bp 10"
done
