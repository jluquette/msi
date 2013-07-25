t <- read.table("stutter_distn.single_mda.flank10.minmapq36_minunits4_minsupp15_maxrefdiff80.txt", sep="\t", header=TRUE)
pcr <- read.table("stutter_distn.pcr_with_binom.flank10.minmapq36_minunits4_minsupp15_maxrefdiff80.txt", sep="\t", header=TRUE)

stdev <- function(x, total) {
  mu <- x/total
  (x*(1 - mu)^2 + (total-x)*(0 - mu)^2) / (total - 1)
}
sem <- function(x, y) stdev(x, y)/sqrt(y)

plotter <- function(unit) {
    x=t[t$unit==unit, 'reflen']
    frac=t[t$unit==unit, 'stutter_reads']
    tot=t[t$unit==unit, 'total_reads']
    y=t[t$unit==unit, 'percent_stutter_reads']
    x2=pcr[pcr$unit==unit, 'reflen']
    frac2=pcr[pcr$unit==unit, 'stutter_reads']
    tot2=pcr[pcr$unit==unit, 'total_reads']
    y2=pcr[pcr$unit==unit, 'percent_stutter_reads']
    par(mar=c(5,6,4,2))
    plot(x, y, type='l', ylab='Fraction of reads attributed to stutter\nError bars: SEM', xlab='Reference STR length', main=sprintf("%snucleotide repeats on chrX and chrY", unit))
    arrows(x0=x, y0=y-sem(frac,tot), y1=y+sem(frac,tot), code=3, angle=90, length=0.025)
    lines(x2, y2, type='l', col='red')
    print(rbind(x2,sem(frac2,tot2)))
    arrows(x0=x2, y0=y2-sem(frac2,tot2), y1=y2+sem(frac2,tot2), col='red', code=3, angle=90, length=0.025)
    legend('topleft', legend=c('PCR', 'MDA single cell'), fill=c('red', 'black'))
}

library(Cairo)

for (unit in c("mono", "di", "tri", "tetra")) {
    CairoPNG(sprintf("%s_stutter.png", unit), width=750, height=450)
    plotter(unit)
}
