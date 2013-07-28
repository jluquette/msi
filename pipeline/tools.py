from cosmos.contrib.ezflow.tool import Tool

class align(Tool):
    inputs = [ 'fastq.gz' ]
    outputs = [ 'bam' ]
    mem_req = 4096 + 10240
    cpu_req = 3
    time_req = 60
    name = 'BWA MEM paired alignment'

    def cmd(self, i, s, p):
        return """{s[bwa_binary]} mem
                    -M
                    -v 2
                    -t 2
                    -R '@RG\tID:{p[readgroup]}\tSM:{p[sample]}\tPL:ILLUMINA'
                    {s[reference_genome]}
                    {i[fastq.gz][0]}
                    {i[fastq.gz][1]}
                    | {s[java_binary]}
                        -Xmx10g
                        -Djava.io.tmpdir={s[tmpdir]}
                        -jar {s[picard_home]}/SortSam.jar
                        SORT_ORDER=coordinate
                        I=/dev/stdin
                        O=$OUT.bam"""


class remdup(Tool):
    inputs = [ 'bam' ]
    outputs = [ 'bam', 'markdup_metrics.txt' ]
    mem_req = 20480
    cpu_req = 1
    time_req = 48*60
    name = 'Remove PCR duplicates'

    def cmd(self, i, s, p):
        inlist = ' '.join('I=' + str(f) for f in i['bam'])
        return """{s[java_binary]}
                    -Xmx20g
                    -Djava.io.tmpdir={s[tmpdir]}
                    -jar {s[picard_home]}/MarkDuplicates.jar
                    VALIDATION_STRINGENCY=LENIENT
                    REMOVE_DUPLICATES=true
                    O=$OUT.bam
                    %s
                    METRICS_FILE=$OUT.markdup_metrics.txt""" % inlist


class index_bam(Tool):
    inputs = [ 'bam' ]
    forward_input = True
    mem_req = 1024
    cpu_req = 1
    time_req = 12*60
    name = "Index BAM"

    def cmd(self, i, s, p):
        return """{s[samtools_binary]} index {i[bam][0]} {i[bam][0]}.bai"""


class sputnik(Tool):
    inputs = [ 'bam' ]
    outputs = [ 'bam' ]
    mem_req = 1024
    cpu_req = 2
    time_req = 12*60
    name = 'sputnik'

    def cmd(self, i, s, p):
        return """{s[sputnik_wrapper]}
                    {i[bam][0]}
                    $OUT.bam
                    {p[chrom]}"""


class msitools(Tool):
    inputs = [ 'bam' ]
    outputs = [ 'bam', 'str_summary.txt' ]
    mem_req = 20480
    cpu_req = 1
    time_req = 12*60
    name = 'msitools'

    def cmd(self, i, s, p):
        return """samtools view -h {i[bam][0]}
                  | {s[msitools_script]}
                    --resource_path {s[resource_path]}
                    --summary $OUT.str_summary.txt
                    --flank_bp 10
                    --chr {p[chrom]}
                  | samtools view -bSo $OUT.bam -"""
