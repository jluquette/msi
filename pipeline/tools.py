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
                    -t {self.cpu_req}
                    -R '@RG\tID:{p[readgroup]}\tSM:{p[sample]}\tPL:ILLUMINA'
                    {s[reference_genome]}
                    {i[fastq.gz][0]}
                    {i[fastq.gz][1]}
                    | {s[java_binary]}
                        -Xmx10g
                        -Djava.io.tmpdir={s[tmpdir]}
                        -jar {s[picard_home]}/SortSam.jar
                        I=/dev/stdin
                        O=/dev/stdout
                        SORT_ORDER=coordinate
                    | samtools view -bSo $OUT.bam -"""


class remdup(Tool):
    inputs = [ 'bam' ]
    outputs = [ 'bam', 'markdup_metrics.txt' ]
    mem_req = 20480
    cpu_req = 1
    time_req = 12*60
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
