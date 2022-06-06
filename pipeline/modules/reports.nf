#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.prev_data = "./input/validation/previous/*.bed"
params.calls = "./results/calls/*-minigraph.*.bed"
params.gene_data = "./input/Mus_musculus.GRCm39.106.chr.gff3"
params.results = "./results"

params.window_lab = 30
params.ranges = "0-50,050-100,100-1000,1000-10000,10000-100000000"

process gene_data {

    input:
        file gff3_data
    
    output:
        file "*.bed"
    
    """
    gff2bed < ${gff3_data} | grep "ID=gene:" | grep ";biotype=protein_coding" | awk -F"\\t" 'BEGIN  {OFS = "\\t"} {print \$1,\$2,\$3,\$4}' > genes.GRCm39.bed
    """
}

process split_ranges {

    input:
        each range // Format -> 50-100
        tuple val(strain), val(type), file(new_data), file(prev_data)

    output:
        tuple val(strain), val(type), val(range), file("*.new.bed"), file("*.prev.bed")
    
    script:

    rangeParts = range.tokenize("-")
    start = rangeParts.get(0)
    end = rangeParts.get(1)

    """
    cat ${new_data} | awk -F"[\\t/]" 'BEGIN {OFS = "\\t"} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 <= ${end}' > "segment.new.bed"
    cat ${prev_data} | awk -F"[\\t/]" 'BEGIN {OFS = "\\t"} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 <= ${end}' > "segment.prev.bed"

    """
}

process intersect_prev_genes{

    input:
        tuple val(strain), val(type), val(range), file(new_data), file(prev_data), file(gene_data)
    
    output:
        tuple val(strain), val(type), val(range), file("stats_*.csv")

    script:

    window = params.window_lab

    """
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${new_data} > new_data
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${prev_data} > prev_data

    bedtools intersect -a new_data -b prev_data -v > novel_elements
    bedtools intersect -a new_data -b ${gene_data} -wa > gene_new
    bedtools intersect -a prev_data -b ${gene_data} -wa > gene_prev

    bedtools intersect -a novel_elements -b ${gene_data} -wa > gene_novel
    
    NEW_COUNT="\$(cat new_data | awk '!x[\$0]++' | wc -l)"
    PREV_COUNT="\$(cat prev_data | awk '!x[\$0]++' | wc -l)"
    GENE_COUNT="\$(cat ${gene_data} | awk '!x[\$0]++' | wc -l)"

    NOVEL_COUNT="\$(cat novel_elements | awk '!x[\$0]++' | wc -l)"
    GENES_INTERSECTED_NEW="\$(cat gene_new | awk '!x[\$0]++' | wc -l)"
    GENES_INTERSECTED_PREV="\$(cat gene_prev | awk '!x[\$0]++' | wc -l)"
    GENES_INTERSECTED_NOVEL="\$(cat gene_novel | awk '!x[\$0]++' | wc -l)"

    echo "STRAIN=${strain}" > gene_stats.data
    echo "TYPE=${type}" >> gene_stats.data
    echo "RANGE=${range}" >> gene_stats.data

    echo "NEW_COUNT=\${NEW_COUNT}" >> gene_stats.data
    echo "PREV_COUNT=\${PREV_COUNT}" >> gene_stats.data
    echo "GENE_COUNT=\${GENE_COUNT}" >> gene_stats.data

    echo "NOVEL_COUNT=\${NOVEL_COUNT}" >> gene_stats.data
    echo "GENES_INTERSECTED_PREV=\${GENES_INTERSECTED_PREV}" >> gene_stats.data
    echo "GENES_INTERSECTED_NEW=\${GENES_INTERSECTED_NEW}" >> gene_stats.data
    echo "GENES_INTERSECTED_NOVEL=\${GENES_INTERSECTED_NOVEL}" >> gene_stats.data

    RND="\$(hexdump -n 16 -v -e '/1 "%02X"' -e '/16 "\\n"' /dev/urandom)"
    # cat gene_stats.data | cut -d'=' -f1 | paste -sd "," > stats_\${RND}.csv
    cat gene_stats.data | cut -d'=' -f2 | paste -sd "," >> stats_\${RND}.csv

    """
}

process collect_stats {

    publishDir file(params.results + '/support-files/analysis'), mode: "copy"

    input:
        file stats
    
    output:
        file "*.csv"
    
    """
    echo "STRAIN,TYPE,RANGE,NEW_COUNT,PREV_COUNT,GENE_COUNT,NOVEL_COUNT,GENES_INTERSECTED_PREV,GENES_INTERSECTED_NEW,GENES_INTERSECTED_NOVEL" > gene_analysis.csv
    cat ${stats} >> gene_analysis.csv

    """

}

def inputTuples(file) { 
    def strain = file.name.tokenize(".").get(0).tokenize("-").get(0)
    def type = file.name.tokenize(".").get(1)
    return tuple(strain, type, file)
}

workflow {

    new_data = Channel.fromPath(params.calls).map{f -> inputTuples(f)}
    prev_data = Channel.fromPath(params.prev_data).map{f -> inputTuples(f)}

    genes_cords = gene_data(file(params.gene_data))

    new_prev = new_data.combine(prev_data, by:[0,1])

    ranges = Channel.of(params.ranges).splitCsv().flatten()
    ranged_new_prev = split_ranges(ranges, new_prev)

    gene_stats = intersect_prev_genes(ranged_new_prev.combine(genes_cords)).map{t -> return t[3]}.collect()
    collect_stats(gene_stats).view()

    
}