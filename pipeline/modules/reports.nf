#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.prev_data = "./input/validation/previous/*.bed"
params.calls = "./results/calls/*-minigraph.*.bed"
params.gene_data = "./input/Mus_musculus.GRCm39.106.chr.gff3"
params.repeats_data = "./input/repeatmasker/*.fasta.out"
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
    cat ${new_data} | awk -F"[\\t/]" 'BEGIN {OFS = "\\t"} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 < ${end}' > "segment.new.bed"
    cat ${prev_data} | awk -F"[\\t/]" 'BEGIN {OFS = "\\t"} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 < ${end}' > "segment.prev.bed"

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
    bedtools intersect -a new_data -b ${gene_data} -u > gene_new
    bedtools intersect -a prev_data -b ${gene_data} -u > gene_prev

    bedtools intersect -a novel_elements -b ${gene_data} -u > gene_novel
    bedtools intersect -a ${gene_data} -b novel_elements -u > genes_intersected
    
    NEW_COUNT="\$(cat new_data | awk '!x[\$0]++' | wc -l)"
    PREV_COUNT="\$(cat prev_data | awk '!x[\$0]++' | wc -l)"
    GENE_COUNT="\$(cat ${gene_data} | awk '!x[\$0]++' | wc -l)"

    NOVEL_COUNT="\$(cat novel_elements | awk '!x[\$0]++' | wc -l)"
    GENE_INTERSECTED_NEW="\$(cat gene_new | awk '!x[\$0]++' | wc -l)"
    GENE_INTERSECTED_PREV="\$(cat gene_prev | awk '!x[\$0]++' | wc -l)"
    GENE_INTERSECTED_NOVEL="\$(cat gene_novel | awk '!x[\$0]++' | wc -l)"
    INTERSECTED_GENES="\$(cat genes_intersected | awk '!x[\$0]++' | wc -l)"
    

    echo "STRAIN=${strain}" > gene_stats.data
    echo "TYPE=${type}" >> gene_stats.data
    echo "RANGE=${range}" >> gene_stats.data

    echo "NEW_COUNT=\${NEW_COUNT}" >> gene_stats.data
    echo "PREV_COUNT=\${PREV_COUNT}" >> gene_stats.data
    echo "GENE_COUNT=\${GENE_COUNT}" >> gene_stats.data

    echo "NOVEL_COUNT=\${NOVEL_COUNT}" >> gene_stats.data
    echo "GENE_INTERSECTED_PREV=\${GENE_INTERSECTED_PREV}" >> gene_stats.data
    echo "GENE_INTERSECTED_NEW=\${GENE_INTERSECTED_NEW}" >> gene_stats.data
    echo "GENE_INTERSECTED_NOVEL=\${GENE_INTERSECTED_NOVEL}" >> gene_stats.data
    echo "INTERSECTED_GENES=\${INTERSECTED_GENES}" >> gene_stats.data

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
    echo "STRAIN,TYPE,RANGE,NEW_COUNT,PREV_COUNT,GENE_COUNT,NOVEL_COUNT,GENE_INTERSECTED_PREVs,GENE_INTERSECTED_NEWs,GENE_INTERSECTED_NOVELs,INTERSECTED_GENES" > gene_analysis.csv
    cat ${stats} >> gene_analysis.csv

    """

}


process repeat_coordinates {

    input:
        tuple val(strain), file(rmasker_out)

    output:
        tuple val(strain), file("*.bed")
    
    script:

    """
    awk -v "OFS=\\t" '{print \$5,\$6,\$7,\$11}' ${rmasker_out} | awk -F"[\\t#]" -v OFS='\\t' '{gsub("chr","",\$3); print \$3,\$4,\$5,\$6}' | awk 'BEGIN {FS='\\t'}; \$1 ~/^[0-9]*\$|^X\$/{print}' > ${strain}.repeats.bed
    """

}

process intersect_repeats {

    input:
        tuple val(strain), val(type), file(calls), file(repeats)

    output:
        tuple val(strain), val(type), file("repeats_stats_*.out.csv"), file("repeats_details_*.out.csv")

    """
    bedtools intersect -a ${repeats} -b  ${calls} -u > repeats_in_calls
    bedtools intersect -a ${calls} -b  ${repeats} -u > calls_in_repeats
    bedtools intersect -a ${calls} -b  ${repeats} -v > calls_out_repeats

    CALLS_COUNT="\$(cat ${calls} | awk '!x[\$0]++' | wc -l)"
    REPEATS_COUNT="\$(cat ${repeats} | awk '!x[\$0]++' | wc -l)"

    REPEATS_IN_CALLS="\$(cat repeats_in_calls | awk '!x[\$0]++' | wc -l)"
    CALLS_IN_REPEATS="\$(cat calls_in_repeats | awk '!x[\$0]++' | wc -l)"
    CALLS_OUT_REPEATS="\$(cat calls_out_repeats | awk '!x[\$0]++' | wc -l)"

    echo "STRAIN=${strain}" > repeats_stats.data
    echo "TYPE=${type}" >> repeats_stats.data
    
    echo "CALLS_COUNT=\${CALLS_COUNT}" >> repeats_stats.data
    echo "REPEATS_COUNT=\${REPEATS_COUNT}" >> repeats_stats.data

    echo "REPEATS_IN_CALLS=\${REPEATS_IN_CALLS}" >> repeats_stats.data
    echo "CALLS_IN_REPEATS=\${CALLS_IN_REPEATS}" >> repeats_stats.data
    echo "CALLS_OUT_REPEATS=\${CALLS_OUT_REPEATS}" >> repeats_stats.data

    RND="\$(hexdump -n 16 -v -e '/1 "%02X"' -e '/16 "\\n"' /dev/urandom)"
    # cat repeats_stats.data | cut -d'=' -f1 | paste -sd "," > repeats_stats_\${RND}.out.csv
    cat repeats_stats.data | cut -d'=' -f2 | paste -sd "," >> repeats_stats_\${RND}.out.csv

    awk '{arr[\$4]+=1}END {for (key in arr) printf("%s\\t%s\\t%s\\t%s\\n","${strain}","${type}",key, arr[key])}' repeats_in_calls | sort -k1,1 | tr '\\t' ','  > repeats_details_\$RND.out.csv

    """
}

process collect_repeat_stats {

    publishDir file(params.results + '/support-files/analysis'), mode: "copy"

    input: 
        file summary
        file detail_repeats

    output:
        file "*.csv"

    """
    echo "STRAIN,TYPE,CALLS_COUNT,REPEATS_COUNT,REPEATS_IN_CALLS,CALLS_IN_REPEATS,CALLS_OUT_REPEATS" > repeats_analysis.csv
    cat ${summary} >> repeats_analysis.csv

    echo "STRAIN,TYPE,CLASS,COUNT" > repeats_detail_analysis.csv
    cat ${detail_repeats} >> repeats_detail_analysis.csv

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
    repeats_data = Channel.fromPath(params.repeats_data).map{f -> inputTuples(f)}.map{t -> return tuple(t[0], t[2])}

    genes_cords = gene_data(file(params.gene_data))
    repeats_cords = repeat_coordinates(repeats_data)

    new_prev = new_data.combine(prev_data, by:[0,1])
    all_new_prev = new_prev.map{t -> return tuple(t[0],t[1],'ALL',t[2],t[3])} // reshape tuple (all ALL range)
    
    ranges = Channel.of(params.ranges).splitCsv().flatten()
    ranged_new_prev = split_ranges(ranges, new_prev)

    gene_stats = intersect_prev_genes(ranged_new_prev.concat(all_new_prev).combine(genes_cords)).map{t -> return t[3]}.collect()
    collect_stats(gene_stats)

    new_repeats = new_data.combine(repeats_cords, by:0)

    repeats_stats = intersect_repeats(new_repeats)
    repeats_stats_summary = repeats_stats.map{t -> return t[2]}.collect()
    repeats_stats_detail = repeats_stats.map{t -> return t[3]}.collect()

    collect_repeat_stats(repeats_stats_summary, repeats_stats_detail).view()
       
    
    
}