#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.tag = "minigraph"
params.lab_data = "./results/mm10/lab/*.bed"
params.prev_data = "./results/mm10/previous/*.bed"
params.calls = "./results/calls/*-${params.tag}.*.bed"
params.results = "./results"

params.window_lab = 30
params.ranges = "0-50,050-100,100-1000,1000-10000,10000-100000000"

process new_lab_prev_stats {

    input:
        tuple val(strain), val(type), val(range), file(lab_data), file(new_data), file(prev_data)

    output:
        tuple val(strain), val(type), val(range), file("new_missed.bed"), file("stats_*.csv")

    script: 

    window = params.window_lab

    """
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${lab_data} > lab_data
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${new_data} > new_data
    awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2-${window},\$3+${window},\$4}' ${prev_data} > prev_data

    bedtools intersect -a lab_data -b new_data -wa > new_intersected
    bedtools intersect -a lab_data -b new_data -v > new_missed.bed

    bedtools intersect -a lab_data -b prev_data -wa > prev_intersected
    bedtools intersect -a lab_data -b prev_data -v > prev_missed.bed

    LAB_COUNT="\$(cat lab_data | awk '!x[\$0]++' | wc -l)"
    NEW_COUNT="\$(cat new_data | awk '!x[\$0]++' | wc -l)"
    PREV_COUNT="\$(cat prev_data | awk '!x[\$0]++' | wc -l)"

    NEW_INTERSECTED="\$(cat new_intersected | awk '!x[\$0]++' | wc -l)"
    NEW_MISSED="\$(cat new_missed.bed | awk '!x[\$0]++' | wc -l)"

    PREV_INTERSECTED="\$(cat prev_intersected | awk '!x[\$0]++' | wc -l)"
    PREV_MISSED="\$(cat prev_missed.bed | awk '!x[\$0]++' | wc -l)"
    
    echo "STRAIN=${strain}" > lab_stats.data
    echo "TYPE=${type}" >> lab_stats.data
    echo "RANGE=${range}" >> lab_stats.data
    echo "LAB_COUNT=\${LAB_COUNT}" >> lab_stats.data
    echo "NEW_COUNT=\${NEW_COUNT}" >> lab_stats.data
    echo "PREV_COUNT=\${PREV_COUNT}" >> lab_stats.data
    echo "NEW_INTERSECTED=\${NEW_INTERSECTED}" >> lab_stats.data
    echo "NEW_MISSED=\${NEW_MISSED}" >> lab_stats.data
    echo "PREV_INTERSECTED=\${PREV_INTERSECTED}" >> lab_stats.data
    echo "PREV_MISSED=\${PREV_MISSED}" >> lab_stats.data

    RND="\$(hexdump -n 16 -v -e '/1 "%02X"' -e '/16 "\\n"' /dev/urandom)"
    # cat lab_stats.data | cut -d'=' -f1 | paste -sd "," > stats_\${RND}.csv
    cat lab_stats.data | cut -d'=' -f2 | paste -sd "," >> stats_\${RND}.csv

    """
}

process split_ranges {

    input:
        each range // Format -> 50-100
        tuple val(strain), val(type), file(lab_data), file(new_data), file(prev_data)

    output:
        tuple val(strain), val(type), val(range), file("*.lab.bed"), file("*.new.bed"), file("*.prev.bed")
    
    script:

    rangeParts = range.tokenize("-")
    start = rangeParts.get(0)
    end = rangeParts.get(1)

    """
    cat ${lab_data} | awk -F"[\\t/]" 'BEGIN {OFS = "\\t"} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 <= ${end}' > "segment.lab.bed"
    cat ${new_data} | awk -F"[\\t/]" 'BEGIN {OFS = "\\t"} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 <= ${end}' > "segment.new.bed"
    cat ${prev_data} | awk -F"[\\t/]" 'BEGIN {OFS = "\\t"} {print \$1,\$2,\$3,\$4<0?\$4*-1:\$4}' | awk '\$4 >= ${start} && \$4 <= ${end}' > "segment.prev.bed"

    """
}

process collect_stats {

    publishDir file(params.results + '/support-files/validation/' + params.tag), mode: "copy"

    input:
        file stats
    
    output:
        file "*.csv"
    
    """
    echo "STRAIN,TYPE,RANGE,LAB_COUNT,NEW_COUNT,PREV_COUNT,NEW_INTERSECTED,NEW_MISSED,PREV_INTERSECTED,PREV_MISSED" > stats_${params.tag}.csv
    cat ${stats} >> stats_${params.tag}.csv
    """

}

process save_missed {

    publishDir file(params.results + '/support-files/validation/' + params.tag), mode: "copy"

    input:
        tuple val(strain), val(type), file(missed_bed)
    
    output:
        tuple val(strain), val(type), file("*.missed.bed")

    """
    cp ${missed_bed} "${strain}.${type}.missed.bed"
    """

}


def inputTuples(file) { 
    def strain = file.name.tokenize(".").get(0).tokenize("-").get(0)
    def type = file.name.tokenize(".").get(1)
    return tuple(strain, type, file)
}

workflow {

    new_data = Channel.fromPath(params.calls).map{f -> inputTuples(f)}
    lab_data = Channel.fromPath(params.lab_data).map{f -> inputTuples(f)}
    prev_data = Channel.fromPath(params.prev_data).map{f -> inputTuples(f)}
    
    ranges = Channel.of(params.ranges).splitCsv().flatten()
    
    lab_new_prev = lab_data.combine(new_data, by:[0,1]).combine(prev_data, by:[0,1])
    all_lab_new_prev = lab_new_prev.map{t -> return tuple(t[0],t[1],'ALL',t[2],t[3],t[4])} // reshape tuple (all ALL range)
    ranged_lab_new_prev = split_ranges(ranges, lab_new_prev)

    all_results = new_lab_prev_stats(all_lab_new_prev.concat(ranged_lab_new_prev))
    
    all_stats = all_results.map{t -> return t[4]}.collect()    
    collect_stats(all_stats)

    missed = all_results.filter{it[2] == "ALL"}.map{t -> return tuple(t[0],t[1],t[3])}
    save_missed(missed)

}

