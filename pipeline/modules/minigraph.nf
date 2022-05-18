#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reference = './input/Mus_musculus.GRCm39.dna.toplevel.fa'
params.chromosomes = './input/chromosomes/*.fasta'
params.results = './results/graph'

process genome_graph {

    publishDir file(params.results + '/support-files/'), mode: "copy"

    cpus 8

    input:
        file reference
        file all_chromosomes

    output:
        file "*.gfa"

    """
    minigraph -xggs -t${task.cpus} ${reference} ${all_chromosomes} > mouse_genomes_graph.gfa
    """

}

process bubble_bed { 

    publishDir file(params.results + '/support-files/'), mode: "copy"
    
    input:
        file genome_graph
        each file chromosomes
    
    output:
        tuple val(strain), file("*.bubbles.bed")

    script:

    strain_name = chromosomes.name.tokenize('.').get(0)

    """
    minigraph -xasm -t${task.cpus} --call ${genome_graph} ${chromosomes} > ${strain_name}.bubbles.bed
    """
}

process call_indels {

    publishDir file(params.results + '/calls/'), mode: "copy"

    input:
        tuple val(strain), file(bubbles_bed)
    
    output:
        tuple val(strain), file("*.bed")
    
    """
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)<\$7)print \$1,\$2,\$3,\$7-(\$3-\$2)}' ${bubbles_bed} > "${strain}-minigraph__INS.bed"
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)>\$7)print \$1,\$2,\$3,\$7-(\$3-\$2)}' ${bubbles_bed} > "${strain}-minigraph__DEL.bed"
    """    
}


workflow {

    reference = file(params.reference)
    chromosomes_ch = Channel.fromPath(params.chromosomes)

    graph = genome_graph(reference, chromosomes_ch.collect())

    bubbles = bubble_bed(graph, chromosomes_ch)

    call_indels(bubbles)
}