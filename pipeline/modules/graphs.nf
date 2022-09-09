#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reference = "./input/Mus_musculus.GRCm39.dna.toplevel.chr19.fa"
params.chromosomes = "./input/chromosomes/*.chr19.fasta"
params.minigraph_bubbles = "./results/support-files/minigraph/*.bubbles.bed"
params.minigraph_gfa = "./results/support-files/minigraph/mouse_genomes_graph.gfa"
params.results = "./results"

maxcpus = Runtime.runtime.availableProcessors()

process genome_graph {

    cpus maxcpus
    publishDir file(params.results + '/support-files/minigraph'), mode: "copy"
    
    input:
        file reference
        file all_chromosomes

    output:
        file "*.gfa"

    """
    echo "Using ${task.cpus} CPUs"
    minigraph -xggs -t${task.cpus} ${reference} ${all_chromosomes} > mouse_genomes_graph.gfa
    """

}

process bubble_bed { 

    cpus 4
    maxForks 3
    publishDir file(params.results + '/support-files/minigraph'), mode: "copy"
    
    input:
        file chromosomes
        each genome_graph
    
    output:
        tuple val(strain), file("*.bubbles.bed")

    script:

    strain = chromosomes.name.tokenize('.').get(0)

    """
    minigraph -xasm -t${task.cpus} --call ${genome_graph} ${chromosomes} > ${strain}.bubbles.bed
    """
}

process call_simple_indels {

    publishDir file(params.results + '/calls/'), mode: "copy"

    input:
        tuple val(strain), file(bubbles_bed)
    
    output:
        tuple val(strain), file("*.bed")
    
    """
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$1==19 && \$2==\$3 && (\$3-\$2)<\$7) print \$1,\$2,\$3,\$7-(\$3-\$2)}' ${bubbles_bed} > "${strain}-minigraph.INS.simple.bed"
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$1==19 && \$6=="*" && (\$3-\$2)>\$7) print \$1,\$2,\$3,\$7-(\$3-\$2)}' ${bubbles_bed} > "${strain}-minigraph.DEL.simple.bed"
    """    
}


workflow {

    // reference = file(params.reference)
    // chromosomes_ch = Channel.fromPath(params.chromosomes)

    // graph = genome_graph(reference, chromosomes_ch.collect())
    // bubbles = bubble_bed(chromosomes_ch, graph)

    graph = file(params.minigraph_gfa)
    bubbles = Channel.fromPath(params.minigraph_bubbles)

    call_simple_indels(bubbles)

}