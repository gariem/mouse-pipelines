#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reference = "./input/Mus_musculus.GRCm39.dna.toplevel.chr19.fa"
params.chromosomes = "./input/chromosomes/*.chr19.fasta"
params.results = "./results/graph"

maxcpus = Runtime.runtime.availableProcessors()

process genome_graph {

    cpus maxcpus
    publishDir file(params.results + '/support-files/'), mode: "copy"
    
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
    publishDir file(params.results + '/support-files/'), mode: "copy"
    
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

process call_indels {

    publishDir file(params.results + '/calls/'), mode: "copy"

    input:
        tuple val(strain), file(bubbles_bed)
    
    output:
        tuple val(strain), file("*.bed")
    
    """
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)<\$7)print \$1,\$2,\$3,\$7-(\$3-\$2)}' ${bubbles_bed} > "${strain}-minigraph.INS.bed"
    awk -F"[\t:]" 'BEGIN {OFS = "\t"} {if(\$6!="."&&(\$3-\$2)>\$7)print \$1,\$2,\$3,\$7-(\$3-\$2)}' ${bubbles_bed} > "${strain}-minigraph.DEL.bed"
    """    
}


workflow {

    reference = file(params.reference)
    Channel.fromPath(params.chromosomes).multiMap{ file ->
        graph: file
        bubbles: file
    }.set{chromosomes}

    all_chromosomes = chromosomes.graph.collect()

    all_chromosomes.view()

    // graph = genome_graph(reference, all_chromosomes)

    // bubbles = bubble_bed(graph, chromosomes.bubbles)

    // call_indels(bubbles)
}