#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reference = "./input/Mus_musculus.GRCm39.dna.toplevel.fa"
params.reads = "./input/reads/long/*/*.gz"
params.results = "./results/graph"

maxcpus = Runtime.runtime.availableProcessors()
usedForks = 2

taskCpus = Math.round(maxcpus / usedForks)


process align_reads {
    
    cpus taskCpus
    maxForks usedForks

    input:
        file reference
        tuple val(strain), file(reads)

    output:
        tuple val(strain), file("*.sorted.bam")

    """
    if zcat ${reads} | head -n 12 | grep -q ccs; then
        PRESET=map-hifi
    else
        PRESET=map-pb
    fi

    minimap2 -R '@RG\tID:${strain}\tSM:${strain}' --MD -Y -t ${taskCpus} -ax map-pb ${reference} ${reads} | samtools view -bS - | samtools sort -T ./tmp -o ${strain}.sorted.bam - 

    """

}

workflow {

    reference = file(params.reference)

    Channel.fromPath(params.reads).map{file ->
        def parent = file.parent.name
        return tuple(parent, file)
    }.groupTuple(by: 0)
    .set{reads}

    align_reads(reference, reads)
}