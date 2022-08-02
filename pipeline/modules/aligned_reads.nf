#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// params.reference = "./input/Mus_musculus.GRCm39.dna.toplevel.chr1.fa"
params.reference = "/media/egarcia/DataBank/mouse/igv_workfiles/Mus_musculus.GRCm39.dna.toplevel.fa"
params.reads = "./input/reads/long/*/*.gz"
params.repeats = "./results/support-files/Mus_musculus.GRCm39.tandemrepeats.bed"
params.aligned_reads = "./results/support-files/alignments/*.bam"
params.results = "./results"

maxcpus = Runtime.runtime.availableProcessors()
usedForks = 2
taskCpus = Math.round(maxcpus / usedForks)

process tandem_repeats { 

    cpus 8
    publishDir file(params.results + '/support-files/'), mode: "copy"

    input: 
        file reference

    output:
        file "Mus_musculus.GRCm39.tandemrepeats.bed"

    """
    python -m pbsvtools.tasks.split_ref_to_chrs ${reference} GRCm39.chr_size.csv
    python -m pbsvtools.tasks.tandem_repeat_finder --nproc 8 ${reference} GRCm39.chr_size.csv Mus_musculus.GRCm39.tandemrepeats.bed
    """
}

process align_reads {
    
    cpus taskCpus
    maxForks usedForks

    publishDir file(params.results + '/support-files/alignments'), mode: "copy"

    input:
        tuple val(strain), file(reads)
        file reference

    output:
        tuple val(strain), file("*.sorted.bam")

    """
    if zcat ${reads} | head -n 12 | grep -q ccs; then
        PRESET=map-hifi
    else
        PRESET=map-pb
    fi

    mkdir tmp
    minimap2 -R '@RG\\tID:${strain}\\tSM:${strain}' --MD -Y -t ${taskCpus} -ax \$PRESET ${reference} ${reads} | samtools view -bS - | samtools sort -T ./tmp -o ${strain}.sorted.bam - 
    """

}

process discover_pbsv {

    cpus taskCpus
    maxForks usedForks

    publishDir file(params.results + '/support-files/pbsv'), mode: "copy"

    input:
        tuple val(strain), file(aligned_reads)
        file repeats

    output:
        tuple val(strain), file("*.svsig.gz")

    """
    pbsv discover --tandem-repeats ${repeats} ${aligned_reads} "${strain}.svsig.gz"
    """
}

process call_pbsv {

    cpus taskCpus
    maxForks usedForks

    publishDir file(params.results + '/support-files/pbsv'), mode: "copy"

    input:
        tuple val(strain), file(sv_signature)
        file reference

    output:
        tuple val(strain), file("*.vcf")

    """
    pbsv call -j ${taskCpus} -t "DEL,INS,INV,DUP" ${reference} ${sv_signature} "${strain}-pbsv.vcf"
    """
}

process filter_pbsv {

     input: 
        tuple val(strain), file(vcf_file)
        val ad_limits
    
    output:
        tuple val(strain), file("*.vcf")

    script:

    outname = vcf_file.name.replace(".vcf", ".filtered.vcf") 
    """
    bcftools view -i'(AD[0:1] - AD[0:0])>=${ad_limits[0]} && (AD[0:1] - AD[0:0])<=${ad_limits[1]}' ${vcf_file} > "${outname}"
    """
}


process bed_files {

    publishDir file(params.results + '/calls/'), mode: "copy"

    input: 
        tuple val(strain), file(vcf_file)
        each type
    
    output:
        tuple val(strain), val(type), file('*.bed')

    script:

    simple_name = vcf_file.name.replace(".vcf","")
    strain = simple_name.tokenize('-').get(0)

    """
    bcftools query -i"SVTYPE='${type}'" -f'%CHROM\\t%POS0\\t%END0\\t%SVLEN\\n' ${vcf_file} | \
            awk -F'\\t' 'BEGIN {OFS = FS} {print \$1,\$2,\$3,\$4}' > "${strain}-pbsv.${type}.bed"
    """
}

workflow {

    reference = file(params.reference)

    // repeats = tandem_repeats(reference)
    repeats = file(params.repeats)

    Channel.fromPath(params.reads).map{file ->
        def parent = file.parent.name
        return tuple(parent, file)
    }.groupTuple(by: 0)
    .set{reads}

    //aligned_reads = align_reads(reads, reference)
    aligned_reads = Channel.fromPath(params.aligned_reads).map {it -> return tuple(it.name.tokenize('.').get(0), it)}

    sv_signatures = discover_pbsv(aligned_reads, repeats)
    vcf = call_pbsv(sv_signatures, reference)

    //filter_pbsv(vcf, [])
    sv_types = Channel.from(['INS', 'DEL', 'INV', 'DUP'])
    bed_files(vcf, sv_types)

}
