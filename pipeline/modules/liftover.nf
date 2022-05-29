#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input_files = "./input/validation/lab/*.bed"
params.remap_api = '/home/egarcia/appdir/bin/remap_api.pl'

//GCF_000001635.20 -> GRCm38 (mm10)
//GCF_000001635.27 -> GRCm39 (mm39)

params.remap_src = 'GCF_000001635.27'
params.remap_dest = 'GCF_000001635.20'

params.alias_src = 'validated'
params.alias_dest = 'validated.mm10'

params.results = "./results/mm10/lab"

process uplift_files {

    maxForks 2
    // publishDir file(params.results + '/' + params.alias_dest), mode: "copy"

    input:
        file input_bed 
    
    output:
        file "*.bed"

    script:
    new_name = input_bed.name.replace(params.alias_src, params.alias_dest)

    """
    perl ${params.remap_api} --mode asm-asm --from ${params.remap_src} --dest ${params.remap_dest} --annotation ${input_bed}  --annot_out ${new_name} 
    """
}

process clean_unmapped {

    publishDir file(params.results), mode: "copy"

    input:
        file input_bed
    
    output:
        file "*.bed"

    script:

    output = input_bed.name.replace(".bed", ".clean.bed")

    """
     grep -v CHR ${input_bed} > ${output}
    """
}


workflow {

    bed_files = Channel.fromPath(params.input_files)
    lifted_over = uplift_files(bed_files)
    clean_unmapped(lifted_over)

}