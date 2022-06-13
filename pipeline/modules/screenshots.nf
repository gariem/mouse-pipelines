#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.call_dir = './results/calls'
params.input_files = './results/support-files/validation/minigraph/DBA_2J*.INS.*.10.bed'
params.igv_workdir = '/media/egarcia/DataBank/mouse/igv_workfiles'
params.results = "./results"

maxcpus = Runtime.runtime.availableProcessors()
usedForks = 4
taskCpus = Math.round(maxcpus / usedForks)

process prepare_screnshot_data {

    cpus taskCpus
    publishDir file(params.igv_workdir), mode: "copy", saveAs: {
                    filename -> filename.split('\\.')[0] + '/' + filename
                }    

    input:
        tuple val(strain), val(type), file(bed_file)
        file igv_workdir
    
    output:
        tuple val(strain), val(type), file("*.bed")
        file("*.session.xml") 
        tuple file("*.bam"), file("*.bam.bai") optional true

    """
    cat ${bed_file} | awk '{print \$1}' | awk '!x[\$0]++' > chromosomes.txt

    while IFS= read -r chr
    do
        if [ ! -f ${igv_workdir}/${strain}/${strain}.\$chr.contigs.bam ]
        then
            samtools view -b -@ ${taskCpus} ${igv_workdir}/${strain}/${strain}.contigs.bam \$chr > ${strain}.\$chr.contigs.bam
            samtools index -@ ${taskCpus} ${strain}.\$chr.contigs.bam
        fi

        if [ ! -f ${igv_workdir}/${strain}/${strain}.\$chr.h1.contigs.bam ]
        then
            MAX_POINT="\$(samtools view -H ${igv_workdir}/${strain}/${strain}.\$chr.contigs.bam | awk '{if(\$1~/@SQ/){print \$2"\\t1\\t"\$3}}' | sed 's/[SL]N://g' | grep -P "^\$chr\\t" | awk '{print \$3}')"
            MID_POINT=\$((\${MAX_POINT}/2))

            samtools view -@ ${taskCpus} -h ${igv_workdir}/${strain}/${strain}.\$chr.contigs.bam \$chr:1-\$MID_POINT -b > ${strain}.h1.\$chr.contigs.bam
            samtools view -@ ${taskCpus} -h ${igv_workdir}/${strain}/${strain}.\$chr.contigs.bam \$chr:\$MID_POINT-\$MAX_POINT -b > ${strain}.h2.\$chr.contigs.bam

            samtools index -@ ${taskCpus} ${strain}.h1.\$chr.contigs.bam
            samtools index -@ ${taskCpus} ${strain}.h2.\$chr.contigs.bam
        fi

        if [ ! -f ${igv_workdir}/${strain}/${strain}.\$chr.ilumina.bam ]
        then
            samtools view -b -@ ${taskCpus} ${igv_workdir}/${strain}/${strain}.ilumina.bam \$chr > ${strain}.\$chr.ilumina.bam
            samtools index -@ ${taskCpus} ${strain}.\$chr.ilumina.bam
        fi

        if [ ! -f ${igv_workdir}/${strain}/${strain}.\$chr.pacbio.bam ]
        then
            samtools view -b -@ ${taskCpus} ${igv_workdir}/${strain}/${strain}.pacbio.bam \$chr > ${strain}.\$chr.pacbio.bam
            samtools index -@ ${taskCpus} ${strain}.\$chr.pacbio.bam
        fi
        
        cp ${igv_workdir}/igv.session.template.xml ${strain}.\$chr.session.xml

        sed -i -e "s/STRAIN/${strain}/g" ${strain}.\$chr.session.xml
        sed -i -e "s/.pacbio.bam/.\$chr.pacbio.bam/g" ${strain}.\$chr.session.xml
        sed -i -e "s/.ilumina.bam/.\$chr.ilumina.bam/g" ${strain}.\$chr.session.xml
        sed -i -e "s/.contigs.bam/.\$chr.contigs.bam/g" ${strain}.\$chr.session.xml

    done < "chromosomes.txt"

    cat ${bed_file} > "${strain}.${type}.ready.bed"

    """
}

process take_screenshots {

    publishDir file(params.results + '/support-files/captures/'), mode: "copy", saveAs: {
                    filename -> filename.split('\\.')[0] + '/' + filename.replace(filename.split('\\.')[0] + ".","")
                }    

    input:
        tuple val(strain), val(type), file(bed_file)
        file igv_workdir
        file call_dir

    output:
        file "*.png"

    script:

    """
    bedToIgv -slop 50 -i ${bed_file} | grep -v snapshotDirectory >> snapshots.tmp
    bedToIgv -slop 500 -i ${bed_file} | grep -v snapshotDirectory >> snapshots.tmp
    bedToIgv -slop 5000 -i ${bed_file} | grep -v snapshotDirectory >> snapshots.tmp

    echo "snapshotDirectory ." >> snapshots.txt

    cat snapshots.tmp | awk '/^goto/ {s=\$0} /^snapshot/ {print s "\\t" \$0}' | sort -t \$'\\t' -k 1,2 | sed 's/\\t/\\n/' | awk  -F"[ :_]" 'BEGIN {prev=\$2; init=0} {if (init=0 || prev != \$2) print "load ${igv_workdir}/${strain}/${strain}."\$2".session.xml \\nload ${bed_file} \\nload ${call_dir}/${strain}-minigraph.DEL.bed \\nload ${call_dir}/${strain}-minigraph.INS.bed"; print; prev = \$2; init=1}' >> snapshots.txt
    
    echo "exit" >> snapshots.txt

    sed -i -e "s/snapshot /snapshot ${strain}.${type}./g" snapshots.txt

    echo "SKIP_VERSION=null,2.12.2" > prefs.properties
    echo "SHOW_SEQUENCE_TRANSLATION=true" > prefs.properties
    echo "IGV.Bounds=129,85,2243,1263" > prefs.properties
    echo "SAM.SHOW_SOFT_CLIPPED=true" > prefs.properties
    echo "IGV.Bounds=0,0,1920,1080" > prefs.properties
    echo "DETAILS_BEHAVIOR=CLICK" >> prefs.properties

    xvfb-run --auto-servernum -s "-screen 0 1920x1080x24" java -Xmx12000m --module-path=/home/mouse/IGV_Linux_2.12.2/lib --module=org.igv/org.broad.igv.ui.Main -b snapshots.txt -o prefs.properties

    """
}

def inputTuples(file) { 
    def strain = file.name.tokenize(".").get(0).tokenize("-").get(0)
    def type = file.name.tokenize(".").get(1)
    return tuple(strain, type, file)
}

workflow { 

    bed_files = Channel.fromPath(params.input_files).map{f -> inputTuples(f)}
    ready_files = prepare_screnshot_data(bed_files, file(params.igv_workdir))[0]
    take_screenshots(ready_files, file(params.igv_workdir), file(params.call_dir))

}

