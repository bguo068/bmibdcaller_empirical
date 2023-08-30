nextflow.enable.dsl=2

params.test = false
params.resdir = "res"

params.min_mac = 20 // assuming homozyogous diploids, eg. min_maf 20 / 2 / 1000 = 0.01
params.mincm = 2.0
params.nchroms = 14

params.hapibd_minoutput = params.mincm
params.hapibd_minseed = params.mincm
params.hapibd_minextend = 1.0
params.hapibd_maxgap = 1000
params.hapibd_minmarkers = 70

params.chrlens = [
    1: 640851,
    2: 947102,
    3: 1067971,
    4: 1200490,
    5: 1343557,
    6: 1418242,
    7: 1445207,
    8: 1472805,
    9: 1541735,
    10: 1687656,
    11: 2038340,
    12: 2271494,
    13: 2925236,
    14: 3291936,
]
def resdir = params.resdir

process CALL_IBD_HAPIBD {
    tag "${args.genome_set_id}_${chrno}_hapibd"
    publishDir "${resdir}/${label}/ibd/hapibd/", \
        mode: "symlink", pattern: "*_hapibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf)
    output:
    tuple val(label), val(chrno), path("*_hapibd.ibd")
    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        chrno: chrno,
        seqlen: args.seqlen,
        minseed: params.hapibd_minseed,
        minoutput: params.hapibd_minoutput,
        maxgap: params.hapibd_maxgap,
        minextend: params.hapibd_minextend,
        minmarkers: params.hapibd_minmarkers,
        minmac: params.min_mac,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_hapibd.py $cmd_options
    """

    stub:
    """touch ${args.genome_set_id}_${chrno}_hapibd.ibd"""
}

workflow {
    ch_input = channel.fromList(1..14).map{chrno->
        def label = "alldom09imp"
        def args = [ 
            r: 0.01/15000, 
            seqlen: params.chrlens[chrno], 
            genome_set_id: 0 
        ]
        def vcf = file("${projectDir}/../../input/impute/pf7_dom_0_9_imp_chr${chrno}.vcf.gz", checkIfExists:true)
        [label, chrno, args, vcf]
    } 
    CALL_IBD_HAPIBD(ch_input)
}