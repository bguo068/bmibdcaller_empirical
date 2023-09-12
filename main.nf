nextflow.enable.dsl=2

// /////////////////////////////////////////////////
// define default parameters for the pipeline
// many of which can be overwritten by command line

params.test = false
params.resdir = "res"
params.recombination_rate = 0.01 / 15000

params.min_mac = 20 // assuming homozyogous diploids, eg. min_maf 20 / 2 / 1000 = 0.01
params.mincm = 2.0
params.nchroms = 14

params.tpbwt_template_opts = 1
params.tpbwt_Lm = params.test ? 152: 100 // optimized
params.tpbwt_Lf = params.mincm
params.tpbwt_use_phase_correction = 0

params.hapibd_minoutput = params.mincm
params.hapibd_minseed = params.mincm
params.hapibd_minextend = 1.0
params.hapibd_maxgap = 1000
params.hapibd_minmarkers = 70

params.refinedibd_length = params.mincm
params.refinedibd_lod = 1.1 // TODO: confirm what does lod mean
params.refinedibd_scale = 0
params.refinedibd_window = 40.0

params.hmmibd_n = 100
params.hmmibd_m = 5

params.isorelate_imiss = 0.3
params.isorelate_vmiss = 0.3
params.isorelate_min_snp = 20 // optimized
params.isorelate_min_mac = 200 // 0.1, which is different from other callers

// params.tsinferibd_max_tmrca = [1000, 3000]

params.filt_ibd_by_ov = false
params.ibdne_mincm = params.mincm
params.ibdne_minregion = 10
params.ibdne_flatmeth = ["none", "keep_hap_1_only", "merge"][0]

params.ifm_transform = ["square", "cube", "none"][0]
params.ifm_ntrials = 1000
params.ifm_mincm = 2.0
params.ifm_mingwcm = 5.0


params.meta = "" // empty for simulation

def resdir = params.resdir

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


// /////////////////////////////////////////////////
// define individual processes

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

process CALL_IBD_TSKIBD {
    tag "${args.genome_set_id}_${chrno}_tskibd"
    publishDir "${resdir}/${label}/ibd/tskibd/", \
        mode: "symlink", pattern: "*_tskibd.ibd"
    input:
    tuple val (label), val(chrno), val(args), path(trees)
    output:
    tuple val(label), val(chrno), path("*_tskibd.ibd")
    script:
    def cmd_options_tskibd = [
        tree: trees,
        chrno: chrno,
        r: args.r,
        seqlen: args.seqlen,
        mincm: params.mincm,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")

    """
    call_tskibd.py ${cmd_options_tskibd}
    """
    stub:
    def prefix = "${args.genome_set_id}_${chrno}"
    """touch  ${prefix}_tskibd.ibd"""
}
process CALL_IBD_REFINEDIBD {
    tag "${args.genome_set_id}_${chrno}_refinedibd"
    publishDir "${resdir}/${label}/ibd/refinedibd/", \
        mode: "symlink", pattern: "*_refinedibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf)
    output:
    tuple val(label), val(chrno), path("*_refinedibd.ibd")
    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        lod: params.refinedibd_lod,
        length: params.refinedibd_length,
        scale: params.refinedibd_scale,
        mem_gb: task.memory.giga,
        minmac: params.min_mac,
        window: params.refinedibd_window,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_refinedibd.py ${cmd_options}
    """
    stub:
    """touch ${args.genome_set_id}_${chrno}_refinedibd.ibd"""
}

process CALL_IBD_TPBWT {
    tag "${args.genome_set_id}_${chrno}_tpbwtibd"
    publishDir "${resdir}/${label}/ibd/tpbwtibd/", \
        mode: "symlink", pattern: "*_tpbwtibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf)
    output:
    tuple val(label), val(chrno), path("*_tpbwtibd.ibd")
    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        template: params.tpbwt_template_opts, 
        use_phase_correction: params.tpbwt_use_phase_correction,
        minmac: params.min_mac,
        Lm: params.tpbwt_Lm,
        Lf: params.tpbwt_Lf,
        mem_gb: task.memory.giga,
        nthreads: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_tpbwt.py ${cmd_options}
    """
    stub:
    """touch ${args.genome_set_id}_${chrno}_tpbwtibd.ibd"""
}

process CALL_IBD_HMMIBD {
    tag "${args.genome_set_id}_${chrno}_hmmibd"
    publishDir "${resdir}/${label}/ibd/hmmibd/", \
        mode: "symlink", pattern: "*_hmmibd.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf)
    output:
    tuple val(label), val(chrno), path("*_hmmibd.ibd")
    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        n: params.hmmibd_n,
        m: params.hmmibd_m,
        mincm: params.mincm,
        minmac: params.min_mac,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    script:
        """
        call_hmmibd.py ${cmd_options}
        """
    stub:
    """touch ${args.genome_set_id}_${chrno}_hmmibd.ibd"""
}

process CALL_IBD_ISORELATE {
    tag "${args.genome_set_id}_${chrno}_isorelate"
    publishDir "${resdir}/${label}/ibd/isorelate/", \
        mode: "symlink", pattern: "*_isorelate.ibd"
    input:
    tuple val(label), val(chrno), val(args), path(vcf)
    output:
    tuple val(label), val(chrno), path("*_isorelate.ibd")
    script:
    def cmd_options = [
        vcf: vcf,
        r: args.r,
        seqlen: args.seqlen,
        chrno: chrno,
        min_snp: params.isorelate_min_snp,
        min_len_bp: Math.round(params.mincm * (0.01/args.r)),
        minmac: params.min_mac,
        imiss: params.isorelate_imiss,
        vmiss: params.isorelate_vmiss,
        cpus: task.cpus,
        genome_set_id: args.genome_set_id,
    ].collect{k, v -> "--${k} ${v}"}.join(" ")
    """
    call_isorelate.py ${cmd_options}
    """
    stub:
    """touch ${args.genome_set_id}_${chrno}_isorelate.ibd"""
}

process PROC_DIST_NE {
    tag "${genome_set_id}"

    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ne_input/", \
        pattern: "*.sh", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ne_input/",\
         pattern: "*.map", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ne_input/",\
         pattern: "*.ibd.gz", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ibddist_ibd/",\
         pattern: "*.ibddist.ibdobj.gz", mode: 'symlink'
    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ibdne_ibd/",\
         pattern: "*.ibdne.ibdobj.gz", mode: 'symlink'

    input:
        tuple val(label),val(ibdcaller), path(ibd_lst), path(ibd_lst_true), path(vcf_lst), val(genome_set_id)
        path(ibdne_jar)
    output:
        tuple val(label),val(ibdcaller),  path("ibdne.jar"), path("*_orig.sh"), \
                path("*_orig.map"), path("*_orig.ibd.gz"), emit: ne_input_orig
        tuple val(label),val(ibdcaller),  path("ibdne.jar"), path("*_rmpeaks.sh"),  \
                path("*_rmpeaks.map"), path("*_rmpeaks.ibd.gz"), emit: ne_input_rmpeaks
        tuple val(label),val(ibdcaller),  path("*.ibddist.ibdobj.gz"), emit: ibddist_ibd_obj
        tuple val(label),val(ibdcaller),  path("*.ibdne.ibdobj.gz"), emit: ibdne_ibd_obj
    script:
    // Whether to pass in true ibd list is determined by params.filt_ibd_by_ov 
    def true_ibd_arg = (params.filt_ibd_by_ov && (ibdcaller!="hmmibd")) ? \
         [ibd_files_true: "${ibd_lst_true}"] : [:]
    def args_local = (true_ibd_arg + [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        vcf_files: "${vcf_lst}", // path is a blank separate list
        genome_set_id: genome_set_id,
        ibdne_mincm: params.ibdne_mincm,
        ibdne_minregion: params.ibdne_minregion,
        ibdne_jar: ibdne_jar,
        ibdne_flatmeth: params.ibdne_flatmeth,
    ]).collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    proc_dist_ne.py ${args_local} 
    """
    stub:
    """
    touch ibdne.jar
    touch ${genome_set_id}{_orig.sh,_orig.map,_orig.ibd.gz}
    touch ${genome_set_id}{_rmpeaks.sh,_rmpeaks.map,_rmpeaks.ibd.gz}
    touch ${genome_set_id}.ibddist.ibdobj.gz
    touch ${genome_set_id}_orig.ibdne.ibdobj.gz
    touch ${genome_set_id}_rmpeaks.ibdne.ibdobj.gz
    """
}

process PROC_INFOMAP {
    tag "${genome_set_id}"

    publishDir "${resdir}/${genome_set_id}_${label}/${ibdcaller}/ifm_input/", \
        pattern: "*.ibdobj.gz", mode: 'symlink'

    input:
        tuple val(label),val(ibdcaller),  path(ibd_lst), path(vcf_lst), val(genome_set_id)

    output:
        tuple val(label),val(ibdcaller),  path("*_orig.ifm.ibdobj.gz"), \
                    emit: ifm_orig_ibd_obj
        tuple val(label),val(ibdcaller),  path("*_rmpeaks.ifm.ibdobj.gz"), \
                    emit: ifm_rmpeaks_ibd_obj
    script:
    def args_local = [
        ibd_files: "${ibd_lst}", // path is a blank separate list
        vcf_files: "${vcf_lst}", // path is a blank separate list
        genome_set_id: genome_set_id,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    proc_infomap.py ${args_local}
    """
    stub:
    """
    touch ${genome_set_id}{_orig.ifm.ibdobj.gz,_rmpeaks.ifm.ibdobj.gz}
    """
}

process RUN_IBDNE {
    tag "${args.genome_set_id}_${are_peaks_removed}"

    publishDir "${resdir}/${args.genome_set_id}_${label}/${ibdcaller}/ne_output/",  mode: 'symlink'

    input:
        tuple val(label), val(ibdcaller), path(ibdne_jar), path(ibdne_sh), path(gmap),\
            path(ibd_gz), val(are_peaks_removed), val(args)
    output:
        tuple val(label), val(ibdcaller), val(are_peaks_removed), path("*.ne")
    script:
    """
    bash ${ibdne_sh}
    """
    stub:
    def src = are_peaks_removed ? "rmpeaks": "orig"
    """
    touch ${args.genome_set_id}_${src}.ne
    """
}

process RUN_INFOMAP {
    tag "${args.genome_set_id}_${are_peaks_removed}"
    publishDir "${resdir}/${args.genome_set_id}_${label}/${ibdcaller}/ifm_output/",  mode: 'symlink'
    input:
        tuple val(label), val(ibdcaller), path(ibd_obj), val(are_peaks_removed), val(args), val(name_map)
    output:
        tuple val(label), val(ibdcaller), val(are_peaks_removed), path("*_member.pq")
    script:
    def cut_mode = are_peaks_removed? 'rmpeaks': 'orig'
    def args_local = [
        ibd_obj: ibd_obj,
        // no using this for making meta file
        // npop: args.npop,
        // nsam: args.nsam,
        name_map: name_map,
        genome_set_id: args.genome_set_id,
        cut_mode: cut_mode,
        ntrials: params.ifm_ntrials,
        transform: params.ifm_transform,
    ].collect{k, v-> "--${k} ${v}"}.join(" ")
    """
    run_infomap.py ${args_local}
    """
    stub:
    def cut_mode = are_peaks_removed? 'rmpeaks': 'orig'
    """
    touch ${args.genome_set_id}_${cut_mode}_member.pq
    """
}


// Pipeline for single-population datasets
workflow WF_SP{
main:

    // *********************** input channel *************************
    ch_sets = Channel.fromList([
        [label: "singlepop_AF-W_Ghana_16_18", genome_set_id: 0],
        [label:"singlepop_AS-SE-E_10_12",     genome_set_id: 1],
    ])
    chr_chrnos = Channel.fromList(1..(params.nchroms))

    ch_in_ibdcall_vcf = ch_sets.combine(chr_chrnos).map{d, chrno ->
        def label = d.label
        def genome_set_id = d.genome_set_id
        def args = [r: params.recombination_rate, seqlen: params.chrlens[chrno], genome_set_id: genome_set_id]
        def vcf = file("${projectDir}/datasets/02_make_dataset/vcf/${label}/chr${chrno}.vcf.gz", checkIfExists: true)
        return [label, chrno, args, vcf]
    }

    // *********************** Call ibd *************************

    CALL_IBD_HAPIBD(ch_in_ibdcall_vcf)
    // CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_REFINEDIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TPBWT(ch_in_ibdcall_vcf)
    CALL_IBD_HMMIBD(ch_in_ibdcall_vcf)
    CALL_IBD_ISORELATE(ch_in_ibdcall_vcf)

    // collect IBD and groupby simulation label and ibdcaller
    ch_out_ibd_grp = CALL_IBD_HAPIBD.out.map{it + ["hapibd"]}
            // .concat( CALL_IBD_TSKIBD.out.map{it + ["tskibd"]})
            .concat(CALL_IBD_REFINEDIBD.out.map{it+ ["refinedibd"]})
            .concat(CALL_IBD_TPBWT.out.map{it + ["tpbwt"]})
            .concat(CALL_IBD_HMMIBD.out.map{it + ["hmmibd"]})
            .concat(CALL_IBD_ISORELATE.out.map{it+["isorelate"]})
            // [label, chrno, ibd, ibdcaller]
        .map{label, chrno, ibd, ibdcaller ->
            def key =  groupKey([label, ibdcaller], params.nchroms)
            def data = [chrno, ibd]
            return [key, data]}
        .groupTuple(by: 0, sort: {a, b -> a[0] <=> b[0]}) //sort by chrno
        .map {key, data_lst ->
            def label = key[0]
            def ibdcaller = key[1]
            def ibd_lst = data_lst.collect{data -> data[1]}
            return ["${label}", ibdcaller, ibd_lst]
        }

    // TODO: will come back to consider remvoing IBD segment overlapping
    // with IBD segments < 1.5 as determined by hmmibd

    ch_true_ibd = CALL_IBD_HMMIBD.out // [label, chrno, ibd]
        .map{label, chrno, ibd-> [groupKey(label, params.nchroms), [chrno, ibd]]}
        .groupTuple(by: 0, sort: {a, b-> a[0]<=>b[0]})  // groupby label, sort by chrno
        .map{key, data_lst -> 
            def label = key.toString()
            def ibd_files_true = data_lst.collect{v -> v[1]}
            [label, ibd_files_true]
        }

    ch_out_ibd_grp = ch_out_ibd_grp
        // combine with true ibd list
        .combine(ch_true_ibd, by: 0)
        .map{label, ibdcaller, ibd_lst, ibd_list_true->
            def vcf_lst = (1..params.nchroms).collect{chrno -> 
                file("${projectDir}/datasets/02_make_dataset/vcf/${label}/chr${chrno}.vcf.gz", checkIfExists: true)
            }
            // fix name collision for tskibd as true and inferrre ibd are the same
            [label, ibdcaller, ibd_lst, ibdcaller=="hmmibd" ? [] : ibd_list_true, vcf_lst]
         }


    ch_grouped_ibd_vcf = ch_out_ibd_grp
        // [label, ibdcaller, ibd_lst, ibd_lst_true, vcf_lst]
        // add genome_set_id
        .combine(ch_sets.map{d->[d.label, d.genome_set_id]}, by: 0)
        // [label, ibdcaller, ibd_lst, ibd_lst_true, vcf_lst, genome_set_id]
        //                                      ^^^^^^^^^^^^^^

    // ********************** Process ibd ***************************
    // NOTE: although the ibd_files_true is passed, but it can still not used 
    // depending the params.filt_ibd_by_ov = True
    PROC_DIST_NE( ch_grouped_ibd_vcf, file("${projectDir}/lib/ibdne.23Apr20.ae9.jar") )
     // out.ne_input_orig
     // out.ne_input_rmpeaks



    // ********************** Run IbdNe ***************************
    ch_in_ibdne = \
            PROC_DIST_NE.out.ne_input_orig.map{it + [false]}  // orig
        .concat( 
            PROC_DIST_NE.out.ne_input_rmpeaks.map {it + [true]}  // rmpeaks
        )
        // label, ibdcaller, jar, sh, map, ibd, are_peaks_removed
        .combine( 
            // [label, args] // args need to only contains genome_set_id
            ch_sets.map{d-> 
                def label = d.label
                def args = [genome_set_id: d.genome_set_id]
                return [label, args]
            },
            by: 0) // add args

    RUN_IBDNE(ch_in_ibdne)

    emit: 

    ch_ibdobj_dist = PROC_DIST_NE.out.ibddist_ibd_obj
    ch_ibdobj_ne   = PROC_DIST_NE.out.ibdne_ibd_obj
    ch_ibdne       = RUN_IBDNE.out
            // label, ibdcaller, are_peaks_removed, ne
}

// Pipeline for multi-population datasets
workflow WF_MP{
    main:

    // *********************** input channel *************************
    ch_sets = Channel.fromList([
        [
            label: "structured", 
            genome_set_id: 100, 
            name_map: file("${projectDir}/datasets/02_make_dataset/vcf/structured/sample_name_map_structured.txt", checkIfExists:true)
        ],
        [
            label:"singlepop_AS-SE-E_10_12", 
            genome_set_id: 1, 
            name_map: file("${projectDir}/datasets/02_make_dataset/vcf/singlepop_AS-SE-E_10_12/sample_name_map_singlepop_AS-SE-E_10_12.txt", checkIfExists:true)
        ],
    ])

    chr_chrnos = Channel.fromList(1..(params.nchroms))

    ch_in_ibdcall_vcf = ch_sets.combine(chr_chrnos).map{d, chrno ->
        def label = d.label
        def genome_set_id = d.genome_set_id
        def args = [r: params.recombination_rate, seqlen: params.chrlens[chrno], genome_set_id: genome_set_id]
        def vcf = file("${projectDir}/datasets/02_make_dataset/vcf/${label}/chr${chrno}.vcf.gz", checkIfExists: true)
        return [label, chrno, args, vcf]
    }

    // CALL_IBD_TSINFERIBD(ch_in_ibdcall_vcf_with_ne)
    CALL_IBD_HAPIBD(ch_in_ibdcall_vcf)
    // CALL_IBD_TSKIBD(ch_in_ibdcall_trees)
    CALL_IBD_REFINEDIBD(ch_in_ibdcall_vcf)
    CALL_IBD_TPBWT(ch_in_ibdcall_vcf)
    CALL_IBD_HMMIBD(ch_in_ibdcall_vcf)
    CALL_IBD_ISORELATE(ch_in_ibdcall_vcf)

    // collect IBD and groupby simulation label and ibdcaller
    ch_out_ibd_grp = CALL_IBD_HAPIBD.out.map{it + ["hapibd"]}
            // .concat( CALL_IBD_TSKIBD.out.map{it + ["tskibd"]})
            .concat(CALL_IBD_REFINEDIBD.out.map{it+ ["refinedibd"]})
            .concat(CALL_IBD_TPBWT.out.map{it + ["tpbwt"]})
            .concat(CALL_IBD_HMMIBD.out.map{it + ["hmmibd"]})
            .concat(CALL_IBD_ISORELATE.out.map{it+["isorelate"]})
            // [label, chrno, ibd, ibdcaller]
        .map{label, chrno, ibd, ibdcaller ->
            def key =  groupKey([label, ibdcaller], params.nchroms)
            def data = [chrno, ibd]
            return [key, data]}
        .groupTuple(by: 0, sort: {a, b -> a[0] <=> b[0]}) //sort by chrno
        .map {key, data_lst ->
            def label = key[0]
            def ibdcaller = key[1]
            def ibd_lst = data_lst.collect{data -> data[1]}
            def vcf_lst = (1..params.nchroms).collect{chrno -> 
                file("${projectDir}/datasets/02_make_dataset/vcf/${label}/chr${chrno}.vcf.gz", checkIfExists: true)
            }
            return [label, ibdcaller, ibd_lst, vcf_lst]
        }

  
    ch_grouped_ibd_vcf = ch_out_ibd_grp
        // [label, ibdcaller, ibd_lst, vcf_lst]
        // add genome_set_id
        .combine(ch_sets.map{d->[d.label, d.genome_set_id]}, by: 0)
        // [label, ibdcaller, ibd_lst, vcf_lst, genome_set_id]

    // ********************** Process ibd ***************************

    PROC_INFOMAP(ch_grouped_ibd_vcf)
    // out.ifm_orig_ibd_obj       // label, ibdcaller, ibd_obj
    // out.ifm_rmpeaks_ibd_obj    // label, ibdcaller, ibd_obj


    // ********************** Run Infomap ***************************
    ch_in_run_infomap = \
        PROC_INFOMAP.out.ifm_orig_ibd_obj.map{it + [false]} // orig
    .concat (
        PROC_INFOMAP.out.ifm_rmpeaks_ibd_obj.map{it + [true]} // rmpeaks
    ).combine(
            // [label, args] // args need to only contains genome_set_id
            ch_sets.map{d-> 
                def label = d.label
                def args = [genome_set_id: d.genome_set_id]
                def name_map = d.name_map
                return [label, args, name_map]
            },
        by: 0)
        // label, ibdcaller, ibd_obj, are_peaks_removed, args, name_map

    RUN_INFOMAP(ch_in_run_infomap)


    emit:

    ch_ifm = RUN_INFOMAP.out
            // label, ibdcaller, are_peaks_removed, member
}

workflow {
    WF_SP()
    WF_MP()
}