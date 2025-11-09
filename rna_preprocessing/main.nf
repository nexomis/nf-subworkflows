include { SPRING_DECOMPRESS } from '../../../modules/process/spring/decompress/main.nf'
include { KALLISTO_QUANT } from '../../../modules/process/kallisto/quant/main.nf'
include { SALMON_QUANT } from '../../../modules/process/salmon/quant/main.nf'
include { checkMeta } from '../utils.nf'

workflow RNA_PREPROCESSING {
    take:
    inputs          // channel: [ meta, reads ]
    reference       // channel: [meta, path(reference)]

    main:

    def expectedMeta = [
        "reference": ["String", "NullObject"]
    ]

    inputs.map { it -> checkMeta(it, expectedMeta) }
    reference.map { it -> checkMeta(it) }


    channel.fromPath(moduleDir + "/multiqc.yml")
    | set {multiqcYml}

    inputs
    | map { it ->  [it[0].reference, it] }
    | combine(reference.map({ it -> [it[0].id, it] }), by:0)
    | map { it -> 
        def reads = it[1]
        def ref = it[2]
        return [reads, ref]
    }
    | set { quantInput }    

    // Branch by method
    quantInput
    | branch { it -> 
        kallisto: it[1][0].method == 'kallisto'
        salmon: it[1][0].method == 'salmon'
    }
    | set { quantBranched }

    // Kallisto quantification
    KALLISTO_QUANT(
        quantBranched.kallisto.map({ it -> it[0] }), 
        quantBranched.kallisto.map({ it -> it[1] })
    )

    // Salmon quantification
    SALMON_QUANT(
        quantBranched.salmon.map({ it -> it[0] }), 
        quantBranched.salmon.map({ it -> it[1] })
    )

    // Prepare MultiQC inputs
    kallistoLogs = KALLISTO_QUANT.out.log

    SALMON_QUANT.out.log
        | concat(SALMON_QUANT.out.flenDist)
        | concat(SALMON_QUANT.out.lib_count)
        | set { salmonLogs }

    ALIGN_MULTIQC(
        kallistoLogs.map({ it -> it[1] }).collect(), 
        salmonLogs.map({ it -> it[1] }).collect(), 
        multiqcYml
    )

    emit:
    kallisto_h5 = KALLISTO_QUANT.out.h5     // channel: [ meta, quant ]
    kallisto_log = kallistoLogs
    salmon_logs = salmonLogs
    salmon_quant = SALMON_QUANT.out.quant   // channel: [ meta, quant ]
    multiqc = ALIGN_MULTIQC.out             // channel: path(html)
}

process ALIGN_MULTIQC {
    container 'multiqc/multiqc:v1.27.1' // version above are bugged

    cpus 1
    memory 8.GB

    input:
    path kallisto_logs, stageAs: 'kallisto_logs/*', arity: '0..*'
    path salmon_logs, stageAs: 'salmon_logs/*', arity: '0..*'
    path conf_yml, arity: 1, stageAs: "multiqc_config.yaml"

    output:
    path("align_multiqc.html", arity: 1)

    script:
    """
    #!/usr/bin/bash

    mkdir -p salmon_logs/aux_info
    mkdir -p salmon_logs/libParams
    
    for file in \$(ls salmon_logs/ | grep meta_info.json); do
        ln -s \$PWD/salmon_logs/\$file \$PWD/salmon_logs/aux_info/\$file
    done
    for file in \$(ls salmon_logs | grep flenDist.txt); do
        ln -s \$PWD/salmon_logs/\$file \$PWD/salmon_logs/libParams/\$file
    done

    multiqc -c $conf_yml --no-data-dir -n align_multiqc.html .
    """

    stub:
    """
    #!/usr/bin/bash

    touch primary_multiqc.html
    """
}
