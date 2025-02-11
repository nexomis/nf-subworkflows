include { SPRING_DECOMPRESS } from '../../../modules/process/spring/decompress/main.nf'
include { KALLISTO_QUANT } from '../../../modules/process/kallisto/quant/main.nf'
include { checkMeta } from '../utils.nf'

workflow RNA_PREPROCESSING {
    take:
    inputs          // channel: [ meta, reads ]
    kallistoIdx    // channel: [meta, path(kallisto_idx)]
    multiqcYml     // channel: path(multiqc_yml)

    main:

    def expectedMeta = [
        "kallisto_idx": ["String", "NullObject"]
    ]


    inputs.map { checkMeta(it, expectedMeta) }
    kallistoIdx.map { checkMeta(it) }

    inputs
    | map { [it[0].kallisto_idx, it] }
    | combine(kallistoIdx.map({ [it[0].id, it] }), by:0)
    | map {
        def reads = it[1]
        def kallisto_idx = it[2]
        return [reads, kallisto_idx]
    }
    | set { kallistoInput }    

    // Kallisto quantification
    KALLISTO_QUANT(kallistoInput.map{ it[0] }, kallistoInput.map{ it[1] })

    // MultiQC report
    KALLISTO_QUANT.out.log
    | map { it[1] }
    | collect
    | set { kallistoLogs }

    ALIGN_MULTIQC(kallistoLogs, multiqcYml)

    emit:
    kallisto_h5 = KALLISTO_QUANT.out.h5    // channel: [ meta, quant ]
    kallisto_log = KALLISTO_QUANT.out.log        // channel: [ meta, log ]
    multiqc = ALIGN_MULTIQC.out                  // channel: path(html)
}

process ALIGN_MULTIQC {
    container 'multiqc/multiqc:v1.21' // version above are bugged

    label 'cpu_x1'
    label 'mem_8G'

    input:
    path kallisto_logs, stageAs: 'kallisto_logs/*', arity: '1..*'
    path conf_yml, arity: 1, stageAs: "multiqc_config.yaml"

    output:
    path("align_multiqc.html", arity: 1)

    script:
    """
    #!/usr/bin/bash

    multiqc -c $conf_yml --no-data-dir -n align_multiqc.html .
    """

    stub:
    """
    #!/usr/bin/bash

    touch primary_multiqc.html
    """
}
