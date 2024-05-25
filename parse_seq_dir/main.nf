workflow PARSE_SEQ_DIR {
  take:
  Path inputDir // Check that inputDir is a Value Channel ?

  main:
  Channel.fromPath("${inputDir}/**.{fastq,fq}{,.gz}")
    | map { file ->
      def basename = file.baseName
      def isPairedEnd = basename =~ /.*[._][Rr](1|2).*/
      def sampleName = (isPairedEnd ? basename.replaceFirst(/[._][Rr](1|2).*$/, "") : basename)
      return tuple(sampleName, file, isPairedEnd)
    }
    | branch {
      paired: it[2]
      single: !it[2]
    }
    | set { fileStreams }

  fileStreams.paired
    | groupTuple(by: [0])
    | map { sampleName, files ->
      files.sort() // Ensure files are sorted by name which typically aligns R1 and R2
      return tuple(sampleName, files)
    }
    | set { pairedInput }

  fileStreams.single
    | map { sampleName, file ->
      return tuple(sampleName, file)
    }
    | set { singleInput }

  emit:
    single: singleInput
    paired: pairedDInput

}