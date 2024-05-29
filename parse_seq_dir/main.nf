/*

Parse sequence directory.

The Inputdir can be either 
(directory)
OR
(directory, sample names) in that case only reads matching sample names are
  treated

A subworkflow without process.
It takes directory as input

It gives list of sequences files with sample name as output:
paired: [sample_name [file R1, file R2]]
single: [sample_name [file R1]]
spring: [sample_name [file spring]]

*/

workflow PARSE_SEQ_DIR {
  take:
  inputDir

  main:

  inputDir.flatMap {
    def files = []
    def list2return = []
    def names = []
    def samplesToKeep = []
    // iterate through it.listFiles()
    if (it instanceof Path) {
      it.listFiles().each { file ->
        names.add(file.getName())
        files.add(file)
      }
    } else {
      it[0].listFiles().each { file ->
        names.add(file.getName())
        files.add(file)
      }
      samplesToKeep = it[1]
    }
    

    def layout = ""
    files.each { file ->
      def name = file.getName()
      def parent = file.getParent()
      layout = "SE"
      if ((match = name =~ /(?i)^(.+?)([\._]r?[12](.*?))?((\.fastq|\.fq)(\.)?(gz|gzip|z|zip|xz|bzip2|b2|bz2)?)$/)) {
        def sampleName = match.group(1)
        def readIndicator = match.group(2)
        def laneIndicator = match.group(3) // maybe confused with readIndicator and is fixed below (weird names such as s_R1_R1)
        def extension = match.group(4)
        def peSampleName = match.group(1)
        def seSampleName = match.group(1)
        if (readIndicator) {
          seSampleName = seSampleName + readIndicator
        }
        if (laneIndicator) {
          if ( mlane = readIndicator =~ /(?i)^(.*)([\._]r?[12])(.*?)$/ ) {
            sampleName = sampleName + mlane.group(1)
            peSampleName = peSampleName + mlane.group(1)
            readIndicator = mlane.group(2)
            laneIndicator = mlane.group(3)
            if (laneIndicator) {
              peSampleName = peSampleName + laneIndicator
            }
          }
        }
        if (readIndicator) {
          def rmate = readIndicator
          if (readIndicator.endsWith("1")) {
            rmate = rmate.replace("1", "2")
          } else {
            rmate = rmate.replace("2", "1")
          }
          def mate_name = sampleName + rmate + extension
          if (laneIndicator) {
            mate_name = sampleName + rmate + laneIndicator + extension
          }
          if (names.contains(mate_name)) {
            layout = "PE"
          }
        }
        if (layout == "PE") {
          sampleName = peSampleName
        } else {
          sampleName = seSampleName
        }
        if (samplesToKeep.size() == 0 || samplesToKeep.contains(sampleName)) {
          list2return.add(tuple(layout, [sampleName, file]))
        }
      } else if ((match = name =~ /(?i)^(.+?)\.spring$/)) {
        def sampleName = match.group(1)
        if (samplesToKeep.size() == 0 || samplesToKeep.contains(sampleName)) {
          list2return.add(tuple("spring", [sampleName, file]))
        }
      }
    }
    return list2return
  }
  | branch {
    single: it[0] == "SE"
      return tuple(it[1][0], [it[1][1]])
    paired: it[0] == "PE"
      return it[1]
    spring: it[0] == "spring"
      return tuple(it[1][0], [it[1][1]])
  }
  | set {allFiles}

  allFiles.paired
  | groupTuple(by: 0)
  | map { it ->
      def sorted = it[1].sort { a, b -> a.name <=> b.name }
      return tuple(it[0], sorted)
  }
  | set {pairedFiles}

  fastqFiles = pairedFiles.concat(allFiles.single)

  emit:
    fastq = fastqFiles
    spring = allFiles.spring

}
