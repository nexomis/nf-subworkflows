
/*

Parse sequence directory.

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
    def names = []
    // iterate through it.listFiles()
    it.listFiles().each { file ->
      names.add(file.getName())
    }
    def layout = ""
    it.listFiles().each { file ->
      def name = file.getName()
      def parent = file.getParent()
      layout = "PE"
      if ((match = name =~ /(?i)^(.+?)([\._]r?[12](.*?))?(\.fastq|\.fq)(\.)?(gz|gzip|z|zip|xz|bzip2|b2|bz2)?$/)) {
        def sampleName = match.group(1)
        def readIndicator = match.group(2)
        def laneIndicator = match.group(3)
        if (readIndicator && laneIndicator){
          readIndicator = readIndicator.replaceAll("${laneIndicator}\$", "").trim()
        }
        if (readIndicator && readIndicator.endsWith("1")) {
          if (names.contains(file.name.replace(readIndicator, readIndicator.replace("1", "2")))) {
          } else {
            sampleName = sampleName + readIndicator
            layout = "SE"
          }
        } else if (readIndicator && readIndicator.endsWith("2")){
          if (names.contains(file.name.replace(readIndicator, readIndicator.replace("2", "1")))) {
          } else {
            sampleName = sampleName + readIndicator
            layout = "SE"
          }
        } else {
          layout = "SE"
        }
        if (laneIndicator) {
          sampleName = sampleName + laneIndicator
        }
        files.add(tuple(layout, [sampleName, file]))
      } else if ((match = name =~ /(?i)^(.+?)\.spring$/)) {
        def sampleName = match.group(1)
        files.add(tuple("spring", [sampleName, file]))
      }
    }
    return files
  }
  | branch {
    single: it[0] == "SE"
      return it[1]
    paired: it[0] == "PE"
      return it[1]
    spring: it[0] == "spring"
      return it[1]
  }
  | set {allFiles}

  allFiles.paired
  | groupTuple(by: 0)
  | map { it ->
      def sorted = it[1].sort { a, b -> a.name <=> b.name }
      return tuple(it[0], sorted)
  }
  | set {pairedFiles}

  emit:
    single = allFiles.single
    paired = pairedFiles
    spring = allFiles.spring

}
