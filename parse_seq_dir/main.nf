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

*/

def listFilesRecurse(_dir, _depth) {
  def files = []
  _dir.listFiles().each { _file -> 
      if (_file.isDirectory()) {
        if (_depth > 0) {
          listFilesRecurse(_file, _depth - 1).each { subfile ->
            files << subfile
          }
        }
      } else {
        files << _file
    }
  }
  return files
}

workflow PARSE_SEQ_DIR {
  take:
  inputDir // [meta, path_dir] where meta can contain depth

  main:

  inputDir.flatMap {
    def dir = it[1]
    def depth = 0
    if (it[0].depth) {
      depth = it[0].depth
    }
    def files = []
    def list2return = []
    def names = []
    def samplesToKeep = []
    // iterate through it.listFiles()
    listFilesRecurse(dir, depth).each { file ->
      names.add(file.getName())
      files.add(file)
    }
    if (it[0].samplesToKeep) {
      samplesToKeep = it[0].samplesToKeep
    }
    
    def layout = ""
    files.each { file ->
      def name = file.getName()
      def parent = file.getParent()
      layout = "SE"
      if ((match = name =~ /(?i)^(.+?)([\._]r?[12](.*?))?((\.fastq|\.fq|\.sfq)(\.)?(gz|gzip|z|zip|xz|bzip2|b2|bz2)?)$/)) {
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
          // Determine read type based on file extension
          def read_type = extension.toLowerCase().contains(".sfq") ? "sfq" : "fastq"
          list2return.add(tuple(layout + "_" + read_type, [sampleName, file]))
        }
      }
    }
    return list2return
  }
  | branch {
    single_fastq: it[0] == "SE_fastq"
      return tuple(["id":it[1][0], "read_type": "SE"], [it[1][1]])
    paired_fastq: it[0] == "PE_fastq"
      return it[1]
    single_sfq: it[0] == "SE_sfq"
      return tuple(["id":it[1][0], "read_type": "sfq"], [it[1][1]])
    paired_sfq: it[0] == "PE_sfq"
      return it[1]
  }
  | set {allFiles}

  // Handle paired PE files
  allFiles.paired_fastq
  | groupTuple(by: 0)
  | map { it ->
      def sorted = it[1].sort { a, b -> a.name <=> b.name }
      return tuple(["id":it[0], "read_type": "PE"], sorted)
  }
  | set {pairedFastqFiles}

  // Handle paired sfq files
  allFiles.paired_sfq
  | groupTuple(by: 0)
  | map { it ->
      def sorted = it[1].sort { a, b -> a.name <=> b.name }
      return tuple(["id":it[0], "read_type": "sfq"], sorted)
  }
  | set {pairedSfqFiles}

  // Combine all fastq files
  fastqFiles = pairedFastqFiles.concat(allFiles.single_fastq)

  // Combine all sfq files
  sfqFiles = pairedSfqFiles.concat(allFiles.single_sfq)

  emit:
    fastq = fastqFiles
    sfq = sfqFiles
}
