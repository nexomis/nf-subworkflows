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

include { listFilesRecurse; parseFilename; detectLayout; groupFilesBySample } from './functions'

workflow PARSE_SEQ_DIR {
  take:
  inputDir // [meta, path_dir] where meta can contain depth, samplesToKeep, and parsingArgs

  main:

  inputDir.flatMap {
    def dir = it[1]
    def meta = it[0]
    def depth = meta.depth ?: 0
    def samplesToKeep = meta.samplesToKeep ?: []
    def parsingArgs = meta.parsingArgs ?: [:]
    
    def files = []
    def parsedFiles = []
    
    // Collect all files
    listFilesRecurse(dir, depth).each { file ->
      files.add(file)
    }
    
    // Parse all filenames using configurable parsing
    files.each { file ->
      def filename = file.getName()
      def parsed = parseFilename(filename, parsingArgs)
      
      if (parsed) {
        // Add the actual file object to the parsed result
        parsed.file = file
        parsedFiles.add(parsed)
      }
    }
    
    // Detect layouts for all samples
    def layouts = detectLayout(parsedFiles)
    
    // Group files by sample and layout
    def groupedFiles = groupFilesBySample(parsedFiles, layouts)
    
    // Filter by samples to keep and return results
    def results = []
    groupedFiles.each { group ->
      if (samplesToKeep.size() == 0 || samplesToKeep.contains(group.sampleName)) {
        if (group.layout == "PE") {
          // For paired-end, return tuple with sample metadata and both files
          // Use original file objects to preserve S3 properties
          def fileList = group.files.collect { it.file }
          results.add(tuple(
            ["id": group.sampleName, "read_type": group.layout],
            fileList
          ))
        } else {
          // For single-end, return each file separately
          // Use original file objects to preserve S3 properties
          group.files.each { parsedFile ->
            results.add(tuple(
              ["id": group.sampleName, "read_type": group.layout],
              [parsedFile.file]
            ))
          }
        }
      }
    }
    
    return results
  }
  | branch {
    fastq: it[0].read_type == "PE" || it[0].read_type == "SE"
      return it
    sfq: it[1][0].name.toLowerCase().contains(".sfq")
      return tuple(["id": it[0].id, "read_type": "sfq"], it[1])
  }
  | set { allFiles }

  emit:
    fastq = allFiles.fastq
    sfq = allFiles.sfq
}
