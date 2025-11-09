/*
 * Utility functions for subworkflows
 */

def listFilesRecurse(_dir, _depth) {
  def files = []
  
  // Handle both local paths and S3 paths using Nextflow's file() function
  def dir = file(_dir)
  
  if (dir.exists() && dir.isDirectory()) {
    dir.listFiles().each { _file -> 
        if (_file.isDirectory()) {
          if (_depth > 0) {
            listFilesRecurse(_file, _depth - 1).each { subfile ->
              files << subfile
            }
          }
        } else {
          files = files + [_file]
        }
    }
  }
  return files
}

/*
 * Parse filename to extract components using configurable regex patterns
 * Returns a map with: sampleName, readNumber, fileType, originalFilename
 */
def parseFilename(filename, parsingArgs = [:]) {
  // Default parsing arguments
  def defaultArgs = [
    tailRegex: /_\d{3}$/,
    readRegex: /_R([12])$/,
    laneRegex: /_S\d+_L\d+$/
  ]
  
  // Merge with provided arguments
  def args = defaultArgs + parsingArgs
  
  def workingName = filename
  def fileType = "fastq"
  def readNumber = 1
  
  // Step 1: Extract file extension and determine file type
  def extensionMatch = workingName =~ /(?i)((\.fastq|\.fq|\.sfq)(\.)?(gz|gzip|z|zip|xz|bzip2|b2|bz2)?)$/
  if (extensionMatch) {
    def extension = extensionMatch.group(1)
    fileType = extension.toLowerCase().contains(".sfq") ? "sfq" : "fastq"
    workingName = workingName.substring(0, extensionMatch.start())
  } else {
    // If no valid sequence file extension is found, return null to filter out the file
    return null
  }
  
  // Step 2: Remove eventual tail (default: "_\d{3}$")
  if (args.tailRegex) {
    workingName = workingName.replaceAll(args.tailRegex, "")
  }
  
  // Step 3: Look for read ID and remove it (default: "_R([12])$")
  if (args.readRegex) {
    def readMatch = workingName =~ args.readRegex
    if (readMatch) {
      readNumber = readMatch.group(1) as Integer
      workingName = workingName.substring(0, readMatch.start())
    }
  }
  
  // Step 4: Remove eventual lane information (default: "_S\d+_L\d+$")
  if (args.laneRegex) {
    workingName = workingName.replaceAll(args.laneRegex, "")
  }
  
  // Step 5: The remaining string is the sample name
  def sampleName = workingName
  
  return [
    sampleName: sampleName,
    readNumber: readNumber,
    fileType: fileType,
    originalFilename: filename
  ]
}

/*
 * Detect if files are paired-end or single-end based on parsed files
 */
def detectLayout(parsedFiles) {
  // Group files by sample name
  def sampleGroups = parsedFiles.groupBy { it -> it.sampleName }
  
  def results = [:]
  sampleGroups.each { sampleName, files ->
    def readNumbers = files.collect { it -> it.readNumber }.unique().sort()
    if (readNumbers.size() == 2 && readNumbers == [1, 2]) {
      results[sampleName] = "PE"
    } else {
      results[sampleName] = "SE"
    }
  }
  
  return results
}

/*
 * Group parsed files by sample and layout
 */
def groupFilesBySample(parsedFiles, layouts) {
  def grouped = []
  
  def sampleGroups = parsedFiles.groupBy { it -> it.sampleName }
  
  sampleGroups.each { sampleName, files ->
    def layout = layouts[sampleName]
    def fileType = files[0].fileType
    
    if (layout == "PE") {
      // Sort files by read number for paired-end
      def sortedFiles = files.sort { it -> it.readNumber }
      grouped.add([
        layout: "PE",
        fileType: fileType,
        sampleName: sampleName,
        files: sortedFiles
      ])
    } else {
      // For single-end, add each file separately
      files.each { file ->
        grouped.add([
          layout: "SE", 
          fileType: fileType,
          sampleName: sampleName,
          files: [file]
        ])
      }
    }
  }
  
  return grouped
}
