inputFiles = Channel.empty()

if( params.inputFileType == "bam" ) {
inputFiles = Channel.fromPath(params.inputFiles).map { file -> [file.getSimpleName(), file.getBaseName(), file, file + ".bai"] }
}
else if( params.inputFileType == "cram" ) {
inputFiles = Channel.fromPath(params.inputFiles).map { file -> [file.getSimpleName(), file.getBaseName(), file, file + ".crai"] }
}
