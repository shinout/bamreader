require("bamreader").create(process.argv[2]).samtools (bam)->
  #console.log bam.rname

require("bamreader").read process.argv[2], (bam)->
  console.log bam.tlen
