console.time "native"
require("bamreader").create(process.argv[2]).fork(
  num: 1
  bam: (bam)->
    bam.qname
    bam.unmapped
    bam.seq
    # bam.qual
    # bam.sam
  finish: ->
    console.timeEnd "native"
    samtools()
)

samtools = ->
  console.time "samtools"
  require("bamreader").create(process.argv[2]).fork_samtools(
    num: 1
    bam: (bam)->
      bam.qname
      bam.unmapped
      bam.seq
      bam.qual
      bam.sam
    finish: ->
      console.timeEnd "samtools"
  )
