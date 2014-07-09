console.time "native"
require("bamreader").create(process.argv[2]).fork(
  num: 1
  bam: (bam)->
    # bam.unmapped
    # bam.seq
    # bam.qual
    # bam.cigar
    # bam.tagstr
    #bam.sam
  finish: ->
    console.timeEnd "native"
    plain()
)

plain = ->
  console.time "plain"
  require("bamreader").create(process.argv[2]).fork(
    num: 1
    bam: (bam)->
    finish: ->
      console.timeEnd "plain"
  )

samtools = ->
  console.time "samtools"
  require("bamreader").create(process.argv[2]).fork_samtools(
    num: 1
    bam: (bam)->
      # bam.qname
      # bam.unmapped
      # bam.seq
      # bam.qual
      # bam.sam
      #console.log bam.tagstr
    finish: ->
      console.timeEnd "samtools"
  )
