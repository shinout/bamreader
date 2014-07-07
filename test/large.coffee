BAMReader = require(__dirname + "/../lib/bamreader.js")
n = 0

BAMReader.create(__dirname + "/large.bam").on "bam", (bam, dOffset, iOffset)->
  n++
  return if n isnt 1
  console.assert bam.qname is 'HS2000-903_160:5:1212:15649:87294'
  console.assert bam.flag is 99
  console.assert bam.rname is "scaffold1"
  console.assert bam.pos is 1
  console.assert bam.mapq is 49
  console.assert bam.cigar.length is "XXSXXM".length * 1024
  console.assert bam.rnext is "scaffold1"
  console.assert bam.pnext is 474
  console.assert bam.tlen is 571
  console.assert bam.length is 100 * 1024
  console.assert bam.qual.length is 100 * 1024
  console.assert bam.multiple is true
  console.assert bam.allmatches is true
  console.assert bam.unmapped is false
  console.assert bam.next_unmapped is false
  console.assert bam.reversed is false
  console.assert bam.next_reversed is true
  console.assert bam.last is false
  console.assert bam.secondary is false
  console.assert bam.lowquality is false
  console.assert bam.unique is true
  console.assert bam.mismatch is 0
  console.assert bam.duplicate is false
  console.assert bam.supplementary is false
  console.assert bam.tagstr is "NM:i:0\tAS:i:87\tXS:i:80"

.on "end", ->
  console.log n
  console.assert n is 101
