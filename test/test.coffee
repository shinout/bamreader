sam = require("fs").readFileSync(__dirname + "/d.sam", "utf-8").split("\n")[0]
BAMReader = require(__dirname + "/../lib/bamreader.js")
Bam = BAMReader.Bam
bam = Bam.createFromSAM(sam)

console.assert bam.qname is 'HS2000-903_160:5:2103:9285:63958'
console.assert bam.flag is 163
console.assert bam.rname is "chr1"
console.assert bam.pos is 12219
console.assert bam.mapq is 0
console.assert bam.cigar is "100M"
console.assert bam.rnext is "="
console.assert bam.pnext is 12431
console.assert bam.tlen is 312
console.assert bam.seq is 'CCTAGGCCAGTAAGTAGTGCTTGTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTACCTGCCGTCTGCTGCCATC'
console.assert bam.qual is '@CCFFFFFHHHHHJEHIGHHIIJHJJJJFJJIJJJJHJJJJEHGIIIJHJGGIJJIJJIHHFDEF>ACCBBDCDBBABB@???@CBDBB@ADDCDDBC@C',
console.assert bam.multiple is true
console.assert bam.allmatches is true
console.assert bam.unmapped is false
console.assert bam.next_unmapped is false
console.assert bam.reversed is false
console.assert bam.next_reversed is true
console.assert bam.last is true
console.assert bam.secondary is false
console.assert bam.lowquality is false
console.assert bam.unique is false
console.assert bam.mismatch is 0
console.assert bam.duplicate is false
console.assert bam.supplementary is false
console.assert bam.tagstr is bam.tagstr_
console.assert bam.tagstr is 'NM:i:0\tAS:i:100\tXS:i:100\tRG:Z:NA12878',
console.assert bam.sam is bam.sam_
console.assert bam.sam is 'HS2000-903_160:5:2103:9285:63958\t163\tchr1\t12219\t0\t100M\t=\t12431\t312\tCCTAGGCCAGTAAGTAGTGCTTGTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCGCTCCAGCAGCTGGACCCCTACCTGCCGTCTGCTGCCATC\t@CCFFFFFHHHHHJEHIGHHIIJHJJJJFJJIJJJJHJJJJEHGIIIJHJGGIJJIJJIHHFDEF>ACCBBDCDBBABB@???@CBDBB@ADDCDDBC@C\tNM:i:0\tAS:i:100\tXS:i:100\tRG:Z:NA12878'

flags = [
  "multiple"
  "allmatches"
  "unmapped"
  "next_unmapped"
  "reversed"
  "next_reversed"
  "first"
  "last"
  "secondary"
  "lowquality"
  "duplicate"
  "supplementary"
]
console.assert bam.flags[flag] is bam[flag] for flag in flags
console.assert bam.start is bam.pos - 1
console.assert bam.length is bam.seq.length
console.assert bam.clipped is false
console.assert bam.soft_clipped is false
console.assert bam.hard_clipped is false
#console.log bam.match_len
console.assert bam.match_len is 100
console.assert bam.left_break is null
console.assert bam.right_break is null
console.assert bam.indel is false
console.assert bam.fully_matched is true
console.assert bam.different_reference is false
console.assert bam.same_strand is false
console.assert bam.has_n is false
#console.log bam.discordant
console.assert bam.discordant is null


n = 0
m = 0
reader = BAMReader.create(__dirname + "/d.bam").on "bam", (bam, dOffset, iOffset)->
  console.log dOffset, iOffset
  if bam.pair
    m++
  n++
.on "end", ->
  console.assert n is m

n2 = 0
m2 = 0
reader = BAMReader.create(__dirname + "/d.bam", native: true).on "bam", (bam)->
  #console.log bam
  if bam.pair
    m2++
  n2++
.on "end", ->
  console.log n2,m2
  console.assert n2 is m2
