BAMReader
==========
BAMReader is a reader to parse .bam files.
Uses samtools If exists, otherwise uses a native parser.

installation
----------------
```bash
$ npm install -g bamreader
```

usage
-------------
```js
var BAMReader = require("bamreader");
var reader = BAMReader.create("/path/to/bamfile.bam");
reader.on("bam", function(bam) {
  // bam: object. see "bam object" section
  console.log(bam.seq, bam.qual);
  console.log(bam.pair); // pair bam object of the bam. To do this, indexing is needed.
});

reader.on("end", function() {
  console.log("all bam alignments have been read.");
});
```

bam object
------------------
- pair          : information of the mate pair. see "bam.pair"
- discordant    : discordant or not. see "bam.discordant"
- qname         : name of the read
- flag          : bitwise flag
- rname         : reference sequence rname
- pos           : 1-based leftmost mapping position
- mapq          : mapping quality
- cigar         : CIGAR string
- rnext         : reference sequence name of the primary alignment of the next read
- next_rname    : the actual reference name of the next read.
- pnext         : position of the primary alignment of the next read
- tlen          : template length
- seq           : segment jsequence
- qual          : ASCII of base quality plus 33
- start         : 0-based leftmost mapping position
- flags         : information of flags
- tagstr        : optional fields as string
- tags          : optional fields as structured data described in "bamdata.tags" section
- multiple      : information from flag: template have multiple segment or not. 
- unmapped      : information from flag: unmapped or not. 
- allmatches    : information from flag: each segment properly aligned or not
- next_unmapped : information from flag: next segment is unmapped or not. 
- reversed      : information from flag: reversely mapped or not.
- next_reversed : information from flag: next segment is reversely mapped or not.
- first         : information from flag: first segment or not.
- last          : information from flag: last segment or not.
- secondary     : information from flag: secondary alignment or not.
- lowquality    : information from flag: not passing quality control or not.
- duplicate     : information from flag: PCR or optic duplicate or not.
- supplementary : information from flag: supplementary alignment or not.
- clipped       : clipped or not.
- soft_clipped  : soft-clipped or not.
- hard_clipped  : hard-clipped or not.
- left_break    : the left position of the breakpoint if clipped.
- right_break   : the right position of the breakpoint if clipped.
- match_len     : length of the matched portion of the read.
- fully_matched : if the whole base of the read is matched or not
- mismatch      : the number of mismatch bases
- unique        : if the read is mapped uniquely or not
- indel         : the mapping result contains insertion/deletion or not.
- different_ref : if both mapped and rname is not rnext or not
- same_strand   : if both mapped and mapped to the same strand or not
- has_n         : the read contains "N" or not
- sam           : SAM formatted string of the bam

bam.pair
-------------------
bamreader can fetch one mate pair of the read using index.

```bash
$ bamreader -c <bamfile>     # creates an index of the bamfile.
```
then you can access bam.pair

bam.pair is also a bam object.

bam.discordant
-------------------
(...)
 
bam.tags
------------------
- (key): tag name (two string character)
- (value): type: data type, value: data

if data type is "B", value[0] is subtype, value.slice(1) is the array of data.
