BAMReader
==========
BAMReader is a reader to parse .bam files.
Uses samtools If exists, otherwise uses a native parser.

installation
----------------
```bash
$ npm install bamreader
```

usage
-------------
```js
var BAMReader = require("bamreader");
var reader = BAMReader.create("/path/to/bamfile.bam");
reader.on("bam", function(bamdata) {
  // bamdata: object. see "bamdata" section
  console.log(bamdata.seq, bamdata.qual);
});

reader.on("sam", function(samline) {
  // samline: sam string
  console.log(samline.split("\t"));
});

reader.on("end", function() {
  console.log("all bam alignments have read.");
});
```
options for BAMReader.create
----------------------------
- native : read bam without samtools (slower than using samtools)

bamdata
------------------
- qname   : name of the read
- flag    : bitwise flag
- rname   : reference sequence rname
- pos     : 1-based leftmost mapping position
- mapq    : mapping quality
- cigar   : CIGAR string
- rnext   : reference sequence name of the primary alignment of the next read
- pnext   : position of the primary alignment of the next read
- tlen    : template length
- seq     : segment jsequence
- qual    : ASCII of base quality plus 33
- start   : 0-based leftmost mapping position
- flags   : information of flags
- tagstr  : optional fields as string
- tags    : optional fields as structured data described in "bamdata.tags" section

bamdata.tags
------------------
- (key): tag name (two string character)
- (value): type: data type, value: data

if data type is "B", value[0] is subtype, value.slice(1) is the array of data.
