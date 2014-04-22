BAMReader
==========
BAMReader is a reader to parse .bam files.

installation
----------------
```bash
$ npm instal bamreader
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
  console.log(samline.split("\t");
});

reader.on("end", function() {
  console.log("all bam alignments have read.");
});
```js

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
- refid   : id of the mapped reference sequence
- nrefid  : id of the reference sequence to which the next read mapped
- start   : 0-based leftmost mapping position
- flags   : information of flags
- tagstr  : optional fields as string
- tags    : optional fields as structured data (for detail, please dump it...)
