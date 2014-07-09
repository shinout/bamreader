defineGetters = (obj, getters)->
  Object.defineProperty(obj, name, get: fn) for name, fn of getters

SEQ_ARR = "=ACMGRSVTWYHKDBN".split("")

CIGAR = module.exports.CIGAR

ARGNAMES = [
  "reader"
  "qname"
  "flag"
  "rname"
  "pos"
  "mapq"
  "cigar"
  "rnext"
  "pnext"
  "tlen"
  "seqbits"
  "length"
  "qualbytes"
  "bytesize"
]


FLAGS = [
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

class Bam
  constructor: (args...)->
    @[v] = args[k] for v,k in ARGNAMES
    @_parseFlag()

  _parseFlag: ->
    @[flagname] = !!(@flag & (0x01 << i)) for flagname,i in FLAGS

  defineGetters @::,
    flags: -> @

    seq: ->
      return @seq_ if @seq_
      seq = []
      for byte in @seqbits
        seq.push SEQ_ARR[byte >>4]
        second = SEQ_ARR[byte & 0x0F]
        seq.push second if second isnt "="
      @seq_ = seq.join("")

    qual: ->
      return @qual_ if @qual_
      @qual_ = (String.fromCharCode(@qualbytes[i]+33) for i in [0...@length]).join("")

    tagstr: ->
      return @tagstr_ if @tagstr_
      @tagstr_ = ([name, tag.type, if Array.isArray tag.value then tag.value.join(",") else tag.value].join(":") for name,tag of @tags).join("\t")

    tags : ->
      return @tags_ if @tags_
      @tags_ = {}
      for tag in @tagstr.split("\t")
        val = tag.split(":")
        tag = val[0]
        type = val[1]
        switch type
          when "i","f" then value = Number val[2]
          when "B"
            value = val[2].split(",")
            subtype = value[0]
            if subtype in ["c","C","s","S","i","I","f"]
              value = (Number v for v in value)
              value[0] = subtype
          else
            value = val[2]
        @tags_[tag] = type: type, value: value
      return @tags_

    next_rname: ->
      if @rnext is "=" then @rname else @rnext
    pair: (bamObj)->
      return null if not @reader or not @reader.dic or not @multiple
      bams = @reader.dic.fetch(@qname)
      for bam in bams
        continue if @secondary or @supplementary or @flag is bam.flag
        if @next_unmapped
          if bam.flags.unmapped
            bam.reader = @reader
            return bam
        else
          next_rname = if @rnext is "=" then @rname else @rnext
          if @pnext is bam.pos and next_rname is bam.rname
            bam.reader = @reader
            return bam
      return null

    start: ->
      if @pos then @pos - 1 else null
    clipped: ->
      if @cigar is "*" then null else !!@cigar.match(/[HS]/)

    soft_clipped: ->
      if @cigar is "*" then null else !!@cigar.match("S")

    hard_clipped: ->
      if @cigar is "*" then null else !!@cigar.match("H")

    match_len: ->
      @CIGAR = new CIGAR(@cigar) unless @CIGAR
      @CIGAR.len()

    left_break: ->
      @CIGAR = new CIGAR(@cigar) unless @CIGAR
      @CIGAR.bpL(@pos)

    right_break: ->
      @CIGAR = new CIGAR(@cigar) unless @CIGAR
      @CIGAR.bpR(@pos)

    indel: ->
      if @cigar is "*" then null else !!@cigar.match(/[ID]/)

    fully_matched: ->
      if @cigar is "*" then null else @cigar is @length+"M"

    different_reference: ->
      if @multiple and not @unmapped and not @next_unmapped then @rnext isnt "=" else null

    same_strand: ->
      if @multiple then @reversed is @next_reversed else null

    has_n: ->
      !!@seq.match("N")

    # only works if the mapper is BWA
    unique: ->
      not @unmapped and @mapq isnt 0

    # only works if the mapper is BWA
    mismatch: ->
      return null if not @tags.NM
      return @tags.NM.value

    discordant: ->
      return null if not @reader or@tlen is 0
      m = @reader.tlen_mean
      sd = @reader.tlen_sd
      return @tlen < m - 2*sd or m + 2*sd < @tlen

    sam: ->
      return @sam_ if @sam_
      @sam_ = [
        @qname
        @flag
        @rname
        @pos
        @mapq
        @cigar or "*"
        if @rnext is @rname and @rname isnt "*" then "=" else @rnext
        @pnext
        @tlen
        @seq
        @qual
        @tagstr
      ].join("\t")

  @createFromSAM = (samline, bamreader)->
    sam = samline.split("\t")
    bam = new Bam(
      bamreader,
      sam[0],              # qname
      Number(sam[1]),      # flag
      sam[2],              # rname
      Number(sam[3]),      # pos
      Number(sam[4]),      # mapq
      sam[5],              # cigar
      sam[6],              # rnext
      Number(sam[7]),      # pnext
      Number(sam[8]),      # tlen
      null,                # seqbits
      sam[9].length,       # length
      null                 # qualbytes
      null                 # bytesize
    )
    bam.seq_  = sam[9]     # seq
    bam.qual_ = sam[10]    # qual
    bam.tagstr_ = sam.slice(11).join("\t")
    bam.sam_ = samline
    return bam

module.exports.Bam = Bam
