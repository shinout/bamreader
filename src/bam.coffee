SEQ_ARR = do ->
  ret = []
  arr = "=ACMGRSVTWYHKDBN".split("")
  for c1,i in arr
    for c2,j in arr
      ret.push if c2 is "=" then c1 else c1 + c2
  return ret

defineGetters = (obj, getters)->
  Object.defineProperty(obj, name, get: fn) for name, fn of getters

CIGAR = module.exports.CIGAR
class BAM
  @createFromSAM = (sam, @reader)->
    # mimics bam object
    d = sam.split("\t")

    # native values
    # "ref_id" and "nref_id" are used for some getter properties
    bam =
      qname : d[0]
      flag  : Number d[1]
      pos   : Number d[3]
      mapq  : Number d[4]
      pnext : Number d[7]
      tlen  : Number d[8]
      rname     : if d[2] is "*" then null else d[2]
      ref_id    : if d[2] is "*" then -1 else d[2]
      rnext     : if d[6] is "*" then null else d[6]
      nref_id   : if d[6] is "*" then -1 else d[2]
      cigar     : if d[5] is "*" then null else d[5]
      seq       : d[9]
      qual      : d[10]
      tagstr_   : d.slice(11).join("\t")
      sam       : sam

    bam.__proto__ = BAM.prototype
    return bam

  constructor: (buf, @reader)->
    @bytesize = buf.readInt32LE(0, true) + 4
    @ref_id   = buf.readInt32LE 4, true
    @pos      = buf.readInt32LE(8, true) + 1
    @mapq     = buf.readUInt8 13, true
    @flag     = buf.readUInt16LE 18, true
    @nref_id  = buf.readInt32LE 24, true
    @pnext    = buf.readInt32LE(28, true) + 1
    @tlen     = buf.readInt32LE 32, true
    #bin      = buf.readUInt16LE 14, true

    l_qname = buf.readUInt8 12, true
    @qname = buf.slice(36, 36 + l_qname - 1).toString("ascii")

    l_cigar = buf.readUInt16LE 16, true
    cursor = 36 + l_qname
    @cigarbytes = buf.slice(cursor, cursor + l_cigar * 4)
    @l_cigar = l_cigar

    l_seq = buf.readInt32LE 20, true
    cursor += l_cigar * 4
    b_seqlen = Math.floor((l_seq+1)/2)
    @seqbytes = buf.slice(cursor, cursor + b_seqlen)
    @length  = l_seq

    cursor += b_seqlen
 
    @qualbytes = buf.slice(cursor, cursor+l_seq)
    cursor += l_seq

    @tagbytes = buf.slice(cursor)

  defineGetters @::,
    ######################
    # FLAG PROPERTIES
    ######################
    multiple      : -> !!(@flag & (0x01))
    allmatches    : -> !!(@flag & (0x02))
    unmapped      : -> !!(@flag & (0x04))
    next_unmapped : -> !!(@flag & (0x08))
    reversed      : -> !!(@flag & (0x10))
    next_reversed : -> !!(@flag & (0x20))
    first         : -> !!(@flag & (0x40))
    last          : -> !!(@flag & (0x80))
    secondary     : -> !!(@flag & (0x100))
    lowquality    : -> !!(@flag & (0x200))
    duplicate     : -> !!(@flag & (0x400))
    supplementary : -> !!(@flag & (0x800))

    ######################
    # BASIC PROPERTIES
    ######################
    rname: -> if @ref_id  is -1 then null else @reader.refs[@ref_id].name
    rnext: -> if @nref_id is -1 then null else @reader.refs[@nref_id].name

    seq: ->
      return @seq_ if @seq_?
      seq = []
      for byte in @seqbytes
        seq.push SEQ_ARR[byte]
      @seq_ = seq.join("")

    qual: ->
      return @qual_ if @qual_?
      qual = []
      for byte in @qualbytes
        qual.push String.fromCharCode(byte + 33)
      @qual_ = qual.join("")

    ######################
    # CIGAR PROPERTIES
    ######################
    CIGAR: ->
      return @CIGAR_ if @CIGAR_?
      @CIGAR_ = new CIGAR(@cigarbytes, @l_cigar)
    cigar: -> @CIGAR.string
    clipped: -> @CIGAR.soft_clipped() or @CIGAR.hard_clipped()
    soft_clipped: -> @CIGAR.soft_clipped()
    hard_clipped: -> @CIGAR.hard_clipped()
    match_len: -> @CIGAR.len()
    left_break: -> @CIGAR.bpL(@pos)
    right_break: -> @CIGAR.bpR(@pos)
    indel: -> @CIGAR.indel()
    fully_matched: -> @CIGAR.fully_matched()

    ######################
    # SAM STRING
    ######################
    sam: ->
      return @sam_ if @sam_
      @sam_ = [
        @qname
        @flag
        if @ref_id  is -1 then "*" else @reader.refs[@ref_id].name
        @pos
        @mapq
        @cigar or "*"
        if @nref_id is -1 then "*" else if @ref_id is @nref_id then "=" else @reader.refs[@nref_id].name
        @pnext
        @tlen
        @seq
        @qual
        @tagstr
      ].join("\t")

    tagstr: ->
      return @tagstr_ if @tagstr_
      @tagstr_ = ([name, tag.type, if Array.isArray tag.value then tag.value.join(",") else tag.value].join(":") for name,tag of @tags).join("\t")

    ######################
    # OPTIONAL PROPERTIES
    ######################
    pair: (bamObj)->
      return null if not @reader or not @multiple
      bams = @reader.find(@qname)
      for bam in bams
        continue if @secondary or @supplementary or @flag is bam.flag
        if @next_unmapped
          if bam.unmapped
            bam.reader = @reader
            return bam
        else if @nref_id isnt -1 and @pnext is bam.pos and @nref_id is bam.ref_id
            bam.reader = @reader
            return bam
      return null

    different_ref: ->
      if @multiple and not @unmapped and not @next_unmapped then @ref_id isnt @nref_id else null

    same_strand: ->
      if @multiple then @reversed is @next_reversed else null

    has_n: ->
      return @has_n_ if @has_n_?
      @has_n_ = false
      for byte in @seqbytes
        return @has_n_ = true if byte >= 0xf0 or (byte & 0x0f) is 0x0f # 15 * 16 (upper) or 15 (lower)
      return @has_n_

    # only works if the mapper is BWA
    unique: ->
      not @unmapped and @mapq isnt 0

    # only works if the mapper is BWA
    mismatch: ->
      return null if not @tags.NM
      return @tags.NM.value

    discordant: ->
      return null if not @reader or @tlen is 0 or not @reader.tlen_mean or not @reader.tlen_sd
      m = @reader.tlen_mean
      sd = @reader.tlen_sd
      return @tlen < m - 2*sd or m + 2*sd < @tlen

    mean_qual: ->
      total = 0
      total += byte for byte in @qualbytes
      return Math.floor(total / @qualbytes.length)

    ######################
    # TAG PROPERTY
    ######################
    tags: ->
      if @tagstr_
        for tag in @tagstr_.split("\t")
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
        return @tags_  if @tags_?

      tags = {}
      cursor = 0
      buflen = @tagbytes.length
      buf = @tagbytes
      loop
        break if cursor-4 >= buflen
        tag = buf.slice(cursor, cursor+2).toString("ascii")
        cursor+=2
        valtype = String.fromCharCode buf[cursor]
        cursor++

        switch valtype
          when "A"
            tags[tag] = type: valtype, value: String.fromCharCode buf[cursor]
            cursor++
          when "c"
            tags[tag] = type: "i", value: buf.readInt8 cursor
            cursor++
          when "C"
            tags[tag] = type: "i", value: buf.readUInt8 cursor
            cursor++
          when "s"
            tags[tag] = type: "i", value: buf.readInt16LE cursor
            cursor+=2
          when "S"
            tags[tag] = type: "i", value: buf.readUInt16LE cursor
            cursor+=2
          when "i"
            tags[tag] = type: "i", value: buf.readInt32LE cursor
            cursor+=4
          when "I"
            tags[tag] = type: "i", value: buf.readUInt32LE cursor
            cursor+=4
          when "f"
            tags[tag] = type: valtype, value: buf.readFloatLE cursor
            cursor+=4
          when "B"
            subtype = String.fromCharCode buf[cursor]
            cursor++
            arrayLen = buf.readInt32LE cursor
            cursor+=4
            switch subtype
              when "c"
                tags[tag] = type: valtype, value: (buf.readInt8 cursor+i, true for i in [0...arrayLen])
                cursor+=arrayLen
              when "C"
                tags[tag] = type: valtype, value: (buf.readUInt8 cursor+i, true for i in [0...arrayLen])
                cursor+=arrayLen
              when "s"
                tags[tag] = type: valtype, value: (buf.readInt16LE cursor+i*2, true for i in [0...arrayLen])
                cursor+=arrayLen*2
              when "S"
                tags[tag] = type: valtype, value: (buf.readUInt16LE cursor+i*2, true for i in [0...arrayLen])
                cursor+=arrayLen*2
              when "i"
                tags[tag] = type: valtype, value: (buf.readInt32LE cursor+i*4, true for i in [0...arrayLen])
                cursor+=arrayLen*4
              when "I"
                tags[tag] = type: valtype, value: (buf.readUInt32LE cursor+i*4, true for i in [0...arrayLen])
                cursor+=arrayLen*4
              when "f"
                tags[tag] = type: valtype, value: (buf.readFloatLE cursor+i*4, true for i in [0...arrayLen])
                cursor+=arrayLen*4
            value.unshift subtype

          when "Z"
            zLen = 0
            zLen++ while buf[cursor+zLen] isnt 0x00
            tags[tag] = type: valtype, value: buf.slice(cursor, cursor+zLen).toString("ascii")
            cursor+=zLen+1
          when "H"
            hLen = 0
            hLen++ while buf[cursor+hLen] isnt 0x00
            tags[tag] = type: valtype, value: buf.slice(cursor, cursor+hLen).toString("hex")
            cursor+=hLen+1
      @tags_ = tags

module.exports.BAM = BAM
