SEQ_ARR = do ->
  ret = []
  arr = "=ACMGRSVTWYHKDBN".split("")
  for c1,i in arr
    for c2,j in arr
      ret.push if c2 is "=" then c1 else c1 + c2
  return ret

QUAL_ARR = do->
  ret = {}
  for i in [0..94]
    ret[i] = String.fromCharCode(i + 33)
  return ret


KNOWN_TAGS = do->
  ret = {}
  for tag in [
    "AM", "AS"
    "BC", "BQ"
    "CC", "CM", "CO", "CP", "CQ", "CS", "CT"
    "E2",
    "FI", "FS", "FZ"
    "H0", "H1", "H2", "HI"
    "IH"
    "LB"
    "MC", "MD", "MQ"
    "NH", "NM"
    "OC", "OP", "OQ"
    "PG", "PQ", "PT", "PU"
    "QT", "Q2"
    "R2", "RG", "RT"
    "SA", "SM"
    "TC"
    "U2", "UQ"
    "XS"
  ]
    k = tag.charCodeAt(0) * 256 + tag.charCodeAt(1)
    ret[k] = tag
  ret

defineGetters = (obj, getters)-> Object.defineProperty(obj, name, get: fn) for name, fn of getters
CIGAR = module.exports.CIGAR
class BAM
  @createFromSAM = (sam, reader)->
    # mimics bam object
    d = sam.split("\t")

    # native values
    # "ref_id" and "nref_id" are used for some getter properties
    bam =
      reader: reader
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
      tagstr    : d.slice(11).join("\t")
      length    : d[9].length
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
      len = @seqbytes.length
      seq = ""
      i = 0
      seq += SEQ_ARR[@seqbytes[i++]] while i < len
      @seq_ = seq

    qual: ->
      return @qual_ if @qual_?
      len = @length
      qual = ""
      i = 0
      qual += QUAL_ARR[@qualbytes[i++]] while i < len
      @qual_ = qual

    ######################
    # CIGAR PROPERTIES
    ######################
    CIGAR: ->
      return @CIGAR_ if @CIGAR_?
      if @cigarbytes
        @CIGAR_ = new CIGAR(@cigarbytes, @l_cigar)
      else
        @CIGAR_ = CIGAR.createFromString(@cigar)
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
      @sam_ =
        @qname + "\t" +
        @flag + "\t" +
        (if @ref_id  is -1 then "*" else @reader.refs[@ref_id].name) + "\t" +
        @pos + "\t" +
        @mapq + "\t" +
        (@cigar or "*") + "\t" +
        (if @nref_id is -1 then "*" else if @ref_id is @nref_id then "=" else @reader.refs[@nref_id].name) + "\t" +
        @pnext + "\t" +
        @tlen + "\t" +
        @seq  + "\t" +
        @qual + "\t" +
        @tagstr

    ######################
    # OPTIONAL PROPERTIES
    ######################
    pair: ->
      return null if not @reader or not @multiple
      return @pair_ if @pair_
      bams = @reader.find(@qname, @d_offset)
      for bam in bams
        continue if @secondary or @supplementary or @flag is bam.flag
        if @next_unmapped
          if bam.unmapped
            bam.reader = @reader
            @pair_ = bam
            return bam
        else if @nref_id isnt -1 and @pnext is bam.pos and @nref_id is bam.ref_id
            bam.reader = @reader
            @pair_ = bam
            return bam
      @pair_ = null
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
      return null if not @reader or @tlen is 0 or not @reader.dic
      m = @reader.tlen_mean
      sd = @reader.tlen_sd
      tlen = Math.abs @tlen
      return tlen < m - 2*sd or m + 2*sd < tlen

    mean_qual: ->
      total = 0
      len = @length
      if @qualbytes
        total += byte for byte in @qualbytes
      else
        i = 0
        while i < len
          total += @qual.charCodeAt(i) - 33
      return Math.floor(total / len)
          

    ######################
    # TAG PROPERTY
    ######################
    tagstr: ->
      return @tagstr_ if @tagstr_
      tagstr = ""
      cursor = 0
      buflen = @tagbytes.length
      buf = @tagbytes
      loop
        break if cursor >= buflen
        tagbyte = buf.readUInt16BE(cursor, true)
        tagname = KNOWN_TAGS[tagbyte] or buf.slice(cursor, cursor+2).toString("ascii")
        switch tagname
          when "NM", "AS", "XS"
            tagstr +=  tagname + ":i:" + buf.readUInt8(cursor+3, true) + "\t"
            cursor += 4
            continue
        cursor+=2
        valtype = String.fromCharCode buf[cursor]
        cursor++
        type = null

        switch valtype
          when "A"
            value = String.fromCharCode buf[cursor]
            cursor++
          when "c"
            value = buf.readInt8 cursor, true
            type = "i"
            cursor++
          when "C"
            value = buf.readUInt8 cursor, true
            type = "i"
            cursor++
          when "s"
            value = buf.readInt16LE cursor,true
            type = "i"
            cursor+=2
          when "S"
            value = buf.readUInt16LE cursor, true
            type = "i"
            cursor+=2
          when "i"
            value = buf.readInt32LE cursor, true
            cursor+=4
          when "I"
            value = buf.readUInt32LE cursor, true
            type = "i"
            cursor+=4
          when "f"
            value = buf.readFloatLE cursor, true
            cursor+=4
          when "B"
            subtype = String.fromCharCode buf[cursor]
            cursor++
            arrayLen = buf.readInt32LE cursor, true
            cursor+=4
            switch subtype
              when "c"
                value = (buf.readInt8 cursor+i, true for i in [0...arrayLen])
                cursor+=arrayLen
              when "C"
                value = (buf.readUInt8 cursor+i, true for i in [0...arrayLen])
                cursor+=arrayLen
              when "s"
                value = (buf.readInt16LE cursor+i*2, true for i in [0...arrayLen])
                cursor+=arrayLen*2
              when "S"
                value = (buf.readUInt16LE cursor+i*2, true for i in [0...arrayLen])
                cursor+=arrayLen*2
              when "i"
                value = (buf.readInt32LE cursor+i*4, true for i in [0...arrayLen])
                cursor+=arrayLen*4
              when "I"
                value = (buf.readUInt32LE cursor+i*4, true for i in [0...arrayLen])
                cursor+=arrayLen*4
              when "f"
                value = (buf.readFloatLE cursor+i*4, true for i in [0...arrayLen])
                cursor+=arrayLen*4
            value.unshift subtype
            value = value.join(",")
          when "Z"
            zLen = 0
            zLen++ while buf[cursor+zLen] isnt 0x00
            value = buf.slice(cursor, cursor+zLen).toString("ascii")
            cursor+=zLen+1
          when "H"
            hLen = 0
            hLen++ while buf[cursor+hLen] isnt 0x00
            value = buf.slice(cursor, cursor+hLen).toString("hex")
            cursor+=hLen+1
          # end of switch
        tagstr += tagname + ":" + (type or valtype) + ":" + value + "\t"
        # end of loop
      @tagstr_ = tagstr.slice(0, -1)

    tags: ->
      return @tags_ if @tags_
      for tag in @tagstr.split("\t")
        val = tag.split(":")
        tag = val[0]
        type = val[1]
        switch type
          when "i","f" then value = Number val[2]
          # when "B"
          #   value = val[2].split(",")
          #   subtype = value[0]
          #   if subtype in ["c","C","s","S","i","I","f"]
          #     value = (Number v for v in value)
          #     value[0] = subtype
          else
            value = val[2]
        @tags_[tag] = type: type, value: value
      return @tags_

module.exports.BAM = BAM
