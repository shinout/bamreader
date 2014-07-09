CIGAR_ARR = "MIDNSHP=X".split("")
nullobj =
  string: null
  soft_clipped: -> null
  hard_clipped: -> null
  bpL: -> null
  bpR: -> null
  len: -> null
  indel: -> null
  fully_matched: -> null


class CIGAR
  constructor: (buf, l_cigar)->
    return nullobj if buf and buf.length is 0
    @_indel = false
    return if buf is null
    @arr = new Array(l_cigar)
    str = ""
    i = 0
    while i < l_cigar
      num  = buf.readUInt32LE(i * 4, true)
      type = CIGAR_ARR[num & 0x0f]
      @_indel = true if type.match(/[ID]/)
      num  = num>>4
      @arr[i] = num: num, type: type
      str += num + type
      i++
    @string = str

  @createFromString = (str)->
    return nullobj if not str or str is "*"
    cigar = new CIGAR()
    arr = []
    cigarr = str.split(/([A-Z=])/)[...-1]
    i = 0
    l_cigar = cigarr.length/2
    while i < l_cigar
      i2 = i*2
      type = cigarr[i2+1]
      cigar._indel = true if type.match(/[ID]/)
      arr.push
        num  : Number cigarr[i2]
        type : type
    cigar.arr = arr
    cigar.string = str

  soft_clipped: -> @arr[0].type is "S" or @arr[@arr.length-1].type is "S"
  hard_clipped: -> @arr[0].type is "H" or @arr[@arr.length-1].type is "H"
  indel:  -> @_indel
  fully_matched: ->
    @arr.length is 1 and @arr.type is "M"

  # leftside breakpoint
  bpL: (pos)-> if @arr[0].type.match(/[SH]/) then pos else null

  # rightside breakpoint
  bpR: (pos)->
    matched = false
    ret = pos
    for info in @arr
      if not matched and info.type is "M"
        matched = true
      if matched
        return ret if info.type.match(/[SH]/)
        ret += info.num
    return null

  len: ->
    ret = 0
    ret += info.num  for info in @arr when info.type.match(/[MS=XI]/)
    return ret

module.exports.CIGAR = CIGAR
