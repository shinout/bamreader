class CIGAR
  constructor: (@str)->
    return null if not str or str is "*"
    @arr = []

    cigarr = @str.split(/([A-Z=])/)[...-1]
    for i in [0...cigarr.length/2]
      i2 = i*2
      @arr.push
        num  : Number cigarr[i2]
        type : cigarr[i2+1]

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
