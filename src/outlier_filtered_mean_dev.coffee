class OutlierFilteredMeanDev
  constructor: (rate)->
    @values = {}
    @rate = if typeof rate is "number" and rate >= 0 and rate < 1 then rate else 0.05
    @n = 0

  add: (v)->
    @values[v] = 0 unless @values[v]
    @values[v]++
    @n++
    
  precalc: ->
    sum = 0
    rate = @rate
    lower_limit = (@n * rate)
    upper_limit = (@n - @n * rate)
    total = 0
    valids = 0
    squared = 0
    for v, n of @values
      sum += n
      continue if sum < lower_limit
      break if sum > upper_limit
      val = v * n
      total += val
      valids += n
      squared += val * v

    return {
      sum    : total
      squared: squared
      n      : valids
    }

    calc: ->
      {total, squared, n } = @precalc()
      mean = sum / n
      return {
        mean : mean
        dev  : squared / n - mean * mean
        n    : n
      }

# ofmd  = new OutlierFilteredMeanDev(0.05)
# rt = require("random-tools")
# i = 0
# 
# while i < 1000000
#   v = Math.round(rt.normalRandom(500, 50))
#   ofmd.add(v)
#   i++
# 
# console.log ofmd.calc()

module.exports.OutlierFilteredMeanDev = OutlierFilteredMeanDev
