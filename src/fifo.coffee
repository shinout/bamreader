class Fifo
  constructor: (@size_limit = INFBUF_CACHE_SIZE)->
    @hash = {}
    @keys = []
    @total_size = 0

  get: (k) -> @hash[k]

  set: (k, v)->
    return if @hash[k]?
    while @total_size > @size_limit
      key_to_del = @keys.shift()
      @total_size -= @hash[key_to_del].length
      delete @hash[key_to_del]

    @keys.push k
    @hash[k] = v
    @total_size += v.length
    return

  clear: ->
    delete @hash[@keys.shift()] while @keys.length
    @total_size = 0
    return

module.exports.Fifo = Fifo
