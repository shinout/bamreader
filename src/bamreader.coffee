fs = require("fs")
inflateBGZF = require("bgzf").inflate
isValidBGZF = require("bgzf").hasValidHeader
INFBUF_CACHE_SIZE = 1024 * 1024 * 400

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


class BAMReader
  constructor: (bamfile, o = {})->
    @bamfile = require("path").resolve(bamfile)
    @cache_size = o.cache_size or INFBUF_CACHE_SIZE
    @infbufs = new Fifo(@cache_size)
    @fd = fs.openSync(@bamfile, "r")

    return if o.from_obj
    @size = fs.statSync(@bamfile).size
    _readHeader_result = @_readHeader() # @header, @refs, @header_offset is set
    throw "couldn't read header" if null is _readHeader_result

  @create: (bamfile)->
    return new BAMReader(bamfile)

  #####################################
  # restore from object(hash)
  #####################################
  @createFromObject: (obj)->
    reader = new BAMReader(obj.bamfile, from_obj: true, cache_size: obj.cache_size)
    for k in ["size", "header", "header_offset", "refs"]
      reader[k] = obj[k]
    return reader
      

  #####################################
  # gets a restorable object(hash)
  #####################################
  toObject: ->
    ret = {}
    ret[k] = @[k] for k in ["size", "header", "header_offset", "refs", "bamfile", "cache_size"]
    return ret

  #####################################
  # creates obj reading bams in order
  #####################################
  createIterator: (o={})->
    o = on_bam: o if typeof o is "function"
    return module.exports.BAMIterator.create(@, o)
    
  #####################################
  # shortcut for iterator
  #####################################
  on : (name, fn)->
    if name is "bam"
      return @createIterator(fn)

   #####################################
  # reads an alignment with the offsets
  #####################################
  read: (i_offset, d_offset)->
    buf = @infbufs.get d_offset
    if buf
      buf = buf.slice(i_offset)
      if buf.length >= 4
        bytesize = buf.readInt32LE(0, true) + 4
        # FIXME: longer bam data cannot be restored
        if buf.length >= bytesize
          return new module.exports.BAM(buf.slice(0, bytesize))
    pitch = 512
    loop
      pitch += pitch
      read_size = Math.min @size - d_offset, pitch
      chunk = new Buffer(read_size)
      fs.readSync @fd, chunk, 0, read_size, d_offset
      [infbuf, i_offsets, d_offsets] = inflateBGZF chunk
      infbuf = infbuf.slice(i_offset)
      if infbuf.length < 4
        throw "couldn't fetch bam" if read_size is @size - d_offset
        continue
      bytesize = infbuf.readInt32LE(0, true) + 4
      if infbuf.length < bytesize
        throw "couldn't fetch bam" if read_size is @size - d_offset
        continue
      bambuf = infbuf.slice(0, bytesize)
      bam = new module.exports.BAM(bambuf, @)
      bam.i_offset = i_offset
      bam.d_offset = d_offset

      # cache
      d_offsets.pop()
      for offset,i in d_offsets
        @infbufs.set d_offset + offset, infbuf.slice(i_offsets[i], i_offsets[i+1])
      return bam

  #####################################
  # splits body section into num parts
  #####################################
  split: (num)->
    pitch = 65535
    num = if typeof num is "number" and num >= 1 then parseInt(num) else 2
    interval = Math.floor((@size - @header_offset)/num)
    positions = []

    k = 0
    # finding accurate position of BGZF
    while k < num
      #start = interval * k + @header_offset - 1
      start = interval * k + @header_offset
      buflen = Math.min(pitch, interval)
      buf = new Buffer(buflen)
      fs.readSync @fd, buf, 0, buflen, start
      cursor = -1
      match = false
      until match or cursor + 16 > buf.length
        cursor++
        continue unless isValidBGZF buf.slice(cursor, cursor+16)
        d_offset = start + cursor
        # checks if the BGZF block contains the start of alignment buffer
        try
          [infbuf] = inflateBGZF buf.slice(cursor)
        catch e
          # invalid format: retry
          infbuf = new Buffer(0)
        if infbuf.length < 24
          if buflen isnt interval
            # retry with much buffer
            k--
            pitch += pitch
          break
        ref_id  = infbuf.readInt32LE 4,  true
        nref_id = infbuf.readInt32LE 24, true
        # if valid inf position
        match = true if (ref_id is -1 or @refs[ref_id]?) and (nref_id is -1 or @refs[nref_id]?)

      positions.push(d_offset) if match
      k++
    if positions.length is 0 and num > 1
      return @split num - 1
    return positions

  #####################################
  # creates child iterators
  #####################################
  fork: (o={})->
    o = on_bam: o if typeof o is "function"
    o.spawn = true
    num = if typeof o.num is "number" and o.num >= 1 then parseInt(o.num) else 2
    positions = @split num
    num = positions.length
    # attach "on_"
    for suffix in ["bam", "end", "finish", "message"]
      o["on_#{suffix}"] = o[suffix] if typeof o[suffix] is "function"
      delete o[suffix]

    positions.push @size
    childs = []
    on_finish  = o.on_finish
    on_message = o.on_message
    delete o.on_finish
    delete o.on_message

    ended_childs = 0
    envs = []
    for n in [0...num]
      _o =
        start : positions[n]
        end   : positions[n+1]
        n     : n
      _o[k] = v for k, v of o
      child = module.exports.BAMIterator.create(@, _o)
      child.on "message", (env)->
        if env.ended
          envs[env.n] = env
        else
          on_message env if typeof on_message is "function"
      child.on "exit", ->
        ended_childs++
        return if ended_childs < num
        on_finish.call(@, envs) if typeof on_finish is "function"
      childs.push child
    return childs

  @fork = (args...)->
    if typeof args[0] is "string"
      throw new Error("argument 2 should be an object.") unless args[1]
      file = args[0]
      o = args[1]
    else
      o = args[0]
      file = o.file
      delete o.file
    return BAMReader.create(file).fork(o)


  # (private) reads header
  _readHeader: ->
    read_size = 32768
    loop
      read_size = Math.min @size, read_size
      buf = new Buffer read_size
      fs.readSync @fd, buf, 0, read_size, 0
      [infbuf, i_offsets, d_offsets] = inflateBGZF buf
      try
        cursor = 0
        refs = {}
        headerLen = infbuf.readInt32LE(4)
        throw new Error("header len") if infbuf.length < headerLen + 16
        headerStr = infbuf.slice(8,headerLen+8).toString("ascii")
        cursor = headerLen + 8
        nRef = infbuf.readInt32LE cursor
        cursor+=4

        blen = infbuf.length

        for i in [0...nRef]
          nameLen = infbuf.readInt32LE cursor
          cursor+=4
          name = infbuf.slice(cursor, cursor+nameLen-1).toString("ascii")
          cursor+=nameLen
          refLen = infbuf.readInt32LE cursor
          cursor+=4
          refs[i] = name: name, len: refLen

        @refs = refs
        @header = headerStr
        loop
          current_i_offset = i_offsets.shift()
          current_d_offset = d_offsets.shift()
          break if cursor <= current_i_offset
        @header_offset = current_d_offset
        break
      catch e
        return null if read_size is @size
        read_size += read_size
    return true

module.exports = BAMReader
