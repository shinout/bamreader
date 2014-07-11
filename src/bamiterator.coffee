fs = require "fs"
inflateBGZF = require("bgzf").inflate
DEFAULT_PITCH = 16384000
BAM = module.exports.BAM
noop = ->

class BAMIterator
  @create = (reader, options = {})->
    return new BAMIterator(reader, options)

  constructor: (@reader, o={})->
    if typeof o.end is "function"
      o.on_end = o.end
      delete o.end
    if typeof o.finish is "function"
      o.on_finish = o.finish
    if typeof o.on_finish is "function"
      o.on_end = o.on_finish if not o.on_end
    if typeof o.bam is "function"
      o.on_bam = o.bam if not o.on_bam

    @offset = if typeof o.start is "number" then o.start else @reader.header_offset
    @end = if typeof o.end is "number" then o.end else @reader.size
    @pitch = if typeof o.pitch is "number" then o.pitch else DEFAULT_PITCH
    @on_bam = if typeof o.on_bam is "function" then o.on_bam else noop
    @on_end = if typeof o.on_end is "function" then o.on_end else noop
    @on_start  = if typeof o.on_start is "function" then o.on_start else noop
    @pause  = if typeof o.pause  is "function" then o.pause  else null
    @env = o.env or o.$ or {} # environment to register variables, especially for child processes
    @paused = false
    @ended  = false
    process.nextTick =>
      @on_start(@env)
      process.nextTick @_read.bind @

  on: (name, fn)->
    switch name
      when "end"
        @on_end = fn
      when "bam"
        @on_bam = fn
    @

  resume: ->
    if @paused
      @paused = false
      if @ended
        process.nextTick @on_end.bind @, @env
      else
        process.nextTick @_read.bind @

  send: (msg)->
    process.send msg if typeof process.send is "function"

  _read: ->
    read_size = Math.min @end - @offset, @pitch
    @ended = @ended or read_size <= 0
    if typeof @pause is "function"
      @paused = @pause(@ended, @env)
    return if @paused
    return @on_end(@env) if @ended
    chunk = new Buffer(read_size)
    fs.readSync @reader.fd, chunk, 0, read_size, @offset
    [infbuf, i_offsets, d_offsets] = inflateBGZF chunk
    infbuf_len = infbuf.length
    if infbuf_len is 0
      @ended = read_size is @end - @offset
      @pitch += @pitch
      return @_read()

    for offset,i in d_offsets
      if i_offsets[i+1]
        @reader.infbufs.set @offset + offset, infbuf.slice(i_offsets[i], i_offsets[i+1])

    i_offset = 0
    current_i_offset = i_offsets.shift()
    current_d_offset = @offset + d_offsets.shift()
    loop
      break if i_offset + 4 > infbuf_len
      bytesize = infbuf.readInt32LE(i_offset, true) + 4
      break if i_offset + bytesize > infbuf_len
      bambuf = infbuf.slice(i_offset, i_offset + bytesize)

      bam = new BAM(bambuf, @reader)
      bam.i_offset = i_offset - current_i_offset
      bam.d_offset = current_d_offset
      @on_bam bam, @env

      i_offset += bytesize

      # updating i_offset, d_offset
      loop
        break if i_offsets[0] is undefined or i_offset < i_offsets[0]
        next_i_offset = i_offsets.shift()
        current_i_offset = next_i_offset
        current_d_offset = @offset + d_offsets.shift()

    @offset = current_d_offset
    setImmediate @_read.bind @

module.exports.BAMIterator = BAMIterator
