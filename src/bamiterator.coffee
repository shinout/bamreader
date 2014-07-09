fs = require "fs"
inflateBGZF = require("bgzf").inflate
DEFAULT_PITCH = 16384
BAM = module.exports.BAM
noop = ->

class BAMIterator
  @create = (reader, options = {})->
    return new BAMIterator(reader, options)

  constructor: (@reader, o={})->
    @offset = if typeof o.start is "number" then o.start else @reader.header_offset
    @end = if typeof o.end is "number" then o.end else @reader.size
    @pitch = if typeof o.pitch is "number" then o.pitch else DEFAULT_PITCH
    @on_bam = if typeof o.on_bam is "function" then o.on_bam else noop
    @on_end = if typeof o.on_end is "function" then o.on_end else noop
    @pause  = if typeof o.pause  is "function" then o.pause  else null
    @env = o.env or o.$ or {} # environment to register variables, especially for child processes
    @paused = false
    @ended  = false
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
    if infbuf.length is 0
      @ended = read_size is @end - @offset
      @pitch += @pitch
      return @_read()
    buf = infbuf
    i_offset = 0
    current_i_offset = i_offsets.shift()
    current_d_offset = @offset + d_offsets.shift()
    loop
      break if buf.length < 4
      bytesize = buf.readInt32LE(0, true) + 4
      break if buf.length < bytesize
      bambuf = buf.slice(0, bytesize)

      bam = new BAM(bambuf, @reader)
      bam.i_offset = i_offset - current_i_offset
      bam.d_offset = current_d_offset
      @on_bam bam, @env

      i_offset += bytesize

      # updating i_offset, d_offset
      loop
        break if i_offsets[0] is undefined or i_offset < i_offsets[0]
        current_i_offset = i_offsets.shift()
        current_d_offset = @offset + d_offsets.shift()
      buf = infbuf.slice(i_offset)
    @offset = current_d_offset
    setImmediate @_read.bind @

module.exports.BAMIterator = BAMIterator
