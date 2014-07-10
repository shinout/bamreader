BAM = module.exports.BAM
cp = require "child_process"
fs = require "fs"
class SAMTools
  constructor: (@reader, o)->
    o = on_bam: o if typeof o is "function"
    o.on_bam = o.bam if o.bam
    if typeof o.end is "function"
      o.on_end = o.end
      delete o.end
    @on_bam = if typeof o.on_bam is "function" then o.on_bam  else ->
    @start = o.start if typeof o.start is "number"
    @end   = o.end   if typeof o.end   is "number"
    @on_end = if typeof o.on_end is "function" then o.on_end else ->
    @env = o.env or o.$ or {} # environment to register variables, especially for child processes
    process.nextTick => @view()

  view: ->
    # 1. if "start" and "end" are given, pipe to stdin
    # 2. otherwise, spawn with file
    if @start? and @end?
      samtools = cp.spawn "samtools", ["view", "-"]
      header_chunk = new Buffer(@reader.header_offset)
      fs.readSync @reader.fd, header_chunk, 0, @reader.header_offset, 0
      samtools.stdin.write header_chunk # send header
      fstream = fs.createReadStream(@reader.bamfile, start: @start, end: @end).pipe samtools.stdin
    else
      samtools = cp.spawn "samtools", ["view", @reader.bamfile]

    rstream = samtools.stdout
    rstream.setEncoding "utf-8"
    _r = ""
    rstream.on "readable", =>
      chunk = rstream.read()
      return if chunk is null
      sams = (_r + chunk).split("\n")
      _r = sams.pop()
      for sam in sams
        bam = BAM.createFromSAM sam, @reader
        @on_bam bam, @env

    rstream.on "end", =>
      @on_end @env

module.exports.SAMTools = SAMTools

BAMReader = module.exports
BAMReader.prototype.samtools = (o)->
  return new SAMTools(@, o)

BAMReader.prototype.fork_samtools = (o={})->
  o = on_bam : o if typeof o is "function"
  o.script = "child_samtools"
  return @fork(o)
