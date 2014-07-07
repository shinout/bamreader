BAMReader = module.exports
cp = require("child_process")
fs = require("fs")
crypto = require("crypto")
ESTIMATE_DELTA = 0.0001
LSIZE = 12

#####################################
# creates dic
#####################################
BAMReader.prototype.createDic = (op = {}, callback)->
  op = num: op if typeof op is "number"
  outfile = @bamfile + ".dic"
  tmpfiles = []
  merged_num = 0
  merging = 0
  finished = false
  $ =
    WFILE_HWM      : 1024*1024*20-1
    MAX_MEMORY_SIZE: 1.2e9
    tmpfile_inc    : 0
    outfile        : outfile
    r_count        : 0
    w_count        : 0
    time           : new Date/1000|0
    debug          : op.debug
    target         : {}

  # spawn children
  @fork
    $: $
    num: op.num
    pitch: 1024 * 1024 * 4

    # calc md5 and store
    bam: (bam, $)->
      binary = bam.d_offset.toString(2) # binary expression of dOffset
      key = crypto.createHash("md5").update(bam.qname).digest().readUInt32BE(0,  true)
      #bam.qname.match(/[0-9]+/g).join("")
      data = new Buffer(LSIZE)
      data.writeUInt32BE(key, 0, true)
      data.writeUInt32BE(parseInt(binary.slice(-32), 2), 4, true) # lower
      upper = if binary.length > 32 then parseInt(binary.slice(0, -32), 2) else 0
      data.writeUInt32BE((upper << 16) + bam.i_offset, 8, true)
      $.target[key] = [] unless $.target[key]?
      $.target[key].push data
      $.r_count++

    # write to tmpfile
    pause: (ended, $)->
      memory = process.memoryUsage()
      console.log [$.n, "R", $.r_count, (new Date/1000|0) - $.time, memory.rss].join("\t") if $.debug
      return false if not ended and memory.rss <= $.MAX_MEMORY_SIZE
      WCHUNK_SIZE = $.WFILE_HWM - 10000
      w_data = ""
      tmpfile = "#{$.outfile}.#{$.n}.#{(++$.tmpfile_inc)}"
      wstream = require("fs").createWriteStream(tmpfile, highWaterMark: $.WFILE_HWM)
      _write = ->
        console.log [$.n, "W", $.w_count, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug
        for key,arr of $.target
          w_data += arr.map((data)-> data.toString("hex") + "\n").join("")
          $.w_count += arr.length
          delete $.target[key]
          if w_data.length > WCHUNK_SIZE
            wstream.write w_data, "utf-8", _write
            w_data = ""
            return
        console.log [$.n, "W", $.w_count, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug
        wstream.end(w_data)
      wstream.on "finish", =>
        @send tmpfile: tmpfile
        @resume()
      _write()
      return true

    # merges tmpfiles
    message: (msg)->
      return unless msg.tmpfile
      tmpfiles.push msg.tmpfile
      return if tmpfiles.length < 2

      files = tmpfiles.join(" ")
      console.log ["M", "M", tmpfiles.length, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug
      merging++
      merge_sort files, ->
        merging--
        on_finish(finished) if merging is 0 and finished
      tmpfiles = []

    end: ($)->
      console.log [$.n, "E", $.w_count, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug

    finish : ($s)->
      finished = $s
      on_finish(finished) if merging is 0

  # merge sort
  merge_sort = (files, cb)->
    new_name = outfile + ".merged" + (++merged_num)
    command = "sort -m #{files} > #{new_name}"
    sort = cp.exec command, ->
      tmpfiles.push new_name
      cp.exec "rm #{files}", cb
    
  # merge first, then binarize
  on_finish = ($s)->
    if tmpfiles.length >= 2

      # merge sort (pipe)
      console.log ["M", "M", tmpfiles.length, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug
      sort = cp.spawn "sort", ["-m"].concat tmpfiles
      binarize($s, sort.stdout)
    else
      binarize $s, fs.createReadStream(tmpfiles[0], highWaterMark: 1024 * 1024 * 10 -1)

  binarize = ($s, rstream)->
    l_count = 0
    rstream.setEncoding "utf-8"
    wstream = fs.createWriteStream(outfile, highWaterMark: 1024 * 1024 * 10 -1)
    # writes header
    idx_header = {tlen_mean: 0, tlen_sd: 0}
    idx_header_str = JSON.stringify(idx_header)
    header_buf = new Buffer(idx_header_str.length + 4)
    header_buf.writeUInt32BE(idx_header_str.length, 0)
    header_buf.write(idx_header_str, 4)
    wstream.write(header_buf)

    # writes body
    remainder = ""
    ended = false
    read_write = ->
      return if ended
      d = rstream.read()
      return rstream.once "readable", read_write if d is null
      str = remainder + d
      lines = str.split("\n")
      remainder = lines.pop()
      buf = new Buffer(LSIZE * lines.length)
      l_count += lines.length
      console.log ["M", "B", l_count, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug and (l_count % 1e5) < 1000
      buf.write(line, i * LSIZE, "hex") for line,i in lines
      wstream.write buf, read_write

    rstream.once "readable", read_write
    rstream.on "end", ->
      ended = true
      wstream.end()

    wstream.on "finish", ->
      console.log ["M", "B", l_count, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug
      cp.exec "rm #{tmpfiles.join(" ")}", ->
        callback($s) if typeof callback is "function"

BAMReader.prototype.find = (qname, debug)->
  @dic = new BAMDic(@) unless @dic
  if @dic is null
    throw new Error(".dic file has not been created. reader.createDic() can make the file.")
  return @dic.fetch(qname, debug)

class BAMDic
  constructor: (@reader)->
    @idxfile = @reader.bamfile + ".dic"
    return null unless fs.existsSync(@idxfile)
    @size = fs.statSync(@idxfile).size
    @fd = fs.openSync @idxfile, "r"

    # read idx header
    _b = new Buffer(4)
    fs.readSync @fd, _b, 0, 4, 0
    headerJSONLen = _b.readUInt32BE(0)
    _b = new Buffer(headerJSONLen)
    fs.readSync @fd, _b, 0, headerJSONLen, 4
    @header = JSON.parse(_b.toString("utf-8"))
    @header_offset = headerJSONLen + 4
    @total_reads = (@size - @header_offset) / LSIZE
    throw "#{@idxfile} is imcomplete bamdic" if @total_reads isnt parseInt(@total_reads)

  fetch: (qname, debug)->
    inputMD5Buf = crypto.createHash("md5").update(qname).digest()
    # estimated values
    inputMD5 = inputMD5Buf.readUInt32BE(0, true)
    estimates =
      center : Math.max(1, Math.floor(inputMD5 / 0xffffffff * @total_reads))
      delta  : Math.floor(@total_reads * ESTIMATE_DELTA)

    if debug
      console.error "#{k}: #{v}" for k, v of name: qname, md5: inputMD5, dicfile: @idxfile
      console.error "#{k}: #{v}" for k, v of estimates

    if delta < 100
      left = 1
      right = @total_reads + 1
    else
      estimates.left  = Math.max(1, estimates.center - estimates.delta)
      estimates.right = Math.min(estimates.center + estimates.delta, @total_reads)
      # LEFT
      lbuf = new Buffer(LSIZE)
      fs.readSync @fd, lbuf, 0, 4, LSIZE * (estimates.left - 1) + @header_offset
      left = if lbuf.readUInt32BE(0, true) < inputMD5 then estimates.left else 1

      # RIGHT
      rbuf = new Buffer(LSIZE)
      fs.readSync @fd, rbuf, 0, 4, LSIZE * (estimates.right - 1) + @header_offset
      right = if inputMD5 < rbuf.readUInt32BE(0, true) then estimates.right else @total_reads + 1

    if debug
      console.error left: left, right: right, reads: @total_reads

    md5 = null

    # binary search
    iterationCount = 0
    until inputMD5 is md5
      iterationCount++
      currentLineNum = Math.floor((left + right)/2)
      buf = new Buffer(LSIZE)
      offset = LSIZE * (currentLineNum-1)
      fs.readSync @fd, buf, 0, LSIZE, offset + @header_offset
      md5 = buf.readUInt32BE(0, true)
      if debug
        console.error "\t#{k}: #{v}" for k, v of input: inputMD5, current: md5, left: left, right: right, lineNum: currentLineNum

      if md5 > inputMD5
        newright = currentLineNum
        break if newright is right
        right = newright
      else
        newleft = currentLineNum
        break if newleft is left
        left = newleft

    console.error "iteration: #{iterationCount}" if debug

    if md5 isnt inputMD5
      return null
    results = [buf]

    # search flanking lines
    for delta in [1, -1]
      num = currentLineNum
      console.error "linenum", num if debug
      loop
        num += delta
        break if num < 1 or @total_reads < num
        console.error "num", num if debug
        buf = new Buffer(LSIZE)
        fs.readSync @fd, buf, 0, LSIZE, LSIZE * (num-1) + @header_offset
        md5 = buf.readUInt32BE(0, true)
        console.error "md5", md5 if debug
        break if md5 isnt inputMD5
        results.push buf

    bams = []
    for buf in results
      lower = buf.readUInt32BE(4, true)
      upper = buf.readUInt16BE(8, true)
      i_offset = buf.readUInt16BE(10, true)
      d_offset = if upper then upper * 0x100000000 + lower else lower
      console.error "#{k}: #{v}" for k, v of d_offset: d_offset, i_offset: i_offset if debug
      bam = @reader.read(i_offset, d_offset)
      bams.push bam if bam and bam.qname is qname
    return bams

module.exports.BAMDic = BAMDic
