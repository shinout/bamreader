require("termcolor").define
BAMReader = module.exports
cp = require("child_process")
fs = require("fs")
crypto = require("crypto")
LINE_SIZE = 11
DIC_SIZE  = 8

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
  tlen_sample_size = if typeof op.tlen_sample_size is "number" then op.tlen_sample_size else 100000
  outlier_rate = if typeof op.outlier_rate is "number" then op.outlier_rate  else 0.02

  $ =
    WFILE_HWM        : 1024*1024*20-1
    MAX_MEMORY_SIZE  : 1.2e9
    tmpfile_inc      : 0
    outfile          : outfile
    r_count          : 0
    w_count          : 0
    # prev_count       : 0
    time             : new Date/1000|0
    debug            : op.debug
    pool             : {}
    pool_count       : 0
    outlier_rate     : outlier_rate
    tlen_sample_size : Math.round tlen_sample_size / op.num
    d_deltas         : [] # broad defbuf lengths
    last_offset      : 0

  # spawn children
  @fork
    $: $
    num: op.num
    nocache: true
    pitch: 8388608

    start: ($)->
      $.tlens = new BAMReader.OutlierFilteredMeanDev($.outlier_rate, $.tlen_sample_size)
      $.d_deltas.push @offset
      $.last_offset = @offset
      
    # calc md5 and store
    bam: (bam, $)->
      binary = bam.d_offset.toString(2) # binary expression of dOffset
      key = crypto.createHash("md5").update(bam.qname).digest().readUInt32BE(0,  true)
      #bam.qname.match(/[0-9]+/g).join("")
      data = new Buffer(LINE_SIZE)
      data.writeUInt32BE(key, 0, true)
      data.writeUInt32BE(parseInt(binary.slice(-32), 2), 4, true) # lower
      upper = if binary.length > 32 then parseInt(binary.slice(0, -32), 2) else 0
      data.writeUInt8(upper, 8, true)
      data.writeUInt16BE(bam.i_offset, 9, true)
      unless $.pool[key]?
        $.pool[key] = []
        $.pool_count++
      $.pool[key].push data

      # mean tlen
      if bam.unmapped is false and bam.next_unmapped is false and bam.same_strand is false and bam.tlen isnt 0
        tlen = Math.abs(bam.tlen)
        $.tlens.add tlen
      $.r_count++

    # write to tmpfile
    pause: ($)->
      $.d_deltas.push(@offset - $.last_offset)
      $.last_offset = @offset
      memory = process.memoryUsage()
      console.log [$.n, "R", $.r_count, (new Date/1000|0) - $.time, memory.rss].join("\t") if $.debug
      return false if memory.rss <= $.MAX_MEMORY_SIZE
      if $.pool_count is 0
        setTimeout =>
          @resume()
        ,1000
        return true
      @write($, @resume.bind @)
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
      @write($, $.exit)
      console.log [$.n, "E", $.w_count, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug
      {sum, squared, n} = $.tlens.precalc()
      $.tlen_sum     = sum
      $.tlen_squared = squared
      $.tlen_n       = n
      delete $.tlens

    props:
      write: ($, cb)->
        #count = $.r_count - $.prev_count
        #$.prev_count = $.r_count
        WCHUNK_SIZE = $.WFILE_HWM - 10000
        w_data = ""
        tmpfile = "#{$.outfile}.#{$.n}.#{(++$.tmpfile_inc)}"
        wstream = require("fs").createWriteStream(tmpfile, highWaterMark: $.WFILE_HWM)
        _write = ->
          console.log [$.n, "W", $.w_count, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug
          for key,arr of $.pool
            for data in arr
              w_data += data.toString("hex") + "\n"
            $.w_count += arr.length
            $.pool_count--
            delete $.pool[key]
            if w_data.length > WCHUNK_SIZE
              wstream.write w_data, "utf-8", _write
              w_data = ""
              return
          $.pool = {}
          console.log [$.n, "W", $.w_count, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug
          wstream.end(w_data)
        wstream.on "finish", =>
          @send tmpfile: tmpfile
          cb()
        _write()

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
    # calculates broad d_offset position
    d_deltas = []
    delta = 0
    for $,i in $s
      s = $.d_deltas.shift()
      d_deltas.push s if i is 0
      for d in $.d_deltas
        delta += d
        if delta > 16777215 # 3byte
          d_deltas.push(delta - d) if d isnt delta
          delta = d
    d_deltas.push(delta) if delta isnt 0

    # calculates tlen statistics information
    tlen_sum     = 0
    tlen_squared = 0
    tlen_n       = 0
    for $ in $s
      tlen_sum     += $.tlen_sum
      tlen_squared += $.tlen_squared
      tlen_n       += $.tlen_n
    tlen_mean = tlen_sum / tlen_n
    tlen_dev  = tlen_squared / tlen_n - tlen_mean * tlen_mean
    tlen_sd   = Math.sqrt tlen_dev
    console.log ["M", "T", Math.round(tlen_mean), (new Date/1000|0) - $.time, "mean"].join("\t") if $.debug
    console.log ["M", "T", Math.round(tlen_sd), (new Date/1000|0) - $.time, "sd"].join("\t") if $.debug

    l_count = 0
    rstream.setEncoding "utf-8"
    wstream = fs.createWriteStream(outfile, highWaterMark: 1024 * 1024 * 10 -1)
    # writes header
    idx_header = {tlen_mean: Math.round(tlen_mean), tlen_sd: Math.round(tlen_sd), tlen_n: tlen_n, outlier_rate: $.outlier_rate, tlen_sample_size: $.tlen_sample_size}
    idx_header_str = JSON.stringify(idx_header)
    header_buf = new Buffer(idx_header_str.length + 4)
    header_buf.writeUInt32BE(idx_header_str.length, 0)
    header_buf.write(idx_header_str, 4)
    wstream.write(header_buf)

    # write d_delta info
    d_delta_len = d_deltas.length
    d_delta_buf = new Buffer( 2 + 3 * d_delta_len)
    d_delta_buf.writeUInt16BE(d_deltas.length, 0)
    offset = 2
    for d_delta in d_deltas
      d_delta_buf.writeUInt16BE(d_delta>>8, offset)
      d_delta_buf.writeUInt8(d_delta&0xff, offset + 2)
      offset += 3
    wstream.write(d_delta_buf)

    # footer info
    three_byte_idx = new Array(256 * 256 * 256)

    # writes body
    remainder = ""
    write_ended = false
    read_write = ->
      return if write_ended
      d = rstream.read()
      return rstream.once "readable", read_write if d is null
      str = remainder + d
      lines = str.split("\n")
      remainder = lines.pop()
      buf = new Buffer(DIC_SIZE * lines.length)
      for line, i in lines
        idx = parseInt(line.slice(0, 6), 16)
        if three_byte_idx[idx]
          three_byte_idx[idx]++
        else
          three_byte_idx[idx] = 1
        buf.write(line.slice(6), i * DIC_SIZE, "hex")
      l_count++
      console.log ["M", "B", l_count, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug and l_count % 10000 is 0
      wstream.write buf, read_write

    rstream.once "readable", read_write
    rstream.on "end", ->
      # write footer
      footer_buf = new Buffer(256 * 256 * 256 * 7)
      offset = 0
      for v, idx in three_byte_idx
        continue unless v
        footer_buf.writeUInt16BE(idx>>8, offset)
        footer_buf.writeUInt8(idx&0xff, offset + 2)
        footer_buf.writeUInt32BE(v, offset + 3)
        offset += 7
      footer_buf.writeUInt32BE(offset, offset)

      wstream.end(footer_buf.slice(0, offset+4))
      write_ended = true

    wstream.on "finish", ->
      console.log ["M", "B", l_count, (new Date/1000|0) - $.time, process.memoryUsage().rss].join("\t") if $.debug
      cp.exec "rm #{tmpfiles.join(" ")}", ->
        callback($s) if typeof callback is "function"

BAMReader.prototype.find = (qname, d_offset_to_filter)->
  if @dic is null
    throw new Error(".dic file has not been created. reader.createDic() can make the file.")
  unless @dic.three_byte_idx
    @dic._read_footer()
  return @dic.fetch(qname, d_offset_to_filter)

class BAMDic
  @create = (reader)->
    idxfile = reader.bamfile + ".dic"
    return null unless fs.existsSync(idxfile)
    return new BAMDic(reader)

  constructor: (@reader)->
    @idxfile = @reader.bamfile + ".dic"
    @size = fs.statSync(@idxfile).size
    @fd = fs.openSync @idxfile, "r"

    # read idx header
    _b = new Buffer(4)
    fs.readSync @fd, _b, 0, 4, 0
    headerJSONLen = _b.readUInt32BE(0)
    _b = new Buffer(headerJSONLen)
    fs.readSync @fd, _b, 0, headerJSONLen, 4
    @header = JSON.parse(_b.toString("utf-8"))
    header_offset = headerJSONLen + 4

    # read d_deltas
    _b = new Buffer(2)
    fs.readSync @fd, _b, 0, 2, header_offset
    d_delta_len = _b.readUInt16BE(0)
    d_delta_buflen = d_delta_len * 3
    _b = new Buffer(d_delta_buflen)
    fs.readSync @fd, _b, 0, d_delta_buflen, header_offset + 2
    cursor = 0
    broad_d_offsets = new Array(d_delta_len + 1)
    pos = 0
    while cursor < d_delta_len
      c3 = cursor * 3
      d_delta = _b.readUInt16BE(c3) * 256 + _b.readUInt8(c3 + 2)
      pos += d_delta
      broad_d_offsets[cursor] = pos
      cursor++
    broad_d_offsets[cursor] = @reader.size
    @broad_d_offsets = broad_d_offsets
    @header_offset = header_offset + 2 + d_delta_buflen

    # calc total reads
    _b = new Buffer(4)
    fs.readSync(@fd, _b, 0, 4, @size - 4)
    @footer_size = _b.readUInt32BE(0)
    @total_reads = (@size - @header_offset - @footer_size - 4) / DIC_SIZE

  _read_footer: ->
    footer_size = @footer_size
    footer = new Buffer(footer_size)
    fs.readSync(@fd, footer, 0, footer_size, @size - footer_size - 4)
    i = 0
    @three_byte_idx = {}
    if (footer_size / 7) < 256 * 256 * 64
      @nums = {}
    total = 0
    while i < footer_size
      idx3byte = footer.readUInt16BE(i) * 256 +  footer.readUInt8(i+2)
      num = footer.readUInt32BE(i+3)
      @three_byte_idx[idx3byte] = total
      @nums[idx3byte] = num if @nums
      total += num
      i+=7
    @three_byte_idx[idx3byte+1] = total # saving the last position

    @bufs = new module.exports.Fifo(1024 * 1024 * 4)

  fetch: (qname, d_offset_to_filter)->
    md5_buf = crypto.createHash("md5").update(qname).digest()
    idx = md5_buf.readUInt16BE(0) * 256 + md5_buf.readUInt8(2)
    obi = md5_buf.readUInt8(3)

    if buf = @bufs.get idx
      read_num = buf.length / DIC_SIZE
    else
      start = @three_byte_idx[idx]
      return null unless start?
      if @nums
        read_num = @nums[idx]
      else
        nx_idx = idx + 1
        while nx_idx <= 16777216 # 256 * 256 * 256
          break if end = @three_byte_idx[nx_idx]
          nx_idx++
        return null unless end?
        read_num = end - start
      buf = new Buffer(DIC_SIZE * read_num)
      fs.readSync @fd, buf, 0, DIC_SIZE * read_num, start * DIC_SIZE + @header_offset
      @bufs.set idx, buf

    if read_num is 1
    # shortcut when only hits one line
      if obi is buf.readUInt8(0, true)
        results = [0]
      else
        return null
    # full scanning
    else if read_num <= 4
      _i = 0
      results = []
      while _i < read_num
        _p = _i * DIC_SIZE
        results.push _i if obi is buf.readUInt8(_p, true)
        _i++
    # binary search
    else
      left = 0
      right = read_num

      md5_o = null
      itr_count = 0
      loop
        itr_count++
        current = Math.floor((left + right)/2)
        md5_o = buf.readUInt8(DIC_SIZE * current, true)
        break if obi is md5_o
        if md5_o > obi
          newright = current
          break if newright is right
          right = newright
        else
          newleft = current
          break if newleft is left
          left = newleft

      return null if md5_o isnt obi

      results = [current]

      # search flanking lines
      for delta in [1, -1]
        num = current
        loop
          num += delta
          break if num < 0 or num >= read_num
          md5_o = buf.readUInt8(DIC_SIZE * num, true)
          break if md5_o isnt obi
          results.push num

    bams = []
    for line_i in results
      _offset = line_i * DIC_SIZE
      lower = buf.readUInt32BE(_offset + 1, true)
      upper = buf.readUInt8(_offset + 5, true)
      d_offset = if upper then upper * 0x100000000 + lower else lower
      continue if d_offset is d_offset_to_filter
      i_offset = buf.readUInt16BE(_offset + 6, true)
      bam = @reader.read(i_offset, d_offset)
      bams.push bam if bam and bam.qname is qname
    return bams

module.exports.BAMDic = BAMDic
