BGZF_ESTIMATED_LEN = 65536
BGZF_HEADER = new Buffer("1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00".split(" ").join(""), "hex")
SEQ_ARR = "=ACMGRSVTWYHKDBN".split("")
CIGAR_ARR = "MIDNSHP=X".split("")

fs = require("fs")
childProcess = require("child_process")
inflateRawSync = require("zlib-raw-sync").inflateRawSync

class BAMReader
  constructor: (@bamfile, @options={})->
    @bamfile.pause() if @bamfile.readable
    try
      @dic = require("bamdic").create(if @bamfile.readable then @options.bamfile else @bamfile)
    catch e
    try
      tlenJSON = require(@options.tlenInfo)
      @tlen_sd   = tlenJSON.purified.sd
      @tlen_mean = tlenJSON.purified.mean
    catch e

    return @ if options.wait

    childProcess.exec "which samtools", (e, stdout,stderr)=>
      return @beginSamtools() if @options.samtools
      if not @options.sam and (e or stderr or @options.native) then @begin() else @beginSamtools(@options.sam)

  @create: (bamfile, options={})->
    return new BAMReader(bamfile, options)

  on: (name, fn)->
    switch name
      when "sam" then @onSam = fn
      when "bam" then @onBam = fn
      when "end" then @onEnd = fn
      when "header" then @onHeader = fn
    return @

  beginSamtools: (isSam)->
    reader = @
    Bam = module.exports.Bam
    samtoolsCmd = @options.samtools or "samtools"
    onBam = @onBam
    onSam = @onSam
    onEnd = @onEnd
    onHeader = @onHeader
    options = @options


    if isSam
      samtools = null
      lines = require("linestream").create(@bamfile)
    else
      file = if @bamfile.readable then "-" else @bamfile
      samtools = childProcess.spawn samtoolsCmd, ["view", "-h", file]
      lines = require("linestream").create(samtools.stdout)

    readingHeader = true
    headerLines = []

    lines.on "end", onEnd if onEnd

    lines.on "data", (samline)->
      if readingHeader
        if samline.charAt(0) is '@'
          headerLines.push samline
          return
        else
          readingHeader = false
          onHeader headerLines.join("\n") if onHeader
          headerLines = null

      onSam samline if onSam
      onBam Bam.createFromSAM samline, reader if onBam

    if @bamfile.readable
      @bamfile.pipe samtools.stdin if samtools
      @bamfile.resume()

  begin: ()->
    reader = @
    onBam = @onBam
    onSam = @onSam
    onEnd = @onEnd
    onHeader = @onHeader
    options = @options

    if @bamfile.readable
      rstream = @bamfile
      rstream.resume()
    else
      rstream = fs.createReadStream(@bamfile, highWaterMark: 1024 * 1024 - 1)

    refs = {}
    bambuf4h = new Buffer(0)
    readingHeader = true
    remainedBuffer = new Buffer(0)

    # read deflated buffers
    _read = (newBuffer)->
      buf = Buffer.concat [remainedBuffer, newBuffer], remainedBuffer.length + newBuffer.length
      # split deflated buffer
      [defBufs, remainedBuffer] = BAMReader.splitDeflatedBuffer(buf)

      for defBuf in defBufs
        bambuf = inflateRawSync defBuf
        # read header
        if readingHeader
          bambuf4h = Buffer.concat [bambuf4h, bambuf]
          try
            {refs, headerStr, bambuf} = BAMReader.readHeaderFromInflatedBuffer bambuf4h, true
            readingHeader = false
            onHeader headerStr if onHeader
            continue if bambuf.length is 0
          catch e
            continue
        # read alignments
        bams = BAMReader.readAlignmentsFromInflatedBuffer bambuf, refs, false, reader
        onBam bam for bam in bams if onBam
        onSam bam.sam for bam in bams if onSam

    rstream.on "data", _read

    rstream.on "end", ()->
      # read remained buffer
      _read(remainedBuffer)
      onEnd() if onEnd

  # read header from a bamfile
  @readHeader = (bamfile)->
    infBuf = new Buffer(0)
    offset = 0
    fd = fs.openSync bamfile, "r"

    loop
      [bufToInflate, next] = BAMReader.getDeflatedBuffer(fd, offset)
      offset = next

      _infBuf = inflateRawSync bufToInflate
      infBuf = Buffer.concat [infBuf, _infBuf]

      try
        headerInfo = BAMReader.readHeaderFromInflatedBuffer(infBuf)
        headerInfo.offset = offset
        headerInfo.fd = fd
        break
      catch e
    return headerInfo

  @splitBody = (bamfile, num, headerInfo)->
    headerInfo = headerInfo or BAMReader.readHeader(bamfile)
    size = (fs.statSync bamfile).size
    offset = headerInfo.offset
    fd = headerInfo.fd or fs.openSync(bamfile, "r")
    interval = Math.floor((size-offset)/num)
    positions = []

    buflen = Math.min(BGZF_ESTIMATED_LEN, interval)

    for k in [0...num]
      # finding accurate position of BGZF
      start = interval * k + offset-1
      buf = new Buffer(buflen)
      fs.readSync fd, buf, 0, buflen, start
      cursor = -1
      match = false
      until match or cursor + 16 > buf.length
        cursor++
        headerCandidate = buf.slice(cursor, cursor+16)
        match = true
        for b,i in BGZF_HEADER
          if b isnt headerCandidate[i]
            match = false
            break
      positions.push(start + cursor) if match
    fs.closeSync(fd)
    return positions: positions, size: size, header: headerInfo

  @getDeflatedBuffer = (fd, offset)->
    defBuf = new Buffer(0)
    k = 0
    loop
      _defBuf = new Buffer(BGZF_ESTIMATED_LEN)
      fs.readSync fd, _defBuf, 0, BGZF_ESTIMATED_LEN, offset + k * BGZF_ESTIMATED_LEN
      for i in [0...16]
        throw new Error("not BGZF (offset=#{offset}, i=#{i})") if _defBuf[i] isnt BGZF_HEADER[i]
      defBuf = Buffer.concat [defBuf, _defBuf]
      delta = defBuf.readUInt16LE(16, true) + 1
      break if defBuf.length >= delta
      k++

    bufToInflate = defBuf.slice(18, delta-8)
    return [bufToInflate, offset + delta]

  @splitDeflatedBuffer = (defBuf)->
    defBufs = []
    loop
      return [defBufs,defBuf] if defBuf.length <= 26
      cdataLen = defBuf.readUInt16LE(16, true)- 25
      return [defBufs,defBuf] if defBuf.length < cdataLen + 26
      # unzip
      defBufs.push defBuf.slice(18, cdataLen + 18)
      defBuf = defBuf.slice(26+cdataLen)

  # reading bam header
  @readHeaderFromInflatedBuffer = (bambuf, ifReturnsBamBuf)->
    refs = {}
    headerLen = bambuf.readInt32LE(4, true)
    throw new Error("header len") if bambuf.length < headerLen + 16
    headerStr = bambuf.slice(8,headerLen+8).toString("ascii")
    cursor = headerLen + 8
    nRef = bambuf.readInt32LE cursor, true
    cursor+=4

    for i in [0...nRef]
      nameLen = bambuf.readInt32LE cursor, true
      cursor+=4
      name = bambuf.slice(cursor, cursor+nameLen-1).toString("ascii")
      cursor+=nameLen
      refLen = bambuf.readInt32LE cursor, true
      cursor+=4
      refs[i] = name: name, len: refLen

    ret = refs: refs, headerStr: headerStr
    ret.bambuf = bambuf.slice(cursor) if ifReturnsBamBuf
    return ret

  # reading bam alignment data
  @readAlignmentsFromInflatedBuffer = (buf, refs, readFirst, reader)->
    Bam = module.exports.Bam
    bams = []
    while buf.length
      cursor = 0
      blockSize = buf.readInt32LE cursor, true

      break if buf.length < blockSize
      cursor+=4

      refId = buf.readInt32LE cursor, true
      rname = if refId is -1 then "*" else refs[refId].name
      cursor+=4

      pos = buf.readInt32LE cursor, true
      cursor+=4

      readNameLen = buf.readUInt8 cursor, true
      cursor++

      mapq = buf.readUInt8 cursor, true
      cursor++

      bin = buf.readUInt16LE cursor, true
      cursor+=2

      cigarLen = buf.readUInt16LE cursor, true
      cursor+=2

      flag = buf.readUInt16LE cursor, true
      cursor+=2

      seqLen = buf.readInt32LE cursor, true
      cursor+=4

      nextRefId = buf.readInt32LE cursor, true
      rnext = if nextRefId is -1 then "*" else refs[nextRefId].name
      cursor+=4

      nextPos = buf.readInt32LE cursor, true
      cursor+=4

      tLen = buf.readInt32LE cursor, true
      cursor+=4

      readName = buf.slice(cursor, cursor+readNameLen-1).toString("ascii")
      cursor+=readNameLen

      cigar = []
      for i in [0...cigarLen]
        num = buf.readUInt32LE(cursor, true)
        char = CIGAR_ARR[num & 0x0f]
        num = num>>4
        cigar.push num + char
        cursor+=4
      cigar = cigar.join("")

      seqLenByte = Math.floor((seqLen+1)/2)

      seqBits = buf.slice(cursor, cursor+seqLenByte)
      seq = []
      for byte in seqBits
        seq.push SEQ_ARR[byte >>4]
        second = SEQ_ARR[byte & 0x0F]
        seq.push second if second isnt "="
      seq = seq.join("")
      cursor+=seqLenByte

      #phredQuals = buf.slice(cursor, cursor+seqLen).toString("hex")
      qual = (String.fromCharCode(buf[cursor+i]+33) for i in [0...seqLen]).join("")
      cursor+=seqLen

      tags = {}
      while true
        break if cursor-4 >= blockSize
        tag = buf.slice(cursor, cursor+2).toString("ascii")
        cursor+=2
        valtype = String.fromCharCode buf[cursor]
        cursor++

        switch valtype
          when "A"
            tags[tag] = type: valtype, value: String.fromCharCode buf[cursor]
            cursor++
          when "c"
            tags[tag] = type: "i", value: buf.readInt8 cursor, true
            cursor++
          when "C"
            tags[tag] = type: "i", value: buf.readUInt8 cursor, true
            cursor++
          when "s"
            tags[tag] = type: "i", value: buf.readInt16LE cursor, true
            cursor+=2
          when "S"
            tags[tag] = type: "i", value: buf.readUInt16LE cursor, true
            cursor+=2
          when "i"
            tags[tag] = type: "i", value: buf.readInt32LE cursor, true
            cursor+=4
          when "I"
            tags[tag] = type: "i", value: buf.readUInt32LE cursor, true
            cursor+=4
          when "f"
            tags[tag] = type: valtype, value: buf.readFloatLE cursor, true
            cursor+=4
          when "B"
            subtype = String.fromCharCode buf[cursor]
            cursor++
            arrayLen = buf.readInt32LE cursor, true
            cursor+=4
            switch subtype
              when "c"
                tags[tag] = type: valtype, value: (buf.readInt8 cursor+i, true for i in [0...arrayLen])
                cursor+=arrayLen
              when "C"
                tags[tag] = type: valtype, value: (buf.readUInt8 cursor+i, true for i in [0...arrayLen])
                cursor+=arrayLen
              when "s"
                tags[tag] = type: valtype, value: (buf.readInt16LE cursor+i*2, true for i in [0...arrayLen])
                cursor+=arrayLen*2
              when "S"
                tags[tag] = type: valtype, value: (buf.readUInt16LE cursor+i*2, true for i in [0...arrayLen])
                cursor+=arrayLen*2
              when "i"
                tags[tag] = type: valtype, value: (buf.readInt32LE cursor+i*4, true for i in [0...arrayLen])
                cursor+=arrayLen*4
              when "I"
                tags[tag] = type: valtype, value: (buf.readUInt32LE cursor+i*4, true for i in [0...arrayLen])
                cursor+=arrayLen*4
              when "f"
                tags[tag] = type: valtype, value: (buf.readFloatLE cursor+i*4, true for i in [0...arrayLen])
                cursor+=arrayLen*4
            value.unshift subtype

          when "Z"
            zLen = 0
            zLen++ while buf[cursor+zLen] isnt 0x00
            tags[tag] = type: valtype, value: buf.slice(cursor, cursor+zLen).toString("ascii")
            cursor+=zLen+1
          when "H"
            hLen = 0
            hLen++ while buf[cursor+hLen] isnt 0x00
            tags[tag] = type: valtype, value: buf.slice(cursor, cursor+hLen).toString("hex")
            cursor+=hLen+1

      buf = buf.slice cursor

      # output
      bams.push new Bam(
        reader or null,
        readName,
        flag,
        rname,
        pos+1,
        mapq,
        cigar,
        rnext,
        nextPos+1,
        tLen,
        seq,
        qual
      )
      bams.tags_ = tags
      return bams[0] if readFirst
    return bams

module.exports = BAMReader
