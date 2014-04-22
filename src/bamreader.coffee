###
# BAMReader by Shin Suzuki(@shinout)
####
createReadStream = require("fs").createReadStream
createInflateRaw = require("zlib").createInflateRaw

SEQ_ARR = "=ACMGRSVTWYHKDBN".split("")
CIGAR_ARR = "MIDNSHP=X".split("")
FLAGS = [
  "multiple"
  "allmatches"
  "unmapped"
  "next_unmapped"
  "reversed"
  "next_reversed"
  "first"
  "last"
  "secondary"
  "lowquality"
  "duplicate"
]

class BAMReader
  constructor: (@bamfile, @options={})->
    process.nextTick ()=> @begin()

  @create: (bamfile, options={})->
    return new BAMReader(bamfile, options)

  on: (name, fn)->
    switch name
      when "sam" then @onSam = fn
      when "bam" then @onBam = fn
      when "end" then @onEnd = fn
      when "header" then @onHeader = fn

  begin: ()->
    onBam = @onBam
    onSam = @onSam
    onEnd = @onEnd
    onHeader = @onHeader
    options = @options

    count = 0
    xi = 0
    refs = {}
    curXi = 1
    lastXi = 0
    inflatedBuffers = {}

    inflateRaw = (i, buffer)->
      engine = createInflateRaw chunkSize: 65535
      nread = 0
      buffers = []

      engine.on "error", (err)->
        engine.removeListener "end"
        engine.removeListener "readable", flow
        console.error err

      engine.on "end", ->
        buf = Buffer.concat buffers, nread
        buffers = []
        if i is 0
          readHeader buf
        else
          inflatedBuffers[i] = buf
        engine.close()

        while inflatedBuffer = inflatedBuffers[curXi]
          readAlignment inflatedBuffer, i
          delete inflatedBuffers[curXi]
          curXi++
        onEnd() if onEnd and curXi is lastXi

      flow = ->
        while null isnt (chunk = engine.read())
          buffers.push chunk
          nread += chunk.length
        engine.once "readable", flow

      engine.end buffer
      flow()

    rstream = createReadStream(@bamfile, highWaterMark: 1024 * 1024 - 1)

    remainedBuffer = new Buffer(0)

    rstream.on "data", (newBuffer)->
      buf = Buffer.concat [remainedBuffer, newBuffer], remainedBuffer.length + newBuffer.length
      loop
        if buf.length <= 26
          remainedBuffer = if buf.length then buf else new Buffer(0)
          break

        cdataLen = buf.readUInt16LE(16)- 25
        if buf.length < cdataLen + 26
          remainedBuffer = buf
          break

        # unzip
        cdataBuffer = buf.slice(18, cdataLen + 18)
        inflateRaw xi++, cdataBuffer
        buf = buf.slice(26+cdataLen)

    rstream.on "end", ()->
      lastXi = xi

    # reading bam header
    readHeader = (bambuf)->
      headerLen = bambuf.readInt32LE(4)
      header = bambuf.slice(8,headerLen+8).toString("ascii")
      cursor = headerLen + 8
      nRef = bambuf.readInt32LE cursor
      cursor+=4

      for i in [0...nRef]
        nameLen = bambuf.readInt32LE cursor
        cursor+=4
        name = bambuf.slice(cursor, cursor+nameLen-1).toString("ascii")
        cursor+=nameLen
        refLen = bambuf.readInt32LE cursor
        cursor+=4
        refs[i] = name: name, len: refLen
      onHeader header if onHeader

    # reading bam alignment data
    readAlignment = (buf, k)->
      itr = 0
      while buf.length
        cursor = 0
        blockSize = buf.readInt32LE cursor

        break if buf.length < blockSize
        cursor+=4
        count++

        refId = buf.readInt32LE cursor
        rname = if refId is -1 then "*" else refs[refId].name
        cursor+=4

        pos = buf.readInt32LE cursor
        cursor+=4

        readNameLen = buf.readUInt8 cursor
        cursor++

        mapq = buf.readUInt8 cursor
        cursor++

        bin = buf.readUInt16LE cursor
        cursor+=2

        cigarLen = buf.readUInt16LE cursor
        cursor+=2

        flag = buf.readUInt16LE cursor
        flags = {}
        flags[flagname] = !!(flag & (0x01 << i)) for flagname,i in FLAGS
        cursor+=2

        seqLen = buf.readInt32LE cursor
        cursor+=4

        nextRefId = buf.readInt32LE cursor
        rnext = if nextRefId is -1 then "*" else refs[nextRefId].name
        cursor+=4

        nextPos = buf.readInt32LE cursor
        cursor+=4

        tLen = buf.readInt32LE cursor
        cursor+=4

        readName = buf.slice(cursor, cursor+readNameLen-1).toString("ascii")
        cursor+=readNameLen

        cigar = []
        for i in [0...cigarLen]
          num = buf.readUInt32LE(cursor, cursor+4)
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
              tags[tag] = type: "i", value: buf.readInt8 cursor
              cursor++
            when "C"
              tags[tag] = type: "i", value: buf.readUInt8 cursor
              cursor++
            when "s"
              tags[tag] = type: "i", value: buf.readInt16LE cursor
              cursor+=2
            when "S"
              tags[tag] = type: "i", value: buf.readUInt16LE cursor
              cursor+=2
            when "i"
              tags[tag] = type: "i", value: buf.readInt32LE cursor
              cursor+=4
            when "I"
              tags[tag] = type: "i", value: buf.readUInt32LE cursor
              cursor+=4
            when "f"
              tags[tag] = type: valtype, value: buf.readFloatLE cursor
              cursor+=4
            when "B"
              subtype = String.fromCharCode buf[cursor]
              cursor++
              arrayLen = buf.readInt32LE cursor
              cursor+=4
              switch subtype
                when "c"
                  tags[tag] = type: valtype + ":#{subtype}", value: (buf.readInt8 cursor+i for i in [0...arrayLen])
                  cursor+=arrayLen
                when "C"
                  tags[tag] = type: valtype + ":#{subtype}", value: (buf.readUInt8 cursor+i for i in [0...arrayLen])
                  cursor+=arrayLen
                when "s"
                  tags[tag] = type: valtype + ":#{subtype}", value: (buf.readInt16LE cursor+i*2 for i in [0...arrayLen])
                  cursor+=arrayLen*2
                when "S"
                  tags[tag] = type: valtype + ":#{subtype}", value: (buf.readUInt16LE cursor+i*2 for i in [0...arrayLen])
                  cursor+=arrayLen*2
                when "i"
                  tags[tag] = type: valtype + ":#{subtype}", value: (buf.readInt32LE cursor+i*4 for i in [0...arrayLen])
                  cursor+=arrayLen*4
                when "I"
                  tags[tag] = type: valtype + ":#{subtype}", value: (buf.readUInt32LE cursor+i*4 for i in [0...arrayLen])
                  cursor+=arrayLen*4
                when "f"
                  tags[tag] = type: valtype + ":#{subtype}", value: (buf.readFloatLE cursor+i*4 for i in [0...arrayLen])
                  cursor+=arrayLen*4
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
        bamline =
          qname   : readName
          flag    : flag
          rname   : rname
          pos     : pos+1
          mapq    : mapq
          cigar   : cigar
          rnext   : rnext
          pnext   : nextPos+1
          tlen    : tLen
          seq     : seq
          qual    : qual
          #phredq : phredQuals
          refid   : refId
          nrefid  : nextRefId
          tags    : tags
          start   : pos
          flags   : flags
          #cigars : cigars # unimplemented
          #count  : count
          tagstr  : ([name, tag.type, tag.value].join(":") for name,tag of tags).join("\t")

        onBam bamline if onBam

        if onSam
          samline = [
            bamline.qname
            bamline.flag
            bamline.rname
            bamline.pos
            bamline.mapq
            bamline.cigar || "*"
            if bamline.rnext is bamline.rname and bamline.rname isnt "*" then "=" else bamline.rnext
            bamline.pnext
            bamline.tlen
            bamline.seq
            bamline.qual
            bamline.tagstr
          ].join("\t")
          onSam samline

module.exports = BAMReader
