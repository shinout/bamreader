class BamStream
  constructor: (options = {})->
    @refs = options.refs or {}
    @readingHeader = if options.refs then false else true #skip reading header when refs is given (mainly for bamdic)
    @remainedDefBuf = new Buffer(0)
    @remainedInfBuf = new Buffer(0)
    @defBufOffset = options.start
    @infBufOffset = 0
    @deltaDefBuf  = 0

  _transform: (chunk, encoding, done)->
    defBuf = Buffer.concat [@remainedDefBuf, chunk], @remainedDefBuf.length + chunk.length
    # split deflated buffer
    [bufsToInflate, @remainedDefBuf] = BAMReader.splitDeflatedBuffer(defBuf)

    for bufToInflate in bufsToInflate
      infBufChunk = inflateRawSync bufToInflate
      infBuf = if @remainedInfBuf.length then Buffer.concat [@remainedInfBuf, infBufChunk] else infBufChunk
      # read header
      if @readingHeader
        try
          headerInfo = BAMReader.readHeaderFromInflatedBuffer infBuf, true
          {refs, headerStr, infBuf} = headerInfo
          @readingHeader = false
          @push headerInfo
          @defBufOffset = @deltaDefBuf
          @deltaDefBuf = 0
        catch e
          @deltaDefBuf += bufToInflate.length + 26
          @remainedInfBuf = infBuf
          continue
      # read alignments
      loop
        if infBuf.length is 0
          @remainedInfBuf = infBuf
          break
        [bam, infBuf] = BAMReader.readAlignmentFromInflatedBuffer infBuf, refs, reader
        if bam is null
          @remainedInfBuf = infBuf
          break
        else
          @push [bam, @defBufOffset, @infBufOffset]
          @infBufOffset = infBufChunk.length - infBuf.length
          onSam bam.sam if onSam
          if @deltaDefBuf.length
            @defBufOffset += @deltaDefBuf
            @deltaDefBuf = 0

      @deltaDefBuf += bufToInflate.length + 26
      if @remainedInfBuf.length is 0
        @defBufOffset += @deltaDefBuf
        @deltaDefBuf = 0
        @infBufOffset = 0

      done()
    
require("util").inherits(BamStream, require("stream").Transform)
