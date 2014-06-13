
/*
 * BAMReader by Shin Suzuki(@shinout)
 */

(function() {
  var BAMReader, BGZF_ESTIMATED_LEN, BGZF_HEADER, CIGAR_ARR, FLAGS, SEQ_ARR, childProcess, createInflateRaw, createReadStream, fs;

  BGZF_ESTIMATED_LEN = 65536;

  BGZF_HEADER = new Buffer("1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00".split(" ").join(""), "hex");

  fs = require("fs");

  createReadStream = fs.createReadStream;

  createInflateRaw = require("zlib").createInflateRaw;

  childProcess = require("child_process");

  require("termcolor").define;

  SEQ_ARR = "=ACMGRSVTWYHKDBN".split("");

  CIGAR_ARR = "MIDNSHP=X".split("");

  FLAGS = ["multiple", "allmatches", "unmapped", "next_unmapped", "reversed", "next_reversed", "first", "last", "secondary", "lowquality", "duplicate", "supplementary"];

  BAMReader = (function() {
    function BAMReader(bamfile, options) {
      this.bamfile = bamfile;
      this.options = options != null ? options : {};
      if (this.bamfile.readable) {
        this.bamfile.pause();
      }
      childProcess.exec("which samtools", (function(_this) {
        return function(e, stdout, stderr) {
          if (_this.options.samtools) {
            _this.beginSamtools(_this.options.sam);
          }
          if (!_this.options.sam && (e || stderr || _this.options["native"])) {
            return _this.begin();
          } else {
            return _this.beginSamtools(_this.options.sam);
          }
        };
      })(this));
    }

    BAMReader.create = function(bamfile, options) {
      if (options == null) {
        options = {};
      }
      return new BAMReader(bamfile, options);
    };

    BAMReader.prototype.on = function(name, fn) {
      switch (name) {
        case "sam":
          return this.onSam = fn;
        case "bam":
          return this.onBam = fn;
        case "end":
          return this.onEnd = fn;
        case "header":
          return this.onHeader = fn;
      }
    };

    BAMReader.prototype.beginSamtools = function(isSam) {
      var file, headerLines, lines, onBam, onEnd, onHeader, onSam, options, readingHeader, samtools, samtoolsCmd;
      samtoolsCmd = this.options.samtools || "samtools";
      onBam = this.onBam;
      onSam = this.onSam;
      onEnd = this.onEnd;
      onHeader = this.onHeader;
      options = this.options;
      if (isSam) {
        samtools = null;
        lines = require("linestream").create(this.bamfile);
      } else {
        file = this.bamfile.readable ? "-" : this.bamfile;
        samtools = childProcess.spawn(samtoolsCmd, ["view", "-h", file]);
        lines = require("linestream").create(samtools.stdout);
      }
      readingHeader = true;
      headerLines = [];
      if (onEnd) {
        lines.on("end", onEnd);
      }
      lines.on("data", function(samline) {
        var bamline, flagname, i, sam, subtype, tag, type, v, val, value, _i, _j, _len, _len1, _ref;
        if (readingHeader) {
          if (samline.charAt(0) === '@') {
            headerLines.push(samline);
            return;
          } else {
            readingHeader = false;
            if (onHeader) {
              onHeader(headerLines.join("\n"));
            }
            headerLines = null;
          }
        }
        if (onSam) {
          onSam(samline);
        }
        if (onBam) {
          sam = samline.split("\t");
          bamline = {
            qname: sam[0],
            flag: Number(sam[1]),
            rname: sam[2],
            pos: Number(sam[3]),
            mapq: Number(sam[4]),
            cigar: sam[5],
            rnext: sam[6],
            pnext: Number(sam[7]) + 1,
            tlen: Number(sam[8]),
            seq: sam[9],
            qual: sam[10],
            tags: {},
            start: Number(sam[3]) - 1,
            flags: {},
            tagstr: sam.slice(11).join("\t")
          };
          for (i = _i = 0, _len = FLAGS.length; _i < _len; i = ++_i) {
            flagname = FLAGS[i];
            bamline.flags[flagname] = !!(bamline.flag & (0x01 << i));
          }
          _ref = sam.slice(11);
          for (_j = 0, _len1 = _ref.length; _j < _len1; _j++) {
            tag = _ref[_j];
            val = tag.split(":");
            tag = val[0];
            type = val[1];
            switch (type) {
              case "i":
              case "f":
                value = Number(val[2]);
                break;
              case "B":
                value = val[2].split(",");
                subtype = value[0];
                if (subtype === "c" || subtype === "C" || subtype === "s" || subtype === "S" || subtype === "i" || subtype === "I" || subtype === "f") {
                  value = (function() {
                    var _k, _len2, _results;
                    _results = [];
                    for (_k = 0, _len2 = value.length; _k < _len2; _k++) {
                      v = value[_k];
                      _results.push(Number(v));
                    }
                    return _results;
                  })();
                  value[0] = subtype;
                }
                break;
              default:
                value = val[2];
            }
            bamline.tags[tag] = {
              type: type,
              value: value
            };
          }
          return onBam(bamline);
        }
      });
      if (this.bamfile.readable) {
        if (samtools) {
          this.bamfile.pipe(samtools.stdin);
        }
        return this.bamfile.resume();
      }
    };

    BAMReader.prototype.begin = function() {
      var bambuf4h, currentIdx, inflatedBuffers, lastIdx, onBam, onEnd, onHeader, onSam, options, readIdx, readInflatedBuffers, readingHeader, refs, remainedBuffer, rstream;
      onBam = this.onBam;
      onSam = this.onSam;
      onEnd = this.onEnd;
      onHeader = this.onHeader;
      options = this.options;
      if (this.bamfile.readable) {
        rstream = this.bamfile;
        rstream.resume();
      } else {
        rstream = createReadStream(this.bamfile, {
          highWaterMark: 1024 * 1024 - 1
        });
      }
      refs = {};
      currentIdx = 0;
      readIdx = 0;
      lastIdx = 0;
      inflatedBuffers = {};
      remainedBuffer = new Buffer(0);
      rstream.on("data", function(newBuffer) {
        var buf, defBuf, defBufs, _i, _len, _ref, _results;
        buf = Buffer.concat([remainedBuffer, newBuffer], remainedBuffer.length + newBuffer.length);
        _ref = BAMReader.splitDeflatedBuffer(buf), defBufs = _ref[0], remainedBuffer = _ref[1];
        _results = [];
        for (_i = 0, _len = defBufs.length; _i < _len; _i++) {
          defBuf = defBufs[_i];
          _results.push((function(k) {
            return BAMReader.inflateRaw(defBuf, function(e, infBuf) {
              inflatedBuffers[k] = infBuf;
              return readInflatedBuffers();
            });
          })(currentIdx++));
        }
        return _results;
      });
      bambuf4h = new Buffer(0);
      readingHeader = true;
      readInflatedBuffers = function() {
        var bambuf, bamline, bams, e, headerStr, _i, _j, _len, _len1, _ref;
        while (bambuf = inflatedBuffers[readIdx]) {
          delete inflatedBuffers[readIdx];
          readIdx++;
          if (readingHeader) {
            bambuf4h = Buffer.concat([bambuf4h, bambuf]);
            try {
              _ref = BAMReader.readHeaderFromInflatedBuffer(bambuf4h, true), refs = _ref.refs, headerStr = _ref.headerStr, bambuf = _ref.bambuf;
              readingHeader = false;
              if (onHeader) {
                onHeader(headerStr);
              }
              if (bambuf.length === 0) {
                return;
              }
            } catch (_error) {
              e = _error;
              continue;
            }
          }
          bams = BAMReader.readAlignmentsFromInflatedBuffer(bambuf, refs);
          if (onBam) {
            for (_i = 0, _len = bams.length; _i < _len; _i++) {
              bamline = bams[_i];
              onBam(bamline);
            }
          }
          if (onSam) {
            for (_j = 0, _len1 = bams.length; _j < _len1; _j++) {
              bamline = bams[_j];
              onSam(BAMReader.bamToSam(bamline));
            }
          }
        }
        if (onEnd && currentIdx === lastIdx) {
          return onEnd();
        }
      };
      return rstream.on("end", function() {
        return lastIdx = currentIdx;
      });
    };

    BAMReader.readHeader = function(bamfile, cb) {
      var fd, getHeader, infBuf, offset;
      infBuf = new Buffer(0);
      offset = 0;
      fd = fs.openSync(bamfile, "r");
      getHeader = function() {
        var bufToInflate, next, _ref;
        _ref = BAMReader.getDeflatedBuffer(fd, offset), bufToInflate = _ref[0], next = _ref[1];
        offset = next;
        return BAMReader.inflateRaw(bufToInflate, function(e, _infBuf) {
          var headerInfo;
          infBuf = Buffer.concat([infBuf, _infBuf]);
          try {
            headerInfo = BAMReader.readHeaderFromInflatedBuffer(infBuf);
            headerInfo.offset = offset;
            headerInfo.fd = fd;
            return cb(null, headerInfo);
          } catch (_error) {
            e = _error;
            return getHeader();
          }
        });
      };
      return getHeader();
    };

    BAMReader.splitBody = function(bamfile, num, headerInfo, cb) {
      var _splitBody;
      _splitBody = function(e, headerInfo) {
        var b, buf, buflen, cursor, fd, headerCandidate, i, interval, k, match, offset, positions, size, start, _i, _j, _len;
        size = (fs.statSync(bamfile)).size;
        offset = headerInfo.offset;
        fd = headerInfo.fd || fs.openSync(bamfile, "r");
        interval = Math.floor((size - offset) / num);
        positions = [];
        buflen = Math.min(BGZF_ESTIMATED_LEN, interval);
        for (k = _i = 0; 0 <= num ? _i < num : _i > num; k = 0 <= num ? ++_i : --_i) {
          start = interval * k + offset - 1;
          buf = new Buffer(buflen);
          fs.readSync(fd, buf, 0, buflen, start);
          cursor = -1;
          match = false;
          while (!(match || cursor + 16 > buf.length)) {
            cursor++;
            headerCandidate = buf.slice(cursor, cursor + 16);
            match = true;
            for (i = _j = 0, _len = BGZF_HEADER.length; _j < _len; i = ++_j) {
              b = BGZF_HEADER[i];
              if (b !== headerCandidate[i]) {
                match = false;
                break;
              }
            }
          }
          if (match) {
            positions.push(start + cursor);
          }
        }
        fs.closeSync(fd);
        return cb(null, {
          positions: positions,
          header: headerInfo,
          size: size
        });
      };
      if (headerInfo) {
        return _splitBody(null, headerInfo);
      } else {
        return BAMReader.readHeader(bamfile, _splitBody);
      }
    };

    BAMReader.getDeflatedBuffer = function(fd, offset) {
      var bufToInflate, defBuf, delta, i, k, _defBuf, _i;
      defBuf = new Buffer(0);
      k = 0;
      while (true) {
        _defBuf = new Buffer(BGZF_ESTIMATED_LEN);
        fs.readSync(fd, _defBuf, 0, BGZF_ESTIMATED_LEN, offset + k * BGZF_ESTIMATED_LEN);
        for (i = _i = 0; _i < 16; i = ++_i) {
          if (_defBuf[i] !== BGZF_HEADER[i]) {
            throw new Error("not BGZF (offset=" + offset + ", i=" + i + ")");
          }
        }
        defBuf = Buffer.concat([defBuf, _defBuf]);
        delta = defBuf.readUInt16LE(16) + 1;
        if (defBuf.length >= delta) {
          break;
        }
        k++;
      }
      bufToInflate = defBuf.slice(18, delta - 8);
      return [bufToInflate, offset + delta];
    };

    BAMReader.splitDeflatedBuffer = function(defBuf) {
      var cdataLen, defBufs;
      defBufs = [];
      while (true) {
        if (defBuf.length <= 26) {
          return [defBufs, defBuf];
        }
        cdataLen = defBuf.readUInt16LE(16) - 25;
        if (defBuf.length < cdataLen + 26) {
          return [defBufs, defBuf];
        }
        defBufs.push(defBuf.slice(18, cdataLen + 18));
        defBuf = defBuf.slice(26 + cdataLen);
      }
    };

    BAMReader.readHeaderFromInflatedBuffer = function(bambuf, ifReturnsBamBuf) {
      var cursor, headerLen, headerStr, i, nRef, name, nameLen, refLen, refs, ret, _i;
      refs = {};
      headerLen = bambuf.readInt32LE(4);
      if (bambuf.length < headerLen + 16) {
        throw new Error("header len");
      }
      headerStr = bambuf.slice(8, headerLen + 8).toString("ascii");
      cursor = headerLen + 8;
      nRef = bambuf.readInt32LE(cursor);
      cursor += 4;
      for (i = _i = 0; 0 <= nRef ? _i < nRef : _i > nRef; i = 0 <= nRef ? ++_i : --_i) {
        nameLen = bambuf.readInt32LE(cursor);
        cursor += 4;
        name = bambuf.slice(cursor, cursor + nameLen - 1).toString("ascii");
        cursor += nameLen;
        refLen = bambuf.readInt32LE(cursor);
        cursor += 4;
        refs[i] = {
          name: name,
          len: refLen
        };
      }
      ret = {
        refs: refs,
        headerStr: headerStr
      };
      if (ifReturnsBamBuf) {
        ret.bambuf = bambuf.slice(cursor);
      }
      return ret;
    };

    BAMReader.readAlignmentsFromInflatedBuffer = function(buf, refs, readFirst) {
      var arrayLen, bams, bin, blockSize, byte, char, cigar, cigarLen, cursor, flag, flagname, flags, hLen, i, mapq, name, nextPos, nextRefId, num, pos, qual, readName, readNameLen, refId, rname, rnext, second, seq, seqBits, seqLen, seqLenByte, subtype, tLen, tag, tags, valtype, zLen, _i, _j, _k, _len, _len1;
      bams = [];
      while (buf.length) {
        cursor = 0;
        blockSize = buf.readInt32LE(cursor);
        if (buf.length < blockSize) {
          break;
        }
        cursor += 4;
        refId = buf.readInt32LE(cursor);
        rname = refId === -1 ? "*" : refs[refId].name;
        cursor += 4;
        pos = buf.readInt32LE(cursor);
        cursor += 4;
        readNameLen = buf.readUInt8(cursor);
        cursor++;
        mapq = buf.readUInt8(cursor);
        cursor++;
        bin = buf.readUInt16LE(cursor);
        cursor += 2;
        cigarLen = buf.readUInt16LE(cursor);
        cursor += 2;
        flag = buf.readUInt16LE(cursor);
        flags = {};
        for (i = _i = 0, _len = FLAGS.length; _i < _len; i = ++_i) {
          flagname = FLAGS[i];
          flags[flagname] = !!(flag & (0x01 << i));
        }
        cursor += 2;
        seqLen = buf.readInt32LE(cursor);
        cursor += 4;
        nextRefId = buf.readInt32LE(cursor);
        rnext = nextRefId === -1 ? "*" : refs[nextRefId].name;
        cursor += 4;
        nextPos = buf.readInt32LE(cursor);
        cursor += 4;
        tLen = buf.readInt32LE(cursor);
        cursor += 4;
        readName = buf.slice(cursor, cursor + readNameLen - 1).toString("ascii");
        cursor += readNameLen;
        cigar = [];
        for (i = _j = 0; 0 <= cigarLen ? _j < cigarLen : _j > cigarLen; i = 0 <= cigarLen ? ++_j : --_j) {
          num = buf.readUInt32LE(cursor, cursor + 4);
          char = CIGAR_ARR[num & 0x0f];
          num = num >> 4;
          cigar.push(num + char);
          cursor += 4;
        }
        cigar = cigar.join("");
        seqLenByte = Math.floor((seqLen + 1) / 2);
        seqBits = buf.slice(cursor, cursor + seqLenByte);
        seq = [];
        for (_k = 0, _len1 = seqBits.length; _k < _len1; _k++) {
          byte = seqBits[_k];
          seq.push(SEQ_ARR[byte >> 4]);
          second = SEQ_ARR[byte & 0x0F];
          if (second !== "=") {
            seq.push(second);
          }
        }
        seq = seq.join("");
        cursor += seqLenByte;
        qual = ((function() {
          var _l, _results;
          _results = [];
          for (i = _l = 0; 0 <= seqLen ? _l < seqLen : _l > seqLen; i = 0 <= seqLen ? ++_l : --_l) {
            _results.push(String.fromCharCode(buf[cursor + i] + 33));
          }
          return _results;
        })()).join("");
        cursor += seqLen;
        tags = {};
        while (true) {
          if (cursor - 4 >= blockSize) {
            break;
          }
          tag = buf.slice(cursor, cursor + 2).toString("ascii");
          cursor += 2;
          valtype = String.fromCharCode(buf[cursor]);
          cursor++;
          switch (valtype) {
            case "A":
              tags[tag] = {
                type: valtype,
                value: String.fromCharCode(buf[cursor])
              };
              cursor++;
              break;
            case "c":
              tags[tag] = {
                type: "i",
                value: buf.readInt8(cursor)
              };
              cursor++;
              break;
            case "C":
              tags[tag] = {
                type: "i",
                value: buf.readUInt8(cursor)
              };
              cursor++;
              break;
            case "s":
              tags[tag] = {
                type: "i",
                value: buf.readInt16LE(cursor)
              };
              cursor += 2;
              break;
            case "S":
              tags[tag] = {
                type: "i",
                value: buf.readUInt16LE(cursor)
              };
              cursor += 2;
              break;
            case "i":
              tags[tag] = {
                type: "i",
                value: buf.readInt32LE(cursor)
              };
              cursor += 4;
              break;
            case "I":
              tags[tag] = {
                type: "i",
                value: buf.readUInt32LE(cursor)
              };
              cursor += 4;
              break;
            case "f":
              tags[tag] = {
                type: valtype,
                value: buf.readFloatLE(cursor)
              };
              cursor += 4;
              break;
            case "B":
              subtype = String.fromCharCode(buf[cursor]);
              cursor++;
              arrayLen = buf.readInt32LE(cursor);
              cursor += 4;
              switch (subtype) {
                case "c":
                  tags[tag] = {
                    type: valtype,
                    value: (function() {
                      var _l, _results;
                      _results = [];
                      for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                        _results.push(buf.readInt8(cursor + i));
                      }
                      return _results;
                    })()
                  };
                  cursor += arrayLen;
                  break;
                case "C":
                  tags[tag] = {
                    type: valtype,
                    value: (function() {
                      var _l, _results;
                      _results = [];
                      for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                        _results.push(buf.readUInt8(cursor + i));
                      }
                      return _results;
                    })()
                  };
                  cursor += arrayLen;
                  break;
                case "s":
                  tags[tag] = {
                    type: valtype,
                    value: (function() {
                      var _l, _results;
                      _results = [];
                      for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                        _results.push(buf.readInt16LE(cursor + i * 2));
                      }
                      return _results;
                    })()
                  };
                  cursor += arrayLen * 2;
                  break;
                case "S":
                  tags[tag] = {
                    type: valtype,
                    value: (function() {
                      var _l, _results;
                      _results = [];
                      for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                        _results.push(buf.readUInt16LE(cursor + i * 2));
                      }
                      return _results;
                    })()
                  };
                  cursor += arrayLen * 2;
                  break;
                case "i":
                  tags[tag] = {
                    type: valtype,
                    value: (function() {
                      var _l, _results;
                      _results = [];
                      for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                        _results.push(buf.readInt32LE(cursor + i * 4));
                      }
                      return _results;
                    })()
                  };
                  cursor += arrayLen * 4;
                  break;
                case "I":
                  tags[tag] = {
                    type: valtype,
                    value: (function() {
                      var _l, _results;
                      _results = [];
                      for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                        _results.push(buf.readUInt32LE(cursor + i * 4));
                      }
                      return _results;
                    })()
                  };
                  cursor += arrayLen * 4;
                  break;
                case "f":
                  tags[tag] = {
                    type: valtype,
                    value: (function() {
                      var _l, _results;
                      _results = [];
                      for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                        _results.push(buf.readFloatLE(cursor + i * 4));
                      }
                      return _results;
                    })()
                  };
                  cursor += arrayLen * 4;
              }
              value.unshift(subtype);
              break;
            case "Z":
              zLen = 0;
              while (buf[cursor + zLen] !== 0x00) {
                zLen++;
              }
              tags[tag] = {
                type: valtype,
                value: buf.slice(cursor, cursor + zLen).toString("ascii")
              };
              cursor += zLen + 1;
              break;
            case "H":
              hLen = 0;
              while (buf[cursor + hLen] !== 0x00) {
                hLen++;
              }
              tags[tag] = {
                type: valtype,
                value: buf.slice(cursor, cursor + hLen).toString("hex")
              };
              cursor += hLen + 1;
          }
        }
        buf = buf.slice(cursor);
        bams.push({
          qname: readName,
          flag: flag,
          rname: rname,
          pos: pos + 1,
          mapq: mapq,
          cigar: cigar,
          rnext: rnext,
          pnext: nextPos + 1,
          tlen: tLen,
          seq: seq,
          qual: qual,
          tags: tags,
          start: pos,
          flags: flags,
          tagstr: ((function() {
            var _results;
            _results = [];
            for (name in tags) {
              tag = tags[name];
              _results.push([name, tag.type, Array.isArray(tag.value) ? tag.value.join(",") : tag.value].join(":"));
            }
            return _results;
          })()).join("\t")
        });
        if (readFirst) {
          return bams[0];
        }
      }
      return bams;
    };

    BAMReader.bamToSam = function(bamline) {
      return [bamline.qname, bamline.flag, bamline.rname, bamline.pos, bamline.mapq, bamline.cigar || "*", bamline.rnext === bamline.rname && bamline.rname !== "*" ? "=" : bamline.rnext, bamline.pnext, bamline.tlen, bamline.seq, bamline.qual, bamline.tagstr].join("\t");
    };

    BAMReader.inflateRaw = function(defBuf, callback) {
      var engine, error, flow, infBufs, nread;
      engine = createInflateRaw({
        chunkSize: 65535
      });
      nread = 0;
      error = false;
      infBufs = [];
      engine.on("error", function(err) {
        error = true;
        return callback(err);
      });
      engine.on("end", function() {
        var infBuf;
        if (error) {
          return;
        }
        infBuf = Buffer.concat(infBufs, nread);
        infBufs = [];
        engine.close();
        return callback(null, infBuf);
      });
      flow = function() {
        var chunk;
        while (null !== (chunk = engine.read())) {
          infBufs.push(chunk);
          nread += chunk.length;
        }
        return engine.once("readable", flow);
      };
      engine.end(defBuf);
      return flow();
    };

    return BAMReader;

  })();

  module.exports = BAMReader;

}).call(this);

(function() {
  var CIGAR;

  CIGAR = (function() {
    function CIGAR(str) {
      var cigarr, i, i2, _i, _ref;
      this.str = str;
      if (!str || str === "*") {
        return null;
      }
      this.arr = [];
      cigarr = this.str.split(/([A-Z=])/).slice(0, -1);
      for (i = _i = 0, _ref = cigarr.length / 2; 0 <= _ref ? _i < _ref : _i > _ref; i = 0 <= _ref ? ++_i : --_i) {
        i2 = i * 2;
        this.arr.push({
          num: Number(cigarr[i2]),
          type: cigarr[i2 + 1]
        });
      }
    }

    CIGAR.prototype.bpL = function(pos) {
      if (this.arr[0].type.match(/[SH]/)) {
        return pos;
      } else {
        return null;
      }
    };

    CIGAR.prototype.bpR = function(pos) {
      var info, matched, ret, _i, _len, _ref;
      matched = false;
      ret = pos;
      _ref = this.arr;
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        info = _ref[_i];
        if (!matched && info.type === "M") {
          matched = true;
        }
        if (matched) {
          if (info.type.match(/[SH]/)) {
            return ret;
          }
          ret += info.num;
        }
      }
      return null;
    };

    CIGAR.prototype.len = function() {
      var info, ret, _i, _len, _ref;
      ret = 0;
      _ref = this.arr;
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        info = _ref[_i];
        if (info.type.match(/[MS=XI]/)) {
          ret += info.num;
        }
      }
      return ret;
    };

    return CIGAR;

  })();

  module.exports.CIGAR = CIGAR;

}).call(this);
