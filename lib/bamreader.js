(function() {
  var BAMReader, BGZF_ESTIMATED_LEN, BGZF_HEADER, CIGAR_ARR, SEQ_ARR, childProcess, fs, inflateRawSync;

  BGZF_ESTIMATED_LEN = 65536;

  BGZF_HEADER = new Buffer("1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00".split(" ").join(""), "hex");

  SEQ_ARR = "=ACMGRSVTWYHKDBN".split("");

  CIGAR_ARR = "MIDNSHP=X".split("");

  fs = require("fs");

  childProcess = require("child_process");

  inflateRawSync = require("zlib-raw-sync").inflateRawSync;

  BAMReader = (function() {
    function BAMReader(bamfile, options) {
      var e, tlenJSON;
      this.bamfile = bamfile;
      this.options = options != null ? options : {};
      if (this.bamfile.readable) {
        this.bamfile.pause();
      }
      try {
        this.dic = require("bamdic").create(this.bamfile.readable ? this.options.bamfile : this.bamfile);
      } catch (_error) {
        e = _error;
      }
      try {
        tlenJSON = require(this.options.tlenInfo);
        this.tlen_sd = tlenJSON.purified.sd;
        this.tlen_mean = tlenJSON.purified.mean;
      } catch (_error) {
        e = _error;
      }
      if (options.wait) {
        return this;
      }
      childProcess.exec("which samtools", (function(_this) {
        return function(e, stdout, stderr) {
          if (_this.options.samtools) {
            return _this.beginSamtools();
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
          this.onSam = fn;
          break;
        case "bam":
          this.onBam = fn;
          break;
        case "end":
          this.onEnd = fn;
          break;
        case "header":
          this.onHeader = fn;
      }
      return this;
    };

    BAMReader.prototype.beginSamtools = function(isSam) {
      var Bam, file, headerLines, lines, onBam, onEnd, onHeader, onSam, options, reader, readingHeader, samtools, samtoolsCmd;
      reader = this;
      Bam = module.exports.Bam;
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
          return onBam(Bam.createFromSAM(samline, reader));
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
      var bambuf4h, onBam, onEnd, onHeader, onSam, options, reader, readingHeader, refs, remainedBuffer, rstream, _read;
      reader = this;
      onBam = this.onBam;
      onSam = this.onSam;
      onEnd = this.onEnd;
      onHeader = this.onHeader;
      options = this.options;
      if (this.bamfile.readable) {
        rstream = this.bamfile;
        rstream.resume();
      } else {
        rstream = fs.createReadStream(this.bamfile, {
          highWaterMark: 1024 * 1024 - 1
        });
      }
      refs = {};
      bambuf4h = new Buffer(0);
      readingHeader = true;
      remainedBuffer = new Buffer(0);
      _read = function(newBuffer) {
        var bam, bambuf, bams, buf, defBuf, defBufs, e, headerStr, _i, _j, _len, _len1, _ref, _ref1, _results;
        buf = Buffer.concat([remainedBuffer, newBuffer], remainedBuffer.length + newBuffer.length);
        _ref = BAMReader.splitDeflatedBuffer(buf), defBufs = _ref[0], remainedBuffer = _ref[1];
        _results = [];
        for (_i = 0, _len = defBufs.length; _i < _len; _i++) {
          defBuf = defBufs[_i];
          bambuf = inflateRawSync(defBuf);
          if (readingHeader) {
            bambuf4h = Buffer.concat([bambuf4h, bambuf]);
            try {
              _ref1 = BAMReader.readHeaderFromInflatedBuffer(bambuf4h, true), refs = _ref1.refs, headerStr = _ref1.headerStr, bambuf = _ref1.bambuf;
              readingHeader = false;
              if (onHeader) {
                onHeader(headerStr);
              }
              if (bambuf.length === 0) {
                continue;
              }
            } catch (_error) {
              e = _error;
              continue;
            }
          }
          bams = BAMReader.readAlignmentsFromInflatedBuffer(bambuf, refs, false, reader);
          if (onBam) {
            for (_j = 0, _len1 = bams.length; _j < _len1; _j++) {
              bam = bams[_j];
              onBam(bam);
            }
          }
          if (onSam) {
            _results.push((function() {
              var _k, _len2, _results1;
              _results1 = [];
              for (_k = 0, _len2 = bams.length; _k < _len2; _k++) {
                bam = bams[_k];
                _results1.push(onSam(bam.sam));
              }
              return _results1;
            })());
          } else {
            _results.push(void 0);
          }
        }
        return _results;
      };
      rstream.on("data", _read);
      return rstream.on("end", function() {
        _read(remainedBuffer);
        if (onEnd) {
          return onEnd();
        }
      });
    };

    BAMReader.readHeader = function(bamfile) {
      var bufToInflate, e, fd, headerInfo, infBuf, next, offset, _infBuf, _ref;
      infBuf = new Buffer(0);
      offset = 0;
      fd = fs.openSync(bamfile, "r");
      while (true) {
        _ref = BAMReader.getDeflatedBuffer(fd, offset), bufToInflate = _ref[0], next = _ref[1];
        offset = next;
        _infBuf = inflateRawSync(bufToInflate);
        infBuf = Buffer.concat([infBuf, _infBuf]);
        try {
          headerInfo = BAMReader.readHeaderFromInflatedBuffer(infBuf);
          headerInfo.offset = offset;
          headerInfo.fd = fd;
          break;
        } catch (_error) {
          e = _error;
        }
      }
      return headerInfo;
    };

    BAMReader.splitBody = function(bamfile, num, headerInfo) {
      var b, buf, buflen, cursor, fd, headerCandidate, i, interval, k, match, offset, positions, size, start, _i, _j, _len;
      headerInfo = headerInfo || BAMReader.readHeader(bamfile);
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
      return {
        positions: positions,
        size: size,
        header: headerInfo
      };
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
        delta = defBuf.readUInt16LE(16, true) + 1;
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
        cdataLen = defBuf.readUInt16LE(16, true) - 25;
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
      headerLen = bambuf.readInt32LE(4, true);
      if (bambuf.length < headerLen + 16) {
        throw new Error("header len");
      }
      headerStr = bambuf.slice(8, headerLen + 8).toString("ascii");
      cursor = headerLen + 8;
      nRef = bambuf.readInt32LE(cursor, true);
      cursor += 4;
      for (i = _i = 0; 0 <= nRef ? _i < nRef : _i > nRef; i = 0 <= nRef ? ++_i : --_i) {
        nameLen = bambuf.readInt32LE(cursor, true);
        cursor += 4;
        name = bambuf.slice(cursor, cursor + nameLen - 1).toString("ascii");
        cursor += nameLen;
        refLen = bambuf.readInt32LE(cursor, true);
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

    BAMReader.readAlignmentsFromInflatedBuffer = function(buf, refs, readFirst, reader) {
      var Bam, arrayLen, bams, bin, blockSize, byte, char, cigar, cigarLen, cursor, flag, hLen, i, mapq, nextPos, nextRefId, num, pos, qual, readName, readNameLen, refId, rname, rnext, second, seq, seqBits, seqLen, seqLenByte, subtype, tLen, tag, tags, valtype, zLen, _i, _j, _len;
      Bam = module.exports.Bam;
      bams = [];
      while (buf.length) {
        cursor = 0;
        blockSize = buf.readInt32LE(cursor, true);
        if (buf.length < blockSize) {
          break;
        }
        cursor += 4;
        refId = buf.readInt32LE(cursor, true);
        rname = refId === -1 ? "*" : refs[refId].name;
        cursor += 4;
        pos = buf.readInt32LE(cursor, true);
        cursor += 4;
        readNameLen = buf.readUInt8(cursor, true);
        cursor++;
        mapq = buf.readUInt8(cursor, true);
        cursor++;
        bin = buf.readUInt16LE(cursor, true);
        cursor += 2;
        cigarLen = buf.readUInt16LE(cursor, true);
        cursor += 2;
        flag = buf.readUInt16LE(cursor, true);
        cursor += 2;
        seqLen = buf.readInt32LE(cursor, true);
        cursor += 4;
        nextRefId = buf.readInt32LE(cursor, true);
        rnext = nextRefId === -1 ? "*" : refs[nextRefId].name;
        cursor += 4;
        nextPos = buf.readInt32LE(cursor, true);
        cursor += 4;
        tLen = buf.readInt32LE(cursor, true);
        cursor += 4;
        readName = buf.slice(cursor, cursor + readNameLen - 1).toString("ascii");
        cursor += readNameLen;
        cigar = [];
        for (i = _i = 0; 0 <= cigarLen ? _i < cigarLen : _i > cigarLen; i = 0 <= cigarLen ? ++_i : --_i) {
          num = buf.readUInt32LE(cursor, true);
          char = CIGAR_ARR[num & 0x0f];
          num = num >> 4;
          cigar.push(num + char);
          cursor += 4;
        }
        cigar = cigar.join("");
        seqLenByte = Math.floor((seqLen + 1) / 2);
        seqBits = buf.slice(cursor, cursor + seqLenByte);
        seq = [];
        for (_j = 0, _len = seqBits.length; _j < _len; _j++) {
          byte = seqBits[_j];
          seq.push(SEQ_ARR[byte >> 4]);
          second = SEQ_ARR[byte & 0x0F];
          if (second !== "=") {
            seq.push(second);
          }
        }
        seq = seq.join("");
        cursor += seqLenByte;
        qual = ((function() {
          var _k, _results;
          _results = [];
          for (i = _k = 0; 0 <= seqLen ? _k < seqLen : _k > seqLen; i = 0 <= seqLen ? ++_k : --_k) {
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
                value: buf.readInt8(cursor, true)
              };
              cursor++;
              break;
            case "C":
              tags[tag] = {
                type: "i",
                value: buf.readUInt8(cursor, true)
              };
              cursor++;
              break;
            case "s":
              tags[tag] = {
                type: "i",
                value: buf.readInt16LE(cursor, true)
              };
              cursor += 2;
              break;
            case "S":
              tags[tag] = {
                type: "i",
                value: buf.readUInt16LE(cursor, true)
              };
              cursor += 2;
              break;
            case "i":
              tags[tag] = {
                type: "i",
                value: buf.readInt32LE(cursor, true)
              };
              cursor += 4;
              break;
            case "I":
              tags[tag] = {
                type: "i",
                value: buf.readUInt32LE(cursor, true)
              };
              cursor += 4;
              break;
            case "f":
              tags[tag] = {
                type: valtype,
                value: buf.readFloatLE(cursor, true)
              };
              cursor += 4;
              break;
            case "B":
              subtype = String.fromCharCode(buf[cursor]);
              cursor++;
              arrayLen = buf.readInt32LE(cursor, true);
              cursor += 4;
              switch (subtype) {
                case "c":
                  tags[tag] = {
                    type: valtype,
                    value: (function() {
                      var _k, _results;
                      _results = [];
                      for (i = _k = 0; 0 <= arrayLen ? _k < arrayLen : _k > arrayLen; i = 0 <= arrayLen ? ++_k : --_k) {
                        _results.push(buf.readInt8(cursor + i, true));
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
                      var _k, _results;
                      _results = [];
                      for (i = _k = 0; 0 <= arrayLen ? _k < arrayLen : _k > arrayLen; i = 0 <= arrayLen ? ++_k : --_k) {
                        _results.push(buf.readUInt8(cursor + i, true));
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
                      var _k, _results;
                      _results = [];
                      for (i = _k = 0; 0 <= arrayLen ? _k < arrayLen : _k > arrayLen; i = 0 <= arrayLen ? ++_k : --_k) {
                        _results.push(buf.readInt16LE(cursor + i * 2, true));
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
                      var _k, _results;
                      _results = [];
                      for (i = _k = 0; 0 <= arrayLen ? _k < arrayLen : _k > arrayLen; i = 0 <= arrayLen ? ++_k : --_k) {
                        _results.push(buf.readUInt16LE(cursor + i * 2, true));
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
                      var _k, _results;
                      _results = [];
                      for (i = _k = 0; 0 <= arrayLen ? _k < arrayLen : _k > arrayLen; i = 0 <= arrayLen ? ++_k : --_k) {
                        _results.push(buf.readInt32LE(cursor + i * 4, true));
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
                      var _k, _results;
                      _results = [];
                      for (i = _k = 0; 0 <= arrayLen ? _k < arrayLen : _k > arrayLen; i = 0 <= arrayLen ? ++_k : --_k) {
                        _results.push(buf.readUInt32LE(cursor + i * 4, true));
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
                      var _k, _results;
                      _results = [];
                      for (i = _k = 0; 0 <= arrayLen ? _k < arrayLen : _k > arrayLen; i = 0 <= arrayLen ? ++_k : --_k) {
                        _results.push(buf.readFloatLE(cursor + i * 4, true));
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
        bams.push(new Bam(reader || null, readName, flag, rname, pos + 1, mapq, cigar, rnext, nextPos + 1, tLen, seq, qual));
        bams.tags_ = tags;
        if (readFirst) {
          return bams[0];
        }
      }
      return bams;
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

(function() {
  var ARGNAMES, Bam, CIGAR, FLAGS, defineGetters,
    __slice = [].slice;

  defineGetters = function(obj, getters) {
    var fn, name, _results;
    _results = [];
    for (name in getters) {
      fn = getters[name];
      _results.push(Object.defineProperty(obj, name, {
        get: fn
      }));
    }
    return _results;
  };

  CIGAR = module.exports.CIGAR;

  ARGNAMES = ["reader", "qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual"];

  FLAGS = ["multiple", "allmatches", "unmapped", "next_unmapped", "reversed", "next_reversed", "first", "last", "secondary", "lowquality", "duplicate", "supplementary"];

  Bam = (function() {
    function Bam() {
      var args, k, v, _i, _len;
      args = 1 <= arguments.length ? __slice.call(arguments, 0) : [];
      for (k = _i = 0, _len = ARGNAMES.length; _i < _len; k = ++_i) {
        v = ARGNAMES[k];
        this[v] = args[k];
      }
      this._parseFlag();
    }

    Bam.prototype._parseFlag = function() {
      var flagname, i, _i, _len, _results;
      _results = [];
      for (i = _i = 0, _len = FLAGS.length; _i < _len; i = ++_i) {
        flagname = FLAGS[i];
        _results.push(this[flagname] = !!(this.flag & (0x01 << i)));
      }
      return _results;
    };

    defineGetters(Bam.prototype, {
      flags: function() {
        return this;
      },
      tagstr: function() {
        var name, tag;
        if (this.tagstr_) {
          return this.tagstr_;
        }
        return this.tagstr_ = ((function() {
          var _ref, _results;
          _ref = this.tags;
          _results = [];
          for (name in _ref) {
            tag = _ref[name];
            _results.push([name, tag.type, Array.isArray(tag.value) ? tag.value.join(",") : tag.value].join(":"));
          }
          return _results;
        }).call(this)).join("\t");
      },
      tags: function() {
        var subtype, tag, type, v, val, value, _i, _len, _ref;
        if (this.tags_) {
          return this.tags_;
        }
        this.tags_ = {};
        _ref = this.tagstr.split("\t");
        for (_i = 0, _len = _ref.length; _i < _len; _i++) {
          tag = _ref[_i];
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
                  var _j, _len1, _results;
                  _results = [];
                  for (_j = 0, _len1 = value.length; _j < _len1; _j++) {
                    v = value[_j];
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
          this.tags_[tag] = {
            type: type,
            value: value
          };
        }
        return this.tags_;
      },
      next_rname: function() {
        if (this.rnext === "=") {
          return this.rname;
        } else {
          return this.rnext;
        }
      },
      pair: function(bamObj) {
        var bam, bams, next_rname, _i, _len;
        if (!this.reader || !this.reader.dic || !this.multiple) {
          return null;
        }
        bams = this.reader.dic.fetch(this.qname);
        for (_i = 0, _len = bams.length; _i < _len; _i++) {
          bam = bams[_i];
          if (this.secondary || this.supplementary || this.flag === bam.flag) {
            continue;
          }
          if (this.next_unmapped) {
            if (bam.flags.unmapped) {
              bam.reader = this.reader;
              return bam;
            }
          } else {
            next_rname = this.rnext === "=" ? this.rname : this.rnext;
            if (this.pnext === bam.pos && next_rname === bam.rname) {
              bam.reader = this.reader;
              return bam;
            }
          }
        }
        return null;
      },
      start: function() {
        if (this.pos) {
          return this.pos - 1;
        } else {
          return null;
        }
      },
      length: function() {
        return this.seq.length;
      },
      clipped: function() {
        if (this.cigar === "*") {
          return null;
        } else {
          return !!this.cigar.match(/[HS]/);
        }
      },
      soft_clipped: function() {
        if (this.cigar === "*") {
          return null;
        } else {
          return !!this.cigar.match("S");
        }
      },
      hard_clipped: function() {
        if (this.cigar === "*") {
          return null;
        } else {
          return !!this.cigar.match("H");
        }
      },
      match_len: function() {
        if (!this.CIGAR) {
          this.CIGAR = new CIGAR(this.cigar);
        }
        return this.CIGAR.len();
      },
      left_break: function() {
        if (!this.CIGAR) {
          this.CIGAR = new CIGAR(this.cigar);
        }
        return this.CIGAR.bpL(this.pos);
      },
      right_break: function() {
        if (!this.CIGAR) {
          this.CIGAR = new CIGAR(this.cigar);
        }
        return this.CIGAR.bpR(this.pos);
      },
      indel: function() {
        if (this.cigar === "*") {
          return null;
        } else {
          return !!this.cigar.match(/[ID]/);
        }
      },
      fully_matched: function() {
        if (this.cigar === "*") {
          return null;
        } else {
          return this.cigar === this.seq.length + "M";
        }
      },
      different_reference: function() {
        if (this.multiple && !this.unmapped && !this.next_unmapped) {
          return this.rnext !== "=";
        } else {
          return null;
        }
      },
      same_strand: function() {
        if (this.multiple) {
          return this.reversed === this.next_reversed;
        } else {
          return null;
        }
      },
      has_n: function() {
        return !!this.seq.match("N");
      },
      unique: function() {
        return !this.unmapped && this.mapq !== 0;
      },
      mismatch: function() {
        if (!this.tags.NM) {
          return null;
        }
        return this.tags.NM.value;
      },
      discordant: function() {
        var m, sd;
        if (!this.reader || this.tlen === 0) {
          return null;
        }
        m = this.reader.tlen_mean;
        sd = this.reader.tlen_sd;
        return this.tlen < m - 2 * sd || m + 2 * sd < this.tlen;
      },
      sam: function() {
        if (this.sam_) {
          return this.sam_;
        }
        return this.sam_ = [this.qname, this.flag, this.rname, this.pos, this.mapq, this.cigar || "*", this.rnext === this.rname && this.rname !== "*" ? "=" : this.rnext, this.pnext, this.tlen, this.seq, this.qual, this.tagstr].join("\t");
      }
    });

    Bam.createFromSAM = function(samline, bamreader) {
      var bam, sam;
      sam = samline.split("\t");
      bam = new Bam(bamreader, sam[0], Number(sam[1]), sam[2], Number(sam[3]), Number(sam[4]), sam[5], sam[6], Number(sam[7]), Number(sam[8]), sam[9], sam[10]);
      bam.tagstr_ = sam.slice(11).join("\t");
      bam.sam_ = samline;
      return bam;
    };

    return Bam;

  })();

  module.exports.Bam = Bam;

}).call(this);
