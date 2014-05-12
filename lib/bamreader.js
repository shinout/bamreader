
/*
 * BAMReader by Shin Suzuki(@shinout)
 */

(function() {
  var BAMReader, CIGAR_ARR, FLAGS, SEQ_ARR, childProcess, createInflateRaw, createReadStream;

  createReadStream = require("fs").createReadStream;

  createInflateRaw = require("zlib").createInflateRaw;

  childProcess = require("child_process");

  require("termcolor").define;

  SEQ_ARR = "=ACMGRSVTWYHKDBN".split("");

  CIGAR_ARR = "MIDNSHP=X".split("");

  FLAGS = ["multiple", "allmatches", "unmapped", "next_unmapped", "reversed", "next_reversed", "first", "last", "secondary", "lowquality", "duplicate"];

  BAMReader = (function() {
    function BAMReader(bamfile, options) {
      this.bamfile = bamfile;
      this.options = options != null ? options : {};
      if (this.bamfile.readable) {
        this.bamfile.pause();
      }
      childProcess.exec("which samtools", (function(_this) {
        return function(e, stdout, stderr) {
          if (e || stderr || _this.options["native"]) {
            return _this.begin();
          } else {
            return _this.beginSamtools();
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

    BAMReader.prototype.beginSamtools = function() {
      var file, headerLines, lines, onBam, onEnd, onHeader, onSam, options, readingHeader, samtools;
      onBam = this.onBam;
      onSam = this.onSam;
      onEnd = this.onEnd;
      onHeader = this.onHeader;
      options = this.options;
      file = this.bamfile.readable ? "-" : this.bamfile;
      samtools = childProcess.spawn("samtools", ["view", "-h", file]);
      lines = require("linestream").create(samtools.stdout);
      readingHeader = true;
      headerLines = [];
      if (onEnd) {
        lines.on("end", onEnd);
      }
      lines.on("data", function(samline) {
        var bamline, flagname, i, sam, subtype, tag, type, v, val, value, _i, _j, _len, _len1, _ref;
        if (readingHeader) {
          if (samline.charAt(0) === '@') {
            return headerLines.push(samline);
          } else {
            readingHeader = false;
            if (onHeader) {
              onHeader(headerLines.join("\n"));
            }
            return headerLines = null;
          }
        } else {
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
              tlen: sam[8],
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
        }
      });
      if (this.bamfile.readable) {
        this.bamfile.pipe(samtools.stdin);
        return this.bamfile.resume();
      }
    };

    BAMReader.prototype.begin = function() {
      var count, curXi, inflateRaw, inflatedBuffers, lastXi, onBam, onEnd, onHeader, onSam, options, readAlignment, readHeader, refs, remainedBuffer, rstream, xi;
      onBam = this.onBam;
      onSam = this.onSam;
      onEnd = this.onEnd;
      onHeader = this.onHeader;
      options = this.options;
      count = 0;
      xi = 0;
      refs = {};
      curXi = 1;
      lastXi = 0;
      inflatedBuffers = {};
      inflateRaw = function(i, buffer) {
        var buffers, engine, flow, nread;
        engine = createInflateRaw({
          chunkSize: 65535
        });
        nread = 0;
        buffers = [];
        engine.on("error", function(err) {
          console.error(err);
          console.error("(lastXi: " + lastXi + ")");
          return console.error("----------------------------------");
        });
        engine.on("end", function() {
          var buf, inflatedBuffer;
          buf = Buffer.concat(buffers, nread);
          buffers = [];
          if (i === 0) {
            readHeader(buf);
          } else {
            inflatedBuffers[i] = buf;
          }
          engine.close();
          while (inflatedBuffer = inflatedBuffers[curXi]) {
            readAlignment(inflatedBuffer, i);
            delete inflatedBuffers[curXi];
            curXi++;
          }
          if (onEnd && curXi === lastXi) {
            return onEnd();
          }
        });
        flow = function() {
          var chunk;
          while (null !== (chunk = engine.read())) {
            buffers.push(chunk);
            nread += chunk.length;
          }
          return engine.once("readable", flow);
        };
        engine.end(buffer);
        return flow();
      };
      if (this.bamfile.readable) {
        rstream = this.bamfile;
        rstream.resume();
      } else {
        rstream = createReadStream(this.bamfile, {
          highWaterMark: 1024 * 1024 * 1024 - 1
        });
      }
      remainedBuffer = new Buffer(0);
      rstream.on("data", function(newBuffer) {
        var buf, cdataBuffer, cdataLen, _results;
        buf = Buffer.concat([remainedBuffer, newBuffer], remainedBuffer.length + newBuffer.length);
        _results = [];
        while (true) {
          if (buf.length <= 26) {
            remainedBuffer = buf.length ? buf : new Buffer(0);
            break;
          }
          cdataLen = buf.readUInt16LE(16) - 25;
          if (buf.length < cdataLen + 26) {
            remainedBuffer = buf;
            break;
          }
          cdataBuffer = buf.slice(18, cdataLen + 18);
          inflateRaw(xi++, cdataBuffer);
          _results.push(buf = buf.slice(26 + cdataLen));
        }
        return _results;
      });
      rstream.on("end", function() {
        return lastXi = xi;
      });
      readHeader = function(bambuf) {
        var cursor, header, headerLen, i, nRef, name, nameLen, refLen, _i;
        headerLen = bambuf.readInt32LE(4);
        header = bambuf.slice(8, headerLen + 8).toString("ascii");
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
        if (onHeader) {
          return onHeader(header);
        }
      };
      return readAlignment = function(buf, k) {
        var arrayLen, bamline, bin, blockSize, byte, char, cigar, cigarLen, cursor, flag, flagname, flags, hLen, i, itr, mapq, name, nextPos, nextRefId, num, pos, qual, readName, readNameLen, refId, rname, rnext, samline, second, seq, seqBits, seqLen, seqLenByte, subtype, tLen, tag, tags, valtype, zLen, _i, _j, _k, _len, _len1, _results;
        itr = 0;
        _results = [];
        while (buf.length) {
          cursor = 0;
          blockSize = buf.readInt32LE(cursor);
          if (buf.length < blockSize) {
            break;
          }
          cursor += 4;
          count++;
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
            var _l, _results1;
            _results1 = [];
            for (i = _l = 0; 0 <= seqLen ? _l < seqLen : _l > seqLen; i = 0 <= seqLen ? ++_l : --_l) {
              _results1.push(String.fromCharCode(buf[cursor + i] + 33));
            }
            return _results1;
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
                        var _l, _results1;
                        _results1 = [];
                        for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                          _results1.push(buf.readInt8(cursor + i));
                        }
                        return _results1;
                      })()
                    };
                    cursor += arrayLen;
                    break;
                  case "C":
                    tags[tag] = {
                      type: valtype,
                      value: (function() {
                        var _l, _results1;
                        _results1 = [];
                        for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                          _results1.push(buf.readUInt8(cursor + i));
                        }
                        return _results1;
                      })()
                    };
                    cursor += arrayLen;
                    break;
                  case "s":
                    tags[tag] = {
                      type: valtype,
                      value: (function() {
                        var _l, _results1;
                        _results1 = [];
                        for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                          _results1.push(buf.readInt16LE(cursor + i * 2));
                        }
                        return _results1;
                      })()
                    };
                    cursor += arrayLen * 2;
                    break;
                  case "S":
                    tags[tag] = {
                      type: valtype,
                      value: (function() {
                        var _l, _results1;
                        _results1 = [];
                        for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                          _results1.push(buf.readUInt16LE(cursor + i * 2));
                        }
                        return _results1;
                      })()
                    };
                    cursor += arrayLen * 2;
                    break;
                  case "i":
                    tags[tag] = {
                      type: valtype,
                      value: (function() {
                        var _l, _results1;
                        _results1 = [];
                        for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                          _results1.push(buf.readInt32LE(cursor + i * 4));
                        }
                        return _results1;
                      })()
                    };
                    cursor += arrayLen * 4;
                    break;
                  case "I":
                    tags[tag] = {
                      type: valtype,
                      value: (function() {
                        var _l, _results1;
                        _results1 = [];
                        for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                          _results1.push(buf.readUInt32LE(cursor + i * 4));
                        }
                        return _results1;
                      })()
                    };
                    cursor += arrayLen * 4;
                    break;
                  case "f":
                    tags[tag] = {
                      type: valtype,
                      value: (function() {
                        var _l, _results1;
                        _results1 = [];
                        for (i = _l = 0; 0 <= arrayLen ? _l < arrayLen : _l > arrayLen; i = 0 <= arrayLen ? ++_l : --_l) {
                          _results1.push(buf.readFloatLE(cursor + i * 4));
                        }
                        return _results1;
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
          bamline = {
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
              var _results1;
              _results1 = [];
              for (name in tags) {
                tag = tags[name];
                _results1.push([name, tag.type, Array.isArray(tag.value) ? tag.value.join(",") : tag.value].join(":"));
              }
              return _results1;
            })()).join("\t")
          };
          if (onBam) {
            onBam(bamline);
          }
          if (onSam) {
            samline = [bamline.qname, bamline.flag, bamline.rname, bamline.pos, bamline.mapq, bamline.cigar || "*", bamline.rnext === bamline.rname && bamline.rname !== "*" ? "=" : bamline.rnext, bamline.pnext, bamline.tlen, bamline.seq, bamline.qual, bamline.tagstr].join("\t");
            _results.push(onSam(samline));
          } else {
            _results.push(void 0);
          }
        }
        return _results;
      };
    };

    return BAMReader;

  })();

  module.exports = BAMReader;

}).call(this);
