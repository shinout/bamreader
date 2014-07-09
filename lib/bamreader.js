(function() {
  var BAMReader, INFBUF_CACHE_SIZE, fs, inflateBGZF, isValidBGZF,
    __slice = [].slice;

  fs = require("fs");

  inflateBGZF = require("bgzf").inflate;

  isValidBGZF = require("bgzf").hasValidHeader;

  INFBUF_CACHE_SIZE = 1024 * 1024 * 400;

  BAMReader = (function() {
    function BAMReader(bamfile, o) {
      var _readHeader_result;
      if (o == null) {
        o = {};
      }
      this.bamfile = require("path").resolve(bamfile);
      this.cache_size = o.cache_size || INFBUF_CACHE_SIZE;
      this.infbufs = new module.exports.Fifo(this.cache_size);
      this.fd = fs.openSync(this.bamfile, "r");
      if (o.from_obj) {
        return;
      }
      this.size = fs.statSync(this.bamfile).size;
      _readHeader_result = this._readHeader();
      if (null === _readHeader_result) {
        throw "couldn't read header";
      }
    }

    BAMReader.create = function(bamfile) {
      return new BAMReader(bamfile);
    };

    BAMReader.prototype.createIterator = function(o) {
      if (o == null) {
        o = {};
      }
      if (typeof o === "function") {
        o = {
          on_bam: o
        };
      }
      return module.exports.BAMIterator.create(this, o);
    };

    BAMReader.prototype.on = function(name, fn) {
      if (name === "bam") {
        return this.createIterator(fn);
      }
    };

    BAMReader.prototype.read = function(i_offset, d_offset) {
      var bam, bambuf, buf, bytesize, chunk, d_offsets, i, i_offsets, infbuf, offset, pitch, read_size, _i, _len, _ref;
      buf = this.infbufs.get(d_offset);
      if (buf) {
        buf = buf.slice(i_offset);
        if (buf.length >= 4) {
          bytesize = buf.readInt32LE(0, true) + 4;
          if (buf.length >= bytesize) {
            return new module.exports.BAM(buf.slice(0, bytesize));
          }
        }
      }
      pitch = 512;
      while (true) {
        pitch += pitch;
        read_size = Math.min(this.size - d_offset, pitch);
        chunk = new Buffer(read_size);
        fs.readSync(this.fd, chunk, 0, read_size, d_offset);
        _ref = inflateBGZF(chunk), infbuf = _ref[0], i_offsets = _ref[1], d_offsets = _ref[2];
        infbuf = infbuf.slice(i_offset);
        if (infbuf.length < 4) {
          if (read_size === this.size - d_offset) {
            throw "couldn't fetch bam";
          }
          continue;
        }
        bytesize = infbuf.readInt32LE(0, true) + 4;
        if (infbuf.length < bytesize) {
          if (read_size === this.size - d_offset) {
            throw "couldn't fetch bam";
          }
          continue;
        }
        bambuf = infbuf.slice(0, bytesize);
        bam = new module.exports.BAM(bambuf, this);
        bam.i_offset = i_offset;
        bam.d_offset = d_offset;
        d_offsets.pop();
        for (i = _i = 0, _len = d_offsets.length; _i < _len; i = ++_i) {
          offset = d_offsets[i];
          this.infbufs.set(d_offset + offset, infbuf.slice(i_offsets[i], i_offsets[i + 1]));
        }
        return bam;
      }
    };

    BAMReader.prototype.split = function(num) {
      var buf, buflen, cursor, d_offset, e, infbuf, interval, k, match, nref_id, pitch, positions, ref_id, start;
      pitch = 65535;
      num = typeof num === "number" && num >= 1 ? parseInt(num) : 2;
      interval = Math.floor((this.size - this.header_offset) / num);
      positions = [];
      k = 0;
      while (k < num) {
        start = interval * k + this.header_offset;
        buflen = Math.min(pitch, interval);
        buf = new Buffer(buflen);
        fs.readSync(this.fd, buf, 0, buflen, start);
        cursor = -1;
        match = false;
        while (!(match || cursor + 16 > buf.length)) {
          cursor++;
          if (!isValidBGZF(buf.slice(cursor, cursor + 16))) {
            continue;
          }
          d_offset = start + cursor;
          try {
            infbuf = inflateBGZF(buf.slice(cursor))[0];
          } catch (_error) {
            e = _error;
            infbuf = new Buffer(0);
          }
          if (infbuf.length < 24) {
            if (buflen !== interval) {
              k--;
              pitch += pitch;
            }
            break;
          }
          ref_id = infbuf.readInt32LE(4, true);
          nref_id = infbuf.readInt32LE(24, true);
          if ((ref_id === -1 || (this.refs[ref_id] != null)) && (nref_id === -1 || (this.refs[nref_id] != null))) {
            match = true;
          }
        }
        if (match) {
          positions.push(d_offset);
        }
        k++;
      }
      if (positions.length === 0 && num > 1) {
        return this.split(num - 1);
      }
      return positions;
    };

    BAMReader.prototype.fork = function(o) {
      var child, childs, ended_childs, envs, k, n, num, on_finish, on_message, positions, script, suffix, v, _i, _j, _len, _ref;
      if (o == null) {
        o = {};
      }
      if (typeof o === "function") {
        o = {
          on_bam: o
        };
      }
      num = typeof o.num === "number" && o.num >= 1 ? parseInt(o.num) : 2;
      positions = this.split(num);
      num = positions.length;
      _ref = ["bam", "end", "finish", "message"];
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        suffix = _ref[_i];
        if (typeof o[suffix] === "function") {
          o["on_" + suffix] = o[suffix];
        }
        delete o[suffix];
      }
      positions.push(this.size);
      childs = [];
      on_finish = o.on_finish;
      on_message = o.on_message;
      delete o.on_finish;
      delete o.on_message;
      if (o.script && o.script.match(/^child[a-z_]+/)) {
        script = o.script;
        delete o.script;
      } else {
        script = "child";
      }
      o._funcs = [];
      for (k in o) {
        v = o[k];
        if (typeof v !== "function") {
          continue;
        }
        o[k] = v.toString();
        o._funcs.push(k);
      }
      o.reader = this.toObject();
      ended_childs = 0;
      envs = [];
      for (n = _j = 0; 0 <= num ? _j < num : _j > num; n = 0 <= num ? ++_j : --_j) {
        child = require("child_process").fork("" + __dirname + "/" + script + ".js");
        child.on("message", function(env) {
          if (env.ended) {
            return envs[env.n] = env;
          } else {
            if (typeof on_message === "function") {
              return on_message(env);
            }
          }
        });
        child.on("exit", function() {
          ended_childs++;
          if (ended_childs < num) {
            return;
          }
          if (typeof on_finish === "function") {
            return on_finish.call(this, envs);
          }
        });
        child.options = {
          start: positions[n],
          end: positions[n + 1],
          n: n
        };
        for (k in o) {
          v = o[k];
          child.options[k] = v;
        }
        childs.push(child);
      }
      process.nextTick((function(_this) {
        return function() {
          var _k, _len1, _results;
          _results = [];
          for (_k = 0, _len1 = childs.length; _k < _len1; _k++) {
            child = childs[_k];
            _results.push(child.send(child.options));
          }
          return _results;
        };
      })(this));
      return childs;
    };

    BAMReader.fork = function() {
      var args, file, o;
      args = 1 <= arguments.length ? __slice.call(arguments, 0) : [];
      if (typeof args[0] === "string") {
        if (!args[1]) {
          throw new Error("argument 2 should be an object.");
        }
        file = args[0];
        o = args[1];
      } else {
        o = args[0];
        file = o.file;
        delete o.file;
      }
      return BAMReader.create(file).fork(o);
    };

    BAMReader.prototype._readHeader = function() {
      var blen, buf, current_d_offset, current_i_offset, cursor, d_offsets, e, headerLen, headerStr, i, i_offsets, infbuf, nRef, name, nameLen, read_size, refLen, refs, _i, _ref;
      read_size = 32768;
      while (true) {
        read_size = Math.min(this.size, read_size);
        buf = new Buffer(read_size);
        fs.readSync(this.fd, buf, 0, read_size, 0);
        _ref = inflateBGZF(buf), infbuf = _ref[0], i_offsets = _ref[1], d_offsets = _ref[2];
        try {
          cursor = 0;
          refs = {};
          headerLen = infbuf.readInt32LE(4);
          if (infbuf.length < headerLen + 16) {
            throw new Error("header len");
          }
          headerStr = infbuf.slice(8, headerLen + 8).toString("ascii");
          cursor = headerLen + 8;
          nRef = infbuf.readInt32LE(cursor);
          cursor += 4;
          blen = infbuf.length;
          for (i = _i = 0; 0 <= nRef ? _i < nRef : _i > nRef; i = 0 <= nRef ? ++_i : --_i) {
            nameLen = infbuf.readInt32LE(cursor);
            cursor += 4;
            name = infbuf.slice(cursor, cursor + nameLen - 1).toString("ascii");
            cursor += nameLen;
            refLen = infbuf.readInt32LE(cursor);
            cursor += 4;
            refs[i] = {
              name: name,
              len: refLen
            };
          }
          this.refs = refs;
          this.header = headerStr;
          while (true) {
            current_i_offset = i_offsets.shift();
            current_d_offset = d_offsets.shift();
            if (cursor <= current_i_offset) {
              break;
            }
          }
          this.header_offset = current_d_offset;
          break;
        } catch (_error) {
          e = _error;
          if (read_size === this.size) {
            return null;
          }
          read_size += read_size;
        }
      }
      return true;
    };

    BAMReader.createFromObject = function(obj) {
      var k, reader, _i, _len, _ref;
      reader = new BAMReader(obj.bamfile, {
        from_obj: true,
        cache_size: obj.cache_size
      });
      _ref = ["size", "header", "header_offset", "refs"];
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        k = _ref[_i];
        reader[k] = obj[k];
      }
      return reader;
    };

    BAMReader.prototype.toObject = function() {
      var k, ret, _i, _len, _ref;
      ret = {};
      _ref = ["size", "header", "header_offset", "refs", "bamfile", "cache_size"];
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        k = _ref[_i];
        ret[k] = this[k];
      }
      return ret;
    };

    return BAMReader;

  })();

  module.exports = BAMReader;

}).call(this);

(function() {
  var Fifo;

  Fifo = (function() {
    function Fifo(size_limit) {
      this.size_limit = size_limit != null ? size_limit : INFBUF_CACHE_SIZE;
      this.hash = {};
      this.keys = [];
      this.total_size = 0;
    }

    Fifo.prototype.get = function(k) {
      return this.hash[k];
    };

    Fifo.prototype.set = function(k, v) {
      var key_to_del;
      if (this.hash[k] != null) {
        return;
      }
      while (this.total_size > this.size_limit) {
        key_to_del = this.keys.shift();
        this.total_size -= this.hash[key_to_del].length;
        delete this.hash[key_to_del];
      }
      this.keys.push(k);
      this.hash[k] = v;
      this.total_size += v.length;
    };

    Fifo.prototype.clear = function() {
      while (this.keys.length) {
        delete this.hash[this.keys.shift()];
      }
      this.total_size = 0;
    };

    return Fifo;

  })();

  module.exports.Fifo = Fifo;

}).call(this);

(function() {
  var CIGAR, CIGAR_ARR, nullobj;

  CIGAR_ARR = "MIDNSHP=X".split("");

  nullobj = {
    string: null,
    soft_clipped: function() {
      return null;
    },
    hard_clipped: function() {
      return null;
    },
    bpL: function() {
      return null;
    },
    bpR: function() {
      return null;
    },
    len: function() {
      return null;
    },
    indel: function() {
      return null;
    },
    fully_matched: function() {
      return null;
    }
  };

  CIGAR = (function() {
    function CIGAR(buf, l_cigar) {
      var i, num, str, type;
      if (buf && buf.length === 0) {
        return nullobj;
      }
      this._indel = false;
      if (buf === null) {
        return;
      }
      this.arr = [];
      str = [];
      i = 0;
      while (i < l_cigar) {
        num = buf.readUInt32LE(i * 4);
        type = CIGAR_ARR[num & 0x0f];
        if (type.match(/[ID]/)) {
          this._indel = true;
        }
        num = num >> 4;
        this.arr.push({
          num: num,
          type: type
        });
        str.push(num, type);
        i++;
      }
      this.string = str.join("");
    }

    CIGAR.createFromString = function(str) {
      var arr, cigar, cigarr, i, i2, l_cigar, type;
      if (!str || str === "*") {
        return null;
      }
      cigar = new CIGAR();
      arr = [];
      cigarr = str.split(/([A-Z=])/).slice(0, -1);
      i = 0;
      l_cigar = cigarr.length / 2;
      while (i < l_cigar) {
        i2 = i * 2;
        type = cigarr[i2 + 1];
        if (type.match(/[ID]/)) {
          cigar._indel = true;
        }
        arr.push({
          num: Number(cigarr[i2]),
          type: type
        });
      }
      cigar.arr = arr;
      return cigar.string = str;
    };

    CIGAR.prototype.soft_clipped = function() {
      return this.arr[0].type === "S" || this.arr[this.arr.length - 1].type === "S";
    };

    CIGAR.prototype.hard_clipped = function() {
      return this.arr[0].type === "H" || this.arr[this.arr.length - 1].type === "H";
    };

    CIGAR.prototype.indel = function() {
      return this._indel;
    };

    CIGAR.prototype.fully_matched = function() {
      return this.arr.length === 1 && this.arr.type === "M";
    };

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
  var BAM, CIGAR, SEQ_ARR, defineGetters;

  SEQ_ARR = (function() {
    var arr, c1, c2, i, j, ret, _i, _j, _len, _len1;
    ret = [];
    arr = "=ACMGRSVTWYHKDBN".split("");
    for (i = _i = 0, _len = arr.length; _i < _len; i = ++_i) {
      c1 = arr[i];
      for (j = _j = 0, _len1 = arr.length; _j < _len1; j = ++_j) {
        c2 = arr[j];
        ret.push(c2 === "=" ? c1 : c1 + c2);
      }
    }
    return ret;
  })();

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

  BAM = (function() {
    BAM.createFromSAM = function(sam, reader) {
      var bam, d;
      this.reader = reader;
      d = sam.split("\t");
      bam = {
        qname: d[0],
        flag: Number(d[1]),
        pos: Number(d[3]),
        mapq: Number(d[4]),
        pnext: Number(d[7]),
        tlen: Number(d[8]),
        rname: d[2] === "*" ? null : d[2],
        ref_id: d[2] === "*" ? -1 : d[2],
        rnext: d[6] === "*" ? null : d[6],
        nref_id: d[6] === "*" ? -1 : d[2],
        cigar: d[5] === "*" ? null : d[5],
        seq: d[9],
        qual: d[10],
        tagstr_: d.slice(11).join("\t"),
        sam: sam
      };
      bam.__proto__ = BAM.prototype;
      return bam;
    };

    function BAM(buf, reader) {
      var b_seqlen, cursor, l_cigar, l_qname, l_seq;
      this.reader = reader;
      this.bytesize = buf.readInt32LE(0, true) + 4;
      this.ref_id = buf.readInt32LE(4, true);
      this.pos = buf.readInt32LE(8, true) + 1;
      this.mapq = buf.readUInt8(13, true);
      this.flag = buf.readUInt16LE(18, true);
      this.nref_id = buf.readInt32LE(24, true);
      this.pnext = buf.readInt32LE(28, true) + 1;
      this.tlen = buf.readInt32LE(32, true);
      l_qname = buf.readUInt8(12, true);
      this.qname = buf.slice(36, 36 + l_qname - 1).toString("ascii");
      l_cigar = buf.readUInt16LE(16, true);
      cursor = 36 + l_qname;
      this.cigarbytes = buf.slice(cursor, cursor + l_cigar * 4);
      this.l_cigar = l_cigar;
      l_seq = buf.readInt32LE(20, true);
      cursor += l_cigar * 4;
      b_seqlen = Math.floor((l_seq + 1) / 2);
      this.seqbytes = buf.slice(cursor, cursor + b_seqlen);
      this.length = l_seq;
      cursor += b_seqlen;
      this.qualbytes = buf.slice(cursor, cursor + l_seq);
      cursor += l_seq;
      this.tagbytes = buf.slice(cursor);
    }

    defineGetters(BAM.prototype, {
      multiple: function() {
        return !!(this.flag & 0x01);
      },
      allmatches: function() {
        return !!(this.flag & 0x02);
      },
      unmapped: function() {
        return !!(this.flag & 0x04);
      },
      next_unmapped: function() {
        return !!(this.flag & 0x08);
      },
      reversed: function() {
        return !!(this.flag & 0x10);
      },
      next_reversed: function() {
        return !!(this.flag & 0x20);
      },
      first: function() {
        return !!(this.flag & 0x40);
      },
      last: function() {
        return !!(this.flag & 0x80);
      },
      secondary: function() {
        return !!(this.flag & 0x100);
      },
      lowquality: function() {
        return !!(this.flag & 0x200);
      },
      duplicate: function() {
        return !!(this.flag & 0x400);
      },
      supplementary: function() {
        return !!(this.flag & 0x800);
      },
      rname: function() {
        if (this.ref_id === -1) {
          return null;
        } else {
          return this.reader.refs[this.ref_id].name;
        }
      },
      rnext: function() {
        if (this.nref_id === -1) {
          return null;
        } else {
          return this.reader.refs[this.nref_id].name;
        }
      },
      seq: function() {
        var byte, seq, _i, _len, _ref;
        if (this.seq_ != null) {
          return this.seq_;
        }
        seq = [];
        _ref = this.seqbytes;
        for (_i = 0, _len = _ref.length; _i < _len; _i++) {
          byte = _ref[_i];
          seq.push(SEQ_ARR[byte]);
        }
        return this.seq_ = seq.join("");
      },
      qual: function() {
        var byte, qual, _i, _len, _ref;
        if (this.qual_ != null) {
          return this.qual_;
        }
        qual = [];
        _ref = this.qualbytes;
        for (_i = 0, _len = _ref.length; _i < _len; _i++) {
          byte = _ref[_i];
          qual.push(String.fromCharCode(byte + 33));
        }
        return this.qual_ = qual.join("");
      },
      CIGAR: function() {
        if (this.CIGAR_ != null) {
          return this.CIGAR_;
        }
        return this.CIGAR_ = new CIGAR(this.cigarbytes, this.l_cigar);
      },
      cigar: function() {
        return this.CIGAR.string;
      },
      clipped: function() {
        return this.CIGAR.soft_clipped() || this.CIGAR.hard_clipped();
      },
      soft_clipped: function() {
        return this.CIGAR.soft_clipped();
      },
      hard_clipped: function() {
        return this.CIGAR.hard_clipped();
      },
      match_len: function() {
        return this.CIGAR.len();
      },
      left_break: function() {
        return this.CIGAR.bpL(this.pos);
      },
      right_break: function() {
        return this.CIGAR.bpR(this.pos);
      },
      indel: function() {
        return this.CIGAR.indel();
      },
      fully_matched: function() {
        return this.CIGAR.fully_matched();
      },
      sam: function() {
        if (this.sam_) {
          return this.sam_;
        }
        return this.sam_ = [this.qname, this.flag, this.ref_id === -1 ? "*" : this.reader.refs[this.ref_id].name, this.pos, this.mapq, this.cigar || "*", this.nref_id === -1 ? "*" : this.ref_id === this.nref_id ? "=" : this.reader.refs[this.nref_id].name, this.pnext, this.tlen, this.seq, this.qual, this.tagstr].join("\t");
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
      pair: function(bamObj) {
        var bam, bams, _i, _len;
        if (!this.reader || !this.multiple) {
          return null;
        }
        bams = this.reader.find(this.qname);
        for (_i = 0, _len = bams.length; _i < _len; _i++) {
          bam = bams[_i];
          if (this.secondary || this.supplementary || this.flag === bam.flag) {
            continue;
          }
          if (this.next_unmapped) {
            if (bam.unmapped) {
              bam.reader = this.reader;
              return bam;
            }
          } else if (this.nref_id !== -1 && this.pnext === bam.pos && this.nref_id === bam.ref_id) {
            bam.reader = this.reader;
            return bam;
          }
        }
        return null;
      },
      different_ref: function() {
        if (this.multiple && !this.unmapped && !this.next_unmapped) {
          return this.ref_id !== this.nref_id;
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
        var byte, _i, _len, _ref;
        if (this.has_n_ != null) {
          return this.has_n_;
        }
        this.has_n_ = false;
        _ref = this.seqbytes;
        for (_i = 0, _len = _ref.length; _i < _len; _i++) {
          byte = _ref[_i];
          if (byte >= 0xf0 || (byte & 0x0f) === 0x0f) {
            return this.has_n_ = true;
          }
        }
        return this.has_n_;
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
        if (!this.reader || this.tlen === 0 || !this.reader.tlen_mean || !this.reader.tlen_sd) {
          return null;
        }
        m = this.reader.tlen_mean;
        sd = this.reader.tlen_sd;
        return this.tlen < m - 2 * sd || m + 2 * sd < this.tlen;
      },
      mean_qual: function() {
        var byte, total, _i, _len, _ref;
        total = 0;
        _ref = this.qualbytes;
        for (_i = 0, _len = _ref.length; _i < _len; _i++) {
          byte = _ref[_i];
          total += byte;
        }
        return Math.floor(total / this.qualbytes.length);
      },
      tags: function() {
        var arrayLen, buf, buflen, cursor, hLen, i, subtype, tag, tags, type, v, val, valtype, value, zLen, _i, _len, _ref;
        if (this.tagstr_) {
          _ref = this.tagstr_.split("\t");
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
          if (this.tags_ != null) {
            return this.tags_;
          }
        }
        tags = {};
        cursor = 0;
        buflen = this.tagbytes.length;
        buf = this.tagbytes;
        while (true) {
          if (cursor - 4 >= buflen) {
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
                      var _j, _results;
                      _results = [];
                      for (i = _j = 0; 0 <= arrayLen ? _j < arrayLen : _j > arrayLen; i = 0 <= arrayLen ? ++_j : --_j) {
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
                      var _j, _results;
                      _results = [];
                      for (i = _j = 0; 0 <= arrayLen ? _j < arrayLen : _j > arrayLen; i = 0 <= arrayLen ? ++_j : --_j) {
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
                      var _j, _results;
                      _results = [];
                      for (i = _j = 0; 0 <= arrayLen ? _j < arrayLen : _j > arrayLen; i = 0 <= arrayLen ? ++_j : --_j) {
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
                      var _j, _results;
                      _results = [];
                      for (i = _j = 0; 0 <= arrayLen ? _j < arrayLen : _j > arrayLen; i = 0 <= arrayLen ? ++_j : --_j) {
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
                      var _j, _results;
                      _results = [];
                      for (i = _j = 0; 0 <= arrayLen ? _j < arrayLen : _j > arrayLen; i = 0 <= arrayLen ? ++_j : --_j) {
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
                      var _j, _results;
                      _results = [];
                      for (i = _j = 0; 0 <= arrayLen ? _j < arrayLen : _j > arrayLen; i = 0 <= arrayLen ? ++_j : --_j) {
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
                      var _j, _results;
                      _results = [];
                      for (i = _j = 0; 0 <= arrayLen ? _j < arrayLen : _j > arrayLen; i = 0 <= arrayLen ? ++_j : --_j) {
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
        return this.tags_ = tags;
      }
    });

    return BAM;

  })();

  module.exports.BAM = BAM;

}).call(this);

(function() {
  var BAM, BAMIterator, DEFAULT_PITCH, fs, inflateBGZF, noop;

  fs = require("fs");

  inflateBGZF = require("bgzf").inflate;

  DEFAULT_PITCH = 16384;

  BAM = module.exports.BAM;

  noop = function() {};

  BAMIterator = (function() {
    BAMIterator.create = function(reader, options) {
      if (options == null) {
        options = {};
      }
      return new BAMIterator(reader, options);
    };

    function BAMIterator(reader, o) {
      this.reader = reader;
      if (o == null) {
        o = {};
      }
      this.offset = typeof o.start === "number" ? o.start : this.reader.header_offset;
      this.end = typeof o.end === "number" ? o.end : this.reader.size;
      this.pitch = typeof o.pitch === "number" ? o.pitch : DEFAULT_PITCH;
      this.on_bam = typeof o.on_bam === "function" ? o.on_bam : noop;
      this.on_end = typeof o.on_end === "function" ? o.on_end : noop;
      this.pause = typeof o.pause === "function" ? o.pause : null;
      this.env = o.env || o.$ || {};
      this.paused = false;
      this.ended = false;
      process.nextTick(this._read.bind(this));
    }

    BAMIterator.prototype.on = function(name, fn) {
      switch (name) {
        case "end":
          this.on_end = fn;
          break;
        case "bam":
          this.on_bam = fn;
      }
      return this;
    };

    BAMIterator.prototype.resume = function() {
      if (this.paused) {
        this.paused = false;
        if (this.ended) {
          return process.nextTick(this.on_end.bind(this, this.env));
        } else {
          return process.nextTick(this._read.bind(this));
        }
      }
    };

    BAMIterator.prototype.send = function(msg) {
      if (typeof process.send === "function") {
        return process.send(msg);
      }
    };

    BAMIterator.prototype._read = function() {
      var bam, bambuf, buf, bytesize, chunk, current_d_offset, current_i_offset, d_offsets, i_offset, i_offsets, infbuf, read_size, _ref;
      read_size = Math.min(this.end - this.offset, this.pitch);
      this.ended = this.ended || read_size <= 0;
      if (typeof this.pause === "function") {
        this.paused = this.pause(this.ended, this.env);
      }
      if (this.paused) {
        return;
      }
      if (this.ended) {
        return this.on_end(this.env);
      }
      chunk = new Buffer(read_size);
      fs.readSync(this.reader.fd, chunk, 0, read_size, this.offset);
      _ref = inflateBGZF(chunk), infbuf = _ref[0], i_offsets = _ref[1], d_offsets = _ref[2];
      if (infbuf.length === 0) {
        this.ended = read_size === this.end - this.offset;
        this.pitch += this.pitch;
        return this._read();
      }
      buf = infbuf;
      i_offset = 0;
      current_i_offset = i_offsets.shift();
      current_d_offset = this.offset + d_offsets.shift();
      while (true) {
        if (buf.length < 4) {
          break;
        }
        bytesize = buf.readInt32LE(0, true) + 4;
        if (buf.length < bytesize) {
          break;
        }
        bambuf = buf.slice(0, bytesize);
        bam = new BAM(bambuf, this.reader);
        bam.i_offset = i_offset - current_i_offset;
        bam.d_offset = current_d_offset;
        this.on_bam(bam, this.env);
        i_offset += bytesize;
        while (true) {
          if (i_offsets[0] === void 0 || i_offset < i_offsets[0]) {
            break;
          }
          current_i_offset = i_offsets.shift();
          current_d_offset = this.offset + d_offsets.shift();
        }
        buf = infbuf.slice(i_offset);
      }
      this.offset = current_d_offset;
      return setImmediate(this._read.bind(this));
    };

    return BAMIterator;

  })();

  module.exports.BAMIterator = BAMIterator;

}).call(this);

(function() {
  var BAMDic, BAMReader, ESTIMATE_DELTA, LSIZE, cp, crypto, fs;

  BAMReader = module.exports;

  cp = require("child_process");

  fs = require("fs");

  crypto = require("crypto");

  ESTIMATE_DELTA = 0.0001;

  LSIZE = 12;

  BAMReader.prototype.createDic = function(op, callback) {
    var $, binarize, finished, merge_sort, merged_num, merging, on_finish, outfile, tmpfiles;
    if (op == null) {
      op = {};
    }
    if (typeof op === "number") {
      op = {
        num: op
      };
    }
    outfile = this.bamfile + ".dic";
    tmpfiles = [];
    merged_num = 0;
    merging = 0;
    finished = false;
    $ = {
      WFILE_HWM: 1024 * 1024 * 20 - 1,
      MAX_MEMORY_SIZE: 1.2e9,
      tmpfile_inc: 0,
      outfile: outfile,
      r_count: 0,
      w_count: 0,
      time: new Date / 1000 | 0,
      debug: op.debug,
      pool: {},
      pool_count: 0
    };
    this.fork({
      $: $,
      num: op.num,
      pitch: 1024 * 1024 * 4,
      bam: function(bam, $) {
        var binary, data, key, upper;
        binary = bam.d_offset.toString(2);
        key = crypto.createHash("md5").update(bam.qname).digest().readUInt32BE(0, true);
        data = new Buffer(LSIZE);
        data.writeUInt32BE(key, 0, true);
        data.writeUInt32BE(parseInt(binary.slice(-32), 2), 4, true);
        upper = binary.length > 32 ? parseInt(binary.slice(0, -32), 2) : 0;
        data.writeUInt32BE((upper << 16) + bam.i_offset, 8, true);
        if ($.pool[key] == null) {
          $.pool[key] = [];
          $.pool_count++;
        }
        $.pool[key].push(data);
        return $.r_count++;
      },
      pause: function(ended, $) {
        var WCHUNK_SIZE, memory, tmpfile, w_data, wstream, _write;
        memory = process.memoryUsage();
        if ($.debug) {
          console.log([$.n, "R", $.r_count, (new Date / 1000 | 0) - $.time, memory.rss].join("\t"));
        }
        if (!ended && memory.rss <= $.MAX_MEMORY_SIZE) {
          return false;
        }
        if ($.pool_count === 0 && !ended) {
          setTimeout((function(_this) {
            return function() {
              return _this.resume();
            };
          })(this), 1000);
          return true;
        }
        WCHUNK_SIZE = $.WFILE_HWM - 10000;
        w_data = "";
        tmpfile = "" + $.outfile + "." + $.n + "." + (++$.tmpfile_inc);
        wstream = require("fs").createWriteStream(tmpfile, {
          highWaterMark: $.WFILE_HWM
        });
        _write = function() {
          var arr, key, _ref;
          if ($.debug) {
            console.log([$.n, "W", $.w_count, (new Date / 1000 | 0) - $.time, process.memoryUsage().rss].join("\t"));
          }
          _ref = $.pool;
          for (key in _ref) {
            arr = _ref[key];
            w_data += arr.map(function(data) {
              return data.toString("hex") + "\n";
            }).join("");
            $.w_count += arr.length;
            $.pool_count--;
            delete $.pool[key];
            if (w_data.length > WCHUNK_SIZE) {
              wstream.write(w_data, "utf-8", _write);
              w_data = "";
              return;
            }
          }
          if ($.debug) {
            console.log([$.n, "W", $.w_count, (new Date / 1000 | 0) - $.time, process.memoryUsage().rss].join("\t"));
          }
          return wstream.end(w_data);
        };
        wstream.on("finish", (function(_this) {
          return function() {
            _this.send({
              tmpfile: tmpfile
            });
            return _this.resume();
          };
        })(this));
        _write();
        return true;
      },
      message: function(msg) {
        var files;
        if (!msg.tmpfile) {
          return;
        }
        tmpfiles.push(msg.tmpfile);
        if (tmpfiles.length < 2) {
          return;
        }
        files = tmpfiles.join(" ");
        if ($.debug) {
          console.log(["M", "M", tmpfiles.length, (new Date / 1000 | 0) - $.time, process.memoryUsage().rss].join("\t"));
        }
        merging++;
        merge_sort(files, function() {
          merging--;
          if (merging === 0 && finished) {
            return on_finish(finished);
          }
        });
        return tmpfiles = [];
      },
      end: function($) {
        if ($.debug) {
          return console.log([$.n, "E", $.w_count, (new Date / 1000 | 0) - $.time, process.memoryUsage().rss].join("\t"));
        }
      },
      finish: function($s) {
        finished = $s;
        if (merging === 0) {
          return on_finish(finished);
        }
      }
    });
    merge_sort = function(files, cb) {
      var command, new_name, sort;
      new_name = outfile + ".merged" + (++merged_num);
      command = "sort -m " + files + " > " + new_name;
      return sort = cp.exec(command, function() {
        tmpfiles.push(new_name);
        return cp.exec("rm " + files, cb);
      });
    };
    on_finish = function($s) {
      var sort;
      if (tmpfiles.length >= 2) {
        if ($.debug) {
          console.log(["M", "M", tmpfiles.length, (new Date / 1000 | 0) - $.time, process.memoryUsage().rss].join("\t"));
        }
        sort = cp.spawn("sort", ["-m"].concat(tmpfiles));
        return binarize($s, sort.stdout);
      } else {
        return binarize($s, fs.createReadStream(tmpfiles[0], {
          highWaterMark: 1024 * 1024 * 10 - 1
        }));
      }
    };
    return binarize = function($s, rstream) {
      var ended, header_buf, idx_header, idx_header_str, l_count, read_write, remainder, wstream;
      l_count = 0;
      rstream.setEncoding("utf-8");
      wstream = fs.createWriteStream(outfile, {
        highWaterMark: 1024 * 1024 * 10 - 1
      });
      idx_header = {
        tlen_mean: 0,
        tlen_sd: 0
      };
      idx_header_str = JSON.stringify(idx_header);
      header_buf = new Buffer(idx_header_str.length + 4);
      header_buf.writeUInt32BE(idx_header_str.length, 0);
      header_buf.write(idx_header_str, 4);
      wstream.write(header_buf);
      remainder = "";
      ended = false;
      read_write = function() {
        var buf, d, i, line, lines, str, _i, _len;
        if (ended) {
          return;
        }
        d = rstream.read();
        if (d === null) {
          return rstream.once("readable", read_write);
        }
        str = remainder + d;
        lines = str.split("\n");
        remainder = lines.pop();
        buf = new Buffer(LSIZE * lines.length);
        l_count += lines.length;
        if ($.debug && (l_count % 1e5) < 1000) {
          console.log(["M", "B", l_count, (new Date / 1000 | 0) - $.time, process.memoryUsage().rss].join("\t"));
        }
        for (i = _i = 0, _len = lines.length; _i < _len; i = ++_i) {
          line = lines[i];
          buf.write(line, i * LSIZE, "hex");
        }
        return wstream.write(buf, read_write);
      };
      rstream.once("readable", read_write);
      rstream.on("end", function() {
        ended = true;
        return wstream.end();
      });
      return wstream.on("finish", function() {
        if ($.debug) {
          console.log(["M", "B", l_count, (new Date / 1000 | 0) - $.time, process.memoryUsage().rss].join("\t"));
        }
        return cp.exec("rm " + (tmpfiles.join(" ")), function() {
          if (typeof callback === "function") {
            return callback($s);
          }
        });
      });
    };
  };

  BAMReader.prototype.find = function(qname, debug) {
    if (!this.dic) {
      this.dic = new BAMDic(this);
    }
    if (this.dic === null) {
      throw new Error(".dic file has not been created. reader.createDic() can make the file.");
    }
    return this.dic.fetch(qname, debug);
  };

  BAMDic = (function() {
    function BAMDic(reader) {
      var headerJSONLen, _b;
      this.reader = reader;
      this.idxfile = this.reader.bamfile + ".dic";
      if (!fs.existsSync(this.idxfile)) {
        return null;
      }
      this.size = fs.statSync(this.idxfile).size;
      this.fd = fs.openSync(this.idxfile, "r");
      _b = new Buffer(4);
      fs.readSync(this.fd, _b, 0, 4, 0);
      headerJSONLen = _b.readUInt32BE(0);
      _b = new Buffer(headerJSONLen);
      fs.readSync(this.fd, _b, 0, headerJSONLen, 4);
      this.header = JSON.parse(_b.toString("utf-8"));
      this.header_offset = headerJSONLen + 4;
      this.total_reads = (this.size - this.header_offset) / LSIZE;
      if (this.total_reads !== parseInt(this.total_reads)) {
        throw "" + this.idxfile + " is imcomplete bamdic";
      }
    }

    BAMDic.prototype.fetch = function(qname, debug) {
      var bam, bams, buf, currentLineNum, d_offset, delta, estimates, i_offset, inputMD5, inputMD5Buf, iterationCount, k, lbuf, left, lower, md5, newleft, newright, num, offset, rbuf, results, right, upper, v, _i, _j, _len, _len1, _ref, _ref1, _ref2, _ref3;
      inputMD5Buf = crypto.createHash("md5").update(qname).digest();
      inputMD5 = inputMD5Buf.readUInt32BE(0, true);
      estimates = {
        center: Math.max(1, Math.floor(inputMD5 / 0xffffffff * this.total_reads)),
        delta: Math.floor(this.total_reads * ESTIMATE_DELTA)
      };
      if (debug) {
        _ref = {
          name: qname,
          md5: inputMD5,
          dicfile: this.idxfile
        };
        for (k in _ref) {
          v = _ref[k];
          console.error("" + k + ": " + v);
        }
        for (k in estimates) {
          v = estimates[k];
          console.error("" + k + ": " + v);
        }
      }
      if (delta < 100) {
        left = 1;
        right = this.total_reads + 1;
      } else {
        estimates.left = Math.max(1, estimates.center - estimates.delta);
        estimates.right = Math.min(estimates.center + estimates.delta, this.total_reads);
        lbuf = new Buffer(LSIZE);
        fs.readSync(this.fd, lbuf, 0, 4, LSIZE * (estimates.left - 1) + this.header_offset);
        left = lbuf.readUInt32BE(0, true) < inputMD5 ? estimates.left : 1;
        rbuf = new Buffer(LSIZE);
        fs.readSync(this.fd, rbuf, 0, 4, LSIZE * (estimates.right - 1) + this.header_offset);
        right = inputMD5 < rbuf.readUInt32BE(0, true) ? estimates.right : this.total_reads + 1;
      }
      if (debug) {
        console.error({
          left: left,
          right: right,
          reads: this.total_reads
        });
      }
      md5 = null;
      iterationCount = 0;
      while (inputMD5 !== md5) {
        iterationCount++;
        currentLineNum = Math.floor((left + right) / 2);
        buf = new Buffer(LSIZE);
        offset = LSIZE * (currentLineNum - 1);
        fs.readSync(this.fd, buf, 0, LSIZE, offset + this.header_offset);
        md5 = buf.readUInt32BE(0, true);
        if (debug) {
          _ref1 = {
            input: inputMD5,
            current: md5,
            left: left,
            right: right,
            lineNum: currentLineNum
          };
          for (k in _ref1) {
            v = _ref1[k];
            console.error("\t" + k + ": " + v);
          }
        }
        if (md5 > inputMD5) {
          newright = currentLineNum;
          if (newright === right) {
            break;
          }
          right = newright;
        } else {
          newleft = currentLineNum;
          if (newleft === left) {
            break;
          }
          left = newleft;
        }
      }
      if (debug) {
        console.error("iteration: " + iterationCount);
      }
      if (md5 !== inputMD5) {
        return null;
      }
      results = [buf];
      _ref2 = [1, -1];
      for (_i = 0, _len = _ref2.length; _i < _len; _i++) {
        delta = _ref2[_i];
        num = currentLineNum;
        if (debug) {
          console.error("linenum", num);
        }
        while (true) {
          num += delta;
          if (num < 1 || this.total_reads < num) {
            break;
          }
          if (debug) {
            console.error("num", num);
          }
          buf = new Buffer(LSIZE);
          fs.readSync(this.fd, buf, 0, LSIZE, LSIZE * (num - 1) + this.header_offset);
          md5 = buf.readUInt32BE(0, true);
          if (debug) {
            console.error("md5", md5);
          }
          if (md5 !== inputMD5) {
            break;
          }
          results.push(buf);
        }
      }
      bams = [];
      for (_j = 0, _len1 = results.length; _j < _len1; _j++) {
        buf = results[_j];
        lower = buf.readUInt32BE(4, true);
        upper = buf.readUInt16BE(8, true);
        i_offset = buf.readUInt16BE(10, true);
        d_offset = upper ? upper * 0x100000000 + lower : lower;
        if (debug) {
          _ref3 = {
            d_offset: d_offset,
            i_offset: i_offset
          };
          for (k in _ref3) {
            v = _ref3[k];
            console.error("" + k + ": " + v);
          }
        }
        bam = this.reader.read(i_offset, d_offset);
        if (bam && bam.qname === qname) {
          bams.push(bam);
        }
      }
      return bams;
    };

    return BAMDic;

  })();

  module.exports.BAMDic = BAMDic;

}).call(this);

(function() {
  var BAM, BAMReader, SAMTools, cp, fs;

  BAM = module.exports.BAM;

  cp = require("child_process");

  fs = require("fs");

  SAMTools = (function() {
    function SAMTools(reader, o) {
      this.reader = reader;
      if (typeof o === "function") {
        o = {
          on_bam: o
        };
      }
      if (o.bam) {
        o.on_bam = o.bam;
      }
      this.on_bam = typeof o.on_bam === "function" ? o.on_bam : function() {};
      if (typeof o.start === "number") {
        this.start = o.start;
      }
      if (typeof o.end === "number") {
        this.end = o.end;
      }
      this.on_end = typeof o.on_end === "function" ? o.on_end : function() {};
      this.env = o.env || o.$ || {};
      process.nextTick((function(_this) {
        return function() {
          return _this.view();
        };
      })(this));
    }

    SAMTools.prototype.view = function() {
      var fstream, header_chunk, rstream, samtools, _r;
      if ((this.start != null) && (this.end != null)) {
        samtools = cp.spawn("samtools", ["view", "-"]);
        header_chunk = new Buffer(this.reader.header_offset);
        fs.readSync(this.reader.fd, header_chunk, 0, this.reader.header_offset, 0);
        samtools.stdin.write(header_chunk);
        fstream = fs.createReadStream(this.reader.bamfile, {
          start: this.start,
          end: this.end
        }).pipe(samtools.stdin);
      } else {
        samtools = cp.spawn("samtools", ["view", this.reader.bamfile]);
      }
      rstream = samtools.stdout;
      rstream.setEncoding("utf-8");
      _r = "";
      rstream.on("readable", (function(_this) {
        return function() {
          var bam, chunk, sam, sams, _i, _len, _results;
          chunk = rstream.read();
          if (chunk === null) {
            return;
          }
          sams = (_r + chunk).split("\n");
          _r = sams.pop();
          _results = [];
          for (_i = 0, _len = sams.length; _i < _len; _i++) {
            sam = sams[_i];
            bam = BAM.createFromSAM(sam, _this.reader);
            _results.push(_this.on_bam(bam, _this.env));
          }
          return _results;
        };
      })(this));
      return rstream.on("end", (function(_this) {
        return function() {
          return _this.on_end(_this.env);
        };
      })(this));
    };

    return SAMTools;

  })();

  module.exports.SAMTools = SAMTools;

  BAMReader = module.exports;

  BAMReader.prototype.samtools = function(o) {
    return new SAMTools(this, o);
  };

  BAMReader.prototype.fork_samtools = function(o) {
    if (o == null) {
      o = {};
    }
    if (typeof o === "function") {
      o = {
        on_bam: o
      };
    }
    o.script = "child_samtools";
    return this.fork(o);
  };

}).call(this);
