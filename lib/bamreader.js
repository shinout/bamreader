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
      this.dic = BAMReader.BAMDic.create(this);
      if (this.dic) {
        this.tlen_mean = this.dic.header.tlen_mean;
        this.tlen_sd = this.dic.header.tlen_sd;
        this.total_reads = this.dic.total_reads;
      } else {
        this.tlen_mean = null;
        this.tlen_sd = null;
        this.total_reads = null;
      }
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
      var bam, bambuf, buf, bytesize, chunk, d_offsets, i, i_offsets, infbuf, len, offset, pitch, read_size, _i, _len, _ref;
      buf = this.infbufs.get(d_offset);
      if (buf) {
        len = buf.length;
        if (i_offset + 4 <= len && (bytesize = buf.readInt32LE(0, true) + 4) <= len) {
          return new module.exports.BAM(buf.slice(i_offset, i_offset + bytesize));
        }
      }
      pitch = 16384;
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
      _ref = ["bam", "start", "end", "finish", "message"];
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

    BAMReader.prototype.iterate = function(o) {
      if (o == null) {
        o = {};
      }
      if (typeof o.num === "number" && o.num >= 2) {
        return this.fork(o);
      } else {
        return this.createIterator(o);
      }
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
  var OutlierFilteredMeanDev;

  OutlierFilteredMeanDev = (function() {
    function OutlierFilteredMeanDev(rate) {
      this.values = {};
      this.rate = typeof rate === "number" && rate >= 0 && rate < 1 ? rate : 0.05;
      this.n = 0;
    }

    OutlierFilteredMeanDev.prototype.add = function(v) {
      if (!this.values[v]) {
        this.values[v] = 0;
      }
      this.values[v]++;
      return this.n++;
    };

    OutlierFilteredMeanDev.prototype.precalc = function() {
      var lower_limit, n, rate, squared, sum, total, upper_limit, v, val, valids, _ref;
      sum = 0;
      rate = this.rate;
      lower_limit = this.n * rate;
      upper_limit = this.n - this.n * rate;
      total = 0;
      valids = 0;
      squared = 0;
      _ref = this.values;
      for (v in _ref) {
        n = _ref[v];
        sum += n;
        if (sum < lower_limit) {
          continue;
        }
        if (sum > upper_limit) {
          break;
        }
        val = v * n;
        total += val;
        valids += n;
        squared += val * v;
      }
      return {
        sum: total,
        squared: squared,
        n: valids
      };
    };

    OutlierFilteredMeanDev.prototype.calc = function() {
      var mean, n, squared, sum, _ref;
      _ref = this.precalc(), sum = _ref.sum, squared = _ref.squared, n = _ref.n;
      mean = sum / n;
      return {
        mean: mean,
        dev: squared / n - mean * mean,
        n: n
      };
    };

    return OutlierFilteredMeanDev;

  })();

  module.exports.OutlierFilteredMeanDev = OutlierFilteredMeanDev;

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
      this.arr = new Array(l_cigar);
      str = "";
      i = 0;
      while (i < l_cigar) {
        num = buf.readUInt32LE(i * 4, true);
        type = CIGAR_ARR[num & 0x0f];
        if (type.match(/[ID]/)) {
          this._indel = true;
        }
        num = num >> 4;
        this.arr[i] = {
          num: num,
          type: type
        };
        str += num + type;
        i++;
      }
      this.string = str;
    }

    CIGAR.createFromString = function(str) {
      var arr, cigar, cigarr, i, i2, l_cigar, type;
      if (!str || str === "*") {
        return nullobj;
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
  var BAM, CIGAR, KNOWN_TAGS, QUAL_ARR, SEQ_ARR, defineGetters;

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

  QUAL_ARR = (function() {
    var i, ret, _i;
    ret = {};
    for (i = _i = 0; _i <= 94; i = ++_i) {
      ret[i] = String.fromCharCode(i + 33);
    }
    return ret;
  })();

  KNOWN_TAGS = (function() {
    var k, ret, tag, _i, _len, _ref;
    ret = {};
    _ref = ["AM", "AS", "BC", "BQ", "CC", "CM", "CO", "CP", "CQ", "CS", "CT", "E2", "FI", "FS", "FZ", "H0", "H1", "H2", "HI", "IH", "LB", "MC", "MD", "MQ", "NH", "NM", "OC", "OP", "OQ", "PG", "PQ", "PT", "PU", "QT", "Q2", "R2", "RG", "RT", "SA", "SM", "TC", "U2", "UQ", "XS"];
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      tag = _ref[_i];
      k = tag.charCodeAt(0) * 256 + tag.charCodeAt(1);
      ret[k] = tag;
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
      d = sam.split("\t");
      bam = {
        reader: reader,
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
        tagstr: d.slice(11).join("\t"),
        length: d[9].length,
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
        var i, len, seq;
        if (this.seq_ != null) {
          return this.seq_;
        }
        len = this.seqbytes.length;
        seq = "";
        i = 0;
        while (i < len) {
          seq += SEQ_ARR[this.seqbytes[i++]];
        }
        return this.seq_ = seq;
      },
      qual: function() {
        var i, len, qual;
        if (this.qual_ != null) {
          return this.qual_;
        }
        len = this.length;
        qual = "";
        i = 0;
        while (i < len) {
          qual += QUAL_ARR[this.qualbytes[i++]];
        }
        return this.qual_ = qual;
      },
      CIGAR: function() {
        if (this.CIGAR_ != null) {
          return this.CIGAR_;
        }
        if (this.cigarbytes) {
          return this.CIGAR_ = new CIGAR(this.cigarbytes, this.l_cigar);
        } else {
          return this.CIGAR_ = CIGAR.createFromString(this.cigar);
        }
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
        return this.sam_ = this.qname + "\t" + this.flag + "\t" + (this.ref_id === -1 ? "*" : this.reader.refs[this.ref_id].name) + "\t" + this.pos + "\t" + this.mapq + "\t" + (this.cigar || "*") + "\t" + (this.nref_id === -1 ? "*" : this.ref_id === this.nref_id ? "=" : this.reader.refs[this.nref_id].name) + "\t" + this.pnext + "\t" + this.tlen + "\t" + this.seq + "\t" + this.qual + "\t" + this.tagstr;
      },
      pair: function() {
        var bam, bams, _i, _len;
        if (!this.reader || !this.multiple) {
          return null;
        }
        if (this.pair_) {
          return this.pair_;
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
              this.pair_ = bam;
              return bam;
            }
          } else if (this.nref_id !== -1 && this.pnext === bam.pos && this.nref_id === bam.ref_id) {
            bam.reader = this.reader;
            this.pair_ = bam;
            return bam;
          }
        }
        this.pair_ = null;
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
        var m, sd, tlen;
        if (!this.reader || this.tlen === 0 || !this.reader.dic) {
          return null;
        }
        m = this.reader.tlen_mean;
        sd = this.reader.tlen_sd;
        tlen = Math.abs(this.tlen);
        return tlen < m - 2 * sd || m + 2 * sd < tlen;
      },
      mean_qual: function() {
        var byte, i, len, total, _i, _len, _ref;
        total = 0;
        len = this.length;
        if (this.qualbytes) {
          _ref = this.qualbytes;
          for (_i = 0, _len = _ref.length; _i < _len; _i++) {
            byte = _ref[_i];
            total += byte;
          }
        } else {
          i = 0;
          while (i < len) {
            total += this.qual.charCodeAt(i) - 33;
          }
        }
        return Math.floor(total / len);
      },
      tagstr: function() {
        var arrayLen, buf, buflen, cursor, hLen, i, subtype, tagbyte, tagname, tagstr, type, valtype, value, zLen;
        if (this.tagstr_) {
          return this.tagstr_;
        }
        tagstr = "";
        cursor = 0;
        buflen = this.tagbytes.length;
        buf = this.tagbytes;
        while (true) {
          if (cursor >= buflen) {
            break;
          }
          tagbyte = buf.readUInt16BE(cursor, true);
          tagname = KNOWN_TAGS[tagbyte] || buf.slice(cursor, cursor + 2).toString("ascii");
          switch (tagname) {
            case "NM":
            case "AS":
            case "XS":
              tagstr += tagname + ":i:" + buf.readUInt8(cursor + 3, true) + "\t";
              cursor += 4;
              continue;
          }
          cursor += 2;
          valtype = String.fromCharCode(buf[cursor]);
          cursor++;
          type = null;
          switch (valtype) {
            case "A":
              value = String.fromCharCode(buf[cursor]);
              cursor++;
              break;
            case "c":
              value = buf.readInt8(cursor, true);
              type = "i";
              cursor++;
              break;
            case "C":
              value = buf.readUInt8(cursor, true);
              type = "i";
              cursor++;
              break;
            case "s":
              value = buf.readInt16LE(cursor, true);
              type = "i";
              cursor += 2;
              break;
            case "S":
              value = buf.readUInt16LE(cursor, true);
              type = "i";
              cursor += 2;
              break;
            case "i":
              value = buf.readInt32LE(cursor, true);
              cursor += 4;
              break;
            case "I":
              value = buf.readUInt32LE(cursor, true);
              type = "i";
              cursor += 4;
              break;
            case "f":
              value = buf.readFloatLE(cursor, true);
              cursor += 4;
              break;
            case "B":
              subtype = String.fromCharCode(buf[cursor]);
              cursor++;
              arrayLen = buf.readInt32LE(cursor, true);
              cursor += 4;
              switch (subtype) {
                case "c":
                  value = (function() {
                    var _i, _results;
                    _results = [];
                    for (i = _i = 0; 0 <= arrayLen ? _i < arrayLen : _i > arrayLen; i = 0 <= arrayLen ? ++_i : --_i) {
                      _results.push(buf.readInt8(cursor + i, true));
                    }
                    return _results;
                  })();
                  cursor += arrayLen;
                  break;
                case "C":
                  value = (function() {
                    var _i, _results;
                    _results = [];
                    for (i = _i = 0; 0 <= arrayLen ? _i < arrayLen : _i > arrayLen; i = 0 <= arrayLen ? ++_i : --_i) {
                      _results.push(buf.readUInt8(cursor + i, true));
                    }
                    return _results;
                  })();
                  cursor += arrayLen;
                  break;
                case "s":
                  value = (function() {
                    var _i, _results;
                    _results = [];
                    for (i = _i = 0; 0 <= arrayLen ? _i < arrayLen : _i > arrayLen; i = 0 <= arrayLen ? ++_i : --_i) {
                      _results.push(buf.readInt16LE(cursor + i * 2, true));
                    }
                    return _results;
                  })();
                  cursor += arrayLen * 2;
                  break;
                case "S":
                  value = (function() {
                    var _i, _results;
                    _results = [];
                    for (i = _i = 0; 0 <= arrayLen ? _i < arrayLen : _i > arrayLen; i = 0 <= arrayLen ? ++_i : --_i) {
                      _results.push(buf.readUInt16LE(cursor + i * 2, true));
                    }
                    return _results;
                  })();
                  cursor += arrayLen * 2;
                  break;
                case "i":
                  value = (function() {
                    var _i, _results;
                    _results = [];
                    for (i = _i = 0; 0 <= arrayLen ? _i < arrayLen : _i > arrayLen; i = 0 <= arrayLen ? ++_i : --_i) {
                      _results.push(buf.readInt32LE(cursor + i * 4, true));
                    }
                    return _results;
                  })();
                  cursor += arrayLen * 4;
                  break;
                case "I":
                  value = (function() {
                    var _i, _results;
                    _results = [];
                    for (i = _i = 0; 0 <= arrayLen ? _i < arrayLen : _i > arrayLen; i = 0 <= arrayLen ? ++_i : --_i) {
                      _results.push(buf.readUInt32LE(cursor + i * 4, true));
                    }
                    return _results;
                  })();
                  cursor += arrayLen * 4;
                  break;
                case "f":
                  value = (function() {
                    var _i, _results;
                    _results = [];
                    for (i = _i = 0; 0 <= arrayLen ? _i < arrayLen : _i > arrayLen; i = 0 <= arrayLen ? ++_i : --_i) {
                      _results.push(buf.readFloatLE(cursor + i * 4, true));
                    }
                    return _results;
                  })();
                  cursor += arrayLen * 4;
              }
              value.unshift(subtype);
              value = value.join(",");
              break;
            case "Z":
              zLen = 0;
              while (buf[cursor + zLen] !== 0x00) {
                zLen++;
              }
              value = buf.slice(cursor, cursor + zLen).toString("ascii");
              cursor += zLen + 1;
              break;
            case "H":
              hLen = 0;
              while (buf[cursor + hLen] !== 0x00) {
                hLen++;
              }
              value = buf.slice(cursor, cursor + hLen).toString("hex");
              cursor += hLen + 1;
          }
          tagstr += tagname + ":" + (type || valtype) + ":" + value + "\t";
        }
        return this.tagstr_ = tagstr.slice(0, -1);
      },
      tags: function() {
        var tag, type, val, value, _i, _len, _ref;
        if (this.tags_) {
          return this.tags_;
        }
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
            default:
              value = val[2];
          }
          this.tags_[tag] = {
            type: type,
            value: value
          };
        }
        return this.tags_;
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

  DEFAULT_PITCH = 16384000;

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
      if (typeof o.end === "function") {
        o.on_end = o.end;
        delete o.end;
      }
      if (typeof o.finish === "function") {
        o.on_finish = o.finish;
      }
      if (typeof o.on_finish === "function") {
        if (!o.on_end) {
          o.on_end = o.on_finish;
        }
      }
      if (typeof o.bam === "function") {
        if (!o.on_bam) {
          o.on_bam = o.bam;
        }
      }
      this.offset = typeof o.start === "number" ? o.start : this.reader.header_offset;
      this.end = typeof o.end === "number" ? o.end : this.reader.size;
      this.pitch = typeof o.pitch === "number" ? o.pitch : DEFAULT_PITCH;
      this.on_bam = typeof o.on_bam === "function" ? o.on_bam : noop;
      this.on_end = typeof o.on_end === "function" ? o.on_end : noop;
      this.on_start = typeof o.on_start === "function" ? o.on_start : noop;
      this.pause = typeof o.pause === "function" ? o.pause : null;
      this.env = o.env || o.$ || {};
      this.paused = false;
      this.ended = false;
      process.nextTick((function(_this) {
        return function() {
          _this.on_start(_this.env);
          return process.nextTick(_this._read.bind(_this));
        };
      })(this));
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
      var bam, bambuf, bytesize, chunk, current_d_offset, current_i_offset, d_offsets, i, i_offset, i_offsets, infbuf, infbuf_len, next_i_offset, offset, read_size, _i, _len, _ref;
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
      infbuf_len = infbuf.length;
      if (infbuf_len === 0) {
        this.ended = read_size === this.end - this.offset;
        this.pitch += this.pitch;
        return this._read();
      }
      for (i = _i = 0, _len = d_offsets.length; _i < _len; i = ++_i) {
        offset = d_offsets[i];
        if (i_offsets[i + 1]) {
          this.reader.infbufs.set(this.offset + offset, infbuf.slice(i_offsets[i], i_offsets[i + 1]));
        }
      }
      i_offset = 0;
      current_i_offset = i_offsets.shift();
      current_d_offset = this.offset + d_offsets.shift();
      while (true) {
        if (i_offset + 4 > infbuf_len) {
          break;
        }
        bytesize = infbuf.readInt32LE(i_offset, true) + 4;
        if (i_offset + bytesize > infbuf_len) {
          break;
        }
        bambuf = infbuf.slice(i_offset, i_offset + bytesize);
        bam = new BAM(bambuf, this.reader);
        bam.i_offset = i_offset - current_i_offset;
        bam.d_offset = current_d_offset;
        this.on_bam(bam, this.env);
        i_offset += bytesize;
        while (true) {
          if (i_offsets[0] === void 0 || i_offset < i_offsets[0]) {
            break;
          }
          next_i_offset = i_offsets.shift();
          current_i_offset = next_i_offset;
          current_d_offset = this.offset + d_offsets.shift();
        }
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
    var $, binarize, finished, merge_sort, merged_num, merging, on_finish, outfile, outlier_rate, tlen_sample_size, tmpfiles;
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
    tlen_sample_size = typeof op.tlen_sample_size === "number" ? op.tlen_sample_size : 100000;
    outlier_rate = typeof op.outlier_rate === "number" ? op.outlier_rate : 0.02;
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
      pool_count: 0,
      outlier_rate: outlier_rate,
      tlen_sample_size: Math.round(tlen_sample_size / op.num)
    };
    this.fork({
      $: $,
      num: op.num,
      pitch: 1024 * 1024 * 4,
      start: function($) {
        return $.tlens = new BAMReader.OutlierFilteredMeanDev($.outlier_rate, $.tlen_sample_size);
      },
      bam: function(bam, $) {
        var binary, data, key, tlen, upper;
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
        if (bam.unmapped === false && bam.next_unmapped === false && bam.same_strand === false && bam.tlen !== 0) {
          tlen = Math.abs(bam.tlen);
          $.tlens.add(tlen);
        }
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
        var n, squared, sum, _ref;
        if ($.debug) {
          console.log([$.n, "E", $.w_count, (new Date / 1000 | 0) - $.time, process.memoryUsage().rss].join("\t"));
        }
        _ref = $.tlens.precalc(), sum = _ref.sum, squared = _ref.squared, n = _ref.n;
        $.tlen_sum = sum;
        $.tlen_squared = squared;
        $.tlen_n = n;
        return delete $.tlens;
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
      var ended, header_buf, idx_header, idx_header_str, l_count, read_write, remainder, tlen_dev, tlen_mean, tlen_n, tlen_sd, tlen_squared, tlen_sum, wstream, _i, _len;
      tlen_sum = 0;
      tlen_squared = 0;
      tlen_n = 0;
      for (_i = 0, _len = $s.length; _i < _len; _i++) {
        $ = $s[_i];
        tlen_sum += $.tlen_sum;
        tlen_squared += $.tlen_squared;
        tlen_n += $.tlen_n;
      }
      tlen_mean = tlen_sum / tlen_n;
      tlen_dev = tlen_squared / tlen_n - tlen_mean * tlen_mean;
      tlen_sd = Math.sqrt(tlen_dev);
      if ($.debug) {
        console.log(["M", "T", Math.round(tlen_mean), (new Date / 1000 | 0) - $.time, "mean"].join("\t"));
      }
      if ($.debug) {
        console.log(["M", "T", Math.round(tlen_sd), (new Date / 1000 | 0) - $.time, "sd"].join("\t"));
      }
      l_count = 0;
      rstream.setEncoding("utf-8");
      wstream = fs.createWriteStream(outfile, {
        highWaterMark: 1024 * 1024 * 10 - 1
      });
      idx_header = {
        tlen_mean: Math.round(tlen_mean),
        tlen_sd: Math.round(tlen_sd),
        tlen_n: tlen_n,
        outlier_rate: $.outlier_rate,
        tlen_sample_size: $.tlen_sample_size
      };
      idx_header_str = JSON.stringify(idx_header);
      header_buf = new Buffer(idx_header_str.length + 4);
      header_buf.writeUInt32BE(idx_header_str.length, 0);
      header_buf.write(idx_header_str, 4);
      wstream.write(header_buf);
      remainder = "";
      ended = false;
      read_write = function() {
        var buf, d, i, line, lines, str, _j, _len1;
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
        for (i = _j = 0, _len1 = lines.length; _j < _len1; i = ++_j) {
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
    if (this.dic === null) {
      throw new Error(".dic file has not been created. reader.createDic() can make the file.");
    }
    return this.dic.fetch(qname, debug);
  };

  BAMDic = (function() {
    BAMDic.create = function(reader) {
      var idxfile;
      idxfile = reader.bamfile + ".dic";
      if (!fs.existsSync(idxfile)) {
        return null;
      }
      return new BAMDic(reader);
    };

    function BAMDic(reader) {
      var headerJSONLen, _b;
      this.reader = reader;
      this.idxfile = this.reader.bamfile + ".dic";
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
      if (typeof o.end === "function") {
        o.on_end = o.end;
        delete o.end;
      }
      if (typeof o.finish === "function") {
        o.on_finish = o.finish;
      }
      if (typeof o.on_finish === "function") {
        if (!o.on_end) {
          o.on_end = o.on_finish;
        }
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
    if (typeof o.num === "number" && o.num >= 2) {
      return this.fork_samtools(o);
    } else {
      return new SAMTools(this, o);
    }
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

(function() {
  var BAMReader, arrayize, parse_query;

  BAMReader = module.exports;

  arrayize = function(v, empty) {
    if (Array.isArray(v)) {
      return v;
    } else if (empty && (v == null)) {
      return [];
    } else {
      return [v];
    }
  };

  parse_query = function(file, setting) {
    var cond_name, condition, conditions, conds, dna, e, file_setting, k, method_name, n_process, on_bam, output, outstream, q, queries, query, reader, samtools, v, v_cond, wstream, _i, _len;
    if (setting == null) {
      setting = {};
    }
    try {
      file_setting = require(file);
    } catch (_error) {
      e = _error;
      console.error("" + file + " : no such file.");
      return;
    }
    if (!file_setting) {
      console.error("" + file + " : invalid format.");
      return;
    }
    setting.__proto__ = file_setting;
    try {
      reader = BAMReader.create(setting.file);
    } catch (_error) {
      e = _error;
      console.error("" + setting.file + " : no such file.");
      return;
    }
    if (typeof setting.query === "function") {
      q = setting.query;
      on_bam = function(bam) {
        if (q(bam)) {
          return output(bam);
        }
      };
    } else if (setting.query != null) {
      queries = arrayize(setting.query);
      conditions = [];
      conds = {};
      for (_i = 0, _len = queries.length; _i < _len; _i++) {
        query = queries[_i];
        condition = {};
        for (k in query) {
          v = query[k];
          if (typeof v === "object") {
            for (cond_name in v) {
              v_cond = v[cond_name];
              if (!conditions[cond_name]) {
                conditions[cond_name] = {};
              }
              conditions[cond_name][k] = v_cond;
            }
          } else {
            if (!condition.equal) {
              condition.equal = {};
            }
            condition.equal[k] = v;
          }
        }
        conditions.push(condition);
      }
      on_bam = function(bam) {
        var _j, _len1, _ref, _ref1, _ref2, _ref3, _ref4, _ref5;
        for (_j = 0, _len1 = conditions.length; _j < _len1; _j++) {
          condition = conditions[_j];
          _ref = condition.equal;
          for (k in _ref) {
            v = _ref[k];
            if (bam[k] !== v) {
              break;
            }
          }
          _ref1 = condition.greater_than;
          for (k in _ref1) {
            v = _ref1[k];
            if (bam[k] <= v) {
              break;
            }
          }
          _ref2 = condition.greater_equal;
          for (k in _ref2) {
            v = _ref2[k];
            if (bam[k] < v) {
              break;
            }
          }
          _ref3 = condition.less_than;
          for (k in _ref3) {
            v = _ref3[k];
            if (bam[k] >= v) {
              break;
            }
          }
          _ref4 = condition.less_equal;
          for (k in _ref4) {
            v = _ref4[k];
            if (bam[k] > v) {
              break;
            }
          }
          _ref5 = condition.values;
          for (k in _ref5) {
            v = _ref5[k];
            if (bam[k] > v) {
              break;
            }
          }
          return output(bam);
        }
      };
    } else {
      console.error("query is required.");
      return;
    }
    outstream = process.stdout;
    if (!setting.output) {
      setting.output = "sam";
    }
    switch (setting.output) {
      case "sam":
        output = function(bam) {
          return outstream.write(bam.sam + "\n");
        };
        break;
      case "bam":
        samtools = require("child_process").spawn("samtools", ["view", "-Sb", "-"]);
        wstream = samtools.stdin;
        samtools.stdout.pipe(outstream);
        output = function(bam) {
          return wstream.write(bam.sam + "\n");
        };
        break;
      case "fastq":
        dna = require("dna");
        output = function(bam) {
          return dna.writeFastq(bam.qname, bam.seq, bam.qual, outstream);
        };
        break;
      case "fastq":
        break;
      default:
        console.log("unknown output type: " + setting.output);
        return;
    }
    method_name = setting["native"] ? "iterate" : "samtools";
    n_process = parseInt(setting.process);
    if (isNaN(n_process || n_process < 1)) {
      n_process = 1;
    }
    return reader[method_name]({
      num: n_process,
      bam: on_bam
    });
  };

}).call(this);
