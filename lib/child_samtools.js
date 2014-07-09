(function() {
  var BAMReader, SAMTools;

  BAMReader = require(__dirname + "/bamreader");

  SAMTools = BAMReader.SAMTools;

  process.on("message", function(msg) {
    var env, k, options, reader, samtools, _i, _j, _len, _len1, _ref, _ref1;
    _ref = msg._funcs;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      k = _ref[_i];
      msg[k] = eval("(" + msg[k] + ")");
    }
    delete msg._funcs;
    reader = BAMReader.createFromObject(msg.reader);
    env = msg.env || msg.$ || {};
    _ref1 = ["start", "end", "n"];
    for (_j = 0, _len1 = _ref1.length; _j < _len1; _j++) {
      k = _ref1[_j];
      env[k] = msg[k];
    }
    Object.defineProperty(env, "exit", {
      get: function() {
        env._exit_will_be_called = true;
        return function() {
          env.ended = true;
          process.send(env);
          return process.exit();
        };
      }
    });
    options = {
      bam: msg.on_bam,
      start: msg.start,
      end: msg.end,
      env: env
    };
    samtools = new SAMTools(reader, options);
    return samtools.on_end = function($) {
      if (typeof msg.on_end === "function") {
        msg.on_end.call(samtools, $);
      }
      if (!$._exit_will_be_called) {
        return $.exit();
      }
    };
  });

}).call(this);
