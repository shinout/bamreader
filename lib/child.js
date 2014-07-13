(function() {
  var BAMIterator, BAMReader, context;

  BAMReader = require(__dirname + "/bamreader.js");

  BAMIterator = BAMReader.BAMIterator;

  context = {
    LINE_SIZE: 11,
    crypto: require("crypto")
  };

  process.on("message", function(msg) {
    var env, itr, k, reader, _i, _len, _ref;
    BAMReader.parseSendable(msg, context);
    reader = BAMReader.createFromObject(msg.reader);
    if (!msg.env) {
      msg.env = msg.$ || {};
    }
    env = msg.env;
    _ref = ["start", "end", "n"];
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      k = _ref[_i];
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
    itr = BAMIterator.create(reader, msg);
    return itr.on_end = function($) {
      if (typeof msg.on_end === "function") {
        msg.on_end.call(itr, $);
      }
      if (!$._exit_will_be_called) {
        return $.exit();
      }
    };
  });

}).call(this);
