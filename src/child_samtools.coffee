BAMReader = require(__dirname + "/bamreader")
SAMTools = BAMReader.SAMTools

process.on "message", (msg)->
  BAMReader.parseSendable msg

  reader = BAMReader.createFromObject msg.reader
  env = msg.env or msg.$ or {}
  env[k] = msg[k] for k in ["start", "end", "n"]
  Object.defineProperty env, "exit",
    get: ->
      env._exit_will_be_called = true
      return ->
        env.ended = true
        process.send(env)
        process.exit()

  options =
    bam    : msg.on_bam
    start  : msg.start
    end    : msg.end
    env    : env

  samtools = new SAMTools(reader, options)
  samtools.on_end = ($)->
    msg.on_end.call(samtools, $) if typeof msg.on_end is "function"
    $.exit() unless $._exit_will_be_called
