path = require('path')
fs   = require('fs')
lib  = path.join(path.dirname(fs.realpathSync(__filename)), '../lib')
BAMReader = require(lib + '/bamreader.js')

main = ->
  try
    ap = require("argparser")
      .nums("p")
      .nonvals("debug", "c", "help")
      .defaults(p: 2)
      .parse()

    return showUsage(true) if ap.opt("help")

    throw message: "bam file is required." if not ap.arg(0)
    bamfile = require("path").resolve(ap.arg(0))
    throw message: "'#{bamfile}': no such file." if not fs.existsSync bamfile

    debug = ap.opt("debug")
    # create index file
    if ap.opt("c")
      num = ap.opt("p")
      reader = BAMReader.create bamfile
      reader.createDic(num: num, debug: debug)
    else
      BAMReader.parse_query(ap.arg(0), ap.opt())

  catch e
    console.error e.message
    showUsage()

showUsage = (out)->
  console[if out then "log" else "error"] """
  [USAGE]
   [query]
  \tbamreader <query file>

   [create dic]
  \tbamreader [memory-size(MB)] -c <bam file> [-p #process] [--debug]
  \tdefault of #process : 2
  \t\t<examples>
  \t\tbamreader -c foo.bam -p 14
  \t\tbamreader 4000 -c bar.bam -p 8 --debug

   [show usage]
  \tbamreader --help
"""

main()
