BAMReader = module.exports
arrayize = (v, empty) ->
  if Array.isArray(v) then v else if (empty and not v?) then [] else [v]

parse_query = (file, setting = {})->
  # 1. initializing
  try
    file_setting = require file
  catch e
    console.error "#{file} : no such file."
    return

  if not file_setting
    console.error "#{file} : invalid format."
    return

  setting.__proto__ = file_setting # arg setting is prior to file setting

  try
    reader = BAMReader.create(setting.file)
  catch e
    console.error "#{setting.file} : no such file."
    return

  # 2. preparing "on_bam" function
  if typeof setting.query is "function"
    q = setting.query
    on_bam = (bam)-> output bam if q bam
  else if setting.query?
    queries = arrayize setting.query
    conditions = []
    conds = {}

    for query in queries
      condition = {}
      for k,v of query
        if typeof v is "object"
          for cond_name, v_cond of v
            conditions[cond_name] = {} unless conditions[cond_name]
            conditions[cond_name][k] = v_cond
        else
          condition.equal = {} unless condition.equal
          condition.equal[k] = v
      conditions.push condition

    on_bam = (bam)->
      for condition in conditions
        for k,v of condition.equal
          break if bam[k] isnt v
        for k,v of condition.greater_than
          break if bam[k] <= v
        for k,v of condition.greater_equal
          break if bam[k] <  v
        for k,v of condition.less_than
          break if bam[k] >= v
        for k,v of condition.less_equal
          break if bam[k] >  v
        for k,v of condition.values
          break if bam[k] >  v
        return output bam

  else
    console.error "query is required."
    return

  # 3. output
  outstream = process.stdout
  setting.output = "sam" unless setting.output
  switch setting.output
    when "sam"
      output = (bam)-> outstream.write bam.sam + "\n"
    when "bam"
      samtools = require("child_process").spawn("samtools", ["view", "-Sb", "-"])
      wstream = samtools.stdin
      samtools.stdout.pipe outstream
      output = (bam)->
        wstream.write bam.sam + "\n"
    when "fastq"
      dna = require "dna"
      output = (bam)->
        dna.writeFastq bam.qname, bam.seq, bam.qual, outstream
    when "fastq"
    else
      console.log "unknown output type: #{setting.output}"
      return
  
  # 4. native or samtools
  method_name = if setting.native then "iterate" else "samtools"

  # 5. run
  n_process = parseInt setting.process
  n_process = 1 if isNaN n_process or n_process < 1

  reader[method_name]
    num: n_process
    bam: on_bam
