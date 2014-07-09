module.exports = (grunt) ->
  grunt.initConfig
    pkg: grunt.file.readJSON "package.json"
    meta:
      shebang: '#!/usr/bin/env node'
    coffee:
      compile:
        files:
          "lib/bamreader.js": [
            "src/bamreader.coffee"
            "src/fifo.coffee"
            "src/CIGAR.coffee"
            "src/bam.coffee"
            "src/bamiterator.coffee"
            "src/bamdic.coffee"
            "src/samtools.coffee"
          ]
          "lib/child.js": "src/child.coffee"
          "lib/child_samtools.js": "src/child_samtools.coffee"
          "bin/_bamreader": "src/exe.coffee"
    concat:
      options:
        banner: '<%= meta.shebang %>\n\n'
      dist:
        src:  "bin/_bamreader"
        dest: "bin/_bamreader"

  grunt.loadNpmTasks "grunt-contrib-concat"
  grunt.loadNpmTasks "grunt-contrib-coffee"
  grunt.registerTask "default", ["coffee", "concat"]
  require("fs").chmodSync "bin/_bamreader", "755"
