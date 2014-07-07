module.exports = (grunt) ->
  grunt.initConfig
    pkg: grunt.file.readJSON "package.json"
    coffee:
      compile:
        files:
          "lib/bamreader.js": ["src/bamreader.coffee", "src/CIGAR.coffee", "src/bam.coffee", "src/bamiterator.coffee", "src/bamdic.coffee"]
          "lib/child.js": ["src/child.coffee"]

  grunt.loadNpmTasks "grunt-contrib-coffee"
  grunt.registerTask "default", ["coffee"]
