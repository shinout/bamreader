module.exports = (grunt) ->
  grunt.initConfig
    pkg: grunt.file.readJSON "package.json"
    coffee:
      compile:
        files:
          "lib/bamreader.js": "src/bamreader.coffee"

  grunt.loadNpmTasks "grunt-contrib-coffee"
  grunt.registerTask "default", ["coffee"]
