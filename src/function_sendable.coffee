# make functions sendable
BAMReader = module.exports
vm = require("vm")

makeSendable = (o)->
  objs = []
  funcs = []
  for k, v of o
    if typeof v is "object" and not Array.isArray v and v not in objs
      objs.push v
      makeSendable v
    else if typeof v is "function"
      o[k] = v.toString()
      funcs.push k
  o._funcs = funcs
 
# parse sendable
parseSendable = (o, scope, context)->
  unless context
    for k, v of global
      scope[k] = v
    scope.BAMReader = BAMReader
    scope.require = require
    scope.fs = require("fs")
    context = vm.createContext(scope)

  if Array.isArray o._funcs
    for k in o._funcs
      #o[k] = eval("(#{o[k]})")
      o[k] = vm.runInContext("(#{o[k]})", context)
    delete o._funcs
  for k, v of o
    parseSendable(v, null, context) if Array.isArray v._funcs

BAMReader.makeSendable  = makeSendable
BAMReader.parseSendable = parseSendable
