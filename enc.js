let yue = require('Yuejs')
let data = yue.readfile(__dirname+'/beisaier.js','utf8')
data = yue.encodeData(data)
yue.writefile(__dirname+'/beisaier-ENC.js',data)