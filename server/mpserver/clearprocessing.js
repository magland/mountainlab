function clearProcessing(config,req,callback) {

	var args=['clear-processing'];

	var spawn=require('child_process').spawn;
	var exe=config.mountainprocess_exe;
	var process=spawn(exe,args);
	var stdout='';
	var stderr='';
	process.stdout.on('data',function(chunk) {
		stdout+=chunk;
	});
	process.stderr.on('data',function(chunk) {
		stderr+=chunk;
	});
	process.on('close',function(code) {
		callback({success:true,exit_code:code,stdout:stdout,stderr:stderr});
	});

	this.close=function() {
		console.log('clearProcessing task closed.');
		process.stdout.pause();
		process.stderr.pause();
		process.kill();
	}
}

module.exports.clearProcessing = clearProcessing;