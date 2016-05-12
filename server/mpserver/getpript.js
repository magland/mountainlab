function getPript(config,req,callback) {

	console.log('getPript');
	console.log(JSON.stringify(req));
	var args=[];
	if (req.prtype=='script') {
		args=['get-script','--id='+req.id];
	}
	else {
		args=['get-process','--id='+req.id];
	}

	var spawn=require('child_process').spawn;
	var exe=config.mountainprocess_exe;
	console.log(JSON.stringify(args));
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
		console.log('getPript task closed.');
		process.stdout.pause();
		process.stderr.pause();
		process.kill();
	}
}

module.exports.getPript = getPript;