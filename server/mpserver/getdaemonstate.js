function getDaemonState(config,req,callback) {

	var args=['daemon-state'];

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
		var state={};
		try {
			state=JSON.parse(stdout);
		}
		catch(err) {
			console.log('Error parsing output of daemon-state.');
		}
		callback({success:true,exit_code:code,stdout:stdout,stderr:stderr,state:state});
	});

	this.close=function() {
		console.log('getDaemonState task closed.');
		process.stdout.pause();
		process.stderr.pause();
		process.kill();
	}
}

module.exports.getDaemonState = getDaemonState;