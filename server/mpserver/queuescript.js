function queueScript(config,req,callback) {
	var m_script=req.script||'';
	var m_scripts=req.scripts||[];
	var m_id='queueScript-'+req.request_id+'-'+make_random_id(10);

	var args=['queue-script'];
	/// TODO make this path relative to source file
	m_scripts.push({fname:'include.js',code:read_text_file('./include.js')}); 

	for (var i in m_scripts) {
		var script=m_scripts[i];
		if ((!safe_file_name(script.fname))||(!ends_with(script.fname,'.js'))) {
			callback({success:false,error:'Not a valid script file name: '+script.fname});
			return;
		}
		if (!script.code) {
			callback({success:false,error:'Code is empty in: '+script.fname});
			return;	
		}
		
		var fname=config.tmp_mpserver_path+'/'+m_id+'-'+script.fname;
		if (!write_text_file(fname,script.code)) {
			console.error('Unable to write file: '+fname);
			callback({success:false,error:'Unable to write file: '+script.fname}); //don't return the entire path, for privacy/security
			return;	
		}
		args.push(fname);
	}

	if (req.detach) {
		args.push('--~detach=1');	
	}

	
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
		console.log('queueScript task closed.');
		process.stdout.pause();
		process.stderr.pause();
		process.kill();
	}
}

function ends_with(str, suffix) {
    return str.indexOf(suffix, str.length - suffix.length) !== -1;
}

function safe_file_name(fname) {
	return /^[a-zA-Z0-9_\-.]*$/.test(fname);
}

function make_random_id(len)
{
    var text = "";
    var possible = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

    for( var i=0; i < len; i++ )
        text += possible.charAt(Math.floor(Math.random() * possible.length));

    return text;
}

function write_text_file(fname,txt) {
	var fs = require('fs');
	try {
		fs.writeFileSync(fname, txt)
		return true;
	}
	catch(err) {
		console.error('Error in write_text_file: '+fname+' '+err.message);
		return false;
	}
}

function read_text_file(fname) {
	var fs = require('fs');
	try {
		var txt=fs.readFileSync(fname, txt)
		return txt;
	}
	catch(err) {
		console.error('Error in read_text_file: '+fname+' '+err.message);
		return false;
	}
}

module.exports.queueScript = queueScript;