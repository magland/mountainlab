var common={};
module.exports=common;

var fs=require('fs');

common.read_text_file=function(path) {
	return fs.readFileSync(path,'utf8');
};

common.write_text_file=function(path,txt) {
	return fs.writeFileSync(path,txt,'utf8');
};

common.remove_temporary_files=function(paths,callback) {
	common.foreach(paths,{},function(ii,path,cb) {
		common.remove_temporary_file(path,cb);
	},callback);
};

common.remove_temporary_file=function(path,callback) {
	console.log ('removing: '+path);
	fs.unlink(path,callback);
};

common.foreach=function(array,opts,step_function,end_function) {
	var num_parallel=opts.num_parallel||1;
	var num_running=0;
	var num_finished=0;
	var already_called_end=false;
	var ii=0;
	next_step();
	function next_step() {
		if (num_finished>=array.length) {
			setTimeout(function() { //important to do it this way so we don't accumulate a call stack
				if (!already_called_end) { 
					already_called_end=true;
					end_function();
				}
			},0);
			return;
		}
		while ((ii<array.length)&&(num_running<num_parallel)) {
			num_running++;
			ii++;
			step_function(ii-1,array[ii-1],function() {
				num_running--;
				num_finished++;
				setTimeout(function() { //important to do it this way so we don't accumulate a call stack
					next_step();
				},0);
			});
		}
	}
};

var mp_exec_process_timers={};
common.mp_exec_process=function(processor_name,inputs,outputs,params,callback) {
	var exe='mountainprocess';
	var args=['exec-process',processor_name];
	args.push('--_parent_pid='+process.pid);
	var all_params=[];
	for (var key in inputs) {
		all_params[key]=inputs[key];
	}
	for (var key in outputs) {
		all_params[key]=outputs[key];
	}
	for (var key in params) {
		all_params[key]=params[key];
	}
	for (var key in all_params) {
		var val=all_params[key];
		if (typeof(val)!='object')
			args.push('--'+key+'='+val);
		else { // a list of values
			for (var j in val)
				args.push('--'+key+'='+val[j]);
		}
	}
	var timer0=new Date();
	common.make_system_call(exe,args,{show_stdout:false,show_stderr:false,num_tries:2},function(aa) {
		if (aa.return_code!=0) {
			console.error('Subprocess '+processor_name+' returned with a non-zero exit code.');
			process.exit(-1);
		}
		if (!(processor_name in mp_exec_process_timers))
			mp_exec_process_timers[processor_name]=0;
		mp_exec_process_timers[processor_name]+=((new Date())-timer0);
		aa.processor_name=processor_name;
		callback(aa);
	});
};

common.print_mp_exec_process_timers=function() {
	console.log('Total mp process times')
	for (var key in mp_exec_process_timers) {
		console.log(key+': '+mp_exec_process_timers[key]/1000+' sec');
	}
}

common.CLParams=function(argv) {
	this.unnamedParameters=[];
	this.namedParameters={};

	var args=argv.slice(2);
	for (var i=0; i<args.length; i++) {
		var arg0=args[i];
		if (arg0.indexOf('--')===0) {
			arg0=arg0.slice(2);
			var ind=arg0.indexOf('=');
			if (ind>=0) {
				var key0=arg0.slice(0,ind);
				var val0=arg0.slice(ind+1);
				val0=convert_to_number_if_appropriate(val0);
				if (!(key0 in this.namedParameters)) {
					this.namedParameters[key0]=val0;
				}
				else {
					var old_val0=this.namedParameters[key0];
					if (typeof(old_val0)=='object')
						old_val0.push(val0);
					else
						this.namedParameters[key0]=[old_val0,val0];
				}
			}
			else {
				//this.namedParameters[arg0]=args[i+1]||'';
				//i++;
				this.namedParameters[arg0]='';
			}
		}
		//else if (arg0.indexOf('-')===0) {
		//	arg0=arg0.slice(1);
		//	this.namedParameters[arg0]='';
		//}
		else {
			this.unnamedParameters.push(arg0);
		}
	}
	function convert_to_number_if_appropriate(val) {
		if (!isNaN(val))
			return Number(val);
		else
			return val;
	}
};

common.read_mda_header=function(path,callback) {
	var exe='mdaconvert';
	var args=[path,'--readheader'];
	common.make_system_call(exe,args,{num_tries:5},function(aa) {
		var ret=JSON.parse(aa.stdout);
		callback(ret);
	});
};

common.lock_counts={};
common.grab_lock=function(lock_name,callback) {
	var max_lock_counts=1;
	if (!common.lock_counts[lock_name])
		common.lock_counts[lock_name]=0;
	if (common.lock_counts[lock_name]<max_lock_counts) {
		common.lock_counts[lock_name]++;
		callback();
	}
	else {
		setTimeout(function() {
			common.grab_lock(lock_name,callback);
		},100);
	}
};
common.release_lock=function(lock_name) {
	common.lock_counts[lock_name]--;
};

/*
//this does not seem to be working as intended
var _all_spawned=[];
function kill_spawned_processes() {
	console.log ('Killing '+_all_spawned.length+' spawned processes');
	for (var i in _all_spawned) {
		console.log ('Terminating process with pid = '+pp.pid);
		var pp=_all_spawned[i];
		process.kill(pp.pid, 'SIGTERM');
	}
}
process.on('SIGINT', function() {
	kill_spawned_processes();
	process.exit();
});
process.on('SIGTERM', function() {
	kill_spawned_processes();
	process.exit();
});
*/

common.make_system_call=function(cmd,args,opts,callback) {
	var num_tries=opts.num_tries||1;
	var child_process=require('child_process');
	console.log ('Calling: '+cmd+' '+args.join(' '));
	var pp=child_process.spawn(cmd,args);
	_all_spawned=[];
	pp.stdout.setEncoding('utf8');
	pp.stderr.setEncoding('utf8');
	var done=false;
	pp.on('close', function(code) {
		return_it(code);
	});
	//pp.on('exit', function(code) {
	//	return_it();
	//});
	pp.on('error',function(err) {
		console.error('Process error: '+cmd+' '+args.join(' '));
		console.error(err);
	});
	var all_stdout='';
	var all_stderr='';
	pp.stdout.on('data',function(data) {
		if (opts.show_stdout) {
			console.log (data);
		}
		all_stdout+=data;
	});
	pp.stderr.on('data',function(data) {
		if (opts.show_stderr) {
			console.log (data);
		}
		all_stderr+=data;
	});
	function return_it(code) {
		if (done) return;
		if (code!=0) {
			if (num_tries>1) {
				console.error ('Re-trying system call: '+cmd+' '+args.join(' '));
				opts.num_tries=opts.num_tries-1; //todo, I really should not manipulate the opts here. very bad idea
				common.make_system_call(cmd,args,opts,callback);
				return;
			}
			else {
				console.log ('***************\n'+all_stderr+'****************\n');
				console.error('Error in system call: '+cmd+' '+args.join(' '));
				process.exit(-1);
			}
		}
  		done=true;
		if (callback) {
			callback({stdout:all_stdout,stderr:all_stderr,return_code:code});
		}
	}
};

common.copy_file=function(src,dst,callback) {
	var rs=fs.createReadStream(src);
	var ws=fs.createWriteStream(dst);
	ws.on("close", function(ex) {
	    callback();
	});
	rs.pipe(ws);
};

common.copy_files=function(src,dst,callback) {
	common.foreach(src,{},function(ii,fname,cb) {
		console.log ('Copying file: '+src[ii]+' --> '+dst[ii]);
		common.copy_file(src[ii],dst[ii],cb);
	},function() {
		callback();
	});
};

common.clone=function(obj) {
	return JSON.parse(JSON.stringify(obj));
}