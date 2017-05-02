var exports = module.exports = {};
var extend=require('extend');
var common=exports;

var fs=require('fs');
var path=require('path');
var child_process=require('child_process');

var config=JSON.parse(fs.readFileSync(__dirname+'/../mountainlab.default.json','utf8'));
try {
	var config_user_json=fs.readFileSync(__dirname+'/../mountainlab.user.json','utf8');
	var config_user=JSON.parse(config_user_json);
	config=extend(true,config,config_user);
}
catch(err) {
	if (config_user_json) {
		console.log(err);
		console.error('Problem parsing mountainlab.user.json');
		process.exit(-1);
	}
}

exports.read_pipelines_from_text_file=function(file_path) {
	var pipelines=[];
	{
		var txt=common.read_text_file(file_path);
		txt=txt.split('\r').join('\n');
		var lines=txt.split('\n');
		for (var i in lines) {
			if (lines[i].trim().slice(0,1)!='#') {
				if (lines[i].trim()) {
					var vals=lines[i].trim().split(' ');
					var absolute_script_path=find_absolute_pipeline_script_path(vals[1],path.dirname(file_path)||'.');
					if (absolute_script_path) {
						if (vals.length>=2) {
							var pip={
								name:vals[0],
								script:vals[1],
								absolute_script_path:absolute_script_path,
								arguments:vals.slice(2)
							};
							var argv=['',''];
							argv=argv.concat(pip.arguments);
							var CLP=new common.CLParams(argv);
							if ('curation' in CLP.namedParameters) {
								pip.absolute_curation_script_path=find_absolute_pipeline_script_path(CLP.namedParameters.curation,path.dirname(file_path)||'.');
							}
							else {
								pip.absolute_curation_script_path='';
							}
							pipelines.push(pip);
						}
						else {
							if (lines[i].trim()) {
								throw 'problem in pipelines file: '+lines[i].trim();
							}
						}
					}
					else {
						console.log ('Warning: Unable to find pipeline script path: '+vals[1]);
					}
				}
			}
		}
	}
	return pipelines;
};

exports.read_datasets_from_text_file=function(file_path) {
	var datasets=[];
	{
		var txt=common.read_text_file(file_path);
		txt=txt.split('\r').join('\n');
		var lines=txt.split('\n');
		for (var i in lines) {
			if (lines[i].trim().slice(0,1)!='#') {
				var vals=lines[i].trim().split(' ');
				if (vals.length>=2) {
					var absolute_folder_path=find_absolute_dataset_folder_path(vals[1],path.dirname(file_path)||'.');
					if (absolute_folder_path) {
						var dataset_params={};
						var params_fname=absolute_folder_path+'/params.json';
						try {
							dataset_params=JSON.parse(common.read_text_file(params_fname));
						}
						catch(err) {
							console.log ('Error parsing json from file: '+params_fname);
							//console.log (common.read_text_file(params_fname));
						}
						datasets.push({
							name:vals[0],
							folder:vals[1],
							absolute_folder_path:absolute_folder_path,
							dataset_params:dataset_params,
							arguments:vals.slice(2)
						});
					}
					else {
						console.log ('Warning: Unable to find dataset folder: '+vals[1]);
					}
				}
				else {
					if (lines[i].trim()) {
						throw 'problem in datasetlist file: '+lines[i].trim();
					}
				}
			}
		}
	}
	//console.log(datasets);
	return datasets;
};

exports.get_interface=function() {
	//for now this is hard-coded for spike sorting
	return {
		parameters:{
			raw:'$dataset_folder$/raw.mda.prv',
			geom:'$dataset_folder$/geom.csv',
			inpath:'$dataset_folder$'
		}
	};
};

exports.find_pipeline=function(pipelines,pipeline_name) {
	for (var i in pipelines) {
		if (pipelines[i].name==pipeline_name) return pipelines[i];
	}
	return null;
};

exports.find_dataset=function(datasets,dsname) {
	for (var i in datasets) {
		if (datasets[i].name==dsname) return datasets[i];
	}
	return null;
};

function find_absolute_dataset_folder_path(folder,text_file_path) {
	var dataset_paths=config.kron.dataset_paths;
	if (text_file_path)
		dataset_paths.push(text_file_path);
	for (var i in dataset_paths) {
		var p=resolve_from_mountainlab(dataset_paths[i]+'/'+folder);
		if (fs.existsSync(p)) {
			return p;
		}
	}
	return null;
}

function find_absolute_pipeline_script_path(script_path,text_file_path) {
	var pipeline_paths=config.kron.pipeline_paths;
	if (text_file_path)
		pipeline_paths.push(text_file_path);
	for (var i in pipeline_paths) {
		var p=resolve_recursively_from_mountainlab(pipeline_paths[i],script_path);
		if (fs.existsSync(p)) {
			return p;
		}
	}
	return null;
}

exports.find_view_program_file=function(program_name) {
	var view_program_paths=config.kron.view_program_paths;
	for (var i in view_program_paths) {
		var p=resolve_from_mountainlab(view_program_paths[i]+'/'+program_name);
		if (fs.existsSync(p)) {
			return p;
		}
	}
	return null;
};

exports.foreach=function(array,opts,step_function,end_function) {
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

function resolve_recursively_from_mountainlab(path,fname) {
	var tmp=false;
	if (path.indexOf('/')===0) tmp=true; //absolute
	if (path.indexOf('.')===0) tmp=true; //prob referring to the working directory. Witold, help!!
	if (!tmp)
		path=__dirname+'/../'+path;
	return find_recursively(path,fname);
}

function find_recursively(path,fname) {
	if (fs.existsSync(path+'/'+fname))
		return path+'/'+fname;
	if (!fs.existsSync(path))
		return '';
	var files=fs.readdirSync(path);
	for (var i in files) {
		var stat = fs.statSync(path+'/'+files[i]);
        if (stat && stat.isDirectory()) {
        	var tmp=find_recursively(path+'/'+files[i],fname);
        	if (tmp) return tmp;
        }
	}
	return '';
}

function resolve_from_mountainlab(path) {
	if (path.indexOf('/')===0) return path; //absolute
	if (path.indexOf('.')===0) return path; //prob referring to the working directory. Witold, help!!
	return __dirname+'/../'+path; //relative
}

exports.CLParams=function(argv) {
	this.unnamedParameters=[];
	this.namedParameters={};

	var args=argv.slice(2);
	for (var i=0; i<args.length; i++) {
		var arg0=args[i];
		if (arg0.indexOf('--')===0) {
			arg0=arg0.slice(2);
			var ind=arg0.indexOf('=');
			if (ind>=0) {
				this.namedParameters[arg0.slice(0,ind)]=arg0.slice(ind+1);
			}
			else {
				//this.namedParameters[arg0]=args[i+1]||'';
				//i++;
				this.namedParameters[arg0]='';
			}
		}
		else if (arg0.indexOf('-')===0) {
			arg0=arg0.slice(1);
			this.namedParameters[arg0]='';
		}
		else {
			this.unnamedParameters.push(arg0);
		}
	}
};

exports.contains_pipeline=function(pipeline_names,pipeline) {
	if (pipeline_names=='all') return true;
	pipeline_names=pipeline_names.split(',');
	for (var i in pipeline_names) {
		if (pipeline_names[i]==pipeline.name)
			return true;
	}
	return false;
};

exports.contains_ds=function(dsnames,ds) {
	if (dsnames=='all') return true;
	dsnames=dsnames.split(',');
	for (var i in dsnames) {
		if (dsnames[i]==ds.name)
			return true;
	}
	return false;
};


exports.mkdir_safe=function(path) {
	try {
		fs.mkdirSync(path);
	}
	catch (err) {

	}
};

exports.read_text_file=function(path) {
	return fs.readFileSync(path,'utf8');
};

exports.write_text_file=function(path,txt) {
	return fs.writeFileSync(path,txt,'utf8');
};

exports.file_exists=function(path) {
	return fs.existsSync(path);
}

exports.copy_file_sync=function(src,dst) {
	if (!fs.existsSync(src)) return;
	var data=fs.readFileSync(src);
	fs.writeFileSync(dst,data);
};

var s_num_system_calls_running=0;
exports.make_system_call=function(cmd,args,callback) {
	s_num_system_calls_running++;
	var pp=child_process.spawn(cmd,args);
	pp.stdout.setEncoding('utf8');
	pp.stderr.setEncoding('utf8');
	var done=false;
	pp.on('close', function(code) {
		if (done) return;
  		done=true;
		if (callback) callback({stdout:all_stdout,stderr:all_stderr});
		s_num_system_calls_running--;
	});
	pp.on('exit', function(code) {
		if (done) return;
  		done=true;
		if (callback) callback({stdout:all_stdout,stderr:all_stderr});
		s_num_system_calls_running--;
	});
	pp.on('error',function(err) {
		console.log ('Process error: '+cmd+' '+args.join(' '));
		console.log (err);
	});
	var all_stdout='';
	var all_stderr='';
	pp.stdout.on('data',function(data) {
		console.log (data);
		all_stdout+=data;
	});
	pp.stderr.on('data',function(data) {
		console.log (data);
		all_stderr+=data;
	});
};

exports.wait_for_system_calls_to_finish=function(callback) {
	setTimeout(check_it,100);
	function check_it() {
		if (s_num_system_calls_running===0) {
			callback();
		}
		else {
			setTimeout(check_it,100);
		}
	}
};

exports.transpose_matrix=function(X) {
	if (X.length===0) return X;
	var Y=[];
	for (var i in X[0]) {
		Y.push([]);
	}
	for (var j in X) {
		for (var i in X[j]) {
			Y[i].push(X[j][i]);
		}
	}
	return Y;
};

exports.read_csv_matrix=function(path) {
	var ret=[];
	var txt=common.read_text_file(path);
	txt=txt.split('\r').join('\n');
	var lines=txt.split('\n');
	for (var i in lines) {
		var vals=lines[i].split(',');
		if (vals.length>0) {
			var row=[];
			for (var k=0; k<vals.length; k++) {
				row.push(Number(vals[k]));
			}
			ret.push(row);
		}
	}
	return common.transpose_matrix(ret); //this is because of a bad decision I made
};

exports.read_csv_vector=function(path) {
	var X=common.read_csv_matrix(path);
	var Y=[];
	for (var i in X) {
		for (var j in X[i]) {
			Y.push(X[i][j]);
		}
	}
	return Y;
};

exports.print_csv_matrix=function(X) {
	var txt='';
	for (var r=0; r<X.length; r++) {
		console.log (X[r].join(','));
	}
};

exports.row_sum=function(X,row) {
	var ret=0;
	for (var i in X[row]) {
		ret=ret+X[row][i];
	}
	return ret;
};
exports.col_sum=function(X,col) {
	var ret=0;
	for (var i in X) {
		ret=ret+X[i][col];
	}
	return ret;
};

exports.topct=function(num) {
	if (isNaN(num)) return '';
	if (num>1) return '>100%';
	if (num<0.1) return Math.round(num*1000)/10+'%';
	else return Math.round(num*1000)/10+'%';
};

exports.clone=function(X) {
	return JSON.parse(JSON.stringify(X));
};