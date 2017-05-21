console.log ('Running prvbucketserver...');

function print_usage() {
	console.log ('Usage: prvbucketserver [data_directory] --listen_port=[port] --server_url=[http://localhost:8000] [--listen_hostname=[host]]');
}

//// requires
var	url=require('url');
var http=require('http');
var fs=require('fs');

CLP=new CLParams(process.argv);
console.log(CLP.namedParameters);
// The listen port comes from the command-line option
var listen_port=CLP.namedParameters['listen_port']||0;
var listen_hostname=CLP.namedParameters['listen_hostname']||'';
var server_url=CLP.namedParameters['server_url']||'';
var prv_exe=__dirname+'/../../bin/prv'
var url_path='/prvbucket';

// The first argument is the data directory -- the base path from which files will be served
var data_directory=CLP.unnamedParameters[0]||'';
if ((!data_directory)||(!listen_port)||(!('server_url' in CLP.namedParameters))) {
	print_usage();
	process.exit(-1);
}

// Exit this program when the user presses ctrl+C
process.on('SIGINT', function() {
    process.exit();
});

// Create the web server!
var SERVER=http.createServer(function (REQ, RESP) {

	//parse the url of the request
	var url_parts = url.parse(REQ.url,true);
	var host=url_parts.host;
	var path=url_parts.pathname;
	var query=url_parts.query;
	console.log(JSON.stringify(url_parts));
	console.log('host='+host);
	
	//by default, the action will be 'download'
	var action=query.a||'download';	
	var info=query.info||'{}';

	if (REQ.method == 'OPTIONS') {
		var headers = {};
		
		//allow cross-domain requests
		/// TODO refine this
		
		headers["Access-Control-Allow-Origin"] = "*";
		headers["Access-Control-Allow-Methods"] = "POST, GET, PUT, DELETE, OPTIONS";
		headers["Access-Control-Allow-Credentials"] = false;
		headers["Access-Control-Max-Age"] = '86400'; // 24 hours
		headers["Access-Control-Allow-Headers"] = "X-Requested-With, X-HTTP-Method-Override, Content-Type, Accept";
		RESP.writeHead(200, headers);
		RESP.end();
	}
	else if (REQ.method=='GET') {
		console.log ('GET: '+REQ.url);

		if (url_path) {
			// Next, if there is a url path, we check to see if that path matches
			if (path.indexOf(url_path)!==0) {
				send_json_response({success:false,error:'Unexpected path: '+path});
				return;
			}
			path=path.slice(url_path.length); // now let's take everything after that path
		}

		if (action=="download") {
			var fname=absolute_data_directory()+"/"+path;
			console.log ('Download: '+fname);
			if (!require('fs').existsSync(fname)) {
				console.log ('Files does not exist: '+fname);
				send_json_response({success:false,error:"File does not exist: "+path});		
				return;	
			}
			var opts={};
			if (query.bytes) {
				//request a subset of bytes
				var vals=query.bytes.split('-');
				if (vals.length!=2) {
					send_json_response({success:false,error:"Error in bytes parameter: "+query.bytes});		
					return;			
				}
				opts.start_byte=Number(vals[0]);
				opts.end_byte=Number(vals[1]);
			}
			else {
				//otherwise get the entire file
				opts.start_byte=0;
				opts.end_byte=get_file_size(fname)-1;
			}
			//serve the file
			serve_file(REQ,fname,RESP,opts);
		}
		else if (action=="stat") {
			var fname=absolute_data_directory()+"/"+path;
			if (!require('fs').existsSync(fname)) {
				send_json_response({success:false,error:"File does not exist: "+path});		
				return;	
			}
			run_process_and_read_stdout(prv_exe,['stat',fname],function(txt) {
				try {
					var obj=JSON.parse(txt);
					send_json_response(obj);
				}
				catch(err) {
					send_json_response({success:false,error:'Problem parsing json response from prv stat: '+txt});
				}
			});
		}
		else if (action=="locate") {
			if ((!query.checksum)||(!query.size)||(!('fcs' in query))) {
				send_json_response({success:false,error:"Invalid query."});	
				return;
			}
			run_process_and_read_stdout(prv_exe,['locate','--search_path='+absolute_data_directory(),'--checksum='+query.checksum,'--size='+query.size,'--fcs='+(query.fcs||'')],function(txt) {
				txt=txt.trim();
				if (!starts_with(txt,absolute_data_directory()+'/')) {
					send_as_text_or_link('<not found>');
					return;
				}
				if (txt) txt=txt.slice(absolute_data_directory().length+1);
				if (server_url) txt=server_url+url_path+'/'+txt;
				send_as_text_or_link(txt);
			});	
		}
		else {
			send_json_response({success:false,error:"Unrecognized action: "+action});
		}
	}
	else if(REQ.method=='POST') {
		send_text_response("Unsuported request method: POST");
	}
	else {
		send_text_response("Unsuported request method.");
	}

	function absolute_data_directory() {
		var ret=data_directory;
		if (ret.indexOf('/')===0) return ret;
		return __dirname+'/../'+ret;
	}
	
	function send_json_response(obj) {
		console.log ('RESPONSE: '+JSON.stringify(obj));
		RESP.writeHead(200, {"Access-Control-Allow-Origin":"*", "Content-Type":"application/json"});
		RESP.end(JSON.stringify(obj));
	}
	function send_text_response(text) {
		console.log ('RESPONSE: '+text);
		RESP.writeHead(200, {"Access-Control-Allow-Origin":"*", "Content-Type":"text/plain"});
		RESP.end(text);
	}
	function send_html_response(text) {
		console.log ('RESPONSE: '+text);
		RESP.writeHead(200, {"Access-Control-Allow-Origin":"*", "Content-Type":"text/html"});
		RESP.end(text);
	}
	function send_as_text_or_link(text) {
		if ((query.return_link=='true')&&(looks_like_it_could_be_a_url(text))) {
			if (text) {
				send_html_response('<html><body><a href="'+text+'">'+text+'</a></body></html>');
			}
			else {
				send_html_response('<html><body>Not found.</html>');	
			}
		}
		else {
			send_text_response(text);
		}
	}
});
if (listen_hostname)
	SERVER.listen(listen_port,listen_hostname);
else
	SERVER.listen(listen_port);
SERVER.timeout=1000*60*60*24; //give it 24 hours!
console.log ('Listening on port: '+listen_port+', host: '+listen_hostname);

function run_process_and_read_stdout(exe,args,callback) {
	console.log ('RUNNING:'+exe+' '+args.join(' '));
	var P;
	try {
		P=require('child_process').spawn(exe,args);
	}
	catch(err) {
		console.log (err);
		console.log ("Problem launching: "+exe+" "+args.join(" "));
		return "";
	}
	var txt='';
	P.stdout.on('data',function(chunk) {
		txt+=chunk;
	});
	P.on('close',function(code) {
		callback(txt);
	});
}

function serve_file(REQ,filename,response,opts) {
	//the number of bytes to read
	var num_bytes_to_read=opts.end_byte-opts.start_byte+1;

	//check if file exists
	fs.exists(filename,function(exists) {
		if (!exists) {
			response.writeHead(404, {"Content-Type": "text/plain"});
			response.write("404 Not Found\n");
			response.end();
			return;
		}

		//write the header for binary data
		response.writeHead(200, {"Access-Control-Allow-Origin":"*", "Content-Type":"application/octet-stream", "Content-Length":num_bytes_to_read});
		
		//keep track of how many bytes we have read
		var num_bytes_read=0;

		//start reading the file asynchronously
		var read_stream=fs.createReadStream(filename,{start:opts.start_byte,end:opts.end_byte});
		var done=false;
		read_stream.on('data',function(chunk) {
			if (done) return;
			if (!response.write(chunk,"binary")) {
				//if return is false, then we need to pause
				//the reading and wait for the drain event on the write stream
				read_stream.pause();
			}
			//increment number of bytes read
			num_bytes_read+=chunk.length;
			if (num_bytes_read==num_bytes_to_read) {
				//we have read everything, hurray!
				console.log ('Read '+num_bytes_read+' bytes from '+filename+' ('+opts.start_byte+','+opts.end_byte+')');
				done=true;
				//end the response
				response.end();
			}
			else if (num_bytes_read>num_bytes_to_read) {
				//This should not happen
				console.log ('Unexpected error. num_bytes_read > num_bytes_to_read: '+num_bytes_read+' '+num_bytes_to_read);
				read_stream.destroy();
				done=true;
				response.end();
			}
		});
		read_stream.on('end',function() {

		});
		read_stream.on('error',function() {
			console.log ('Error in read stream: '+filename);
		});
		response.on('drain',function() {
			//see above where the read stream was paused
			read_stream.resume();
		});
		REQ.socket.on('close',function() {
			read_stream.destroy(); //stop reading the file because the request has been closed
			done=true;
			response.end();
		});
	});
}

function looks_like_it_could_be_a_file_path(txt) {
	if (txt.indexOf(' ')>=0) return false;
	if (txt.indexOf('http://')==0) return false;
	if (txt.indexOf('https://')==0) return false;
	return true;
}

function looks_like_it_could_be_a_url(txt) {
	if (txt.indexOf('http://')==0) return true;
	if (txt.indexOf('https://')==0) return true;
	return false;
}

function is_valid_checksum(str) {
	if (str.length!=40) return false;
	return /^[a-zA-Z0-9]+$/.test(str);
}

function make_random_id(len)
{
    var text = "";
    var possible = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

    for( var i=0; i < len; i++ )
        text += possible.charAt(Math.floor(Math.random() * possible.length));

    return text;
}

function remove_file(fname) {
	try {
		fs.unlinkSync(fname);
		return true;
	}
	catch(err) {
		return false;
	}
}

function rename_file(fname,fname_new) {
	try {
		fs.renameSync(fname,fname_new);
		return true;
	}
	catch(err) {
		return false;
	}	
}

function get_file_size(fname) {
	try {
		var s=fs.statSync(fname);
		return s.size;
	}
	catch(err) {
		return 0;
	}
}

function compute_file_checksum(fname,callback) {
	run_process_and_read_stdout(prv_exe,['sha1sum',fname],function(txt) {
		callback(txt.trim());
	});
}

function starts_with(str,substr) {
	return (str.slice(0,substr.length)==substr);
}

function CLParams(argv) {
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
				this.namedParameters[arg0]=args[i+1]||'';
				i++;
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
