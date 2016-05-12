//// requires
var	url=require('url');
var http=require('http');
var fs=require('fs');

///polyfill
if (!String.prototype.endsWith) {
  String.prototype.endsWith = function(searchString, position) {
      var subjectString = this.toString();
      if (typeof position !== 'number' || !isFinite(position) || Math.floor(position) !== position || position > subjectString.length) {
        position = subjectString.length;
      }
      position -= searchString.length;
      var lastIndex = subjectString.indexOf(searchString, position);
      return lastIndex !== -1 && lastIndex === position;
  };
}

//// tasks
var actions={};
actions['queueScript']=require('./queuescript.js').queueScript;
actions['getDaemonState']=require('./getdaemonstate.js').getDaemonState;
actions['clearProcessing']=require('./clearprocessing.js').clearProcessing;
actions['getPript']=require('./getpript').getPript;

//// configuration
var config_fname='config.json';
try {
	var config = JSON.parse(fs.readFileSync(config_fname, 'utf8'));
}
catch(err) {
	console.log('Problem parsing config file. Copy config.example.json to config.json and then modify');
	return;
}

//make the configuration paths absolute!
for (var key in config) {
	if (config[key].indexOf('/')>=0) {
		config[key]=require('path').resolve(__dirname,config[key]);	
	}
}
console.log(JSON.stringify(config));

//// setup
mkdir_if_needed(config.tmp_mpserver_path);
var X=new MPManager();
var last_request_id=0;

//// The manager of all the tasks
function MPManager() {
	// Request is sent here
	this.handleRequest=function(req,callback) {
		var T=initialize_task(req,callback);
		if (!T) return;
		m_tasks[req.request_id]=T;
	}
	// Called when the connection is disconnected by the client
	this.closeRequest=function(request_id) {close_request(request_id);}
	
	// The collection of active tasks to be managed
	var m_tasks={};

	function close_request(request_id) {
		/// if we find the task, let's close it!
		if (m_tasks[request_id]) {
			m_tasks[request_id].close();
			delete m_tasks[request_id];
		}
	}

	function initialize_task(req,callback) {
		if (req.action in actions) {
			return new actions[req.action](config,req,function(resp) {
				callback(resp);
			});
		}
		else {
			console.log('Unrecognized action: '+req.action);
			callback({success:false,error:'Unrecognized action: '+req.action});
			return null;
		}
	}
}

http.createServer(function (REQ, RESP) {
	var url_parts = url.parse(REQ.url,true);	
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
	else if(REQ.method=='GET') {
		var path=url_parts.pathname;
		var query=url_parts.query;
		var suffixes=['html','js','css','json','txt','jpg','png','gif'];
		var suf=path.split('.').pop();
		if (suffixes.indexOf(suf)>=0) {
			var filename = require('path').join(config.html_path, path);
			console.log('SERVING FILE: '+filename);
			serve_file(filename,RESP);
		}
		else if (query.action=='getDaemonState') {
			X.handleRequest({action:query.action},function(resp) {
				send_json_response(resp);
			});
		}
		else {
			send_json_response({success:false,error:'GET not supported for this action'});
		}
	}
	else if(REQ.method=='POST') {
		receive_json_post(REQ,function(req) {
			var request_id=last_request_id+1; last_request_id=request_id;
			req.request_id=request_id;
			REQ.on('close',function() {
				X.closeRequest(request_id);
			});
			X.handleRequest(req,function(resp) {
				send_json_response(resp);	
			});
		});
	}

	function serve_file(filename,response) {
		fs.exists(filename,function(exists) {
			if (!exists) {
				response.writeHead(404, {"Content-Type": "text/plain"});
				response.write("404 Not Found\n");
				response.end();
				return;
			}

			fs.readFile(filename, "binary", function(err, file) {
				if(err) {        
					response.writeHead(500, {"Content-Type": "text/plain"});
					response.write(err + "\n");
					response.end();
					return;
				}

				response.writeHead(200);
				response.write(file, "binary");
				response.end();
			});
		});
	}

	function receive_json_post(REQ,callback) {
		var body='';
		REQ.on('data',function(data) {
			body+=data;
		});
		
		REQ.on('end',function() {
			callback(JSON.parse(body));;
		});
	}
	
	function send_json_response(obj) {
		RESP.writeHead(200, {"Access-Control-Allow-Origin":"*", "Content-Type":"application/json"});
		RESP.end(JSON.stringify(obj));
	}
	function send_text_response(text) {
		RESP.writeHead(200, {"Access-Control-Allow-Origin":"*", "Content-Type":"text/plain"});
		RESP.end(text);
	}
}).listen(config.listen_port);
console.log ('Listening on port '+config.listen_port);

function mkdir_if_needed(path) {
	var fs=require('fs');
	if (!fs.existsSync(path)){
    	fs.mkdirSync(path);
	}
}

