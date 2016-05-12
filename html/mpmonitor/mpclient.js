function MPClient(url) {
	this.getDaemonState=function(callback) {
        $.ajax({
            type:'POST',
            url:url,
            dataType:'json',
            data:JSON.stringify({action:'getDaemonState'}),
            success:function(resp) {
                if (resp.success) {
                    callback(resp);
                }
                else {
                    console.log('Problem getting daemon state: '+resp.error);
                    callback(resp);;
                }
            },
            error:function() {
            	console.log('Failed to get daemon state: '+url);
            	callback({success:false,error:'Failed to get daemon state. Most likely we are unable to connect to server.'});
            }
        });
	};
	this.clearProcessing=function() {
    	$.ajax({
			type:'POST',
			url:url,
			dataType:'json',
			data:JSON.stringify({action:'clearProcessing'}),
			success:on_success,
			error:on_error
		});
    	function on_success() {alert('Processing cleared');}
    	function on_error() {alert('Error clearing processing');}
    }

    this.queueScript=function(scripts,opts,callback) {
		$.ajax({
			type:'POST',
			url:url,
			dataType:'json',
			data:JSON.stringify({action:'queueScript',scripts:scripts,detach:opts.detach}),
			success:callback
		});
	}
	this.getScript=function(id,callback) {
		$.ajax({
			type:'POST',
			url:url,
			dataType:'json',
			data:JSON.stringify({action:'getPript',id:id,prtype:'script'}),
			success:function(resp) {
				try {
					resp.script=JSON.parse(resp.stdout);
				}
				catch(parse_err) {
					resp.script={};
				}
				callback(resp);
			},
			error:function() {
				console.log('Error in POST of getScript.');
			}
		});	
	}
	this.getProcess=function(id,callback) {
		$.ajax({
			type:'POST',
			url:url,
			dataType:'json',
			data:JSON.stringify({action:'getPript',id:id,prtype:'process'}),
			success:function(resp) {
				try {
					resp.process=JSON.parse(resp.stdout);
				}
				catch(parse_err) {
					resp.process={};
				}
				callback(resp);
			},
			error:function() {
				console.log('Error in POST of getProcess.');
			}
		});	
	}
}