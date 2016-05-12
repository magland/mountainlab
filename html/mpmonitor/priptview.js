function PriptView() {
	this.setClient=function(client) {m_client=client; debug_init();};
	this.loadScript=load_script;
	this.loadProcess=load_process;
	this.div=function() {return m_div;};
	this.setGeom=setGeom;

	var m_div=$('#templates #priptview').clone();
	m_div.find('textarea').attr('readonly','readonly');
	
	setGeom(800,800);
	function setGeom(x,y,W,H) {
		absgeom(m_div,x,y,W,H);
		m_div.find('textarea').css({width:W-30,height:200});
	}

	m_div.css('overflow','scroll');

	var m_client=null;
	var m_pript={};

	refresh();

	function debug_init() {
		//debugging
		load_script('ZZ0K5EE8ML0ZD4AX7F64');
	}
	

	function load_script(key) {
		if (!m_client) return;
		m_pript={};
		refresh();
		m_client.getScript(key,function(resp) {
    		if (resp.success) {
    			m_pript=resp.script;
        		refresh();
        	}
        	else {
        		console.log('Error loading script: '+resp.error);
        		m_pript={error:'ERR: '+resp.error};
        		refresh();
        	}
    	});
	}

	function load_process(key) {
		if (!m_client) return;
		m_pript={};
		refresh();
		m_client.getProcess(key,function(resp) {
    		if (resp.success) {
    			m_pript=resp.process;
        		refresh();
        	}
        	else {
        		console.log('Error loading process: '+resp.error);
        		m_pript={error:'ERR: '+resp.error};
        		refresh();
        	}
    	});
	}

	function reload() {
		var id=m_pript.id||'';
		var prtype=m_pript.prtype||'';
		m_pript={};
		refresh();
		if (prtype=='script')
			load_script(id);
		else
			load_process(id);
	}

	function refresh() {
		if (!m_pript.id) {
			m_div.css('visibility','hidden');
		}
		else {
			var status='not running';
			if (m_pript.is_running) status='running';
			else if (m_pript.is_finished) {
				status='finished';
				if (m_pript.error) status='finished with error';
			}

			var a0=$('<a href=# title="Click to refresh">'+m_pript.id+'</a>');
			a0.click(reload);
			m_div.css('visibility','visible');
			m_div.find('#pript_id').empty();
			m_div.find('#pript_id').append(a0);


			m_div.find('#type').html(m_pript.prtype);
			if (m_pript.prtype=='script') {
				var str=get_script_file_name_list(m_pript);
				m_div.find('#name').html(str);
				
			}
			else {
				m_div.find('#name').html(m_pript.processor_name||'<>');
			}
			m_div.find('#status').html(status);
			m_div.find('#raw').val('RAW:\n'+JSON.stringify(m_pript));
			m_div.find('#stdout').val('STDOUT:\n'+get_stdout(m_pript));
		}
	}
	function get_stdout(P) {
		var rr=P.runtime_results||{};
		return rr.stdout||'';
	}
	function get_script_file_name_list(S) {
		var paths=S.script_paths||[];
		var ret=[];
		for (var i=0; i<paths.length; i++) {
			ret.push(get_file_name(paths[i]));
		}
		return ret.join(', ');
	}
	function get_file_name(path) {
		var ind=path.lastIndexOf('/');
		if (ind<0) return path;
		return path.slice(ind+1);
	}
	function absgeom(div,x,y,w,h) {
		div.css({
			position:'absolute',
			width:w,
			height:h,
			left:x,
			top:y
		});
	}
}