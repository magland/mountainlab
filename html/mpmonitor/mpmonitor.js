function MPMonitor() {

	this.setClient=function(client) {m_client=client; m_dsv.setClient(client); m_pript_view.setClient(client);};
	this.div=function() {return m_div;};

	var m_dsv=new DaemonStateView();
	var m_pript_view=new PriptView();

	this.loadState=m_dsv.loadState;
	this.setGeom=setGeom;

	var m_div=$('<div></div>');
	m_div.append(m_dsv.div());
	m_div.append(m_pript_view.div());
	
	setGeom(0,0,800,800);
	function setGeom(x,y,W,H) {
		var topmargin=150;
		absgeom(m_div,x,y,W,H);
		m_dsv.setGeom(0,0,W*2/3,H);
		m_pript_view.setGeom(W*2/3,topmargin,W/3,H-topmargin);
	}

	m_div.find('#scripts').css('overflow','scroll');
	m_div.find('#processes').css('overflow','scroll');

	var m_client=null;

	m_dsv.onScriptItemClicked(function(evt,key) {
		m_pript_view.loadScript(key);
	});
	m_dsv.onProcessItemClicked(function(evt,key) {
		m_pript_view.loadProcess(key);
	});

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