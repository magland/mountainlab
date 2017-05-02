function JObjectManager() {
	this.objectById(id) {
		if (id in m_objects) {
			return m_objects[id];
		}
		else return null;
	}
	var m_objects={};
}
var GlobalJObjectManager=new JObjectManager;

function JObject() {
	this.id=function() {return m_id;}
	this.connect=function(signal_name,receiver_object,signal_or_slot_name) {connect(signal_name,receiver_object,signal_or_slot_name);}
	this.emit=function(signal_name) {emit(signal_name);}
	this.invoke_signal_or_slot=function(signal_or_slot_name) {invoke_signal_or_slot(signal_or_slot_name);}

	var m_id=JUtil.make_random_id();
	var m_signals={};
	GlobalJObjectManager[m_id]=this;

	function connect(signal_name,receiver_object,signal_or_slot_name) {
		if (!(signal_name in m_signals)) {
			m_signals[signal_name]={
				connections:[]
			};
		}
		var SS=m_signals[signal_name];
		var con={
			receiver:receiver_object.id(),
			signal_or_slot_name:signal_or_slot_name
		}
		SS.connections.push(con);
	}
	function emit(signal_name) {
		if (signal_name in m_signals) {
			var SS=m_signals[signal_name];
			for (var i=0; i<SS.length; i++) {
				var receiver_id=SS[i].receiver;
				var receiver_obj=GlobalJObjectManager.objectById(receiver_id);
				if (!receiver_obj) {
					//should clean it up here since this object no longer exists
					return;
				}
				receiver_obj.invoke_signal_or_slot(SS[i].signal_or_slot_name)
			}
		}
	}
	function invok
}

function BasicProperty() {
	this.get=function() {
		return m_value;
	}
	this.set=function(val) {
		if (m_value===val) return;
		m_value=val;
		m_action.trigger();
	}
	this.onChanged(handler) {m_handler.add()}
	
	var m_value;
	var m_action;

}

function MVContext {

	/////////////////////////////////////////////////
    this.clusterMerge=function() {return d.m_cluster_merge;}
    void setClusterMerge(const ClusterMerge& CM);

    /////////////////////////////////////////////////
    QJsonObject clusterAttributes(int num) const;
    QList<int> clusterAttributesKeys() const;
    void setClusterAttributes(int num, const QJsonObject& obj);
    QSet<QString> clusterTags(int num) const; //part of attributes
    QList<QString> clusterTagsList(int num) const;
    QSet<QString> allClusterTags() const;
    void setClusterTags(int num, const QSet<QString>& tags); //part of attributes

    /////////////////////////////////////////////////
    QJsonObject clusterPairAttributes(const ClusterPair& pair) const;
    QList<ClusterPair> clusterPairAttributesKeys() const;
    void setClusterPairAttributes(const ClusterPair& pair, const QJsonObject& obj);
    QSet<QString> clusterPairTags(const ClusterPair& pair) const; //part of attributes
    QList<QString> clusterPairTagsList(const ClusterPair& pair) const;
    QSet<QString> allClusterPairTags() const;
    void setClusterPairTags(const ClusterPair& pair, const QSet<QString>& tags); //part of attributes

    /////////////////////////////////////////////////
    ClusterVisibilityRule clusterVisibilityRule() const;
    QList<int> visibleClusters(int Kmax) const;
    QList<int> visibleClustersIncludingMerges(int Kmax) const;
    bool clusterIsVisible(int k) const;
    void setClusterVisibilityRule(const ClusterVisibilityRule& rule);

    /////////////////////////////////////////////////
    MVEvent currentEvent() const;
    int currentCluster() const;
    QList<int> selectedClusters() const;
    QList<int> selectedClustersIncludingMerges() const;
    double currentTimepoint() const;
    MVRange currentTimeRange() const;
    void setCurrentEvent(const MVEvent& evt);
    void setCurrentCluster(int k);
    void setSelectedClusters(const QList<int>& ks);
    void setCurrentTimepoint(double tp);
    void setCurrentTimeRange(const MVRange& range);
    void clickCluster(int k, Qt::KeyboardModifiers modifiers);

    /////////////////////////////////////////////////
    QSet<ClusterPair> selectedClusterPairs() const;
    void setSelectedClusterPairs(const QSet<ClusterPair>& pairs);
    void clickClusterPair(const ClusterPair& pair, Qt::KeyboardModifiers modifiers);

    /////////////////////////////////////////////////
    QColor clusterColor(int k) const;
    QColor channelColor(int m) const;
    QColor color(QString name, QColor default_color = Qt::black) const;
    QMap<QString, QColor> colors() const;
    QList<QColor> channelColors() const;
    QList<QColor> clusterColors() const;
    void setClusterColors(const QList<QColor>& colors);
    void setChannelColors(const QList<QColor>& colors);
    void setColors(const QMap<QString, QColor>& colors);

    /////////////////////////////////////////////////
    DiskReadMda currentTimeseries() const;
    QString currentTimeseriesName() const;
    QStringList timeseriesNames() const;
    void addTimeseries(QString name, DiskReadMda32 timeseries);
    void setCurrentTimeseriesName(QString name);

    /////////////////////////////////////////////////
    DiskReadMda firings();
    void setFirings(const DiskReadMda& F);
    MVEventFilter eventFilter();
    void setEventFilter(const MVEventFilter& EF);

    /////////////////////////////////////////////////
    // these should be set once at beginning
    double sampleRate() const;
    void setSampleRate(double sample_rate);
    QString mlProxyUrl() const;
    void setMLProxyUrl(QString url);

    /////////////////////////////////////////////////
    QVariant option(QString name, QVariant default_val = QVariant());
    void setOption(QString name, QVariant value);
    void onOptionChanged(QString name, const QObject* receiver, const char* member, Qt::ConnectionType type = Qt::DirectConnection);

    /////////////////////////////////////////////////
    void copySettingsFrom(MVContext* other);	
}
