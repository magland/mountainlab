#include "sstimeserieswidget_prev.h"
#include <QList>
#include "extractclipsdialog.h"
#include "ssabstractview_prev.h"
#include "sstimeseriesview_prev.h"
#include "sscontroller.h"
#include <QSplitter>
#include <QDebug>
#include <QVBoxLayout>
#include <QLabel>
#include <QMenuBar>
#include <QMenu>
#include <QMessageBox>
#include <QProcess>
#include <QInputDialog>
#include <QFileInfo>
#include <QApplication>
#include <QDir>
#include <QTime>
#include <QFileDialog>
#include <QObjectList>
#include <QTextBrowser>
#include "cvcombowidget.h"
#include "cvcommon.h"
#include "mlutils.h"

void do_wait(int ms) {
	QTime timer; timer.start();
	while (timer.elapsed()<ms) {
		qApp->processEvents();
	}
}

class SSTimeSeriesWidgetPrivate {
public:
	SSTimeSeriesWidget *q;

	QList<SSAbstractView *> m_views;
	SSAbstractView *m_current_view;
	QSplitter m_splitter;
	QLabel m_info;
	double m_sampling_frequency; //Hz
	QString m_widget_type;
	ExtractClipsDialog m_extract_clips_dlg;
	QMenuBar *m_menubar;
    Mda m_clip_data;

	void set_info(QString txt);
	void connect_view(SSAbstractView *V,bool current_timepoint_only);
	void update_info(SSAbstractView *V);
	void tell_current_view();
	void update_menu_bar();
};

SSTimeSeriesWidget::SSTimeSeriesWidget(QWidget *parent) : QWidget(parent) {

	d=new SSTimeSeriesWidgetPrivate;
	d->q=this;

	d->m_current_view=0;
	d->m_sampling_frequency=0; //eg 20000
	d->m_widget_type="raw";

	d->m_splitter.setStyleSheet("QSplitter {background: gray;}");
	d->m_splitter.setOrientation(Qt::Vertical);

	d->m_info.setFixedHeight(30);
	d->m_info.setFont(QFont("Arial",10));
	d->m_info.setMargin(20);

	QWidget *CW=new QWidget();
	QVBoxLayout *CL=new QVBoxLayout();
	CL->setMargin(0); CL->setContentsMargins(4,4,4,4);
	CL->setSpacing(0);
	CW->setLayout(CL);

	CL->addWidget(&d->m_splitter);
	CL->addWidget(&d->m_info);

	QVBoxLayout *vlayout=new QVBoxLayout;
	vlayout->setContentsMargins(0,0,0,0);
	vlayout->setSpacing(0);
	setLayout(vlayout);

	setAttribute(Qt::WA_DeleteOnClose);

	//menu
	{
		QMenuBar *menubar=new QMenuBar(this);
		d->m_menubar=menubar;
		vlayout->addWidget(menubar);
		{ //Tools
			QMenu *menu=new QMenu("Tools",menubar);
			menubar->addMenu(menu);
			{
				QAction *A=new QAction("Extract clips",menu);
				menu->addAction(A);
				A->setObjectName("extract_clips");
				A->setShortcut(QKeySequence("Ctrl+E"));
				connect(A,SIGNAL(triggered()),this,SLOT(slot_extract_clips()));
			}
			{
				QAction *A=new QAction("Extract comparison clips",menu);
				menu->addAction(A);
				A->setObjectName("extract_comparison_clips");
				connect(A,SIGNAL(triggered()),this,SLOT(slot_extract_comparison_clips()));
			}
			{
				QAction *A=new QAction("View clusters",menu);
				menu->addAction(A);
				A->setObjectName("view_clusters");
				A->setShortcut(QKeySequence("Ctrl+C"));
				connect(A,SIGNAL(triggered()),this,SLOT(slot_view_clusters()));
			}
		}
		{ //Instructions
			QMenu *menu=new QMenu("Instructions",menubar);
			menubar->addMenu(menu);
			{
				QAction *A=new QAction("Navigation Instructions",menu);
				menu->addAction(A);
				A->setObjectName("navigation_instructions");
				A->setShortcut(QKeySequence("Ctrl+I"));
				connect(A,SIGNAL(triggered()),this,SLOT(slot_navigation_instructions()));
			}
		}
	}

	vlayout->addWidget(CW);

}

SSTimeSeriesWidget::~SSTimeSeriesWidget() {
	qDeleteAll(d->m_views);
	delete d;
}

void SSTimeSeriesWidget::addView(SSAbstractView *V) {

	if (d->m_views.isEmpty()) V->setTimelineVisible(true);
	else V->setTimelineVisible(false);

	d->m_views.append(V);
	d->m_splitter.addWidget(V);
	d->connect_view(V,false);
	V->initialize();

	if (!d->m_current_view) {
		d->m_current_view=V;
		d->update_menu_bar();
	}
    d->tell_current_view();
}

void SSTimeSeriesWidget::setClipData(const Mda &X)
{
    d->m_clip_data=X;
}

void SSTimeSeriesWidget::hideMenu()
{
    d->m_menubar->hide();
}

SSAbstractView *SSTimeSeriesWidget::view(int index)
{
    if (index>=d->m_views.count()) return 0;
    if (index<0) return 0;
    return d->m_views[index];
}

void SSTimeSeriesWidgetPrivate::set_info(QString txt) {
	m_info.setText(txt);
}
void SSTimeSeriesWidgetPrivate::connect_view(SSAbstractView *V,bool current_timepoint_only) {
	q->connect(V,SIGNAL(currentXChanged()),q,SLOT(slot_current_x_changed()));
	if (!current_timepoint_only) {
		q->connect(V,SIGNAL(currentChannelChanged()),q,SLOT(slot_current_channel_changed()));
		q->connect(V,SIGNAL(xRangeChanged()),q,SLOT(slot_x_range_changed()));
		q->connect(V,SIGNAL(replotNeeded()),q,SLOT(slot_replot_needed()));
		q->connect(V,SIGNAL(selectionRangeChanged()),q,SLOT(slot_selection_range_changed()));
		q->connect(V,SIGNAL(clicked()),q,SLOT(slot_view_clicked()));
	}
}
void SSTimeSeriesWidget::slot_current_x_changed() {
	SSAbstractView *V=(SSAbstractView *)sender();
	d->update_info(V);

	foreach (SSAbstractView *V2,d->m_views) {
		if (V!=V2) {
			if (d->m_views.indexOf(V)>=0) V2->disableSignals();
			V2->setCurrentTimepoint(V->currentTimepoint());
			if (d->m_views.indexOf(V)>=0) V2->enableSignals();
		}
	}
}
void SSTimeSeriesWidget::slot_current_channel_changed() {
	SSAbstractView *V=(SSAbstractView *)sender();
	d->update_info(V);
}
void SSTimeSeriesWidget::slot_x_range_changed() {
	SSAbstractView *V=(SSAbstractView *)sender();
	d->update_info(V);

	foreach (SSAbstractView *V2,d->m_views) {
		if (V!=V2) {
			V2->disableSignals();
			V2->setXRange(V->xRange());
			V2->enableSignals();
		}
	}
}

void SSTimeSeriesWidget::slot_replot_needed()
{
	/*
	SSAbstractView *V=(SSAbstractView *)sender();

	foreach (SSAbstractView *V2,d->m_views) {
		if (V!=V2) {
			V2->disableSignals();
			V2->setVerticalZoomFactor(V->verticalZoomFactor());
			V2->enableSignals();
		}
	}
	*/
}

void SSTimeSeriesWidget::slot_selection_range_changed()
{
	SSAbstractView *V=(SSAbstractView *)sender();
	d->update_info(V);

	foreach (SSAbstractView *V2,d->m_views) {
		if (V!=V2) {
			V2->disableSignals();
			V2->setSelectionRange(V->selectionRange());
			V2->enableSignals();
		}
	}
}

void SSTimeSeriesWidget::slot_view_clicked()
{
	SSAbstractView *W=(SSAbstractView *)sender();
	if (d->m_current_view!=W) {
		d->m_current_view=W;
		d->tell_current_view();
		d->update_menu_bar();
	}
}

int get_num_clips_from_process_output(QString output) {
	qDebug()  << output;
	int ind1=output.indexOf("Number of clips: ");
	if (ind1<0) return 0;
	output=output.mid(ind1+QString("Number of clips: ").count());
	int ind2=output.indexOf("\n");
	if (ind2<0) return 0;
	output=output.mid(0,ind2).trimmed();
	return output.toInt();
}

void SSTimeSeriesWidget::slot_extract_clips()
{
	SSAbstractView *V=d->m_current_view;
	if (!V) {
		QMessageBox::critical(this,"Extract clips","No view selected.");
		return;
	}
	if (V->viewType()!="SSTimeSeriesView") {
		QMessageBox::critical(this,"Extract clips","You must select a time series view.");
		return;
	}

	ExtractClipsDialog *dlg=&d->m_extract_clips_dlg;
	dlg->setComparisonMode(false);
	dlg->setModal(true);
	if (dlg->exec()!=QDialog::Accepted) {
		return;
	}
	QString labels_txt=dlg->param("labels").toString();

	if (labels_txt.isEmpty()) return;
	if (labels_txt.contains(" ")) {
		QMessageBox::critical(this,"Extract clips","Spaces are not allowed");
		return;
	}

	int clipsize=dlg->param("clipsize").toInt();
	if ((clipsize<=3)||(clipsize>=1000)) {
		QMessageBox::critical(this,"Extract clips","Invalid clip size");
		return;
	}

	bool fixed_clipsize=dlg->param("fixed_clipsize").toBool();

	QString save_to_path;
	if (dlg->param("save_to_disk").toBool()) {
		QString filter="*.mda";
		save_to_path=QFileDialog::getSaveFileName(0,"Save clips to file","",filter,&filter);
		if (save_to_path.isEmpty()) return;
		if (QFileInfo(save_to_path).suffix()!="mda") save_to_path+=".mda";
	}

	SSTimeSeriesView *VV=(SSTimeSeriesView *)V;
	DiskArrayModel *DD=VV->data();
	QString data_path=DD->path();
	SSLabelsModel *TL=VV->getLabels();
    MemoryMda TL0;
	if (TL) TL0=TL->getTimepointsLabels(0,TL->getMaxTimepoint());
    else qWarning() << "TL is null in slot_extract_clips";
	QString TL_path=ssTempPath()+"/spikespy-tmp-extract-clips-TL.mda";
	removeOnClose(TL_path);
	TL0.write(TL_path.toLatin1().data());
	QString out_path=QString(ssTempPath()+"/spikespy-tmp-extract-clips-out-%1.mda").arg(qAbs(qrand()));
	removeOnClose(out_path);
	QString out_path_TM=QFileInfo(out_path).path()+"/"+QFileInfo(out_path).completeBaseName()+".TM.mda";
	removeOnClose(out_path_TM);
	QString out_path_TL=QFileInfo(out_path).path()+"/"+QFileInfo(out_path).completeBaseName()+".TL.mda";
	removeOnClose(out_path_TL);

	QStringList args; args << data_path << TL_path << out_path << QString("--clipsize=%1").arg(clipsize) << "--labels="+labels_txt;
    if (fixed_clipsize) args << "--fixed-clipsize";
	QString exe=cfp(qApp->applicationDirPath()+"/extractclips");
	if (!QFile::exists(exe)) {
		QString exe2=cfp(qApp->applicationDirPath()+"/../src/spikespy/bin/extractclips");
		if (QFile::exists(exe2)) exe=exe2;
    }
	qDebug()  << exe << args;
	QProcess process;
	process.start(exe,args);
	process.waitForFinished();
	process.waitForReadyRead();
	QString output=QString(process.readAll());
	int num_clips=get_num_clips_from_process_output(output);
    do_wait(100); //may be necessary for files to finish writing (probably not)

	if (save_to_path.isEmpty()) {
        DiskReadMdaOld out_TM(out_path_TM);
		SSController CC;
		SSTimeSeriesWidget *WWW=(SSTimeSeriesWidget *)CC.createTimeSeriesWidget();
		SSTimeSeriesView *VVV=(SSTimeSeriesView *)CC.createTimeSeriesView();
		VVV->setData((DiskArrayModel *)CC.loadArray(out_path),true); //added true on 11/24/15
		VVV->setProperty("data_path",out_path); //a hack!
		VVV->setProperty("TM_data_path",out_path_TM); //a hack!
		VVV->setTimepointMapping(out_TM);
		//VVV->setConnectZeros(false);
		VVV->setClipMode(true);
		VVV->setProperty("fixed_clipsize",fixed_clipsize);
		connect(VVV,SIGNAL(requestCenterOnCursor()),this,SLOT(slot_center_on_cursor()));
		WWW->addView(VVV);
		d->connect_view(VVV,true);
		QString title=QString("%1 - Labels: %2; %3 clips").arg(V->title()).arg(labels_txt).arg(num_clips);
		WWW->topLevelWidget()->setWindowTitle(title);
	}
	else {
		QFile::remove(save_to_path);
		if (!QFile::rename(out_path,save_to_path)) {
			QMessageBox::critical(0,"Save clips to file",
								  QString("Unable to save clips to file. %1->%2.").arg(out_path).arg(save_to_path)
								  );
		}
		else {
			QMessageBox::critical(0,"Save clips to file","Clips saved to file.");
		}
	}
}

void SSTimeSeriesWidget::slot_extract_comparison_clips()
{
	SSAbstractView *V=d->m_current_view;
	SSAbstractView *V2=0;
	if (!V) {
		QMessageBox::critical(this,"Extract clips","No view selected.");
		return;
	}
	if (V->viewType()!="SSTimeSeriesView") {
		QMessageBox::critical(this,"Extract clips","You must select a time series view.");
		return;
	}
	for (int i=0; i<d->m_views.count(); i++) {
		if (d->m_views[i]==V) {
			if (i+1<d->m_views.count()) V2=d->m_views[i+1];
		}
	}
	if (!V2) {
		QMessageBox::critical(this,"Extract clips","The view below the selected view must also be a time series view.");
		return;
	}
	if (V2->viewType()!="SSTimeSeriesView") {
		QMessageBox::critical(this,"Extract clips","The view below the selected view must be a time series view.");
		return;
	}

	ExtractClipsDialog *dlg=&d->m_extract_clips_dlg;
	dlg->setComparisonMode(true);
	dlg->setModal(true);
	if (dlg->exec()!=QDialog::Accepted) {
		return;
	}
	QString labels_txt=dlg->param("labels").toString();

	if (labels_txt.isEmpty()) return;
	if (labels_txt.contains(" ")) {
		QMessageBox::critical(this,"Extract clips","Spaces are not allowed");
		return;
	}
	if (labels_txt.split(",").count()!=2) {
		QMessageBox::critical(this,"Extract clips","Exactly two label numbers required!");
		return;
	}

	int clipsize=dlg->param("clipsize").toInt();
	if ((clipsize<=3)||(clipsize>=1000)) {
		QMessageBox::critical(this,"Extract clips","Invalid clip size");
		return;
	}

	QString save_to_path;
	if (dlg->param("save_to_disk").toBool()) {
		QString filter="*.mda";
		save_to_path=QFileDialog::getSaveFileName(0,"Save clips to file","",filter,&filter);
		if (save_to_path.isEmpty()) return;
		if (QFileInfo(save_to_path).suffix()!="mda") save_to_path+=".mda";
	}

	SSTimeSeriesView *VV1=(SSTimeSeriesView *)V;
	DiskArrayModel *DD1=VV1->data();
	QString data_path1=DD1->path();
	SSLabelsModel *TL1=VV1->getLabels();
	if (!TL1) {
		QMessageBox::critical(this,"Extract clips","No labels found for first view");
		return;
	}
	MemoryMda TL1_0=TL1->getTimepointsLabels(0,TL1->getMaxTimepoint());
	QString TL1_path=ssTempPath()+"/spikespy-tmp-extract-clips-TL1.mda";
	TL1_0.write(TL1_path.toLatin1().data());
	QString out_path1=QString(ssTempPath()+"/spikespy-tmp-extract-clips-out1-%1.mda").arg(qAbs(qrand()));
	removeOnClose(out_path1);

	SSTimeSeriesView *VV2=(SSTimeSeriesView *)V2;
	DiskArrayModel *DD2=VV2->data();
	QString data_path2=DD2->path();
	SSLabelsModel *TL2=VV2->getLabels();
	if (!TL2) {
		QMessageBox::critical(this,"Extract clips","No labels found for second view");
		return;
	}
	MemoryMda TL2_0=TL2->getTimepointsLabels(0,TL2->getMaxTimepoint());
	QString TL2_path=ssTempPath()+"/spikespy-tmp-extract-clips-TL2.mda";
	removeOnClose(TL2_path);
	TL2_0.write(TL2_path.toLatin1().data());
	QString out_path2=QString(ssTempPath()+"/spikespy-tmp-extract-clips-out2-%1.mda").arg(qAbs(qrand()));

	QString out_path_TM=QFileInfo(out_path1).path()+"/"+QFileInfo(out_path1).completeBaseName()+".TM.mda";

	QStringList args; args << data_path1 << data_path2 << TL1_path << TL1_path << out_path1 << out_path2 << QString("--clipsize=%1").arg(clipsize) << "--labels="+labels_txt;
	QString exe=qApp->applicationDirPath()+"/extractclips2";
	if (!QFile::exists(exe)) {
		QString exe2=qApp->applicationDirPath()+"/../src/spikespy/bin/extractclips2";
		if (QFile::exists(exe2)) exe=exe2;
	}
	qDebug()  << exe << args;
	QProcess process;
	process.start(exe,args);
	process.waitForFinished();
	process.waitForReadyRead();
	QString output=QString(process.readAll());
	int num_clips=get_num_clips_from_process_output(output);
	do_wait(100); //may be necessary for files to finish writing (probably not)

	if (save_to_path.isEmpty()) {
        DiskReadMdaOld out_TM(out_path_TM);
		SSController CC;
		SSTimeSeriesWidget *WWW=(SSTimeSeriesWidget *)CC.createTimeSeriesWidget();

		SSTimeSeriesView *VVV1=(SSTimeSeriesView *)CC.createTimeSeriesView();
		VVV1->setData((DiskArrayModel *)CC.loadArray(out_path1),true);
		VVV1->setTimepointMapping(out_TM);
		//VVV1->setConnectZeros(false);
		VVV1->setClipMode(true);
		WWW->addView(VVV1);
		d->connect_view(VVV1,true);

		SSTimeSeriesView *VVV2=(SSTimeSeriesView *)CC.createTimeSeriesView();
		VVV2->setData((DiskArrayModel *)CC.loadArray(out_path2),true); //added on 11/24/15
		VVV2->setTimepointMapping(out_TM);
		//VVV2->setConnectZeros(false);
		VVV2->setClipMode(true);
		WWW->addView(VVV2);
		d->connect_view(VVV2,true);

		QString title=QString("%1 - %2 - Labels: %3; %4 clips").arg(VV1->title()).arg(VV2->title()).arg(labels_txt).arg(num_clips);
		WWW->topLevelWidget()->setWindowTitle(title);
	}
	else {
		QFile::remove(save_to_path);
		if (!QFile::rename(out_path1,save_to_path)) {
			QMessageBox::critical(0,"Save clips to file",
								  QString("Unable to save clips to file. %1->%2.").arg(out_path1).arg(save_to_path)
								  );
		}
		else {
			QMessageBox::critical(0,"Save clips to file","Clips saved to file.");
		}
	}
}

void SSTimeSeriesWidget::slot_center_on_cursor()
{
	d->m_current_view->centerOnCursor();
}

void SSTimeSeriesWidget::slot_view_clusters()
{

    if (!d->m_current_view) return;
    CVComboWidget *W=new CVComboWidget();
    W->setAttribute(Qt::WA_DeleteOnClose);
    QString clip_data_path=d->m_current_view->property("data_path").toString();
    QString clip_TM_data_path=d->m_current_view->property("TM_data_path").toString();

    if (!clip_data_path.isEmpty()) {
        if (clip_data_path.isEmpty()) {
            qWarning() << "clip_data_path is empty.";
            return;
        }
        DiskReadMdaOld X(clip_data_path); //a hack!
        W->setClips(X);
        DiskReadMdaOld TM(clip_TM_data_path);
        W->setTimepointMapping(TM);
    }
    else {
        DiskReadMdaOld X(d->m_clip_data);
        W->setClips(X);
        DiskReadMdaOld TM(clip_TM_data_path);
    }
    W->show();
    W->resize(500,500);
    W->setWindowTitle(this->windowTitle());
	d->connect_view(W->timeSeriesView(),true);
}

void SSTimeSeriesWidget::slot_navigation_instructions()
{

	QString html;
	html+="<h3>SPIKESPY</h3>\n";
	html+="<ul>\n";
	html+="<li>Mousewheel scroll: zoom in, zoom out</li>\n";
	html+="<li>+/- keys also zoom in, zoom out</li>\n";
	html+="<li>up/down arrow keys do vertical zoom (on all views together)</li>\n";
	html+="<li>Click and drag to pan</li>\n";
	html+="<li>\"0\" to zoom back to full view</li>\n";
	html+="<li>Left/right arrow keys scroll left/right by one page</li>\n";
	html+="<li>Ctrl + left/right arrow keys move cursor left/right</li>\n";
	html+="<li>Right-click and drag to select a time interval</li>\n";
	html+="<li>ENTER to zoom in to the selected time interval</li>\n";
	//html+="<li>F (when clicked on a view) flips the channel ordering</li>\n";
	html+="</ul>\n";

	QTextBrowser *W=new QTextBrowser;
	W->setAttribute(Qt::WA_DeleteOnClose);
	W->setHtml(html);
	W->show();
	W->resize(800,800);
	W->move(this->geometry().topLeft()+QPoint(100,100));
}

void SSTimeSeriesWidgetPrivate::update_info(SSAbstractView *V) {
	QString str="";
	Vec2 SR=V->selectionRange();
	if (SR.x>=0) {
		str=QString("Duration = %1").arg(V->getTimeLabelForX(SR.y-SR.x));
	}
	else {
		double x0=V->currentX();
		double t0=V->currentTimepoint();
		if ((x0>=0)&&(t0>=0)) {
			str=QString("Time = %1 (%2)").arg(V->getTimeLabelForX(x0)).arg((long)t0+1);

			//str+=QString("; Value = %1  ").arg(V->currentValue());
		}
	}
	double tmp_mb_1=((int)jbytesallocated())*1.0/1000000;
	double tmp_mb_2=((int)jnumbytesread())*1.0/1000000;
	str+=QString("       %1/%2 MB used/read").arg((int)tmp_mb_1).arg((int)tmp_mb_2);

	set_info(str);
}

void SSTimeSeriesWidgetPrivate::tell_current_view()
{
	for (int i=0; i<m_views.count(); i++) {
		if (m_views[i]==m_current_view) {
			m_views[i]->setActivated(true);
		}
		else {
			m_views[i]->setActivated(false);
		}
	}
}

void set_action_enabled(QObject *par,QString name,bool val) {
	QObjectList children=par->children();
	foreach (QObject *obj,children) {
		QAction *A=obj->findChild<QAction *>(name);
		if (A) {
			A->setEnabled(val);
		}
		set_action_enabled(obj,name,val);
	}
}

void SSTimeSeriesWidgetPrivate::update_menu_bar()
{
	if (!m_current_view) return;
	bool clip_mode=false;
	bool time_series_view=false;
	bool fixed_clipsize=m_current_view->property("fixed_clipsize").toBool();
	if (m_current_view->viewType()=="SSTimeSeriesView") {
		time_series_view=true;
		clip_mode=((SSTimeSeriesView *)m_current_view)->clipMode();
	}
	if (clip_mode) {
		set_action_enabled(m_menubar,"extract_clips",false);
		set_action_enabled(m_menubar,"extract_comparison_clips",false);
		set_action_enabled(m_menubar,"view_clusters",fixed_clipsize);
	}
	else {
		set_action_enabled(m_menubar,"extract_clips",time_series_view);
		set_action_enabled(m_menubar,"extract_comparison_clips",time_series_view);
		set_action_enabled(m_menubar,"view_clusters",false);
	}
}
